#!/usr/bin/perl -w

###The standard methyratio.py script that comes with the BSmap program is pretty memory and speed inefficient, mainly due to the fact that it takes a non-sorted input,
###and must load the entire genome into memory. This alternate version processes a sorted bsmap .bam file on the fly which improves memory and speed, also referecne sequence is included in the
###BSmap bam output with the -R option.

use strict;
use Getopt::Long;


# How to execute this script - Getopt options straight from Jixian Zhai
my $USAGE = "\nUSAGE: methratio_alt.pl 
                                   --input /BSmap/output/file.bam
				   --output_file    /ratio/file/location.ratio
                                   ";
                                 

my $input;
my $output;
my $log_file;



my @coverage;
my @conversion_array;#matches,total,rc_matches,rc_total,strand,context
my @methylation_report;

my $line;

my $chromosome="new";
my $low_position=0;

my @data_point;
my $BS_conversion="CG:T:GC:A";
my ($cr, $pos, $cigar, $seq, $strand, $insert);
my $ref;
my ($match, $convert, $rc_match, $rc_convert,$search,$rc_search);
my @matches;
my @asymmetric_matches;
my $match_data;
my $context;
my @seq_array;
my $seq_pos;
my $export_pos;
my $shift_amount;

my $repeat="FALSE";
my ($meth_percent,$meth_adjust);
my $effective_CT_counts;
my $meth_ratio;

my $read_count=0;
my $alignment_count=0;
my $cytosine_count=0;
my $average_coverage=0;

my $debug;
my ($f_match,$r_match);
  
GetOptions('input=s'=>\$input,
	   'output_file=s'=>\$output
	   );
die $USAGE unless (defined($input) && defined($output));

if ($output=~/(.+)\.\w+/){
    $log_file=$1.".log";
}

open LOG, ">$log_file" or die "Cannot open $log_file\n";
print LOG "#methratio_alt.pl log file\n";
print LOG "Input: $input\n";
print LOG "Output: $output\n";
print LOG "Start time: ".localtime(time)."\n";

open ALIGNMENTS, "samtools view $input |" or die "Cannot open pipe to $input\n";


open OUTPUT, ">$output" or die "Cannot open $output\n";
print OUTPUT "chr\tpos\tstrand\tcontext\tratio\teff_CT_count\tC_count\tCT_count\trev_G_count\trev_GA_count\n";

while(!eof(ALIGNMENTS)){
    if($repeat eq "TRUE"){
	$repeat="FALSE";
    }else{
	$line=<ALIGNMENTS>;
	$read_count++;
	if ($read_count%1000000==0){
	    print "\t\t$read_count reads processed\n";
	}
    }
    @data_point=split /\t/, $line;
    if (!($data_point[1] & 4)){
	$alignment_count++;
	($cr, $pos, $cigar, $seq, $strand, $insert) = ($data_point[2], $data_point[3], $data_point[5], $data_point[9], '', $data_point[8]);
	if($chromosome eq "new"){
	    $chromosome=$cr;
	}
	if($cr eq $chromosome){
	    if($data_point[12]=~/XR:Z:\w\w([A-Z]+)/){
		$ref=$1;
	    }
	    #$ref=$data_point[12];
	    #$ref=substr($ref,2,length($ref)-4); #remove soft clipped bp
	    $strand=$data_point[15];
	    if(($pos>$low_position) && (length(scalar(@conversion_array))>0)){
		$shift_amount=($pos-$low_position)-1;
		for $seq_pos (0..$shift_amount){
		    if (defined($conversion_array[0])){
			@methylation_report=@{shift @conversion_array};
			#print "woo\n";
			#print "@methylation_report\n";
			if(defined($methylation_report[4])){
			    if(!defined($methylation_report[2])){
				$methylation_report[2]=0;
			    }
			    
			    if(!defined($methylation_report[3])){
				$meth_adjust=-1;
				$methylation_report[3]=0;
			    }else{
				$meth_adjust=$methylation_report[2]/$methylation_report[3];
			    };
			    
			    if(!defined($methylation_report[0])){
				$methylation_report[0]=0;
			    }
			    
			    if(!defined($methylation_report[1])){
				$meth_percent=-1.000;
				$methylation_report[1]=0;
			    }else{
				if($meth_adjust>0){
				    $effective_CT_counts=$methylation_report[1]*$meth_adjust;
				}else{
				    $effective_CT_counts=$methylation_report[1];
				}
				$meth_percent=$methylation_report[0]/$effective_CT_counts;
				$meth_percent*=1000;
				$meth_percent=int($meth_percent);
				$meth_percent/=1000;
				$export_pos=$seq_pos+$low_position;
				$cytosine_count++;
				$average_coverage+=$effective_CT_counts;
				print OUTPUT "$chromosome\t$export_pos\t$methylation_report[4]\t$methylation_report[5]\t$meth_percent\t$effective_CT_counts\t$methylation_report[0]\t$methylation_report[1]\t$methylation_report[2]\t$methylation_report[3]\n";
			    };
    			};
		    }else{
			shift @conversion_array;
		    }
		};
		$low_position=$pos;
	    };
	    $coverage[$pos-$low_position]+=1;
	    if($strand=~/ZS:Z:(\S{2})/){
		$strand=$1;
	    }else{
		print 'missing strand information "ZS:Z:xx"';
	    }
	    if($data_point[1] & 16){
		$strand="-";
		($match, $convert,$rc_match, $rc_convert)=split /:/, "G:A:C:T";
	    }else{
		$strand="+";
		($match, $convert,$rc_match, $rc_convert)=split /:/, "C:T:G:A";
	    }
	    
	    
	    
	    #print "$match,$convert,$rc_match,$rc_convert\n";
	    while ($ref=~/CG/gi){
		$f_match=$-[0];
		$r_match=$f_match+1;
		if($strand eq "-"){
		    push @matches, $r_match;###@+ stores the location of the beginning of each match
		}else{
		    push @matches, $f_match;###@- stores the location of the beginning of each match
		}
	    }
	    @seq_array = split //, $seq;
	    if(scalar(@matches)>0){
		foreach $seq_pos (@matches){
		    #print "$seq_pos\n";
		    $conversion_array[($pos+$seq_pos)-$low_position][4]=$strand;
		    $conversion_array[($pos+$seq_pos)-$low_position][5]="CG";
		    if ($seq_array[$seq_pos] eq $convert){
			$conversion_array[($pos+$seq_pos)-$low_position][1]+=1;
			$debug=scalar(@conversion_array);
			#print "woohoo\t$debug\n";
		    }elsif($seq_array[$seq_pos] eq $match){
			$conversion_array[($pos+$seq_pos)-$low_position][0]+=1;
			$conversion_array[($pos+$seq_pos)-$low_position][1]+=1;
		    }
		}
	    }
	    @matches=();
	    
	    ###Look for assymetric methylation
	    if($strand eq "+"){
		while ($ref=~/(C[A|T|C]\w)/gi){
		    $context=$1;
		    $f_match=$-[0];
		    pos($ref)=(pos($ref)-2);
		    push @matches, "$f_match\t$context";###@- stores the location of the beginning of each match
		}
	    }
	    if($strand eq "-"){
		while ($ref=~/(\w[A|T|G]G)/gi){
		    $context=$1;
		    $f_match=$-[0];
		    $r_match=$f_match+2;
		    pos($ref)=(pos($ref)-2);
		    $context=reverse($context);
		    $context=~tr/ATCG/TAGC/;
		    push @matches, "$r_match\t$context";###@+ stores the location of the beginning of each match
		}
	    }
	    @seq_array = split //, $seq;
	    if(scalar(@matches)>0){
		foreach $match_data (@matches){
		    ($seq_pos,$context)=split /\t/, $match_data;
		    #print "$seq_pos\n";
		    $conversion_array[($pos+$seq_pos)-$low_position][4]=$strand;
		    $conversion_array[($pos+$seq_pos)-$low_position][5]=$context;
		    if ($seq_array[$seq_pos] eq $convert){
			$conversion_array[($pos+$seq_pos)-$low_position][1]+=1;
			$debug=scalar(@conversion_array);
			#print "woohoo\t$debug\n";
		    }elsif($seq_array[$seq_pos] eq $match){
			$conversion_array[($pos+$seq_pos)-$low_position][0]+=1;
			$conversion_array[($pos+$seq_pos)-$low_position][1]+=1;
		    }
		}
	    }
	    @matches=();
	    
	    
	    
	    
	    
	    ###reverse complement
	    
	    while ($ref=~/CG/gi){
		$f_match=$-[0];
		$r_match=$f_match+1;
		if($strand eq "-"){
		    push @matches, $f_match;###@+ stores the location of the beginning of each match
		}else{
		    push @matches, $r_match;###@- stores the location of the beginning of each match
		}
	    }
	    if(scalar(@matches)>0){
		foreach $seq_pos (@matches){
		    if ($strand eq "+"){
			$conversion_array[($pos+$seq_pos)-$low_position][4]="-";
		    }else{
			$conversion_array[($pos+$seq_pos)-$low_position][4]="+";
		    }
		    $conversion_array[($pos+$seq_pos)-$low_position][5]="CG";
		    if ($seq_array[$seq_pos] eq $rc_convert){
			$conversion_array[($pos+$seq_pos)-$low_position][3]+=1;
		    }elsif($seq_array[$seq_pos] eq $rc_match){
			$conversion_array[($pos+$seq_pos)-$low_position][2]+=1;
			$conversion_array[($pos+$seq_pos)-$low_position][3]+=1;
		    }
		}
	    }
	    @matches=();
	    
	    
	    ###Look for assymetric methylation
	    if($strand eq "-"){
		while ($ref=~/(C[A|T|C]\w)/gi){
		    $context=$1;
		    $f_match=$-[0];
		    pos($ref)=(pos($ref)-2);
		    push @matches, "$f_match\t$context";###@- stores the location of the beginning of each match
		}
	    }
	    if($strand eq "+"){
		while ($ref=~/(\w[A|T|G]G)/gi){
		    $context=$1;
		    $f_match=$-[0];
		    $r_match=$f_match+2;
		    pos($ref)=(pos($ref)-2);
		    $context=reverse($context);
		    $context=~tr/ATCG/TAGC/;
		    push @matches, "$r_match\t$context";###@+ stores the location of the beginning of each match
		}
	    }
	    @seq_array = split //, $seq;
	    if(scalar(@matches)>0){
		foreach $match_data (@matches){
		    ($seq_pos,$context)=split /\t/, $match_data;
		    #print "$seq_pos\n";
		    if ($strand eq "+"){
			$conversion_array[($pos+$seq_pos)-$low_position][4]="-";
		    }else{
			$conversion_array[($pos+$seq_pos)-$low_position][4]="+";
		    }
		    $conversion_array[($pos+$seq_pos)-$low_position][5]=$context;
		    if ($seq_array[$seq_pos] eq $convert){
			$conversion_array[($pos+$seq_pos)-$low_position][3]+=1;
			$debug=scalar(@conversion_array);
			#print "woohoo\t$debug\n";
		    }elsif($seq_array[$seq_pos] eq $match){
			$conversion_array[($pos+$seq_pos)-$low_position][2]+=1;
			$conversion_array[($pos+$seq_pos)-$low_position][3]+=1;
		    }
		}
	    }
	    @matches=();
	    #print "$low_position\n";
	    #<STDIN>;
	}else{ #new chromosome
	    for $seq_pos (0..$#conversion_array){
		if (defined($conversion_array[0])){
		    @methylation_report=@{shift @conversion_array};
		    #print "woo\n";
		    #print "@methylation_report\n";
		    if(defined($methylation_report[4])){
			if(!defined($methylation_report[2])){
			    $methylation_report[2]=0;
			}
			
			if(!defined($methylation_report[3])){
			    $meth_adjust=-1;
			    $methylation_report[3]=0;
			}else{
			    $meth_adjust=$methylation_report[2]/$methylation_report[3];
			};
			
			if(!defined($methylation_report[0])){
			    $methylation_report[0]=0;
			}
			
			if(!defined($methylation_report[1])){
			    $meth_percent=-1.000;
			    $methylation_report[1]=0;
			}else{
			    if($meth_adjust>0){
				$effective_CT_counts=$methylation_report[1]*$meth_adjust;
			    }else{
				$effective_CT_counts=$methylation_report[1];
			    }
			    $meth_percent=$methylation_report[0]/$effective_CT_counts;
			};
			
			#if($meth_adjust>0){
			#    $meth_percent*=$meth_adjust;
			#};
			$meth_percent*=1000;
			$meth_percent=int($meth_percent);
			$meth_percent/=1000;
			$export_pos=$seq_pos+$low_position;
			$cytosine_count++;
			$average_coverage+=$effective_CT_counts;
			print OUTPUT "$chromosome\t$export_pos\t$methylation_report[4]\t$methylation_report[5]\t$meth_percent\t$effective_CT_counts\t$methylation_report[0]\t$methylation_report[1]\t$methylation_report[2]\t$methylation_report[3]\n";
		    };
		}else{
		    shift @conversion_array;
		}
	    };
	    $low_position=0;
	    @conversion_array=();
	    @coverage=();
	    $chromosome=$cr;
	    $repeat="TRUE";
	}
    }
}

close ALIGNMENTS;
close OUTPUT;
print LOG "Finish: ".localtime(time)."\n";
print LOG "Total reads: $read_count\n";
print LOG "Successful alignments: $alignment_count\n";
print LOG "Cytosines covered: $cytosine_count\n";
$average_coverage/=$cytosine_count;
print LOG "Average coverage: $average_coverage\n";
close LOG;



system("gzip $output");