#!/usr/bin/perl -w
#Chris Hale - UCLA
use strict;
use Getopt::Long;

#use lib '~/project-mcdb/PerlModules/Statistics-Basic/lib';
#use Statistics::Basic qw(:all);
###RRBS-parse_correlation_methylation_region.pl
###Take the output from RRBS-correlation_methylation_region.pl and create bed file of regions


my @working_array;
my $USAGE = "\nUSAGE: RRBS-autocorrelation_methylation_state.pl 
                                   --correlation_file /input/directory/file.txt
				   --output_bed /input/directory/file.bed
				   --max_gap [max number of gaps (correlation<cutff) allowed DEFAULT:1]
				   --correlation_cutoff [minimum Rsquared value required to extend region DEFAULT:0.1]
                                   ";

my ($correlation_file,$output_file);
my $max_gap=1;
my $correlation_cutoff=0.1;
my ($line,$chromosome,$position);

my @current_line;
my @trailing_array;
my @leading_array;

my $trailing_count=0;
my $leading_count=0;
my $line_count;

my ($left_end,$right_end);

my $line_position;
my $correlation_value;
my $cytosine_count;

my $current_gap=0;

my $sample_number=0;


GetOptions('correlation_file=s'=>\$correlation_file,
	   'output_bed=s'=>\$output_file,
	   'max_gap=i'=>\$max_gap,
	   'correlation_cutoff=f'=>\$correlation_cutoff
	   );
die $USAGE unless (defined($correlation_file) && defined($output_file));

if ($correlation_file=~/\.gz$/){
    open CORRELATION, "gunzip -c $correlation_file |" or die "Cannot open pipe to $correlation_file\n";
}else{
    open CORRELATION, "$correlation_file" or die "Cannot open $correlation_file\n";
}


open OUTPUT, ">$output_file" or die "Cannot open $output_file";


while ($line=<CORRELATION>){
    chomp($line);
    @current_line=split /\t/, $line;
    
    $sample_number++;
    
    $chromosome=shift @current_line;
    $position=shift @current_line;
    
    $left_end=$position;
    $right_end=$position;
 
    $line_count=$#current_line;
    
    $correlation_value=shift @current_line;
    while ($correlation_value!~/C/i){
        unshift @trailing_array, $correlation_value;
        $correlation_value=shift @current_line;
    };
    if (scalar(@current_line)>0){
        while (scalar(@current_line)>0){
            $correlation_value=shift @current_line;
            push @leading_array, $correlation_value;
        }
    }
    
    $current_gap=0;
    if (scalar(@trailing_array)>0){
        while (($current_gap<=$max_gap) && (scalar(@trailing_array)>0)){
            $correlation_value=shift @trailing_array;
            ($line_position,$correlation_value)=split /:/, $correlation_value;
            if ($correlation_value>$correlation_cutoff){
                $trailing_count+=(1+$current_gap);
                $left_end=$line_position;
                $current_gap=0;
            }else{
                $current_gap++;
            }
        }
    };
    
    $current_gap=0;
    if (scalar(@leading_array)>0){
	while (($current_gap<$max_gap) && (scalar(@leading_array)>0)){
	    $correlation_value=shift @leading_array;
	    ($line_position,$correlation_value)=split /:/, $correlation_value;
	    if ($correlation_value>$correlation_cutoff){
            $leading_count+=(1+$current_gap);
            $right_end=$line_position;
            $current_gap=0;
	    }else{
            $current_gap++;
	    }
	}
    };
    
    $current_gap=0;
    @trailing_array=();
    @leading_array=();
    @current_line=();
    
    $cytosine_count=($trailing_count+$leading_count+1);
    
    print OUTPUT "$chromosome\t$left_end\t$right_end\traw_region_$sample_number\t$cytosine_count\t.\n";
    
    $trailing_count=0;
    $leading_count=0;
}




