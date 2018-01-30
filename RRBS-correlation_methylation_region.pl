#!/usr/bin/perl -w
use strict;
use Getopt::Long;
##Chris Hale
##UCLA

##Find the correlation of "local" cytosine methylation to the methylation state of every cytosine. Essentially, create a sliding window across the methylome that takes each cytosine and finds the correlation of nearby methylation to the methylation of that cytosine and generates output that can be parsed into regions later. Input is a large matrix of cytosine methylation values by cytosine position for an array of individuals.

use lib '~/project-mcdb/PerlModules/Statistics-Basic/lib';
use Statistics::Basic qw(:all);



my @working_array;
my $USAGE = "\nUSAGE: RRBS-autocorrelation_methylation_state.pl 
                                   --matrix_file /input/directory/file.txt.gz
				   --output /input/directory/file.txt
				   --max_region_size [maximum size in bp of a region DEFAULT:10000]
                                   ";

my ($matrix_file,$output_file);
my $max_size=10000;
my ($line,$chromosome,$position);
my ($new_chromosome,$new_position);
my $holding_line="empty";
my $current_distance;
my $max_size_reached="FALSE";

my @correlation_array;
my @vector_array;
my @current_line;
my @compare_array;

my $compare_line;

my @vector_1;
my @vector_2;
my @export_line;

my $correlation_value;
my $sample_number;
my $cytosine_number;
my $distance;
my $cytosine_count;

#my $autocorrelation_ref=\@autocorrelation_array;

GetOptions('matrix_file=s'=>\$matrix_file,
	   'output=s'=>\$output_file,
	   'max_region_size=s'=>\$max_size
	   );
die $USAGE unless (defined($matrix_file) && defined($output_file) && defined($max_size));

open MATRIX, "gunzip -c $matrix_file |" or die "Cannot open pipe to $matrix_file\n";
open OUTPUT, ">$output_file" or die "Cannot open $output_file";

#Header
$line=<MATRIX>;

#First cytosine
$line=<MATRIX>;
chomp($line);
@current_line=split /\t/, $line;
$chromosome=$current_line[0];
$position=$current_line[1];

push @vector_array, [@current_line];



while (!eof(MATRIX)){
    if($holding_line eq "empty"){
        $line=<MATRIX>;
    }else{
        $line=$holding_line;
        $holding_line="empty";
    }
    chomp($line);

    @current_line=split /\t/, $line;
    if ($chromosome eq "empty"){
        $chromosome=$current_line[0];
        $position=$current_line[1];
        push @vector_array, [@current_line];
    };
    
    $new_chromosome=$current_line[0];
    $new_position=$current_line[1];

    if (($chromosome eq $new_chromosome) && (($new_position-$position)<$max_size) && (($new_position-$position)!=0)){
        push @vector_array, [@current_line];
        $distance=$new_position-$position;
    }elsif(($new_position-$position)!=0){
        if ($chromosome ne $new_chromosome){
            print "$new_chromosome\n";
        }
        $holding_line=$line;
        push @{$correlation_array[0]}, "C" ;
        $cytosine_count=1;
        @current_line=@{shift @vector_array};
        foreach $compare_line (@vector_array){
            @compare_array=@$compare_line;
            for $sample_number (2..$#current_line){
                if (($current_line[$sample_number] ne "NA") && ($compare_array[$sample_number] ne "NA")){
                    push @vector_1, $current_line[$sample_number];
                    push @vector_2, $compare_array[$sample_number];
                }
            }
            $distance=$compare_array[1]-$current_line[1];
            $correlation_value=correlation([@vector_1],[@vector_2]);
            #print "$correlation_value\n";
            if (!defined($correlation_value)){
                $correlation_value=0;
            }#else{
            #    $correlation_value*=$correlation_value;
            #}
            push @{$correlation_array[0]}, "$compare_array[1]:$correlation_value";
            push @{$correlation_array[$cytosine_count]}, "$current_line[1]:$correlation_value";
            $cytosine_count++;
            @vector_1=();
            @vector_2=();
        }
        @export_line=@{shift @correlation_array};
        print OUTPUT join("\t",$chromosome,$position,@export_line);
        print OUTPUT "\n";
        if (scalar(@vector_array)>0){
            #print "Hoo!\n";
            $compare_line=$vector_array[0];
            $chromosome=$compare_line->[0];
            $position=$compare_line->[1];            
        }else{
            #print "Bingo!\n";
            #<STDIN>;
            $chromosome="empty";
            $position="empty";
        }
    }    
}




