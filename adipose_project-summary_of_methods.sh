# methylation region pipeline:

# prior to region calling the following steps were performed to get cytosine methylation ratios

trim_galore -q -1 --phred64 -o ./adipose_project/ --rrbs $fastq_file ; 

bsmap -a $trimmed_fastq -d hg19_unmasked_only_whole_csomes.fa -o $trimmed_fastq"_4mis_uniq_w_reference.sam" -R -p 8 -w 1 -r 0 -v 4 -D C-CGG ; 

perl methratio_alt_v3.pl --input $alignment_bam --output_file $ratio_file ;

# The resultant methylation ratios for each cytosine were combined into a master table for all samples which was then filtered to only include
# cytosines covered in >=150 samples. This table had the following format (with a header row).
# chromosome\tposition\tsample1ratio\tsample2ratio...
# Name of this matrix file is "combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries.txt.gz"


perl RRBS-correlation_methylation_region.pl \
--matrix_file combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries.txt.gz \
--output combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS.txt \
--max_region_size 500 ;

perl RRBS-parse_correlation_methylation_region.pl \
--correlation_file combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS.txt \
--output_bed combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_2_cutoff_2e-1.bed \
--max_gap 0 \
--correlation_cutoff 0.447213595 ;

awk '$5>3' combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_0_cutoff_2e-1.bed > combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_0_cutoff_2e-1_min_cyto_3.bed ; 

sort -k1n,1 -k2n,2 -k3n,3 combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_0_cutoff_2e-1_min_cyto_3.bed > combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_0_cutoff_2e-1_min_cyto_3-sorted.bed ; 

bedtools merge -i combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_0_cutoff_2e-1_min_cyto_3-sorted.bed > combined_methylation_table--AT1-AT52_covered_in_150_seed_libraries--CORRELATION_REGIONS_max_size_1000_max_gap_0_cutoff_2e-1_min_cyto_3-MERGED.bed ; 

