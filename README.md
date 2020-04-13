
Bumpfinder protocol (*assuming_replicate_data n =2)


Requirements: R, UNIX, Samtools, bedtools, bowtie

1. The first step is to align FASTQ file to “gencode” transcripts FASTA file:

bowtie -S gencode.v19.pc_transcripts.fa fastq_name output_sam

2. Sam to bam conversion in SAMTOOLS

samtools view -Sb output_sam > output_bam

3. Sort the bam file

samtools sort output_bam output_bam_sort

4. BAM to BED conversion

bedtools bamtobed -i output_bam_sort > output_bam_sort_bed

5. In the folder ‘alignment_file_preparation: follow the steps in README_2.txt

#### DATA preparation code : ALIGNMENT FILE TO COUNT FILES #####

#####REQUIRES: UNIX, BEDTOOLS

BED files are obrained from conversion of alignment (BAM) file to BED file and awk proccessing to keep only coordinates, length, strand (since transcript alignment, this is always '+' and ENST IDS)

Count files are obtained (overlap with longest_transcript_select.bed) through use of BEDTOOLS:

EXAMPLE:

bedtools coverage -a ../human_data_windows/longest_transcript_select_100_windows.bed -b ribo_final_1.bed | cut -f1,2,3,4 > ribo_final_1.bed_counts.txt


Note:

longest_transcript_select_100_windows.bed- is obtained from analysis in folder human_data_windows: Read /human_data_windows/README_1.txt...*This is one time step

Ribo_final_1.bed is example of a bed file obtained from step 4.

Output is “output_bed_counts”


6. Next step is to create read and difference profiles per transcript per window (in that order).

Follow steps in window_overlap/README_3.txt


#### Processing code : MAKE GENE PROFILES #####

#####REQUIRES: UNIX, BEDTOOLS, R


For a comparison condition (Example, for MD55, minus/plus IFN - 2 replicates (m1,m2 - minusIFN, p1,p2- plusIFN)), a data_md55.txt is made from adding reads, from all conditions, that overlap with 'longest_transcript_select.bed '




ids file (Eg, ids_mdd5.txt ) was created by cut -f1 data_md55.txt | sort -u > ids_md55.txt

'RUN R script profiler_md55.R'


Similarily;

data_tyr.txt for tyrosine depletion
data_md55.txt for infy induction (48hrs) in md55 cellsmore 
data_48h for ifny induction (48hrs) in 12T cells
data_24h for ifny induction (24hrs) in 12T cells
data_108T for ifny induction (48hrs) in 108T cells

and profiler_tyr.R, profiler_md55.R, profilier_48hr.R, profiler_24hr.R and profiler_108T.R respectively,

OUTPUT files are GENE_DETAILS* and PROF_*


The immportant output is GENE_DETAILS*


7. Calling Bumps from GENE_DETAILS: processing folder

Make a seperate folder for every data: For example, put GENE_DETAILS_md55.txt in the folder: md55. 


#### Processing code : IDENTIFY bumps (or here peaks) #####

#####REQUIRES: UNIX, BEDTOOLS, R

In the folders corresponding to particular datasets (eg, md55)

On 'GENE_DETAILS*", run 'script_1.R'

Later run bedtools.sh

Run ‘script_2.R”

copy PLUS_ONLY_AWAY_DETAILED.xls to plus folder..…
copy MINUS_ONLY_AWAY_DETAILED.xls to minus folder


In both these folders (make sure you have gencode.v19.pc_translations.fa and CDS.bed) :

(a) On unix:

cut -f1 PLUS_ONLY_AWAY_DETAILED.xls | sed 's/_/\t/g' | grep "ENST" > PLUS_ONLY_AWAY_DETAILED.bed 

cut -f1 MINUS_ONLY_AWAY_DETAILED.xls | sed 's/_/\t/g' | grep "ENST" > PLUS_ONLY_AWAY_DETAILED.bed 

(b) 

run; script_3.R

(c) On UNIX CONSOLE:

awk '{print $1 "\t" $7 -29 "\t" $7 +29}' bumps_cds_atleast30.txt | grep "ENST" > bumps_30each.bed

(d) Run PERL script

perl sequence_extractor.pl | sort -u > bumps_30_each.txt

(e) On unix console;

cut -f5 bumps_30_each.txt  | sed 's/A/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > A.txt
cut -f5 bumps_30_each.txt  | sed 's/C/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > C.txt
cut -f5 bumps_30_each.txt  | sed 's/D/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > D.txt
cut -f5 bumps_30_each.txt  | sed 's/E/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > E.txt
cut -f5 bumps_30_each.txt  | sed 's/F/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > F.txt
cut -f5 bumps_30_each.txt  | sed 's/M/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > M.txt
cut -f5 bumps_30_each.txt  | sed 's/N/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > N.txt
cut -f5 bumps_30_each.txt  | sed 's/G/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > G.txt
cut -f5 bumps_30_each.txt  | sed 's/H/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > H.txt
cut -f5 bumps_30_each.txt  | sed 's/I/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > I.txt
cut -f5 bumps_30_each.txt  | sed 's/K/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > K.txt
cut -f5 bumps_30_each.txt  | sed 's/L/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > L.txt
cut -f5 bumps_30_each.txt  | sed 's/P/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > P.txt
cut -f5 bumps_30_each.txt  | sed 's/Q/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > Q.txt
cut -f5 bumps_30_each.txt  | sed 's/R/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > R.txt
cut -f5 bumps_30_each.txt  | sed 's/S/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > S.txt
cut -f5 bumps_30_each.txt  | sed 's/T/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > T.txt
cut -f5 bumps_30_each.txt  | sed 's/V/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > V.txt
cut -f5 bumps_30_each.txt  | sed 's/Y/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > Y.txt
cut -f5 bumps_30_each.txt  | sed 's/W/1\t/g' |  sed 's/[A-Z]/0\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,49,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58  > W.txt

(f) Run plot.R





# End of Flowchart






