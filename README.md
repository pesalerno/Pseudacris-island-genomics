# Pseudacris/Xantusia-island-genomics
Following workflow is for processing raw data from RADseq libraries. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage and single-end. The higher coverage reads will be used for generating a reference genome. Thus, this pipeline actually runs Stacks three separate times:
- First: denovo_map for paired-end reads (Step 1.4)
- Second: ref_map for single-end reads (Step 3)
- Third: rxstacks for fixing genotyping errors with population statistics (Step 4). 

##Step 1: de-multiplexing

######1.1. De-multiplex each library individually
De-multiplexing was done with program [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) individually for each library within its directory, and renamed with sample names within Stacks ([here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/barcodes-1933.txt) is an example barcodes file). 

- Commands for process_radtags for the two Paired-end libraries were:

		process_radtags -P -p ./PE-lib-1610/ -b barcodes-1610.txt -i gzfastq - \
		o ./processed-1610/ -e sbfI -c -q -r -D
	(example for library #1610 for *Xantusia*)
- Commands for process_radtags for all the other single-end libraries were:

		process_radtags -p ./SR-lib5-1994/ -b barcodes-1994.txt -i gzfastq \ 
		-o ./processed-1994/ -e sbfI -c -q -r -D
	(example for library #1994 for shared *Xantusia* and *Pseudacris*)


######1.3. Merge all libraries into single directory

Copy all renamed libraries for all individuals into their "species" directories (up to this point we had one library where ***Xantusia*** and ***Pseudacris*** were mixed together - library #1994)


##Step 2: purge PCR duplicates of PE reads


I ran the open source perl script [purge_PCR_duplicates.pl](https://github.com/claudiuskerth/scripts_for_RAD/blob/master/purge_PCR_duplicates.pl) by Claudius Kerth. It needs the perl module [Parallel::ForkManager](http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm) since it is set up for running parallelized, which Dan Sloan installed for me on the server's root.

For the program to run, files need to be unzipped (.fq) and end with either
\"fq_1\" for the SE file or \"fq_2\" for the PE file. Example: XXX.fq_1 and XXX.fq_2. Use rename script to rename all files to add -1 termination.

For usage, I typed:

	perl purge_PCR_duplicates.pl > logfile

while having all files to be purged within the same directory as the script. The script does not overwrite original files, but outputs new files that have been purged (and adds "purged" to file names).

The reads per individual pre and post purging of PCR duplicates, and the percent purged, were:

Individual | initial reads | porst-purge reads | percent purged |
------------ | ------------- | ------------ | ------------- | 
Pr_SCI_03	|	6,703,487	|	1,871,919	|	72	|
Pr_SCI_04	|	2,159,994	|	589,841	|	73	|
Pr_SCI_05	|	2,960,879	|	818,600	|	72	|
Pr_SCI_06	|	4,413,774	|	1,193,235	|	73	|
Pr_SRI_02	|	7,572,368	|	2,094,594	|	72	|
Pr_SRI_03	|	1,298,586	|	353,847	|	73	|
Pr_SRI_04	|	2,636,033	|	720,639	|	73	|
Pr_SRI_05	|	2,145,094	|	591,698	|	72	|
Xr_SBI_03	|	1,727,677	|	987,327	|	43	|
Xr_SBI_04	|	2,582,799	|	1,475,626	|	43	|
Xr_SCL_27	|	1,368,256	|	769,389	|	44	|
Xr_SCL_28	|	1,646,492	|	922,189	|	44	|
Xr_SNI_27	|	1,109,304	|	632,239	|	43	|
Xr_SNI_28	|	3,094,797	|	1,766,134	|	43	|
Xv_JTS_03	|	1,650,289	|	938,941	|	43	|
Xv_JTS_04	|	3,754,638	|	2,124,663	|	43	|



##Step 3: Generate reference contigs from PE reads


I initally genotyped and created reference contigs with STACKS, but now I am re-doing the analysis with pyrad to account for the high divergence across islands, particularly for *Xantusia riversiana/vigilis*. 


#####3.1. Merge paired-reads
Because in RADseq reads for the same loci can be of different lengths, and most will be overlapping segments (as in, R1 and R2 will overlap) then we will merge the reads following flow-cell information so that they are processed together when making the stacks (greatly reduces computational time, and also prevents a messy analysis). We have to merge the reads with the program [PEAR](https://github.com/xflouris/PEAR).

=>The first thing we have to do is unzip the reads if they are gzipped


From the Manual:

-----------
***How to use:** PEAR  can  robustly  assemble most of the data sets with default parameters. The basic command to run PEAR is:*

		./pear -f forward_read.fastq -r reverse_read.fastq -o output_prefix
	
*The forward_read file usually has "R1" in the name, and the reverse_read file usually has "R2" in the name.*

***How to cite:** Zhang, J., K. Kobert, T. Flouri, A. Stamatakis. PEAR: A fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30(5): 614-620, 2014.*

I renamed files with:

	rename -1.fq_1 _R1.fq *fq*
	
Then I used the following for loop script from Deren Eaton (within folder with sequences):

	for gfile in *.fq;
    do pear -f $gfile \
            -r ${gfile/_R1.fq/_R2.fq} \
            -o ${gfile/_R1.fq/} \
            -n 33 \
            -t 33 \
            -q 10 \
            -j 20  >> pear.log 2>&1;
	done

Then I transferred only the *assembled* on ran step#3, changing line # 11 (data type) to ***merged***:

	
And also changed the path to the files so that it only picks files that have *assembled* in the name. 




\

\

\

\

\



====> **RE-WRITE** THIS SECTION FOLLOWING RECENT CHANGES IN PIPELINE







Setting up experiments (permutations) in Stacks for testing which parameters combinations are better for retrieving a higher number of "good" loci that have decent coverage within and among populations. Because we have decently high coverage, I'm varying the -m parameter (# reads required to form a stack) from 3-7, skipping even numbers (just to have a bigger range initially). Also, I'm not at all using -n 0 or 1, since they are biologically unrealistic given these datasets, and likely to eliminate and oversplit loci among individuals/populations with higher divergence. 

Permutations | -m | -M | -n | --max_locus_stacks 
------------ | ------------- | ------------ | ------------- | ------------ |
a | 3 | 2 | 2 | 3 | 
b | 5 | 2 | 2 | 3 |
c | 7 | 2 | 2 | 3 | 
d | 3 | 3 | 2 | 3 |
e | 3 | 4 | 2 | 3 |
f | 3 | 5 | 2 | 3 |
g | 3 | 2 | 3 | 3 |
h | 3 | 2 | 4 | 3 |
i | 3 | 2 | 5 | 3 |
j | 3 | 2 | 2 | 4 |
k | 3 | 2 | 2 | 5 |

----------------------------------------------
----------------------------------------------


##Step 4: Genotyping of SR reads based on PE contigs

######4.1. Check for among-library repeats.

Libraries #1612 and #1835 are essentially duplicates, with few exceptions. Library #1834 is mostly unique with a few that are repeats from #1612. So for ***Xantusia***, I'm renaming all sequence files and adding a -1.fq.gz to library #1 (1612), -2.fq.gz to library 2 (1834), -3.fq.gz to library 3 (1835), -4.fq.gz to library 4 (1994), and -5.fq.gz to library 5 (1995). Renaming all files at once using the following code:

	rename .fq.gz -1.fq.gz *.fq.gz
	
	


------------------------------------


######4.2. Merge fasta files for library duplicates

After being renamed, move all files back to the SR-denovo-prelim folder and there I merge the fasta files. I merge following these guidelines (from this [source](http://www.researchgate.net/post/How_do_I_merge_several_multisequence-fasta_files_to_create_one_tree_for_subsequent_Unifrac_analysis)):

*To merge several files use the SHELL, go to your folder where the files are and use the cat command. E.g. to merge seqfile001.fasta, seqfile002.fasta and seqfile003.fasta type*

	cat seqfile001.fasta seqfile002.fasta seqfile003.fasta > seqcombined.fasta


*or if you have more files use*

	cat *.fasta > seqcombined.fasta


The duplicated files are sorted into a separate folder before merging, just to keep track of what's being merged. Then the post-merged files are sorted back into the general directory containing all sequences. Total number of files before merging duplicates from different ***Xantusia*** library preps was 187, and after merging duplicate individuals we now have 142 files for denovo_map input. Total number of files before merging duplicates from different ***Pseudacris*** library preps was 180, and after merging duplicate individuals we now have 132 files for denovo_map input. 

------------------------------------
#####How many reads on average for each species?? 
![read counts](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/read-counts-SR.png)

/

Histogram of reads per individual for ***Pseudacris***:
![PSeu-reads](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/Pseudacris-reads.png)

/

Histogram of reads per individual for ***Xantusia***:
![Xari-reads](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/Xantusia-reads.png)

/

######4.3. Genotyping with pyrad
***===>folders ready to go for genotyping of all single-end illumina reads.***


/

/

/


In this step, all forward reads are used, including the ones from the initial Paired-end libaries. 

The code used for running denovo_map for only the SR reads was:

	denovo_map.pl -m 3 -M 2 -n 1 -T 16 -b 1 -t -S -o ./denovo-1/ \



------------------------------------


######4. running program populations for exporting SNP matrix

I'm running populations with the filters for keeping 50% of SNPs that are present in each populations and with two settings for numbers of populations kept. Here is the popmap for [*Xantusia*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/popmap_Xari.txt) and the popmap for [*Pseudacris*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/popmap-Pseu.txt). 

The script for populations for ***Xantusia*** and for **Pseudacris** was:

	populations -b 1 -P ./denovo-1 -M ./popmap_Xari.txt  -t 36 -p 6 -r 0.5 --write_random_snp --structure --genepop --vcf

I ran this twice, first with 6/7 populations (-p 6), and second with all populations (-p 7).I tried to run the filter of minor allele frequency (**--min_maf**) and it kept failing, maybe I need to upgrade version of Stacks, or maybe I don't know how to set it! This is why I saved it in **--vcf** format so that it can be viewed and filtered in [GATK](https://www.broadinstitute.org/gatk/) easily (or other variant calling format software).

The total number of SNPs were:


- ***Xantusia*** p=6 => 7739 snps
	
- ***Xantusia*** p=7 => 1319 snps 

- ***Pseudacris*** p=6 => 36921 snps

- ***Pseudacris*** p=7 => 36921 snps (did I make a mistake in the previous one? this one's correct according to script). it's just odd they're identical.... 


Just out of curiosity, running ***Pseudacris*** again with a much more stringent setting of r=0.8 to see how many SNPs are still kept. 

Based on number of SNPs, it's already suggesting that ***Xantusia*** has a much higher divergence among populations.... 

-> When doing -p7 -r 0.7 for ***Xantusia*** you get ZERO SNPs.... which means that there is either a huge amount of missing data, or that the mainland population is too divergent for this constraint.... 




----------------------------------------------
----------------------------------------------
----------------------------------------------
----------------------------------------------



#

#

##Step 2: create reference genome with PE reads

######2.1. Run denovo_map program for paired-end reads

Run in Stacks the script [denovo_map](http://creskolab.uoregon.edu/stacks/comp/denovo_map.php) with libraries with higher depth of coverage and paired-end reads for creating the reference genome. I'm trying to be highly conservative about the flag parameters. 

	> denovo_map.pl -m 3 -M 1 -n 0 -T 16 -b 1 -t -S -o ./denovo-1/ \


From Hohenlohe et al. 2013: *"We grouped the forward and reverse reads from all individuals in these populations into a separate file for each RAD locus, using the STACKS program sort_read_pairs.pl"*

More from Hohenlohe et al. 2013: *"We created a catalog of RAD tag loci using cstacks and matched indi- viduals against the catalog using sstacks. We populated and indexed a MYSQL database of loci using load_ radtags.pl and index_radtags.pl and then exported the data using export_sql.pl. Finally, we grouped the for- ward and reverse reads from each individual corre- sponding to each RAD locus using sort_read_pairs.pl."*

More from Hohenlohe et al. 2013: *"We assembled the reads in each file separately to produce a set of RAD contigs (Fig. 2b), using both VELVET (Zerbino & Bir- ney 2008) and CAP3 (Huang & Madan 1999) assembly software."*


#

#

#

#

#

#

#

#

#




---> May want to run [exec_velvet.pl](http://catchenlab.life.illinois.edu/stacks/comp/exec_velvet.php) to generate collated fasta file for reference genome.

Use the output consensus sequence file from denovo_map (*"catalogs.tags.tsv"*) and transform to fasta format for input into [bwa](http://bio-bwa.sourceforge.net/bwa.shtml), using Kelly's R script (need to modify once I run it):

	
	# Import the data and check the structure
	tags<-read.table('batch_1.catalog.tags.tsv', header=FALSE)
	# all the sequences are in $V9

	unique(tags$V2) # sample ID is meaningless here
	length(unique(tags$V3)) #417153 -- each row has a unique locus ID

	# verify that all of the tags are consensus:
	consensus.tags<-subset(tags, V6=='consensus')
	length(consensus.tags[,1]) #417153
	length(tags[,1]) #417153

	# each sequence needs a fasta header to uniquely identify it
	fa.id<-paste('>', tags$V3, '_pseudoreference_pe_concatenated_without_rev_complement', sep='')

	fa<-cbind(fa.id, as.character(tags$V9))
	write.table(fa, file='D_variabilis_denovo_psuedoreference.fa', quote=FALSE, sep='\n', row.names=FALSE,
            col.names=FALSE)
 

######2.2. Index reference genome and align in bwa

Create reference genome from fasta consensus sequences:

	> bwa index paired-end.fasta

Then align all of the de-multiplexed files (single-end 100bp reads) to reference genome:

	> bwa mem -t 6 paired-end.fasta all-reads-demultiplexed.fq > align-allRADs.sam
	
Transform .sam file to .bam file for visualization in IGV:

	> module load samtools
	> samtools view -b -S -o 454-align.bam 454-align.sam
	> samtools sort 454-align.bam 454-align.sorted
	> samtools index 454-align.sorted.bam

Visualize in IGV

Use either .sam or .bam alignment for input in ref_map.pl

##Step 3: Map to new reference genome in Stacks

Use ref_map.pl for mapping the raw read files to the reference genome we created using the in-depth paired-end reads plus the single-end reads of all individuals.

	> ref_map.pl -o /path/for/output -n 2 -m 2 -T 12 -O popmap.txt -b 1 -S -s ./all-sequences-here\

##Step 4: Filter data with haplotype corrections in Stacks

Here we essentially re-run the stacks pipeline (post-demultiplexing) to make corrections based on population haplotype statistics (reduced probability of obtainind fake haplotypes created by paralog stacks and PCR duplicate errors)

	> rxstacks -b 1 -P ./input-stacksoutput/ -o /new-output-path/ --prune_haplo --lnl_filter --lnl_dist 

 