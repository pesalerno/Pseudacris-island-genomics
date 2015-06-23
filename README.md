# Pseudacris-island-genomics
Following workflow is for processing raw data from RADseq libraries. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage and single-end. The higher coverage reads will be used for generating a reference genome. Thus, this pipeline actually runs Stacks three separate times:
- First: denovo_map for paired-end reads (Step 1.4)
- Second: ref_map for single-end reads (Step 3)
- Third: rxstacks for fixing genotyping errors with population statistics (Step 4). 

## Step 1: de-multiplexing

######1.1. De-multiplex each library individually
De-multiplexing was done with program [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) individually for each library within its directory, and renamed with sample names within Stacks ([here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/barcodes-1933.txt) is an example barcodes file). 

- Commands for process_radtags for the two Paired-end libraries were:

		process_radtags -P -p ./PE-lib-1610/ -b barcodes-1610.txt -i gzfastq - \
		o ./processed-1610/ -e sbfI -c -q -r -D
	(example for library #1610 for *Xanthusia*)
- Commands for process_radtags for all the other single-end libraries were:

		process_radtags -p ./SR-lib5-1994/ -b barcodes-1994.txt -i gzfastq \ 
		-o ./processed-1994/ -e sbfI -c -q -r -D
	(example for library #1994 for shared *Xanthusia* and *Pseudacris*)


######1.3. Merge all libraries into single directory

Copy all renamed libraries for all individuals into their "species" directories (up to this point we had one library where *Xanthusia* and *Pseudacris* were mixed together - library #1994)


######2. purge PCR duplicates from within each file

I can use [clone_filter](http://catchenlab.life.illinois.edu/stacks/comp/clone_filter.php) from Stacks to purge PCR duplicates. The difference between the stacks script and the purge_PCR_duplicates script (below) is that the second one retains quality data and the first doesn't. 

----------------------------------

I need to run the open source perl script [purge_PCR_duplicates.pl](https://github.com/claudiuskerth/scripts_for_RAD/blob/master/purge_PCR_duplicates.pl) by Claudius Kerth. It needs the perl module [Parallel::ForkManager](http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm) since it is set up for running parallelized, which Dan Sloan installed for me on the server's root.

For the program to run, files need to be unzipped (.fq) and end with either
\"fq_1\" for the SE file or \"fq_2\" for the PE file. Example: XXX.fq_1 and XXX.fq_2. Use above script to rename all files to add -1 termination.

For usage, type:

	./purge_PCR_duplicates.pl [options] > logfile

Having all files to be purged within the same directory as the script.


----------------------------------------------
----------------------------------------------


###Intermediate step: get a quick SNP matrix only with SR sequences for Chris' report

######1. check to see which files are repeats in different libraries, or else when moving them to do denovo_map they will be rewritten/lost.

Libraries #1612 and #1835 are essentially duplicates, with few exceptions. Library #1834 is mostly unique with a few that are repeats from #1612. So, I'm renaming all sequence files and adding a -1.fq.gz to library #1 (1612), -2.fq.gz to library 2 (1834), -3.fq.gz to library 3 (1835), -4.fq.gz to library 4 (1994), and -5.fq.gz to library 5 (1995). Renaming all files at once using the following code:

	rename .fq.gz -1.fq.gz *.fq.gz
	
	


------------------------------------


######2. Merge fasta files for library duplicates

After being renamed, move all files back to the SR-denovo-prelim folder and there I merge the fasta files. I merge following these guidelines (from this [source](http://www.researchgate.net/post/How_do_I_merge_several_multisequence-fasta_files_to_create_one_tree_for_subsequent_Unifrac_analysis)):

*To merge several files use the SHELL, go to your folder where the files are and use the cat command. E.g. to merge seqfile001.fasta, seqfile002.fasta and seqfile003.fasta type*

	cat seqfile001.fasta seqfile002.fasta seqfile003.fasta > seqcombined.fasta


*or if you have more files use*

	cat *.fasta > seqcombined.fasta


The duplicated files are sorted into a separate folder before merging, just to keep track of what's being merged. Then the post-merged files are sorted back into the general directory containing all sequences. Total number of files before merging duplicates from different ***Xanthusia*** library preps was 187, and after merging duplicate individuals we now have 142 files for denovo_map input. Total number of files before merging duplicates from different ***Pseudacris*** library preps was 180, and after merging duplicate individuals we now have 132 files for denovo_map input. 

------------------------------------


######3. Running denovo_map for SR reads

The code used for running denovo_map for only the SR reads was:

	denovo_map.pl -m 3 -M 2 -n 1 -T 16 -b 1 -t -S -o ./denovo-1/ \



------------------------------------


######4. running program populations for exporting SNP matrix

I'm running populations with the filters for keeping SNPs that are present in all populations for ***Xanthusia***, and in xx/xx populations for ***Pseudacris*** (to avoid losing too many SNPs with mainland species). Here is the popmap for [*Xanthusia*]() and the popmap for [*Pseudacris*](). 

The script for populations for ***Xanthusia*** was:

	populations 

The script for populations for ***Pseudacris*** was:

	populations 


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

 