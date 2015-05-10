# Pseudacris-island-genomics
Following workflow is for processing raw data from RADseq libraries. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage and single-end. The higher coverage reads will be used for generating a reference genome. Thus, this pipeline actually runs Stacks three separate times:
- First: denovo_map for paired-end reads (Step 1.4)
- Second: ref_map for single-end reads (Step 3)
- Third: rxstacks for fixing genotyping errors with population statistics (Step 4). 

## Step 1: de-multiplexing

######1.1. De-multiplex each library individually
If too many libraries, find a way to automate with a for loop to run within each directory at the same time.

######1.2. Re-name files for adding individual IDs, populations, etc.  
First change file terminations to add the extra "-" required in python script. Then run python script [rename_barcodes_fixed.py](https://github.com/pesalerno/Chimanta-genomics/blob/master/rename_barcodes_fixed.py) within each demultiplexing directory. 

######1.3. Merge all libraries into single directory

Copy all renamed libraries for all individuals into their "species" directories (at this point we have *Xanthusia* and *Pseudacris* mixed in the same libraries)
Make a ls > list-of-files.txt file for generating input for [denovo_map](http://creskolab.uoregon.edu/stacks/comp/denovo_map.php) 

######1.4. Run denovo_map for more in-depth sequences.

Run in Stacks the script denovo_map.pl with libraries with higher depth of coverage and paired-end reads for creating the reference genome:

	> denovo_map.pl enter-rest-of-script


##Step 2: create reference genome in bwa

######2.1. Transform consensus sequence to fasta format.
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

 