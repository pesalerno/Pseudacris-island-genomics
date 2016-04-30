# ***Pseudacris/Xantusia***-island-genomics
Following workflow is for processing raw data from several RADseq libraries from species *Pseudacris regilla* and *Xantusia riversiana*. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage and single-end. The higher coverage reads will be used for generating reference contigs. Then all of the forward reads (including initial higher coverage individuals) are mapped onto reference contigs for genotyping. The main pipelineused for genotyping is pyrad. 



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


I ran the open source perl script [purge_PCR_duplicates.pl](https://github.com/claudiuskerth/scripts_for_RAD/blob/master/purge_PCR_duplicates.pl) by Claudius Kerth. It needs the perl module [Parallel::ForkManager](http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm) since it is set up for running parallelized.

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


######3.1. Merge paired-reads
Because in RADseq reads for the same loci can be of different lengths, and most will be overlapping segments (as in, R1 and R2 will overlap) then we will merge the reads following flow-cell information so that they are processed together when making the stacks (greatly reduces computational time, and also prevents a messy analysis). We have to merge the reads with the program [PEAR](https://github.com/xflouris/PEAR).

=>The first thing we have to do is unzip the reads if they are gzipped


From the Manual:

-----------
-----------

***How to use:** PEAR  can  robustly  assemble most of the data sets with default parameters. The basic command to run PEAR is:*

		./pear -f forward_read.fastq -r reverse_read.fastq -o output_prefix
	
*The forward_read file usually has "R1" in the name, and the reverse_read file usually has "R2" in the name.*

***How to cite:** Zhang, J., K. Kobert, T. Flouri, A. Stamatakis. PEAR: A fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics 30(5): 614-620, 2014.*

-----------
-----------

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

Then I transferred only the *assembled* to the ***'/edits/'*** folder.

######3.2. Run pyrad: *within-sample clustering* (step 3)

I started the pyrad pipeline on step#3, making sure that line # 11 (data type) is set to ***merged***:

	merged       ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)



See [*Pseudacris* parameters](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-params-PE.txt) and [*Xantusia* parameters](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-params-PE.txt) file for the rest of the specifications of parameters files for the paired end reference contigs. 

pyrad was called as follows:

	pyrad -p Pr-params-d.txt -s 3

correns .sh file was set up as follows:

	#!/bin/bash
	#
	#SBATCH --time=100:00:00
	#SBATCH --job-name=pyrad-Xr-step3-t1
	#SBATCH --mail-user=patriciasalerno@gmail.com
	#SBATCH --error=stderr-pyrad-Xr-step3-t1
	#SBATCH --output=stdout-pyrad-Xr-step3-t1

	/opt/software/Python-2.7.10/python
	pyrad -p /home/salerno/Xantusia/pyrad/Xr-params-t1.txt -s 3

The full output for the within-sample clustering can be found here for [*Xantusia*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-within-clusters-s3.txt) and for [*Pseudacris*](). The summary of coverage per loci within sample is:


taxa	|	total	|	mean-depth	|	std-dev	|	total_d>9	|	mean_d>9	|	std-dev_d>9
----------------- | ------------- | ------------ |------------ | ------------- |------------ | ------------- |
Xr_SBI_03	|	417666	|	1.499	|	1.38	|	1394	|	14.249	|	8.826
Xr_SBI_04	|	581256	|	1.555	|	1.734	|	4692	|	13.598	|	8.596
Xr_SCL_27	|	112121	|	3.059	|	3.52	|	4680	|	13.343	|	9.422
Xr_SCL_28	|	108005	|	3.682	|	4.244	|	9479	|	13.951	|	6.378
Xr_SNI_27	|	274638	|	1.463	|	1.168	|	439	|	15.59	|	9.477
Xr_SNI_28	|	599975	|	1.689	|	2.169	|	9527	|	13.911	|	7.323
Xv_JTS_03	|	303511	|	1.791	|	2.087	|	3701	|	13.357	|	8.52
Xv_JTS_04	|	174084	|	4.967	|	6.667	|	29066	|	16.529	|	8.921


######3.2. Run pyrad: *error rate and heterozygosity estimation* (step 4):

	pyrad -p Pr-params-d.txt -s 4

The result of this step for both species is:

Individual | heterozygosity | error	|
------------ | ------------- | ------------ |
Xr_SCL_27-PE	|	0.00798693	|	0.00233735	|
Xr_SCL_28-PE	|	0.00524442	|	0.00178935	|
Xr_SNI_27-PE	|	0.03099332	|	0.01102511	|
Xv_JTS_04-PE	|	0.00605539	|	0.00147974	|
Xv_JTS_03-PE	|	0.01107941	|	0.00316072	|
Xr_SBI_03-PE	|	0.01819725	|	0.00687464	|
Xr_SBI_04-PE	|	0.00845829	|	0.00365785	|
Xr_SNI_28-PE	|	0.00587565	|	0.00232468	|


######3.4. Run pyrad: *within-sample consensus sequences* (step5):
This analysis uses the error rate and heterozygosity estimations/corrections from step 4. For example, for ***Xantusia***:

	pyrad -p Xr-params-t1.txt -s 5
	
	     ------------------------------------------------------------
      pyRAD : RADseq for phylogenetics & introgression analyses
     ------------------------------------------------------------


	step 5: creating consensus seqs for 8 samples, using H=0.01174 E=0.00408

And for ***Pseudacris***:

	pyrad -p Pr-params... -s 5
	
**From the [manual](http://nbviewer.jupyter.org/gist/dereneaton/af9548ea0e94bff99aa0/pyRAD_v.3.0.ipynb):** *"create consensus sequences. Using the mean error rate and heterozygosity estimated in step 4, this step creates consensus sequences for each cluster. Those which have less than the minimum coverage, more than the maximum number of undetermined sites, or more than the maximum number of heterozygous sites, or more than the allowed number of alleles, are discarded. In diploid data if two alleles are present the phase of heterozygous sites are retained in the consensus sequences"*

######3.5. Run pyrad: *among-sample clustering* (step6)

**From the manual**: *"Consensus sequences are clustered across samples using the same settings as in step 3. If heterozygous, one allele is randomly sampled and used in this step, although both alleles are retained in the final data set."*



\

\

\



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
######=>How many reads on average for each species?? 


We counted reads for each individual using the unzipped files and with the following script:

	echo -e 'SAMPLE_ID_FULL\tNUM_READS'
	for file in ~/path/to/denovo-map/*.fq ##edit your path!!
	do
	echo -n $(basename $file .fq)$'\t'
	cat $file | grep '^@.*' | wc -l
	done

![read counts](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/read-counts-SR.png)

/

Histogram of reads per individual for ***Pseudacris***:
![PSeu-reads](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/Pseudacris-reads.png)

/

Histogram of reads per individual for ***Xantusia***:
![Xari-reads](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/Xantusia-reads.png)

/

######4.3. Genotyping with pyrad
***===>folders ready to go for genotyping of all single-end illumina reads.*** In this step, all forward reads are used, including the ones from the initial Paired-end libaries. So they all need to be merged into one directory.


/

/

/


 

