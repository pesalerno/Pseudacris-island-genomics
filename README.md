# ***Pseudacris/Xantusia***-island-genomics
Following workflow is for processing raw data from several RADseq libraries from species *Pseudacris regilla* and *Xantusia riversiana*. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage and single-end. The higher coverage reads were used for generating "cleaner reads" with higher coverage and filtered by PCR duplicates. 


The rest of the libraries were single read, so they were genotyped using denovo_map.pl together with the PCR-duplicate filtered reads. Following is a step-by-step of the workflow, from raw data to many of the final analyses. 



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


##Step 2: CLEAN and prepare paired-end libraries for denovo_map.pl



####2.1. Purge PCR duplicates of PE reads


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



####2.2. Merge matching PE reads


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

####2.3. Estimating coverage of high-depth libraries using pyrad

*Within-sample clustering* in pyrad (step 3)

I started the pyrad pipeline on step#3, making sure that line # 11 (data type) is set to ***merged***:

	merged       ## 11. Datatype: rad,gbs,pairgbs,pairddrad,(others:see docs)(all)



See [*Pseudacris* parameters](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-params-PE.txt) and [*Xantusia* parameters](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-params-PE.txt) file for the rest of the specifications of parameter files for the paired end  contigs. 

pyrad was called as follows:

	pyrad -p Pr-params-d.txt -s 3


The full output for the within-sample clustering can be found here for [*Xantusia*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-within-clusters-s3.txt) and for [*Pseudacris*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-within-clusters-s3.txt). The summary of coverage per loci within sample is:


taxa	|	total	|	mean-depth	|	std-dev	|	total_d>9	|	mean_d>9	|	std-dev_d>9
----------------- | ------------- | ------------ |------------ | ------------- |------------ | ------------- |
Xr_SBI_03-PE	|	417666	|	1.499	|	1.38	|	1394	|	14.249	|	8.826
Xr_SBI_04-PE	|	581256	|	1.555	|	1.734	|	4692	|	13.598	|	8.596
Xr_SCL_27-PE	|	112121	|	3.059	|	3.52	|	4680	|	13.343	|	9.422
Xr_SCL_28-PE	|	108005	|	3.682	|	4.244	|	9479	|	13.951	|	6.378
Xr_SNI_27-PE	|	274638	|	1.463	|	1.168	|	439	|	15.59	|	9.477
Xr_SNI_28-PE	|	599975	|	1.689	|	2.169	|	9527	|	13.911	|	7.323
Xv_JTS_03-PE	|	303511	|	1.791	|	2.087	|	3701	|	13.357	|	8.52
Xv_JTS_04-PE	|	174084	|	4.967	|	6.667	|	29066	|	16.529	|	8.921
Pr_SCI_03-PE	|	193208	|	13.856	|	18.895	|	99915	|	22.382	|	23.092
Pr_SCI_04-PE	|	109581	|	9.339	|	12.522	|	39358	|	17.557	|	17.87
Pr_SCI_05-PE	|	130115	|	10.226	|	13.523	|	51999	|	18.417	|	18.323
Pr_SCI_06-PE	|	163965	|	11.35	|	16.678	|	72843	|	19.553	|	22.283
Pr_SRI_02-PE	|	193637	|	15.231	|	21.232	|	106953	|	23.762	|	25.451
Pr_SRI_03-PE	|	76988	|	8.387	|	10.075	|	24117	|	16.538	|	14.595
Pr_SRI_04-PE	|	121641	|	9.903	|	12.471	|	47253	|	17.998	|	16.814
Pr_SRI_05-PE	|	107048	|	9.322	|	13.387	|	38560	|	17.459	|	19.553


######Before moving on to denovo_map.pl, merge fasta files for library duplicates

After being renamed, move all files back to the SR-denovo-prelim folder and there I merge the fasta files. I merge following these guidelines (from this [source](http://www.researchgate.net/post/How_do_I_merge_several_multisequence-fasta_files_to_create_one_tree_for_subsequent_Unifrac_analysis)):

*To merge several files use the SHELL, go to your folder where the files are and use the cat command. E.g. to merge seqfile001.fasta, seqfile002.fasta and seqfile003.fasta type*

	cat seqfile001.fasta seqfile002.fasta seqfile003.fasta > seqcombined.fasta


*or if you have more files use*

	cat *.fasta > seqcombined.fasta


The duplicated files were sorted into a separate folder before merging, just to keep track of what's being merged. Then the post-merged files were sorted back into the general directory containing all sequences. Total number of files before merging duplicates from different ***Xantusia*** library preps was 187, and after merging duplicate individuals we now have 142 files for denovo_map input. Total number of files before merging duplicates from different ***Pseudacris*** library preps was 180, and after merging duplicate individuals we now have 132 files for denovo_map input. 

------------------------------------
######4.1.3. How many reads on average for each species?? 


We counted reads for each individual using the unzipped files and with the following script:

	echo -e 'SAMPLE_ID_FULL\tNUM_READS'
	for file in ~/path/to/denovo-map/*.fq ##edit your path!!
	do
	echo -n $(basename $file .fq)$'\t'
	cat $file | grep '^@.*' | wc -l
	done


#Step 6: *de novo* genotyping in STACKS
 


The following is the workflow/code for the denovo_map pipeline for generating a SNP matrix in .stru and genepop file formats for downstream analyses. This was done after several permutations of parameters -m -M and -n and picked based on stability of loci retrieved and population Fsts (both outputs in stacks).


####Final Code used for ***denovo_map.pl*** in Stacks:

Here, we selected a conservative combination of parameters, and used the same for both datasets. 

	
	denovo_map.pl -m 3 -M 2 -n 2 -T 16 -b 1 -t -S -o ./denovo/ -s ./



####Code used for ***populations*** in Stacks: 

Here, we used two "levels" of filters as output matrices to import into plink. The first with stringent filters:

	populations -b 1 -P ./input-sequences -M ./popmap-Pseu.txt -t 36 -p 1 -r 0.5 --write_random_snp --structure --plink --vcf --genepop --fstats

This keeps a single (random) SNP per read that is present in at least one population at a rate of 50% or higher. After this minimal filtering in populations, we filtered a few different ways in ***plink***, and the results of the permutations can be found [here](). After picking the optimal filters for each dataset, the final results of these data filters and outputs/stats can be found here for [*Pseudacris*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pseudacris-data-filters-results.md) and for [*Xantusia*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xantusia-data-filters-results.md).


***NOTE***: OUTLIER Xr_SNI_03 in Xantusia NOT a result of missing data, since it's present after these stringent filters. It seems likely that it's contamination, since in the PCA it seems closer to Santa Barbara, but it's always correctly assigned in the DAPC.... 





#downstream analyses
###a. population stats and structure
We used [R code](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/adegenet-Xantusia-NEW.R) with several packages to estimate/evaluate patterns of population structure and to estimate population-levels statistics. 

