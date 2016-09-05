# ***Pseudacris/Xantusia***-island-genomics
Following workflow is for processing raw data from several RADseq libraries from species *Pseudacris regilla* and *Xantusia riversiana*. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage and single-end. The higher coverage reads were used for generating reference contigs. Then all of the forward reads (including initial higher coverage individuals) were mapped onto reference contigs for genotyping. The main pipelines used in this workflow are STACKS (demultiplexing, denovo mapping, and reference mapping) and pyRAD (denovo mapping and building reference contigs). 

Following is a step-by-step of the workflow, from raw data to many of the final analyses. 



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


##Step 2: Generate PE reference contigs
In order to generate references contigs with the PE reads, we first need to purge PCR duplicates and then merge the PE reads. 


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

####3.3. Runing pyrad steps 3–7

######3.3.1. *within-sample clustering* (step 3)

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

***NOTE:*** we actually set the clustering threshold differently for *Xantusia* (N>0.90) and for *Pseudacris* (N>0.93) to account for the difference in depth of seuqencing for the two libraries. 


######3.3.2. *error rate and heterozygosity estimation* (step 4):

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
Pr_SRI_03-PE	|	0.00511534	|	0.00127006	|
Pr_SRI_05-PE	|	0.00462439	|	0.00121035	|
Pr_SCI_04-PE	|	0.00491929	|	0.00126408	|
Pr_SRI_04-PE	|	0.0046619	|	0.00124301	|
Pr_SCI_05-PE	|	0.0050002	|	0.00127653	|
Pr_SCI_06-PE	|	0.00490685	|	0.00114311	|
Pr_SCI_03-PE	|	0.00516595	|	0.00118411	|
Pr_SRI_02-PE	|	0.00466262	|	0.00117116	|

######3.3.3. *within-sample consensus sequences* (step5):
This analysis uses the error rate and heterozygosity estimations/corrections from step 4. For example, for ***Xantusia***:

	pyrad -p Xr-params-t1.txt -s 5
	
	     ------------------------------------------------------------
      pyRAD : RADseq for phylogenetics & introgression analyses
     ------------------------------------------------------------


	step 5: creating consensus seqs for 8 samples, using H=0.01174 E=0.00408

And for ***Pseudacris***:

	pyrad -p Pr-params... -s 5
	
	     ------------------------------------------------------------
      pyRAD : RADseq for phylogenetics & introgression analyses
     ------------------------------------------------------------


	step 5: creating consensus seqs for 8 samples, using H=0.00488 E=0.00122
	
**From the [manual](http://nbviewer.jupyter.org/gist/dereneaton/af9548ea0e94bff99aa0/pyRAD_v.3.0.ipynb):** *"create consensus sequences. Using the mean error rate and heterozygosity estimated in step 4, this step creates consensus sequences for each cluster. Those which have less than the minimum coverage, more than the maximum number of undetermined sites, or more than the maximum number of heterozygous sites, or more than the allowed number of alleles, are discarded. In diploid data if two alleles are present the phase of heterozygous sites are retained in the consensus sequences"*

The results from this step are:

Individual | nloci | f1loci	 | f2loci | nsites | npoly | poly |
------------ | ------------- | ------------ |------------ | ------------- | ------------ |------------ | 
Xr_SBI_04-PE	|	581256	|	4692	|	4026	|	667573	|	1115	|	0.0016702
Xr_SBI_03-PE	|	417666	|	1394	|	1002	|	156595	|	616	|	0.0039337
Xr_SNI_28-PE	|	599975	|	9527	|	8685	|	1526160	|	1663	|	0.0010897
Xv_JTS_03-PE	|	303511	|	3701	|	3193	|	538514	|	1807	|	0.0033555
Xr_SNI_27-PE	|	274638	|	439	|	232	|	33780	|	306	|	0.0090586
Xr_SCL_28-PE	|	108005	|	9479	|	8907	|	1625061	|	1972	|	0.0012135
Xr_SCL_27-PE	|	112121	|	4680	|	4220	|	747718	|	1099	|	0.0014698
Xv_JTS_04-PE	|	174084	|	29066	|	27573	|	5216308	|	14625	|	0.0028037
Pr_SCI_03-PEq	|	193208	|	99881	|	92219	|	19120806	|	22398	|	0.0011714
Pr_SRI_02-PEq	|	193637	|	106908	|	99310	|	20631482	|	19323	|	0.0009366
Pr_SCI_05-PEq	|	130115	|	51993	|	47994	|	9955953	|	9521	|	0.0009563
Pr_SCI_06-PEq	|	163965	|	72826	|	67459	|	14022674	|	14215	|	0.0010137
Pr_SRI_04-PEq	|	121641	|	47250	|	43775	|	9092347	|	6878	|	0.0007565
Pr_SCI_04-PEq	|	109581	|	39350	|	36340	|	7479041	|	6557	|	0.0008767
Pr_SRI_03-PEq	|	76988	|	24117	|	22158	|	4585386	|	3087	|	0.0006732
Pr_SRI_05-PEq	|	107048	|	38552	|	35732	|	7412961	|	5491	|	0.0007407

Where:

	## nloci = number of loci
    ## f1loci = number of loci with >N depth coverage
    ## f2loci = number of loci with >N depth and passed paralog filter
    ## nsites = number of sites across f loci
    ## npoly = number of polymorphic sites in nsites
    ## poly = frequency of polymorphic sites


######3.3.4. *among-sample clustering* (step6)

**From the manual**: *"Consensus sequences are clustered across samples using the same settings as in step 3. If heterozygous, one allele is randomly sampled and used in this step, although both alleles are retained in the final data set."*

######3.3.5. *alignment and paralog filters* (step7)



The full output from step 7 post-genotyping and filtering can be found here for [*Pseudacris*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-step7-output.txt) and [*Xantusia*](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-step7-output.txt).



####3.4 Generate reference genome with bwa/samtools/picard

From a [gatk help blog](http://gatkforums.broadinstitute.org/wdl/discussion/2798/howto-prepare-a-reference-for-use-with-bwa-and-gatk), the requirements are bwa, samtools, and picard, which I need to install on the cluster. 

The steps to generating the input reads for ref_map in STACKS is:

1. Generate bwa index with fasta contigs
2. Align sequences to index
3. Transform align to sam input format



->Before anything, we need to transform the /loci output from pyrad into .fasta for generating the reference genome. Becca Tarvin helped me creating this program in python to transform the files:

	import sys

	def convert_to_fasta(file,file_out): ## def defines function, followed by 	"function_name(args):"
		## load in file
		with open(file,'r') as f:  # indented once
			lines = f.readlines()
	
	## initiate arbitrary locus number
	locus_num=0
	
	## open new file
	
	# 	prefix=file.split('.')[0]
	# 	file_out = prefix+'.fasta'
		with open(file_out,'w') as new_f:
			for line in lines:
				if line[0] == '>':
					new_f.write(line.split()[0]+'_'+str(locus_num)+'\n')
					new_f.write(line.split()[1]+'\n')
				elif line[0] == '/':
					locus_num += 1
					new_f.write('\n')

	## convention in python is to call functions from main() clause
	def main():	
	# 	print sys.argv  ## prints what you typed to command line
		file = sys.argv[1]
		file_out = sys.argv[2]
		convert_to_fasta(file,file_out)
	
	## don't forget to execute main()
	main()



Using the output fasta file we just generated with the previous python script, we can now generate the bwa index:

	./bwa index -a bwtsw Pr-output-c.fasta

where the flag -a is to specify that we want to handle a large genome. 

Then we aligned the sequences using multiple threads (in this example, 4 CPUs). 

	./bwa aln -t 4 Pr-output-c.fasta ~/scratch_dir/Pseudacris/SR-denovo-prelim/demultiplexed-reads/Pr_1 > align-bwa-Pr-1

**NOTE:**The alignment needs to be done for **EACH** sequence file that will go into ref_map.pl. Write a recursive script or a find/replace text file of all commands in sequence. 
 
We also aligned the sequences changing the default of -n 0.04 to -n 0.08 – *the default being of 4% threshold of seq. diff.* – to see effect of a lower threshold:

	./bwa aln -t 4 -n 0.08 Xr-input-bwa.fasta ~/Xantusia/SR-denovo-prelim/demultiplexed-sequences/edits/Xr_SNI_26.fq.gz.edit > Xr_SNI_26.sai

Then you transform the file which is currenty in *.sai* format to *'.sam'* format using the ***samse*** command (***sampe*** for paired-end reads):

	./bwa samse Pr-output-c.fasta align-bwa-Pr-2 /home/salerno/Pseudacris/SR-denovo-prelim/demultiplexed-reads/edits/Pr_SCA_29-4.fq.gz.edit > align-bwa-2.sam







----------------------------------------------
----------------------------------------------


##Step 4: STACKS ref_map with SR reads based on PE reference contigs


####4.1. Prep files for genotyping and check for raw read numbers
######4.1.1. Check for among-library repeats.

Libraries #1612 and #1835 are essentially duplicates, with few exceptions. Library #1834 is mostly unique with a few that are repeats from #1612. So for ***Xantusia***, I'm renaming all sequence files and adding a -1.fq.gz to library #1 (1612), -2.fq.gz to library 2 (1834), -3.fq.gz to library 3 (1835), -4.fq.gz to library 4 (1994), and -5.fq.gz to library 5 (1995). Renaming all files at once using the following code:

	rename .fq.gz -1.fq.gz *.fq.gz
	
	


------------------------------------


######4.1.2. Merge fasta files for library duplicates

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

![read counts](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/read-counts-SR.png)

/

Histogram of reads per individual for ***Pseudacris***:
![PSeu-reads](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/Pseudacris-reads.png)

/

Histogram of reads per individual for ***Xantusia***:
![Xari-reads](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/figures-reads/Xantusia-reads.png)


####4.2 Prep input files for ref_map.pl 
Generate list of sequences in sam format to be used for the input for ref_map.pl. In this step, all forward reads are used, including the ones from the initial Paired-end libraries. So they all need to be merged into one directory. 

####4.3. Run ref_map.pl in STACKS


-> We ran [**ref_map.pl**](http://catchenlab.life.illinois.edu/stacks/comp/ref_map.php) using two settings for *-m* (3 and 5) combined with the two settings in *bwa* (0.04% and 0.08% sequence difference), for a total of four reference mapping analyses for each species. The general code for running **ref_map.pl** was:


	mkdir ./Pr-ref-m3-n04/
	ref_map.pl -m 3 -T 15 -B Pr_refmap_radtags -b 1 --create_db -D "ref aligned Pseudacris populations, align 0.04, m 3" -o ./Pr-ref-m3-n04/ -O ./Pr-popmap.txt -X "populations:--fstats" \
		-s ./samples/Pr_SCA_29-4.sam \
		###plus all other sequences
		
--> in this example, we used *-m 3* in STACKS and *0.04%* in bwa.
 

We ran populations for all eight ref_map.pl analyses using this general code: 

	populations -P ./Xr-refmap-n08-m5/ -M ./popmap_Xari.txt -b 4 -t 16 -r 0.4 -p 5 -m 5 -k \
	--plink --phylip --vcf --genepop --structure --fasta --fstats --write_random_snp 
 
The results obtained for each analysis are:

#####a. *Pseudacris regilla*:

a.1. ***n 0.04, m=3:*** 1044 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n04-m3.log), [Structure output](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n04-m3.stru) and [populations Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n04-m3-FST.tsv).

a.2. ***n 0.04, m=5:*** 557 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n04-m5.log), [Structure output](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n04-m5.stru) and [populations Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n04-m5-FST.tsv).

a.3. ***n 0.08, m=3:*** 817 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n08-m3.log), [Structure output](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n08-m3.stru) and [populations Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n08-m3-FST.tsv).

a.4. ***n 0.08, m=5:*** 445 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n08-m5.log), [Structure output](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n08-m5.stru) and [populations Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Pseudacris/Pr-n08-m5-FST.tsv).

#####b. *Xantusia riversiana*:

a.1. ***n 0.04, m=3:*** 441 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n04-m3.log), [Structure output](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n04-m3.stru) and [populations Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n04-m3-FST.tsv).

a.2. ***n 0.04, m=5:*** 180 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n04-m5.log). No loci resulted from running populations (none passed filters).

a.3. ***n 0.08, m=3:*** 441 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n08-m3.log), [Structure output](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n08-m3.stru) and [populations Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n08-m3-FST.tsv).

a.4. ***n 0.08, m=5:*** 180 loci; see full outputs for [ref_map.log](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/ref-map-results/Xantusia/Xr-n08-m5.log). No loci resulted from running populations (none passed filters).


 

 




##Step 5: *de novo* genotyping in pyRAD

We performed denovo genotyping both in stacks and pyrad, using the longer reads only as forward reads (as in, without generating any references with the initial 16 individual higher coverage paired-end data). 

####pyRAD denovo genotyping

For pyrad, stats per step, starting from step 2 (since demultiplexing was done in stacks) were:

######***Pseudacris***: all 140 individuals were used for this analysis


**Step3**: *"within-sample clustering of 140 samples at 
	        '.90' similarity. Running 2 parallel jobs
	 		with up to 6 threads per job."*

To see the number of clusters and coverage per individual see [attached](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/pyrad-denovo-ALL/Pseudacris/s3.clusters.txt).  

**Step4**:

To see the results for the error and heterozygosity estimates, see [attached](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/pyrad-denovo-ALL/Pseudacris/Pi_E_estimate.txt). 

**Step5**: *"creating consensus seqs for 140 samples, using H=0.02735 E=0.00652"*

To see full resuts of the consensus sequences, see [attached](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/pyrad-denovo-ALL/Pseudacris/s5.consens.txt).


**Step6**: *clustering across 140 samples at '.90' similarity* 

		step 6: clustering across 140 samples at '.90' similarity 

	Reading file /home/salerno/Pseudacris/SR-denovo-prelim/demultiplexed-reads/clust.90/cat.haplos_ 100%
	202357374 nt in 2069325 seqs, min 80, max 148, avg 98
	Counting unique k-mers 100%
	Clustering 100%
	Sorting clusters 100%
	Writing clusters 100%
	Clusters: 196531 Size min 1, max 572, avg 10.5
	Singletons: 43820, 2.1% of seqs, 22.3% of clusters
	
######***Xantusia***: all 149 individuals were used for this analysis


**Step3**: *"within-sample clustering of 141 samples at 
	        '.90' similarity. Running 2 parallel jobs
	 	with up to 6 threads per job. If needed, 
		adjust to avoid CPU and MEM limits"*

To see the number of clusters and coverage per individual see [attached](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/pyrad-denovo-ALL/Xantusia/s3.clusters.txt).  

**Step4**:

To see the results for the error and heterozygosity estimates, see [attached](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/pyrad-denovo-ALL/Xantusia/Pi_E_estimate.txt). 

**Step5**: *"creating consensus seqs for 149 samples, using H=0.01605 E=0.00518"*

To see full resuts of the consensus sequences, see [attached](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/pyrad-denovo-ALL/Xantusia/s5.consens.txt).


**Step6**: *clustering across 149 samples at '.90' similarity* 

		step 6: clustering across 149 samples at '.90' similarity 

	Reading file /home/salerno/Xantusia/SR-denovo-prelim/demultiplexed-sequences/clust.90/cat.haplos_ 100%
	324798892 nt in 3316590 seqs, min 80, max 148, avg 98
	Counting unique k-mers 100%
	Clustering 100%
	Sorting clusters 100%
	Writing clusters 100%
	Clusters: 411718 Size min 1, max 179, avg 8.1
	Singletons: 292073, 8.8% of seqs, 70.9% of clusters



/

/

/

##Step 6: *de novo* genotyping in STACKS
 
We performed **denovo_map.pl** using the following general code: 

	mkdir denovo-1
	denovo_map.pl -m 3 -M 2 -n 1 -T 16 -b 1 -t -S -o ./denovo-1/ \
	-s ./Pr_SCA_29-4.fq.gz \

And ran the program **populations** using the following general code:

	populations -b 1 -P ./denovo-1/PARALELL-populations-p7-r5 -M ./popmap-Pseu-b.txt  -t 36 -p 6 -r 0.5 --fstats --write_random_snp --structure --genepop --vcf
	
The results for each time populations was run is: 

####*1. Pseudacris regilla*

a. p7r7: 68 SNPS, [structure](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Pseudacris/denovo-p7r7.stru) and [Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Pseudacris/denovo-p7r7-FST.tsv) from STACKS. 

b. p7r5: 2,589 SNPs, [structure](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Pseudacris/denovo-p7r5.stru) and [Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Pseudacris/denovo-p7r5-FST.tsv) from STACKS.

c. p6r5:  8,942 SNPs, [structure](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Pseudacris/denovo-p6r5.stru) and [Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Pseudacris/denovo-p6r5.FST.tsv) from STACKS.



####*2. Xantusia riversiana/vigilis*

a. p7r7: O SNPS, [populations error file](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Xantusia/stderr-pops-Xr-p7-r7) from STACKS. 

b. p7r5: 1039 SNPs, [structure](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Xantusia/denovo-p7r5.stru) and [Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Xantusia/denovo-p7r5-FST.tsv) from STACKS.

c. p6r5:  7461 SNPs, [structure](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Xantusia/denovo-p6r5.stru) and [Fst](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Stacks-denovo-results/Xantusia/denovo-p6r5-FST.tsv) from STACKS.

d. n=2, p6r5: 


####*2. Xantusia riversiana* (no *X. vigilis*)

a. n=1, p6r5:  7461 SNPs,[structure]() and [Fst]() from STACKS.

b. n=2, p6r5: 


/

/


##Step 7: plink filters and statistics for STACKS and ipyRAD outputs
We obtained plink input files from STACKS (in populations)for the matrices to be used for downstream analyses.

For ***Pseudacris regilla*** we used:

1. denovo, p=7, r=0.5, n=2
2. ref_map, n=0.04, m=3
3. ref_map, n=0.08, m=3

For ***Xantusia riversiana*** we used:

1. denovo with *X. vigilis*, p=6, r=0.5
3. denovo without *X. vigilis*, p=6, r=0.5

Plink was installed on local computer. The commands used for plink were the same for all files:

First, we make a list of the loci/individuals that are filtered due to 

	./plink --ped filename.plink.ped --map filename.plink.map --maf 0.01
	

	
	./plink --ped filename.plink.ped --map filename.plink.map
	
	 --mind 0.5 --missing --hardy --freq 
	
Linkage disequilibrium pruning:

	plink --file data --indep-pairwise 50 5 0.5 
	
	###this creates files:
	>plink.prune.in
	>plink.prune.out

They are simple lists of SNP IDs, which can then be used as argument for --extract or --exclude commands in plink.


Results from first filters are:

Pseudacris regilla:

1. denovo, p=7, r=0.5 - 1773/2589 SNPs filtered - 24/132 individuals filtered

