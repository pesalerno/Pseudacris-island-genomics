***Pseudacris/Xantusia***-island-genomics
=====
Following workflow is for processing raw data from several RADseq libraries from species *Pseudacris regilla* and *Xantusia riversiana*. Two "sets" of libraries were made, one with higher depth of coverage and paired-end reads, and another with lower coverage single-end longer reads. The higher coverage reads were used for generating "cleaner reads" with higher coverage and filtered by PCR duplicates. The rest of the libraries were single read, so they were genotyped using denovo_map.pl together with the PCR-duplicate filtered reads. Following is a step-by-step of the workflow, from raw data to many of the final analyses. All of this workflow is intellectual property and **copyright of Patricia E. Salerno**, and is freely available for usage upon citation. 

*--> In Review - Journal of Biogeography*



Step 1: de-multiplexing
----
De-multiplexing was done with program [process_radtags](http://creskolab.uoregon.edu/stacks/comp/process_radtags.php) individually for each library within its directory, and renamed with sample names within Stacks ([here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/barcodes-1933.txt) is an example barcodes file). 

- Commands for process_radtags for the two Paired-end libraries were:

		process_radtags -P -p ./PE-lib-1610/ -b barcodes-1610.txt -i gzfastq - \
		o ./processed-1610/ -e sbfI -c -q -r -D
	(example for library #1610 for *Xantusia*)
- Commands for process_radtags for all the other single-end libraries were:

		process_radtags -p ./SR-lib5-1994/ -b barcodes-1994.txt -i gzfastq \ 
		-o ./processed-1994/ -e sbfI -c -q -r -D
	(example for library #1994 for shared *Xantusia* and *Pseudacris*)



Step 2: prepare paired-end libraries for denovo_map.pl
----


**2.1. Purge PCR duplicates from PE reads**


I ran the open source perl script [purge_PCR_duplicates.pl](https://github.com/claudiuskerth/scripts_for_RAD/blob/master/purge_PCR_duplicates.pl) by Claudius Kerth. It needs the perl module [Parallel::ForkManager](http://search.cpan.org/~dlux/Parallel-ForkManager-0.7.5/ForkManager.pm) since it is set up for running parallelized.

For the program to run, files need to be unzipped (.fq) and end with either
`fq_1` for the SE file or `fq_2` for the PE file. Example: `XXX.fq_1` and `XXX.fq_2`. Used `rename` script to rename all files to add `-1` or `-2` termination.

For usage, I typed:

	perl purge_PCR_duplicates.pl > logfile

while having all files to be purged within the same directory as the script. The script does not overwrite original files, but outputs new files that have been purged (and adds "purged" to file names).

The reads per individual pre and post purging of PCR duplicates, and the percent purged, were:

Individual | initial reads | porst-purge reads | percent purged |
------------ | ------------- | ------------ | ------------- | 
Pr-SCI-03	|	6,703,487	|	1,871,919	|	72	|
Pr-SCI-04	|	2,159,994	|	589,841	|	73	|
Pr-SCI-05	|	2,960,879	|	818,600	|	72	|
Pr-SCI-06	|	4,413,774	|	1,193,235	|	73	|
Pr-SRI-02	|	7,572,368	|	2,094,594	|	72	|
Pr-SRI-03	|	1,298,586	|	353,847	|	73	|
Pr-SRI-04	|	2,636,033	|	720,639	|	73	|
Pr-SRI-05	|	2,145,094	|	591,698	|	72	|
Xr-SBI-03	|	1,727,677	|	987,327	|	43	|
Xr-SBI-04	|	2,582,799	|	1,475,626	|	43	|
Xr-SCL-27	|	1,368,256	|	769,389	|	44	|
Xr-SCL-28	|	1,646,492	|	922,189	|	44	|
Xr-SNI-27	|	1,109,304	|	632,239	|	43	|
Xr-SNI-28	|	3,094,797	|	1,766,134	|	43	|
Xv-JTS-03	|	1,650,289	|	938,941	|	43	|
Xv-JTS-04	|	3,754,638	|	2,124,663	|	43	|



2.2. Merge matching PE reads
-----

Because in traditional RADseq, sequencing reads that belong to the same loci can be of different lengths, and many will have overlapping segments in the paired reads (as in, R1 and R2 will overlap) then we will merge the reads using flow-cell information so that they are processed together when making the stacks (greatly reduces computational time, and also prevents a messy analysis). We have to merge the reads with the program [PEAR](https://sco.h-its.org/exelixis/web/software/pear/).

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


Step 3: prepare single-end libraries for denovo_map.pl
---

**First, we merged fasta files for library duplicates**

Some libraries were re-sequenced, so after being renamed the fasta files were merged. ([source](http://www.researchgate.net/post/How_do_I_merge_several_multisequence-fasta_files_to_create_one_tree_for_subsequent_Unifrac_analysis)):

*To merge several files use the SHELL, go to your folder where the files are and use the cat command. E.g. to merge seqfile001.fasta, seqfile002.fasta and seqfile003.fasta type*

	cat seqfile001.fasta seqfile002.fasta seqfile003.fasta > seqcombined.fasta


*or if you have more files use*

	cat *.fasta > seqcombined.fasta


Total number of files before merging duplicates from different ***Xantusia*** library preps was 187, and after merging duplicate individuals we now have 142 files for denovo_map input. Total number of files before merging duplicates from different ***Pseudacris*** library preps was 180, and after merging duplicate individuals now had 132 files for denovo_map input. 

------------------------------------
**Then, we estimated reads per individual/species**


We counted reads for each individual using the unzipped files and with the following script:

	echo -e 'SAMPLE_ID_FULL\tNUM_READS'
	for file in ~/path/to/denovo-map/*.fq 
	do
	echo -n $(basename $file .fq)$'\t'
	cat $file | grep '^@.*' | wc -l
	done

The total of reads per individual can be found [here](). 

Step 6: *de novo* genotyping in STACKS
---
 


We performed various permutations of key parameters based on the supposed quality of data and divergence among individuals in each dataset, in order to assess which combination of parameters most likely reduces the chances of over- or under-merging loci.


**1.Permutations for *Pseudacris***

Permutations | -m | -M | -n | 
------------ | ------------- | ------------ | ------------ |
a | 3 | 2 | 2 | 
b | 5 | 2 | 2 | 
c | 7 | 2 | 2 | 
d | 3 | 2 | 3 | 
e | 5 | 2 | 3 | 
f | 7 | 2 | 3 |
g | 3 | 4 | 4 | 
h | 5 | 4 | 4 | 
i | 7 | 4 | 4 |
j | 3 | 4 | 5 | 
k | 5 | 4 | 5 | 
l | 7 | 4 | 5 |
m | 3 | 6 | 6 | 
n | 5 | 6 | 6 | 
o | 7 | 6 | 6 |
p | 3 | 6 | 7 | 
q | 5 | 6 | 7 | 
r | 7 | 6 | 7 |

**1.Permutations for *Xantusia***

Permutations | -m | -M | -n | 
------------ | ------------- | ------------ | ------------ |
a | 3 | 4 | 4 | 
b | 5 | 4 | 4 | 
c | 7 | 4 | 4 | 
d | 3 | 4 | 5 | 
e | 5 | 4 | 5 | 
f | 7 | 4 | 5 |
g | 3 | 6 | 6 | 
h | 5 | 6 | 6 | 
i | 7 | 6 | 6 |
j | 3 | 6 | 7 | 
k | 5 | 6 | 7 | 
l | 7 | 6 | 7 |
m | 3 | 8 | 8 | 
n | 5 | 8 | 8 | 
o | 7 | 8 | 8 |
p | 3 | 8 | 9 | 
q | 5 | 8 | 9 | 
r | 7 | 8 | 9 |
s | 3 | 10 | 10 | 
t | 5 | 10 | 10 | 
u | 7 | 10 | 10 |
v | 3 | 10 | 11 | 
w | 5 | 10 | 11 | 
x | 7 | 10 | 11 |


--------------



> ***Xantusia*** **parameter decisions**: Based on the above parameter permutations and corresponding results ([See here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/XANTUSIA_Param_tests_Graphs.pdf)), we chose to keep two stacks runs for downstream filters and analyses: `Xa357` and `Xa567`. 
> 
> -------------------------------
> 
> ***Pseudacris*** **parameter decisions**: Based on the above parameter permutations and corresponding results ([See here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/PSEUDACRIS_Param_tests_Graphs.pdf)), we chose to keep two stacks runs for downstream filters and analyses: `Pr323` and `Pr345`. 

--------------

**CODE FOR `DENOVO_MAP` RUNS**

The general code used in `denovo_map.pl` (executed within an `.sh` file) was: 


	denovo_map.pl -M 3 -n 2 -o ./denovo-test --popmap ./popmap-Pr.txt \
		--samples ./raw

--------------
> the **COVERAGE RESULTS** for the final `denovo_map` runs can be found here: 
> 
> - [`Xa367`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Xa_367_COVERAGE.txt)
> 
> - [`Xa567`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Xa_567-COVERAGE.txt)
> 
> - [`Pr323`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Pr_323_COVERAGE.txt)
> 
> - [`Pr345`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Pr_345_COVERAGE.txt)

exporting the initial SNP matrices
---------

We expoted the SNP matrix with minimal filter and arguments in order to obtain a .vcf file to filter in `vcftools` and `plink`. 


	mkdir pops-final-Xa567

	populations -P ./denovo-Xa567 --popmap ./popmap-Xa.txt -O ./pops-final-Xa567/ --write_random_snp -r 0.1 -t 8 --vcf --plink --structure


## SNP matrix filtering

**1.Filtering loci with too much missing data:**

		/Users/patriciasalerno/bash-programs/vcftools_0.1.13/bin/vcftools --vcf Xa367.snps.vcf --max-missing 0.75 --recode --out Xa367-b

> `Xa357`: After filtering, kept 16123 out of a possible 76247 Sites
> 
> `Xa567`: After filtering, kept 7234 out of a possible 63102 Sites
>
> `Pr323`: After filtering, kept X out of a possible X Sites
> 
> `Pr345`: After filtering, kept X out of a possible X Sites

**2.Filtering by minor allele frequency**

		/Users/patriciasalerno/bash-programs/vcftools_0.1.13/bin/vcftools --vcf Xa367-b.recode.vcf --maf 0.02 --recode --out Xa367-c
	
> `Xa357`: After filtering, kept 5084 out of a possible 16123 Sites
> 
> `Xa567`: After filtering, kept 2440 out of a possible 7234 Sites
>
> `Pr323`: After filtering, kept X out of a possible X Sites
> 
> `Pr345`: After filtering, kept X out of a possible X Sites

**3.Filtering by position:** We saw the number of times base #85-96 were found in a given SNP list using the following code: 

		cat loci-rows.txt | awk '/_90/ {count++} END {print count}'


> *Xantusia*: We decided to not eliminate any of the loci towards end of sequence due to a lack of incremental SNPs (potential error) towards end of sequence. 



**4.Filtering by individuals with too much missing data** (using `plink` and CEDIA cluster)

First, we exported the VCF matrix as a `.ped` file: 

	/Users/patriciasalerno/bash-programs/vcftools_0.1.13/bin/vcftools --vcf Xa367-c.recode.vcf --plink --out Xa365-c
	
Then we estimated individuals that had more than 50% missing data using the CEDIA cluster:

	/home/patricia.salerno/programs/plink --file Xa365-c --mind 0.5 --recode --out Xa365-d.ped  --noweb

> *Xantusia*: The list of removed individuals can be seen [here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Xantusia-removed-inds.txt). 
> *Pseudacris*: The list of removed individuals can be seen [here](link).


Re-filtering in **populations** with a whitelist of loci and individuals that passed filters
------
	
For downstream analyses, we used the final list of retained individuals and loci for each matrix as a `whitelist` in the program `populations` to ibtain final matrices and also basi population stats. The whitelist requires a file that only has the locis ID and excludes the SNP position ID. Thus, only the first string before the underscore needs to be kept. The whitelist file format is ordered as a simple text file containing one catalog locus per line: 

		3
		7
		521
		11
		46

We used the ***.map*** output from the last ***plink*** filter in Text Wrangler, and generated the populations whitelist using find and replace arguments using **grep**:


	search for \d\t(\d*)_\d*\t\d\t\d*$
	replace with \1

Based the **.irem** file obtained in *plink* we removed from the popmap (to use in populations input) the individuals that did not pass the 50% missing data filter. Now we can run populations again using the whitelist of loci and the updated popmap file for loci and individuals to retain based on the plink filters. 

	populations -b 1 -P ./ -M ./popmap.txt  -p 1 -r 0.5 -W Pr-whitelist --write_random_snp --structure --plink --vcf --genepop --fstats --phylip

	
















------------------


------------------


-------------------


---------------------
#

 #
 
#

#

#

#

#

 
 
 
 This same filtering scheme was done for each island individually for each taxon. The resulting final matrices can be found here: 
 
 *Pseudacris*: 
 - [Mainland only](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-03-22-6c-MNLND.stru)
 - [Santa Rosa Island](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-03-22-6c-SRI.stru)
 - [Santa Cruz Island](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-03-22-6c-SCI.stru)
 - [Santa Catalina Island](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-03-22-6c-SCA.stru)
 - [Santarosae paleoisland](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr-03-22-6c-ROSAE.stru)
 
 *Xantusia*:
 - [X. vigilis](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-03-22-2c-MNLND.stru) (mainland congener)
 - [Santa Barbara Island](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-03-22-2c-SBI.stru)
 - [San Clemente Island](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-03-22-2c-SCL.stru)
 - [San Nicolas Island](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-03-22-2c-SNI.stru)
	


***NOTE***: We found an odd outlier, Xr_SNI_03 in the *Xantusia* dataset that, after many iterations between downstream analyses and filters, the individual did not seem to fall as outlier as a result of missing data, so we attributed lab contamination to it. It seems likely that it's contamination, since in the PCA it seems closer to Santa Barbara, but it's always correctly assigned to San Nicolas in the DAPC.



downstream analyses
===

a. population stats and structure
---

We used [this R code](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/adegenet-Xantusia-NEW.R) that uses several packages to estimate/evaluate patterns of population structure and to estimate population-levels statistics. 

For obtaining Pi (nucleotide diversity) estimates, I re-ran ***populations*** in stacks using a whitelist of the loci that remain post-filtering in ***plink***, as such:

**First I generated whitelist using find and replace commands with grep in TextWrangler:**


	find: 		\_\d\d\n
	replace:	\n 
which fixes input files for SNP names, since they can't contain SNP position (55609_56), the regular Stacks output format, but only the actual SNP ID (55609). 

**Then, I re-ran populations using the whitelist to obtain per-population pi stats:** 

	populations -b 1 -P ./input-sequences -M ./popmap-Pseu.txt -t 36 -p 1 -r 0.5 -W whitelist-SNPs --write_random_snp --structure --plink --vcf --genepop --fstats
