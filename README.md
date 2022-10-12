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
> - [`Xa367`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Xa-367_coverage_B.txt)
> 
> - [`Xa567`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Xa-567_coverage_B.txt)
> 
> - [`Pr323`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Pr-323_coverage_B.txt)
> 
> - [`Pr345`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/denovo_results/Pr-345_coverage_B.txt)

exporting the initial SNP matrices
---------

We exported the SNP matrix with minimal filter and arguments in order to obtain a `.vcf` file to process in `vcftools` and filter in `plink`. 


	populations -P ./denovo-Xa567 --popmap ./popmap-Xa.txt -O ./pops-final-Xa567/ --write_random_snp -r 0.1 -t 8 --vcf --plink --structure


When then transformed the `.vcf` files into `.ped` and `.map` files using `vcftools` and the following code: 

	/Users/patriciasalerno/bash-programs/vcftools_0.1.13/bin/vcftools --vcf Xr367-populations.snps.vcf --plink --out new-Xr367
	

## SNP matrix filtering

Using the CEDIA server, we filtered the matrices in three general steps, in order to find the whitelist of loci and individuals that were retained after filters. We initially permuted the filters in order to see the ideal stringency for retaining the most amount of individuals and loci. 

**1.Filtering loci with too much missing data:**

		/home/patricia.salerno/programs/plink --file new_Pr323.plink --geno 0.3 --recode --out new_Pr323-b --noweb

> `Xa367`: Filtered sites that were present in less than 80% inds. 
> 
> 		After filtering, kept 2002 out of a possible 70607 Sites
> 
> `Xa567`: Filtered sites that were present in less than 70% inds. 
> 
> 		After filtering, kept 4591 out of a possible 56310 Sites
>
> `Pr323`: Filtered sites that were present in less than 70% inds. 
> 
> 		After filtering, kept 8676 out of a possible 147330 Sites
> 
> `Pr345`: Filtered sites that were present in less than 70% inds. 
> 
> 		After filtering, kept 6615 out of a possible 131659 Sites


**2.Filtering by individuals with too much missing data**. For all analyses we filtered out individuals that had more than 50% missing data. 

	/home/patricia.salerno/programs/plink --file new_Pr323-b --mind 0.5 --recode --out new_Pr323-c --noweb

> *Xantusia*: The list of removed individuals can be seen here for [Xr367](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/new-Xr367-c.irem) and [Xr567](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/new_new_Xr-567-c.irem). 
> 
> *Pseudacris*: The list of removed individuals can be seen here for [Pr323](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/new_Pr323-c.irem) and [Pr345](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/NEWPr-345-c.irem).


**3.Filtering by minor allele frequency**

		/home/patricia.salerno/programs/plink --file new_Pr323-b --maf 0.02 --recode --out new_Pr323-d --noweb
	
> `Xa357`: Due to low number of alleles, we used maf <0.01
> 		
> 		After frequency and genotyping pruning, there are 1111 SNPs
> 
> `Xa567`: used maf < 0.02
> 
> 		After frequency and genotyping pruning, there are 1473 SNPs
>
> `Pr323`: used maf < 0.02
> 
> 		After frequency and genotyping pruning, there are 4417 SNPs
> 
> `Pr345`: used maf < 0.02
> 
> 		After frequency and genotyping pruning, there are 3251 SNPs


We obtained the list of loci retained, to be used as the `whitelist` for the FINAL `populations` runs, using plink: 

	/home/patricia.salerno/programs/plink --file new_new_Xr-567-d --write-snplist --noweb


**4.Filtering by position:** We saw the number of times base #85-96 were found in a given SNP list using the following code: 

		cat loci-rows.txt | awk '/_90/ {count++} END {print count}'


>  We decided to not eliminate any of the loci towards end of sequence in either dataset due to a lack of incremental SNPs (potential error) towards end of sequence. 


> **FINAL MATRICES:** Based on the above permutations of filters (locus missingness, individual missingness, minor allele frequency), and with he intention of keeping as many islands and individuals as possible, we kept the following matrices for downstream analyses: 
> 
> 	- [Pr345](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Pr345.gen): maf 0.02, max_missing 0.7, mind 0.5; XXX SNPs, XX inds
> 
> 	- [Xr567-islands](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-noSNI03-islands.gen): x, x, x, XXX SNPs, XX inds
> 
> 	- [Xr567-allpops](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/Xr-noSNI03-allpops.gen): x, x, x, XXX SNPs, XX inds

	
>***NOTE***: We found an odd outlier, Xr_SNI_03 in the *Xantusia* dataset that, after many iterations between downstream analyses and filters, the individual did not seem to fall as outlier as a result of missing data, so we attributed lab contamination to it. It also had outlier heterozygosity estimates (as estimated using --het in `plink`), suggesting that it was likely the result of laboratory contamination.



Re-filtering in **populations** with a whitelist of loci and individuals that passed filters
------
	
For downstream analyses, we used the final list of retained individuals and loci for each matrix as a `whitelist` in the program `populations` to obtain final matrices and also basic population stats. The whitelist requires a file that has the locus ID on the left column and the SNP position on the right. For example:

	2070	6
	4893	6
	6115	6
	9158	6
	19214	6
	19427	6
	20960	6
	27494	6
	28614	6



Based the **.irem** file obtained in *plink* we removed from the popmap (to use in populations input) the individuals that did not pass the 50% missing data filter. 

Finally, we ran populations again using the `whitelist` of loci and the updated popmap file for loci and individuals to retain based on the plink filters. 

	populations -b 1 -P ./ -M ./popmap.txt  -p 1 -r 0.5 -W Pr-whitelist --structure --plink --vcf --genepop --fstats --treemix 



>--------------
>
***Xantusia*** whitelist loci can be found here for [`Xa357`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/whitelist_Xr-367) and [`Xa567`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/whitelist_Xr-567). 
> 
***Pseudacris*** whitelist can be found here for [`Pr323`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/whitelist-Pr-323) and [`Pr345`](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/whitelist-Pr345). 
> 
> ---------------



 
Estimating individual heterozygosities using `plink`
------

After filtering the matrix for linked loci, we estimated the individual heterozygosities using the `--het` flag in `plink`: 

	/home/patricia.salerno/programs/plink --file final-Pr345 --het --out Pr345-het --noweb
	
	
After this, we generated a distribution graph using ggplot2 in R and using the following code: 

	```R 
	insert ggplot code here




Estimating admixed trees using `TreeMix`
------

We used TreeMix using x y z parameters... etc. 

	treemix insert code here


Estimating population structure using `adegenet`
----
We used `adegenet` in `R` to estimate population structure with the use of PCA and DAPC functions. The code used for this section can be found [here](https://github.com/pesalerno/Pseudacris-island-genomics/blob/master/CIpopgen_adegenet-final-code.R).

