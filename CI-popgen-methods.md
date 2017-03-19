
#denovo genotyping for *Pseudacris* and *Xantusia*

The following is the workflow/code for the denovo_map pipeline for generating a SNP matrix in .stru and genepop file formats for downstream analyses. This was done after several permutations of parameters -m -M and -n and picked based on stability of loci retrieved and population Fsts (both outputs in stacks).


####Final Code used for ***denovo_map.pl*** in Stacks:

Here, we selected a conservative combination of parameters, and used the same for both datasets. 

	
	denovo_map.pl -m 3 -M 2 -n 2 -T 16 -b 1 -t -S -o ./denovo/ -s ./



####Code used for ***populations*** in Stacks: 

Here, we used two "levels" of filters as output matrices to import into plink. The first with stringent filters:

	populations -b 1 -P ./denovo-1-n2 -M ./popmap-Pseu-b.txt -t 36 -p 7 -r 0.5 --write_random_snp --structure --plink --vcf --genepop --fstats

For Xantusia, using ***-p 6*** (total number of populations). Thus, this filter keeps SNPs present in all populations at a rate or 50%. 

The second level/type of filter was essentially to keep all SNPs present in at least one population (-p 1) but present in at least 50% of individuals in that population (-r 0.5). So there should be less bias related to WHICH SNPs are kept (as in, less ascertainment bias that is related to phylogenetic signal).

	

I used the STRUCTURE output format from ***populations*** and transformed it to .ped in PGDSpider (the .ped output from Stacks does not save the ind. names!!)




***NOTE***: OUTLIER Xr_SNI_03 in Xantusia NOT a result of missing data, since it's present after these stringent filters. It seems likely that it's contamination, since in the PCA it seems closer to Santa Barbara, but it's always correctly assigned in the DAPC.... 



		
		
