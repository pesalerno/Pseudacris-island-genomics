Pseudacris output from populations with no filters, then plink filters in order

1. filter loci with less than 60% sequenced

		Options in effect:
		--ped Pr.all.short-z.ped
		--map Pr.all.short-z.map
		--geno 0.35  ##results directly below were done with 0.4, but rest of results are accurate
		--out Pr-03.18
		--recode
		--noweb
		Before frequency and genotyping pruning, there are 37597 SNPs
		132 founders and 0 non-founders found
		Total genotyping rate in remaining individuals is 0.388049
		35841 SNPs failed missingness test ( GENO > 0.4 )##at first I did this but then decided to filter more SNPs on the first step so that I would retain more of the poorly genotyped Legacy (mainland) populations later on.
		0 SNPs failed frequency test ( MAF < 0 )
		After frequency and genotyping pruning, there are 1756 SNPs
		

2. filter individuals that have less than 50% data 

		Options in effect:
		--ped Pr-03.18.ped
		--map Pr-03.18.map
		--mind 0.5
		--out Pr-03.18-b
		--recode
		--noweb
		Before frequency and genotyping pruning, there are 1756 SNPs
		132 founders and 0 non-founders found
		Writing list of removed individuals to [ Pr-03.18-b.irem ]
		38 of 132 individuals removed for low genotyping ( MIND > 0.5 )
		Total genotyping rate in remaining individuals is 0.804131
		0 SNPs failed missingness test ( GENO > 1 )
		0 SNPs failed frequency test ( MAF < 0 )
		After frequency and genotyping pruning, there are 1756 SNPs
		

3. filter loci genotyped in less than 50% individuals and MAF < 0.02 in remaining individuals

		Options in effect:
		--ped Pr-03.18-b.ped
		--map Pr-03.18-b.map
		--geno 0.5
		--maf 0.02
		--out Pr-03.18-c
		--recode
		--noweb
		Before frequency and genotyping pruning, there are 1756 SNPs
		94 founders and 0 non-founders found
		Total genotyping rate in remaining individuals is 0.804131
		0 SNPs failed missingness test ( GENO > 0.5 )
		68 SNPs failed frequency test ( MAF < 0.02 )
		After frequency and genotyping pruning, there are 1688 SNPs