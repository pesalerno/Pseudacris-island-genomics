output from populations with no filters, then plink filters in order

1. filter loci with less than 60% sequenced

		Options in effect:
		--ped Xr.ALL.new-z.ped
		--map Xr.SBI.map
		--geno 0.4
		--out Xr.denovo-a
		--recode
		--noweb
	
		Before frequency and genotyping pruning, there are 83445 SNPs
		141 founders and 0 non-founders found
		Total genotyping rate in remaining individuals is 0.301829
		80300 SNPs failed missingness test ( GENO > 0.4 )
		0 SNPs failed frequency test ( MAF < 0 )
		After frequency and genotyping pruning, there are 3145 SNPs

2. filter individuals that have less than 50% data 

		Options in effect:
		--ped Xr.denovo-a.ped
		--map Xr.denovo-a.map
		--mind 0.5
		--out Xr.denovo-b
		--recode
		--noweb
		Before frequency and genotyping pruning, there are 3145 SNPs
		141 founders and 0 non-founders found
		Writing list of removed individuals to [ Xr.denovo-b.irem ]
		49 of 141 individuals removed for low genotyping ( MIND > 0.5 )
		Total genotyping rate in remaining individuals is 0.815812
		0 SNPs failed missingness test ( GENO > 1 )
		0 SNPs failed frequency test ( MAF < 0 )
		After frequency and genotyping pruning, there are 3145 SNPs
	


3. filter loci genotyped in less than 50% individuals and MAF < 0.02 in remaining individuals 

		Options in effect:
		--ped Xr.denovo-b.ped
		--map Xr.denovo-b.map
		--geno 0.5
		--maf 0.02
		--out Xr.denovo-c
		--recode
		--noweb
		Before frequency and genotyping pruning, there are 3145 SNPs
		92 founders and 0 non-founders found
		Total genotyping rate in remaining individuals is 0.815812
		0 SNPs failed missingness test ( GENO > 0.5 )
		2027 SNPs failed frequency test ( MAF < 0.02 )
		After frequency and genotyping pruning, there are 1118 SNPs


Final genotyping rate:

	Total genotyping rate in remaining individuals is 0.816938

