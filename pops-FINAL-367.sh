#!/bin/bash
#SBATCH --time=10000:00:00
#SBATCH --job-name=pops-Xa-367
#SBATCH --mail-type=ALL
#SBATCH --mail-user=psalerno@colostate.edu
#SBATCH --error=stderr-pops-Xa-367
#SBATCH --output=stdout-pops-Xa-367



mkdir NEW-denovo-Xr367/final-pops


/home/patricia.salerno/programs/bin/populations populations -P ./NEW-denovo-Xr367 -M ./Xr_popmap_purged.txt -O ./NEW-denovo-Xr367/final-pops  -p 1 -r 0.5 -W Pr-whitelist --write_random_snp --structure --plink --vcf --genepop --fstats --phylip --treemix