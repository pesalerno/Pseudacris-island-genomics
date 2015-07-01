###Channel Island demography, phylogeography, conservation genomics

1. Searching the literature for important papers regarding phylogenetic, phylogeographic, and demographic inference using RADseq-type datasets and what the issues are with different approaches. 

######Issues with missing data, type of data, pipeline for cleaning
-[This paper](http://arxiv.org/abs/1312.6439) illustrates comparisons, both with actual datasets and with simulations, of RADseq and UCE approaches. It finds RADseq is as appropriate and accurate as UCEs when using 5000SNPs. More than 5K doesn not improve estimations. However, longer reads of up to 500bp can improve estimations. Locus length has no effect on demographic parameters, only phylogenetic inference. They used custom python scripts for generating input files for programs found [here](https://github.com/mgharvey/misc_python), used STRUCTURE for population assignments, [BUCKy](http://www.stat.wisc.edu/~ane/bucky/) forestimating topological concordance across loci, and MRBayes for posterior distributions of gene tress to import into BUCKy. They inferred demographic parameters using [∂a∂i](https://code.google.com/p/dadi/). 

-[This paper](http://sysbio.oxfordjournals.org/content/early/2014/07/27/sysbio.syu046.short) simulates RADseq data to evaluate the effect of missing data on phylogenetic and phylogeographic inference with this type of data. Decisions made on levels of missing data impacts the type of loci sampled. When stringent filters are applied, not surprinsingly, loci with higher mutation rates are exlucede. "intuitive appeals about being conservative by removing loci may be misguided"!! <-- important... read further.... 

-[This paper](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291/abstract) does a really thorough evaluation of different issues when it comes to pipeline filters specifically for stacks denovo genotyping, and what the effects are in population genetic inferences. The idea is to find parameter combinations of assembly parameters that minimize errors, maximize the retrieval of informative loci. The *"model system"* is a plant from high altitude mountains in Mexico.... so very relevant for patch island-like distributions. Even though they say every dataset is different and thus will require different filters, they recommend in general to use replicates of dataset when using low-coverage datasets for more transparent reporting of genotyping error with RADseq data and the effects thereof. 


######Species delimitation with RADseq and/or SNP datasets
-In [this paper](http://sysbio.oxfordjournals.org/content/63/4/534.short), Leache and company take on species delimitation with RADseq approaches (they were working on this while I was at the Albuquerque herp conference in 2013). **READ IT** and add more....

-I should also look into the program [DISSECT](http://bioinformatics.oxfordjournals.org/content/31/7/991.short), which is essentially a Bayesian discovery method for species delimitation (multispecies coalescent), but ***"without restricting the parameter space by requiring a guide tree and/or prior assignment of individuals to clusters of species"***. It's implemented within BEAST and it's very similar in terms of implementation to a *BEAST analysis. You can find the packages for download [here](http://www.indriid.com/software.html).



-This other [paper](http://onlinelibrary.wiley.com/doi/10.1111/mec.12228/abstract?deniedAccessCustomisedMessage=&userIsAuthenticated=false) deals with phylogeographic inference of marine invertebrates using RADseq data. I can't access the article off-campus, so I'll get more details later.











 


