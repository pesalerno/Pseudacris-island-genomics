###Channel Island demography, phylogeography, conservation genomics

1. Searching the literature for important papers regarding phylogenetic, phylogeographic, and demographic inference using RADseq-type datasets and what the issues are with different approaches. 

######Issues with missing data, type of data, pipeline for cleaning
-[This paper](http://arxiv.org/abs/1312.6439) illustrates comparisons, both with actual datasets and with simulations, of RADseq and UCE approaches. It finds RADseq is as appropriate and accurate as UCEs when using 5000SNPs. More than 5K doesn not improve estimations. However, longer reads of up to 500bp can improve estimations. Locus length has no effect on demographic parameters, only phylogenetic inference. They used custom python scripts for generating input files for programs found [here](https://github.com/mgharvey/misc_python), used STRUCTURE for population assignments, [BUCKy](http://www.stat.wisc.edu/~ane/bucky/) for estimating topological concordance across loci, and MrBayes for posterior distributions of gene tress to import into BUCKy. They inferred demographic parameters using [∂a∂i](https://code.google.com/p/dadi/). 

-[This paper](http://sysbio.oxfordjournals.org/content/early/2014/07/27/sysbio.syu046.short) simulates RADseq data to evaluate the effect of missing data on phylogenetic and phylogeographic inference with this type of data. Decisions made on levels of missing data impacts the type of loci sampled. When stringent filters are applied, not surprinsingly, loci with higher mutation rates are exlucede. "intuitive appeals about being conservative by removing loci may be misguided"!! <-- important... read further.... 

-[This paper](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291/abstract) does a really thorough evaluation of different issues when it comes to pipeline filters specifically for stacks denovo genotyping, and what the effects are in population genetic inferences. The idea is to find combinations of assembly parameters that minimize errors and maximize the retrieval of informative loci. Even though they say every dataset is different and thus will require different filters, they recommend in general to use replicates of dataset when using low-coverage datasets for more transparent reporting of genotyping error with RADseq data and the effects thereof. 

Main take homes: 

1. Having replicate individuals sequenced is super useful to get an estimate of error and of which parameters are better at retrieving "Real" loci and SNPs (excluding errors), while keeping the most amount of informative loci, and this is praticularly useful for non-model organisms, since not having a reference locus precludes strategies with sliding windows and mapping loci to references (like Catchen et al. 2013 did). *"Replicates can be used not only to estimate error rates, but also to optimize the denovo assembly of RADseq data"*.
2. In general, setting -m (minimum number of raw reads required to form a stack) too high can lead to locus/allele dropout large enough to cause incorrect inferences of individual differentiation (this should be especially true if coverage is low).
3. In general, setting -M (#mismatches allowed between stacks) too high causes differences between replicates, because it merges paralogs in a single locus forming nonsensical loci.
4. In general, parameters -m and -n (#mismatches allowed between loci) contributed the most to the variance in amount of data and missing loci (n drops a ton of loci because ).
5. In general, it seems like the default n=0 is completely incorrect, biologically speaking, for ANY of the RADseq datasets we've been getting. As in, you are expected to have some differences among individuals within the loci (SNPs!!!), and n=0 *"will result in loci represented independently across individuals that are actually alleles of the same locus. This is important for population studies where monomorphic or fixed loci may exist in different individuals. If n>0, the consensus sequence from each locus is used to attemps to merge loci, which can increase the probability of assembling real loci and decreasing error rate in genotyping"*. 



-In terms of pipelines, [PyRAD](http://dereneaton.com/software/pyrad/) seems to be the best one for trying to get a dataset that is more suited for phylogeographic and phylogenetic inference. Main difference with STACKS (and improvement over it) is the way that loci are clustered, because STACKS doesn't allow for indels, so one indel can cause complete shift in reading frame and discard obviously orthologous reads. This should make the algorithm perform better with more distantly related samples. Also, barcodes can be of varying length, which saves a lot of effort when demultiplexing. [Here](http://nbviewer.ipython.org/gist/dereneaton/af9548ea0e94bff99aa0/pyRAD_v.3.0.ipynb) is a super helpful tutorial for general PyRAD stuff, including output files and such.


######Species delimitation with RADseq and/or SNP datasets
-In [this paper](http://sysbio.oxfordjournals.org/content/63/4/534.short), Leache and company take on species delimitation with RADseq approaches (they were working on this while I was at the Albuquerque herp conference in 2013). **READ IT** and add more....

-I should also look into the program [DISSECT](http://bioinformatics.oxfordjournals.org/content/31/7/991.short), which is essentially a Bayesian discovery method for species delimitation (multispecies coalescent), but ***"without restricting the parameter space by requiring a guide tree and/or prior assignment of individuals to clusters of species"***. It's implemented within BEAST and it's very similar in terms of implementation to a *BEAST analysis. You can find the packages for download [here](http://www.indriid.com/software.html).


######General papers with a phylogeographic approach to RADseq data
-This other [paper](http://onlinelibrary.wiley.com/doi/10.1111/mec.12228/abstract?deniedAccessCustomisedMessage=&userIsAuthenticated=false) deals with phylogeographic inference of marine invertebrates using RADseq data. I can't access the article off-campus, so I'll get more details later.

-[This paper](http://onlinelibrary.wiley.com/store/10.1002/ece3.1366/asset/ece31366.pdf?v=1&t=ibpb3esg&s=304673628a522c9b16d3b39f494d51a1dd70c762) looks at phylogeographic history of caddisflies, and makes explicit tests of phylogeographic hypotheses using ABC approaches. They also use COI and "surprisingly" find that RADseq is more informative for phylogeography.... no shit Sherlock. They used STACKS for processing and genotyping, and then used the [stacks2fasta.pl](https://github.com/evoeco/radtools) to convert tsv files produced by STACKS into other formats (Arlequin, fasta, mega, nexus, etc) with the following usage:

	stacks2fasta.pl INPUT [OPTIONS] > OUTPUT

Could be super useful!

They also have another script, for obtaining reliable estimates of heterozygosity via a bootstrapping permutation test,called [arpsampler](https://github.com/evoeco/radtools) from the same place. They used [DYIABC](http://bioinformatics.oxfordjournals.org/content/early/2014/01/02/bioinformatics.btt763) for testing the biogeographic scenarios/hypotheses, but their explanations on what they did are fairly complicated, I'll need to read more carefully to understand! 

-I also need to re-read the "modern classic" [Emerson 2010](http://www.pnas.org/content/107/37/16196.full) paper, even though it's a bit outdated in terms of methods, approaches, programs, etc, it's still useful to look at more carefully.... 



######Landscape analyses

Summarize ideas between the Wang et al 2012 and the Bradburd et al 2013 papers. How are these two approaches different? Are the programs themselves easily compared, in terms of ease of analyses? 

#####*"Flexible and scalable genotyping-by-sequencing strategies for population studies"* - BMC Genomics










 


