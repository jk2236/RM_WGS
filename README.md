# Record-matching of STR profiles with fragmentary genomic SNP data

### Description
This repository contains tools for performing record matching of STR profiles with fragmentary genomic SNP profiles, as described in [Kim and Rosenberg](https://doi.org/10.1101/2022.09.01.505545). The record matching pipeline starts with a set of reference files containing SNP-STR genotypes on reference individuals, a set of files containing STR profiles only on test individuals, a set of files containing SNP profiles only on test individuals, and a set of genetic map files. It outputs a match-score matrix for all STR-SNP pairs, with rows indicating STR profiles and columns indicating SNP profiles. It enables hypothesis tests for a variety of relatedness hypotheses based on SNP and STR profiles. For a demonstration of the record matching pipeline, please see [Example](https://github.com/jk2236/RecordMatching/tree/main/examples). The repository also includes scripts to generate simulated pedigrees and fragmentary genomic SNP data. 

### Dependencies
* [BEAGLE 4.1](https://faculty.washington.edu/browning/beagle/b4_1.html) 
* Java version 8 - required by BEAGLE. See [BEAGLE manual](https://faculty.washington.edu/browning/beagle/beagle_4.1_21Jan17.pdf) for details.
* [VCFtools](https://github.com/vcftools/vcftools). `conda` version is also available [here](https://anaconda.org/bioconda/vcftools).

### Dataset
* A phased reference SNP-STR haplotype panel of [Saini et al.](https://doi.org/10.1038/s41467-018-06694-0) from the 1000 Genomes Project phase 3 containing 2,504 individuals can be downloaded from [here](https://gymreklab.com/2018/03/05/snpstr_imputation.html). Processed data containing 1-Mb SNP windows extending 500 kb in each direction from each CODIS locus midpoint can be found under [data/1KGP](https://github.com/jk2236/RM_WGS/tree/main/data/1KGP).
* HGDP SNP-STR data containing 872 individuals can be downloaded from [here](https://rosenberglab.stanford.edu/data/edgeEtAl2017/unphased_all_vcf.zip).
* Human genetic maps. HapMap GrCh36 and GrCh37 genetic maps in PLINK format. Can be downloaded from [BEAGLE page](http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/).
* All genotypes in the reference panel must be **non-missing and phased** for BEAGLE imputation. 

### Reference
Kim J, Rosenberg NA (2022). Record-matching of STR profiles with fragmentary genomic SNP data. *bioRxiv*, 2022.09.01.505545. [10.1101/2022.09.01.505545](https://doi.org/10.1101/2022.09.01.505545).

Kim J, Edge MD, Algee-Hewitt BFB, Li JZ, Rosenberg NA (2018). Statistical detection of relatives typed with disjoint forensic and biomedical loci. *Cell*, 175(3):848-858.e6. [10.1016/j.cell.2018.09.008](https://doi.org/10.1016/j.cell.2018.09.008).

Edge MD, Algee-Hewitt BFB, Pemberton TJ, Li JA, Rosenberg NA (2017). Linkage disequilibrium matches forensic genetic records to disjoint genomic marker sets. *PNAS*, 114(22):5671-5676. [10.1073/pnas.1619944114](https://doi.org/10.1073/pnas.1619944114).
