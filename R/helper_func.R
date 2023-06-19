ind.to.gt <- function(index) {
    # Genotype index to genotype. The output is sorted by a1 first and then a2.
    # For each genotype, a1 >= a2.
    # It's in the order or (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), (3,0), ...
    a1 <- floor((sqrt(8 * index - 7) -1) / 2)
    a2 <- index - a1 * (a1 + 1) / 2 - 1
    return(data.frame('a1' = a1, 'a2' = a2))
}


gt.to.ind <- function(gt){
    # Genotype to genotpe index.
    # Input is a vector (or data frame) of length two (or dim 1x2).
    # Output is an index (scalar) that correspond to the input genotype.
    gt <- sort(gt, decreasing=TRUE, na.last=TRUE)
    ind <- gt[1] * (gt[1] + 1) / 2 + 1 + gt[2]
    return(ind)
}


read.al.freq <- function(filename) {
    # This function reads files containing allele frequencies and return
    # allele frequencies in training set allele frequency files.
    # Input:
    #   filename: filename containing allele frequencies from training set.
    #             the extension must be 'frq'.
    #   correct: T/F, default = T
    #            correcting for zero allele frequencies
    # Output:
    #   result: list with the following items:
    #       $al.freq: allele frequency vector
    #       $min.af: minimum nonzero allele frequency value BEFORE correction

    # stopifnot(strsplit(filename, "\\.")[[1]][2] == 'frq') # check input format
    al.freq.tab <- read.table(file=filename, head=FALSE, as.is=TRUE, skip=1)
    freq.char <- as.character(al.freq.tab[-(1:4)])
    al.freq <- as.numeric(matrix(unlist(strsplit(freq.char, ":")), nrow = 2)[2,])
    n.al <- al.freq.tab[1, 3]
    stopifnot(n.al == length(al.freq)) # check to see if we extracted everything

    smallest.af <- min(al.freq[al.freq > 0])
    n.zeros <- sum(al.freq == 0)
    zero.indices <- which(al.freq == 0)
    al.freq[al.freq == 0] <- (smallest.af ^ 2) / 1000
    al.freq <- al.freq / sum(al.freq)

    result <- list('al.freq'=al.freq, 'min.af'=smallest.af)
    return(result)
}


read.gt <- function(filename) {
    # This function reads files containing genotypes from BEAGLE.
    # Note that the allele indexing starts from "0" from BEAGLE.
    # Input:
    #   filename: name of the file containing genotype data
    #   phased: Whether genotype in the filename has been phased or not.
    #           Default is "T". If phased, separator is "|", and if not,
    #           the separator is "/".
    # Output:
    #   Genotypes of individuals in the file.
    #   row: individuals,
    #   columns: allele 1 and allele 2, NOT sorted

    # stopifnot(strsplit(filename, "\\.")[[1]][2] == 'GT') # check input format

    gt.str <- read.table(filename, head=TRUE, as.is=TRUE, na.strings=c("NA","-9"))

    if (length(which(strsplit(as.character(gt.str[3]), "")[[1]]=='|')) != 0) {
        split.char <- '|'
    } else {
        split.char <- '/'
    }
    gt.data <- matrix(as.numeric(unlist(sapply(gt.str[-c(1,2)], strsplit,
                                               split=split.char, fixed=TRUE))),
                      ncol=2, byrow=TRUE)
    return(gt.data)
}



read.gp <- function(filename) {
    # This function reads files containing genotype probabilities from BEAGLE.
    # Input:
    #   filename: name of the file containing genotype probability data
    # Output:
    #   Genotypes probabilities of individuals in the file.
    #   row: individuals,
    #   columns: genotype probabilities of each genotype
    #            the ordering of genotype probabilities follows the ones in
    #            the genotype file.
    # Note: The first and the second element of gp.str are chromosom number
    #       and position, respectively.

    # stopifnot(strsplit(filename, "\\.")[[1]][2] == 'GP') # check input format
    gp.str <- read.table(filename, head=TRUE, as.is=TRUE, na.strings = c("NA","-9"))
    gp.data <- matrix(as.numeric(unlist(sapply(gp.str[-c(1,2)], strsplit,
                                               split=",", fixed=TRUE))),
                      nrow=dim(gp.str)[2] - 2, byrow=TRUE)
    rownames(gp.data) <- names(gp.str)[-c(1,2)]
    return(gp.data)
}



get.ind.id <- function(vcf.f) {
    # This function takes vcf file and returns individual IDs in the file
    vcf.data <- split.vcf.file(vcf.f)
    return(strsplit(vcf.data$header, split='\t')[[1]][-c(1:9)])
}



get.snp.id <- function(vcf.f) {
    # This function takes vcf file and returns SNP IDs in the file
    vcf.data <- split.vcf.file(vcf.f)
    return(vcf.data$data[,3])
}




