cott.to.jacq <- function(cott.vec) {
    ## This function maps Cotterman identity coefficients (C_0, C_1, C_2) to
    #  9 Jacquard coefficients.
    if (length(cott.vec) != 3) {
        stop('The relatedness vector has to be length 3.')
    }
    if (sum(cott.vec) != 1) {
        stop('The Cotterman identity coefficients must sum to 1.')
    }
    if (any(cott.vec > 1) | any(cott.vec < 0)) {
        stop('The Cotterman identity coefficients must be in [0, 1].')
    }

    del.vec <- rep(0, 9)
    del.vec[9] <- cott.vec[1]
    del.vec[8] <- cott.vec[2]
    del.vec[7] <- cott.vec[3]

    return(del.vec)
}


#' Compute match score matrix based on a single STR locus.
#'
#' insert description.
#'
#' @param cott.vec A vector of Cotterman identity coefficients representing
#'                 the test relatedness hypothesis.
#' @param marker The name of a STR locus under consideration.
#' @param str.f A name (including full path) of a file containing
#'              *known (NOT imputed)* STR genotypes at the STR locus `marker`.
#'
#' @return A match score matrix based on a single STR locus. An `(i,j)` entry of
#' the match score matrix is a match score between individual `i` with a STR *observed* genotype
#' and individual `j` with an *observed* SNP profile based on a single STR locus `marker`
#' and a test relatedness hypothesis `cott.vec`.
#'
#' @note
#' `str.f` should be a `GT` file format of `vcftools`. The first column should be
#' a chromosome and the second column should be a marker genomic position.
#' The sample data starts from the third column. The first row should be a header
#' and the second row should be data. The genotype can be either a phased or unphased format.
#'
#' @export
comp.match.mat <- function(cott.vec, marker, str.f) {

    check.setup()
    base.dir <- get("base.dir", envir=.RMEnv)

    del.vec <- cott.to.jacq(cott.vec)
    mat.save.dir <- file.path(base.dir, 'result', sep='')
    if(!dir.exists(mat.save.dir)) {
        dir.create(mat.save.dir)
    }

    mat.save.dir<- file.path(mat.save.dir, gsub('\\.', '', paste(cott.vec, collapse="_")))
    if(!dir.exists(mat.save.dir)) {
        dir.create(mat.save.dir)
    }


    # Allele frequency from training set
    al.freq.filename <- file.path(base.dir, 'ref_alfrq',
                                  paste('ref_', marker, '.frq', sep=''))

    # STR genotype probability from imputation (BEAGLE) of B
    gp.filename <- file.path(base.dir, 'imputed_str',
                             paste('imp_str_', marker, '.GP.FORMAT', sep=''))

    # Matching genotype (known) for A
    gt.filename <- str.f


    ## TRUE genotypes of A, i.e. known STR genotype, NOT from imputation.
    ## The missing STR genotype is represented with NA
    val.gen.str.true <- read.gt(gt.filename)
    val.str.true.ind <- apply(val.gen.str.true, 1, gt.to.ind)
    n.ppl <- length(val.str.true.ind)

    ## Load allele frequencies from training set of A
    al.freq <- read.al.freq(al.freq.filename)
    al.freq.vec <- al.freq$al.freq
    smallest.af <- al.freq$min.af
    n.al <- length(al.freq.vec)
    n.gen <- n.al * (n.al + 1) / 2

    ## Construct HW genotype frequency matrix of A based on
    #  the training set allele frequencies of A
    outer.mat <- outer(al.freq.vec, al.freq.vec)
    hw.gf.mat <- 2 * outer.mat - diag(diag(outer.mat))
    hw.gf.vec <- hw.gf.mat[!lower.tri(hw.gf.mat)] #GP in same order as BEAGLE output
    hw.gf.min <- apply(cbind(hw.gf.vec, rep(0.0005, length(hw.gf.vec))),
                       1, min) #minimum values to be taken for imputation probabilities

    unobserved.mat <- hw.gf.mat
    unobserved.mat[hw.gf.mat >= (smallest.af^2)/10] <- 0
    unobserved.mat[hw.gf.mat < (smallest.af^2)/10] <- 2
    unobserved.vec <- unobserved.mat[!lower.tri(unobserved.mat)]
    #vector with "2" everywhere with a genotype involving
    #an allele unobserved in the training data.

    ## Step 1: Load BEAGLE imputed genotype probabilities from validation set,
    #  P(R_Bl | S_Bl) and set imputation probabilities for anything
    #  unobserved in training set to be 0, which will make the algorithm
    #  reset them to the small value in gf.min
    val.gp <- read.gp(gp.filename)
    val.gp[ , unobserved.vec == 2] <- 0

    ## Step 2: Set imputation probabilities of genotypes unobserved in
    #  the training set to be min(training set GF, 0.005) and normalize
    #  imputatin probability vector to sum to 1.
    hw.gf.min.mat <- matrix(rep(hw.gf.min, n.ppl), nrow=n.ppl, byrow=TRUE)
    scale.mat <- diag((1 - rowSums(hw.gf.min.mat * (val.gp == 0)))
                      / rowSums(val.gp))
    val.gp.scale <- scale.mat %*% val.gp
    val.gp.scale[val.gp == 0] <- hw.gf.min.mat[val.gp == 0]

    ## Compute P(R_Al | S_Bl, Del)
    ibd.term <- const.ibd.mat(al.freq.vec, del.vec)
    bgl.term <- val.gp.scale
    match.prob.mat <- bgl.term %*% ibd.term
    rownames(match.prob.mat) <- paste('S', 'B', 1:n.ppl, sep='.')

    if (any(is.na(match.prob.mat))) {
        print(paste("marker:", marker))
        stop("something went wrong in processing const.match.prob.mat function.")
    }

    ## Compute log-likelihood ratio matrix
    #  conditional.llmat: ln[P(R_Al = i | S_Bl = j, Del)]
    #                     in entry (i,j), the log probability of observing
    #                     STR genotype R_Al=i on SNP haplotype genotype S_Bl=j
    #  unconditional.llmat: ln[P(R_Al = i)]
    #                       in entry (i,j), the log probability of observing
    #                       STR genotype i, regardless of SNP
    #  lrmat: the difference of these two matrices,
    #         the log-likelihood of match vs. non-match hypotheses.

    conditional.llmat <- log(apply(match.prob.mat, 1, "[", i=val.str.true.ind))
    train.gf <- matrix(rep(hw.gf.vec, n.ppl), nrow=n.ppl, byrow=TRUE)
    unconditional.llmat <- log(apply(train.gf, 1, "[", i=val.str.true.ind))
    lrmat <- conditional.llmat - unconditional.llmat

    lrmat[is.na(lrmat)] <- 0

    ## Save llr matrix
    save.name <- file.path(mat.save.dir, paste('llr_', marker, '.RData', sep=''))
    save(lrmat, file=save.name)

    return(lrmat)
}


true.greater.than.false <- function(mat){
    trues <- diag(mat)
    falses <- c(as.numeric(mat[lower.tri(mat, diag = FALSE)]),
                as.numeric(mat[upper.tri(mat, diag = FALSE)]))
    maxfalse <- max(falses)
    return(trues > maxfalse)
}




#' Compute match accuracies.
#'
#' Compute match accuracies given a match score matrix under four schemes:
#' one-to-one (Hungarian algorithm), one-to-many:pick STR, one-to-many: pick SNP, one match.
#'
#' @param mat A match score matrix
#'
#' @return A data frame containing accuracies under four matching schemes.
#'         `one_to_one`, `SNPquery`, `STRquery`, `needle_in_haystack`.
#'
#' @details
#'
#' @export
comp.match.acc <- function(mat) {
    if (dim(mat)[1] != dim(mat)[2]) {
        stop("The input match score matrix must be a square matrix.")
    }

    #Remove anyone with missing data at all loci.
    mat <- mat[(rowSums(mat^2) != 0),(rowSums(mat^2) != 0)]

    #  One-to-one matching
    matched.LSAP <- clue::solve_LSAP(mat - min(mat) + 1, maximum=TRUE)
    hungarian.acc <- mean(matched.LSAP == 1:length(matched.LSAP))

    #  One-to-many matching
    pickSTR.acc <- mean(apply(mat, 2, which.max) == 1:length(mat[,1]))
    pickSNP.acc <- mean(apply(mat, 1, which.max) == 1:length(mat[,1]))

    #  Needle-in-hayatack matching
    onematch.acc <- mean(true.greater.than.false(mat))

    return(data.frame('one_to_one'=hungarian.acc,
                      'SNPquery'=pickSTR.acc,
                      'STRquery'=pickSNP.acc,
                      'needle_in_haystack'=onematch.acc))
}





