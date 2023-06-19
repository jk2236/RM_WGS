#' Impute STR genotypes from surrounding SNPs and extract imputed STR genotypes and genotype probabilites.
#'
#' This function imputes genotypes of individuals at a STR `marker` locus from its
#' surrounding SNPs `snp.f` using a reference panel `ref.f` and a genetic map `map.f`.
#'
#' @param snp.f A name (including full path) of a file containing SNPs.
#'              Must be a standard vcf format.
#' @param ref.f A name (including full path) of a file containing SNP-STR
#'              haplotypes to be used as a reference panel for imputation.
#'              *Note that the reference panel must be phased and have no
#'              missing values.*
#' @param map.f A name (including full path) of a file genetic map in plink format.
#' @param marker Name of a STR locus to be imputed.
#' @param nthreads (positive integer) Number of threads for BEAGLE, default=1.
#' @param niterations (positive integer) Number of phasing iterations in BEAGLE, default=10.
#' @param maxlr (number ≥ 1) The maximum likelihood ratio at a genotype, default=1000000.
#' @param lowmem (true/false) Whether a memory efficient algorithm should be used, default=false.
#' @param window (positive integer) The number of markers to include in each sliding window, default=50000.
#' @param overlap (positive integer) The number of markers of overlap between sliding windows, default=3000.
#' @param cluster (non-negative number) The maximum cM distance between individual markers that are
#'                combined into an aggregate marker when imputing ungenotyped markers, default=0.005.
#' @param ne (integer) The effective population size when imputing ungenotyped markers, default=1000000.
#' @param err (nonnegative number) The allele miscall rate, default=0.0001.
#' @param seed (integer) The seed for the random number generator, default=-99999.
#' @param modelscale (positive number) the model scale parameter when sampling haplotypes
#'                   for unrelated individuals, default=0.8.
#'
#' @return None
#'
#' @details
#' This function imputes genotypes of individuals at a STR `marker` locus from its
#' surrounding SNPs `snp.f` using a reference panel `ref.f` and a genetic map `map.f`.
#' *The reference panel must be phased and have no missing values.*
#' The imputed and phased SNP-STR vcf file is stored at `(base.dir)/imputed_str` where
#' the `base.dir` was set by running \code{\link{setup}} in the beginning of the pipeline.
#' From the imputed SNP-STR vcf, the function extracts imputed genotypes and genotype probabilities
#' at the STR `marker` locus and save it in the `(base.dir)/imputed_str` folder.
#' The output file name containing imputed genotypes probabilities is `imp_str_(marker).GP.FORMAT`.
#'
#' @seealso
#' \href{https://faculty.washington.edu/browning/beagle/beagle_4.1_21Jan17.pdf}{BEAGLE manual}.
#' \href{https://github.com/vcftools/vcftools}{vcftools manual}.
#'
#' @export
impute.str <- function(snp.f, ref.f, map.f, marker,
                       nthreads=1, niterations=10,
                       maxlr=1000000, lowmem='false', window=50000,
                       overlap=3000, cluster=0.005, ne=1000000,
                       err=0.0001, seed=-99999, modelscale=0.8) {

    check.setup()
    base.dir <- get("base.dir", envir=.RMEnv)
    bgl.jar <- get("bgl.jar", envir=.RMEnv)

    imputed.save.dir <- file.path(base.dir, 'imputed_str')
    if (!dir.exists(imputed.save.dir)) {
        dir.create(imputed.save.dir)
    }

    out.pre <- file.path(imputed.save.dir,
                         paste(basename(tools::file_path_sans_ext(snp.f)),
                               '_imp', sep=''))

    cat(paste("\n---- Imputing STR", marker, "from SNPs ----\n\n"))

    beag.str <- paste("java -Xmx2g -jar ", bgl.jar,
                      " gt=", snp.f,
                      " out=", out.pre,
                      " ref=", ref.f,
                      " map=", map.f,
                      " gprobs=true impute=true",
                      " nthreads=", as.character(nthreads),
                      " niterations=", as.character(niterations),
                      " maxlr=", as.character(maxlr),
                      " lowmem=", as.character(lowmem),
                      " window=", as.character(window),
                      " overlap=", as.character(overlap),
                      " cluster=", as.character(cluster),
                      " ne=", as.character(ne),
                      " err=", as.character(err),
                      " seed=", as.character(seed),
                      " modelscale=", as.character(modelscale),
                      sep = "")

    system(beag.str)

    curr.dir <- getwd()
    setwd(imputed.save.dir)
    system("gunzip *.gz")
    setwd(curr.dir)


    out.pre2 <- file.path(imputed.save.dir, paste('imp_str_', marker, sep=''))

    vcftools.str.1 <- paste("vcftools --vcf ", paste(out.pre, '.vcf', sep=''),
                            " --snp ", marker, " --recode --recode-INFO-all",
                            " --out ", out.pre, sep = "")
    system(vcftools.str.1)

    vcftools.str.2 <- paste("vcftools --vcf ",
                            paste(out.pre, '.recode.vcf', sep=''),
                            " --extract-FORMAT-info GP",
                            " --out ", out.pre2, sep = "")
    system(vcftools.str.2)

    cat("\n")
    cat(strrep("=", 70))
    cat("\n    The imputed vcf file has been saved to:\n")
    cat(strrep(" ", 5), paste0(out.pre, '.vcf'), "\n")

    cat("\n    The imputed STR genotype probabilities have been saved to:\n")
    cat(strrep(" ", 5), paste0(out.pre2, '.GP.FORMAT'), "\n")

    cat(strrep("=", 70))
    cat("\n")

}










