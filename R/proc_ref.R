#' Phase reference vcf file.
#'
#' This function phase reference vcf file (`ref.f`) with `BEAGLE`.
#'
#' @param snp.f A name (including full path) of a file containing SNPs.
#'              Must be a standard vcf format.
#' @param ref.f A name (including full path) of a file containing SNP-STR
#'              haplotypes to be used as a reference panel for imputation.
#'              Note that the reference panel must be phased and have no
#'              missing values.
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
#' The phased reference file is to be used as a reference panel for `BEAGLE` imputation.
#' The output file is saved under `ref_phased` folder within the base folder.
#' The output file name is the same as `ref.f` but with `_phs` postfix.
#'
#' @seealso
#' \href{https://faculty.washington.edu/browning/beagle/beagle_4.1_21Jan17.pdf}{BEAGLE manual}.
#'
#' @export
phase.ref <- function(ref.f,
                      nthreads=1, niterations=10,
                      lowmem='false', window=50000, overlap=3000,
                      impute='true', cluster=0.005, gprobs='false',
                      ne=1000000, err=0.0001, seed=-99999, modelscale=0.8) {

    check.setup()
    base.dir <- get("base.dir", envir=.RMEnv)
    bgl.jar <- get("bgl.jar", envir=.RMEnv)

    phased.save.dir <- file.path(base.dir, 'ref_phased')
    if (!dir.exists(phased.save.dir)) {
        dir.create(phased.save.dir)
    }

    out.pre <- file.path(phased.save.dir,
                         paste(basename(tools::file_path_sans_ext(ref.f)),
                               '_phs', sep=''))

    cat("\n---- Phasing reference panel ----\n\n")

    beag.str <- paste("java -Xmx2g -jar ", bgl.jar,
                      " gt=", ref.f,
                      " out=", out.pre,
                      " nthreads=", as.character(nthreads),
                      " niterations=", as.character(niterations),
                      " lowmem=", as.character(lowmem),
                      " window=", as.character(window),
                      " overlap=", as.character(overlap),
                      " impute=", as.character(impute),
                      " cluster=", as.character(cluster),
                      " gprobs=", as.character(gprobs),
                      " ne=", as.character(ne),
                      " err=", as.character(err),
                      " seed=", as.character(seed),
                      " modelscale=", as.character(modelscale),
                      sep = "")
    system(beag.str)

    curr.dir <- getwd()
    setwd(phased.save.dir)
    system("gunzip *.gz")
    setwd(curr.dir)

    cat(strrep("=", 70))
    cat("\n    The phased reference file has been saved to:\n")
    cat(strrep(" ", 5), paste0(out.pre, '.vcf'), "\n")
    cat(strrep("=", 70))
    cat("\n")

}



#' Compute reference STR allele frequencies from a vcf file.
#'
#' This function extract STR allele frequences from reference vcf file (`ref.f`) with `vcftools`.
#'
#' @param ref.f A name (including full path) of a file containing SNP-STR
#'              haplotypes to be used as a reference panel for imputation.
#'              Note that the reference panel must be phased and have no
#'              missing values.
#' @param marker STR marker name of which allele frequency is computed.
#'
#' @return None
#'
#' @details
#' The STR allele frequency is necessary for computing match scores.
#' The output file is saved under `ref_alfrq` folder within the base folder.
#' The output file name is `ref_(marker).frq` format.
#'
#' @seealso
#' \href{https://github.com/vcftools/vcftools}{vcftools manual}.
#'
#' @export
ref.al.freq <- function(ref.f, marker) {
    check.setup()
    base.dir <- get("base.dir", envir=.RMEnv)
    vcf.exe <- get("vcf.exe", envir=.RMEnv)

    al.save.dir <- file.path(base.dir, 'ref_alfrq')
    if (!dir.exists(al.save.dir)) {
        dir.create(al.save.dir)
    }

    out.pre <- file.path(al.save.dir, paste('ref_', marker, sep=''))

    cat("\n---- Computing STR allele frequency in the reference panel ----\n\n")
    vcftools.str <- paste(vcf.exe, " --vcf ", ref.f,
                          " --snp ", marker,
                          " --freq --out ", out.pre, sep = "")
    system(vcftools.str)

    cat("\n")
    cat(strrep("=", 70))
    cat("\n", strrep(" ", 3), "The allele frequency of STR marker", marker,
        "\n", strrep(" ", 3), "in the reference panel has been saved to: \n")
    cat(strrep(" ", 5), paste0(out.pre, '.frq'), "\n")
    cat(strrep("=", 70))
    cat("\n")
}














