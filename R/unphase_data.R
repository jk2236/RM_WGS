# This script is for unphasing data

rm(list=ls())

library(parallel)
library(foreach)
library(doParallel)

base.dir <- './'
vcf.all.dir <- paste(base.dir, 'phased/', sep='') # dir containing phased data
save.dir <- paste(base.dir, 'unphased/', sep='') # dir to save unphased data

if(!dir.exists(save.dir)){
    dir.create(save.dir)
}

marker.file <- './marker_positions.txt'
marker.info <- read.table(marker.file, header=T, as.is=T)

stump.w.str <- '_halfwindow500000WithCODIS.vcf' 

rand.seed <- 1
set.seed(rand.seed)


## ==========================================
## Unphase original file with STR
## ==========================================
n.core <- detectCores()
cl <- makeCluster(n.core)
registerDoParallel(cl)
result <- foreach(i=1:dim(marker.info)[1], .verbose=F) %dopar% {
    source('./R/gen_ped/helper_gen_par_child_vcf.R')
    
    m.name <- marker.info$str[i]
    print('=============================')
    print(paste('Processing', m.name))
    print('=============================')
    
    tmp.vcf.name <- paste(vcf.all.dir, m.name, stump.w.str, sep='')
    tmp.save.name <- paste(save.dir, m.name, '_halfwindow500000WithCODIS', sep='')
    
    data.phased <- split.vcf.file(tmp.vcf.name)
    meta.str.phased <- data.phased[[1]]
    header.str.phased <- data.phased[[2]]
    data.fixed <- data.phased[[3]][, 1:9]
    header.str <- strsplit(header.str.phased, split='\t')[[1]]
    ind.lab <- header.str[-c(1:9)]
    header.str.fixed <- header.str[1:9]
    
    # write meta header file for save 
    meta.file <- paste(save.dir, 'tmp.meta.str.', m.name, '.txt', sep='')
    con <- file(meta.file)
    writeLines(meta.str.phased, con)
    close(con)
    
    # unphase genotype
    data.gt <- data.phased[[3]][, 10:length(header.str)] #rows: locus, cols: individuals
    n.loc <- dim(data.gt)[1]
    
    unphased.gen.mat <- matrix(NA, nrow=dim(data.gt)[1], ncol=dim(data.gt)[2])
    
    for (i in 1:dim(data.gt)[2]) {
        
        if(i %% 100 == 0) {
            print(paste('processing', i, 'out of', dim(data.gt)[2]))
        }
        
        tmp.gen.mat <- matrix(unlist(strsplit(data.gt[,i], split='[|:]'),
                                     use.names=F), ncol=3, byrow=T)
        tmp.gen.mat.rand <- t(apply(tmp.gen.mat[,1:2], 1, sample, size=2))
        unphased.gen.mat[,i] <- paste(paste(tmp.gen.mat.rand[,1], 
                                            tmp.gen.mat.rand[,2], sep='/'),
                                      tmp.gen.mat[,3], sep=':')
    }
    
    stopifnot(!any(is.na(unphased.gen.mat)))
    
    write.ped.data(meta.file=meta.file, 
                   vcf.header=header.str, 
                   gt.data=cbind(data.fixed, unphased.gen.mat), 
                   ped.lab='', phased='n', 
                   tmp.save.name)
    
    rm(data.phased, data.gt)  
}

stopCluster(cl)


## ==========================================
## Write SNP only vcf for each STR
## ==========================================

for (i in 1:dim(marker.info)[1]) {
    m.name <- marker.info$str[i]
    print(m.name)
    vcf.name <- paste(save.dir, m.name, stump.w.str, sep='')

    # save current STR name to be excluded
    exclude.file <- paste(save.dir, 'exclude_STR.txt', sep='')
    write.table(m.name, exclude.file, quote=F, row.names=F, col.names=F)

    vcf.name <- paste(save.dir, m.name, stump.w.str, sep='')
    snp.only.vcf <- paste(save.dir, m.name, '_halfwindow500000SNPsOnly', sep='')

    # write files with SNP only
    tmp.vcftools.str <- paste("vcftools --vcf ", vcf.name,
                              " --exclude ", exclude.file,
                              " --recode --out ", snp.only.vcf, sep = "")
    system(tmp.vcftools.str)
    system(paste('mv ', snp.only.vcf, '.recode.vcf ', snp.only.vcf, '.vcf', sep=''))
}







