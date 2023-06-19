## This script creates fragmentary genomic SNP data based on 
#  the random genomic coverage model. 
#
# n : number of SNPs in genome in Gymrek's data
# p : coverage
# f : fraction of SNPs in regions of interest collectively 
#     i.e., # SNPs in 15 MB in Gymrek's data / n
# 
# 1. Total  number of SNPs sequenced is floor(np)
# 2. Choose X ~ Binomial(floor(np), f) to assign number of SNPs
#    sequenced in relevent region.
# 3. Choose the X SNPs that are sequenced in whole genome.
# 4. Based on # of snps around each CODIS, distribute X SNPs 
#    across 15 regions based on multinomial distribution. 
#    Note that snps were drawn for each individuals in the test set. 

rm(list=ls())

library(parallel)
library(foreach)
library(doParallel)
source('./R/helper_func.R')

rand.seed <- 922
set.seed(rand.seed)

n.core <- detectCores()

base.dir <- '.'
data.dir <- './data/1KGP'
save.dir <- file.path(base.dir, 'snps')

if(!dir.exists(save.dir)) {
    dir.create(save.dir)
}

# ============================================================
#                   Read CODIS loci info
# ============================================================
marker.file <- file.path(data.dir, 'marker_positions.txt')
marker.info <- read.table(marker.file, header=T, as.is=T)
n.codis <- dim(marker.info)[1]
win.size <- 1e6 # 1Mbp window size around each CODIS locus

# ============================================================
#              Read snp info around each CODIS
# ============================================================
snp.info.f <- file.path(data.dir, 'snp_win.txt')
snp.info.comb <- read.table(snp.info.f, header=T, as.is=T)
n.snp.by.codis <- as.data.frame(table(snp.info.comb$codis)[marker.info$str])
colnames(n.snp.by.codis) <- c('codis', 'freq')

# ============================================================
#            Read individual IDs in Data
# ============================================================
id.f <- file.path(data.dir, 'ind_id_all.txt')
ind.id.all <- read.table(id.f, header=F, as.is=T)
n.ind.all <- dim(ind.id.all)[1]
n.ind.test <- n.ind.all * 0.25

# ============================================================
#       Set base parameters
# ============================================================

n.snp.ref <- 27185239 # number of snps in data 
n.snp.codis <- sum(n.snp.by.codis$freq) # number of snps within windows (1MB per each codis)
# n.str.ref <- 445725   # number of strs in data
stopifnot(n.snp.codis == dim(snp.info.comb)[1])

f <- n.snp.codis / n.snp.ref #fraction of SNPs in regions of interest collectively

# ============================================================
#       Concatenate all SNPs from different chromosoms.
#       This SNPs are the ones within 1Mb window of each CODIS
# ============================================================
snp.info.dir <- file.path(data.dir, 'snp_only')
snp.all.list <- list()
snp.ind.list <- list()

counter <- 0

for (m.ind in 1:n.codis) {
    m <- marker.info$str[m.ind]
    tmp.f <- file.path(snp.info.dir, paste0(m, '_loci_info.txt'))
    tmp.data <- read.table(tmp.f, header=F, as.is=T)
    snp.all.list[[m.ind]] <- data.frame(rsid=tmp.data$V1, 
                                        codis=rep(m, dim(tmp.data)[1]),
                                        chr=rep(marker.info$chr[m.ind], dim(tmp.data)[1]),
                                        stringsAsFactors=F)
    snp.ind.list[[m.ind]] <- data.frame(start.ind=counter + 1, 
                                        end.ind=counter + dim(tmp.data)[1],
                                        stringsAsFactors=F)
    counter <- counter + dim(tmp.data)[1]
}

snp.all <- do.call('rbind', snp.all.list)
snp.ind <- do.call('rbind', snp.ind.list)

stopifnot(dim(snp.all)[1] == snp.ind$end.ind[n.codis])
rm(m.ind, m)

n.snp.win <- dim(snp.all)[1] #number of snps within 1Mb window of CODIS

write.table(snp.all, file.path(snp.info.dir, 'snp_win.txt'),
            quote=F, row.names=F, col.names=T, sep='\t')

# ============================================================
#   Construct probability vector for multinomial sampling
# ============================================================
# In the same order as marker.info
snp.prob.vec <- n.snp.by.codis$freq / n.snp.codis

# ============================================================
#               Sample Random Sets of SNPs
# ============================================================
# coverage of snps vector
cov.vec <- sort(c(seq(0.001, 0.009, by=0.001),
                  seq(0.01, 0.1, by=0.01)), decreasing=T)

n.iter <- 100

for (p in cov.vec) {
    print(paste('snp coverage =', p))
    tmp.save.dir <- file.path(save.dir, paste0('cov_', gsub('\\.', '_', p)))
    if (!dir.exists(tmp.save.dir)) {
        dir.create(tmp.save.dir)
    }
    
    n.snps <- floor(n.snp.ref * p)
    
    # Draw X, the number of snps within 1MB CODIS window, aggregated, 
    # from X ~ binomial(n.snps, f)
    X.vec <- rbinom(n.iter, n.snps, f)
    write.table(X.vec, file.path(tmp.save.dir, 'X.txt'), 
                quote=F, row.names=F, col.names=T)

    cl <- makeCluster(n.core)
    registerDoParallel(cl)
    foreach (iter.id=1:n.iter) %dopar% {
        # print(iter.id)
        set.seed(iter.id)
        tmp.save.dir.2 <- file.path(tmp.save.dir, paste0('iter_', iter.id))
        if (!dir.exists(tmp.save.dir.2)) {
            dir.create(tmp.save.dir.2)
        }
        
        X <- X.vec[iter.id]
        
        # generate every iteration iter.id.
        # row: individuals
        # col: codis loci
        n.snp.sample.mat <- rmultinom(n=n.ind.test, size=X, prob=snp.prob.vec)
        stopifnot(all(colSums(n.snp.sample.mat) == X))
        write.table(n.snp.sample.mat, 
                    file.path(tmp.save.dir, paste0('multinom_', iter.id, '.txt')),
                    quote=F, col.names=F, row.names=F)
        
        ## Process for each CODIS locus
        for (m.ind in 1:n.codis) {
            m <- marker.info$str[m.ind]
            codis.snps <- snp.all[snp.all$codis==m, 1]
            n.codis.snps <- length(codis.snps)
            n.snp.sample.vec <- n.snp.sample.mat[m.ind, ]
            
            ## Sample snps for each individual
            sample.snp.ind.list <- list()
            sample.snp.list <- list()
            
            for (ppl.ind in 1:n.ind.test) {
                # print(ppl.ind)
                tmp.snp.ind <- sort(sample.int(n=n.codis.snps,
                                               size=n.snp.sample.vec[ppl.ind], 
                                               replace=F))
                tmp.snp <- codis.snps[tmp.snp.ind]
                stopifnot(length(tmp.snp) != 0)
                sample.snp.ind.list[[ppl.ind]] <- tmp.snp.ind
                sample.snp.list[[ppl.ind]] <- tmp.snp
                rm(tmp.snp.ind, tmp.snp)
            }
            
            snp.ind.save.f <- file.path(tmp.save.dir.2, paste0(m, '_ind.RData'))
            save(sample.snp.ind.list, file=snp.ind.save.f)
            
            snp.data.save.f <- file.path(tmp.save.dir.2, paste0(m, '.RData'))
            save(sample.snp.list, file=snp.data.save.f)
            
            rm(m, m.ind, codis.snps, n.codis.snps, n.snp.sample.vec, 
               sample.snp.ind.list, snp.ind.save.f,
               sample.snp.list, snp.data.save.f)
        }
        
        rm(n.snp.sample.mat)
    }
    
    stopCluster(cl)

}


save(rand.seed, n.snp.ref, f, snp.all, snp.ind, n.snp.win,
     cov.vec, n.iter, n.snp.by.codis,
     file=file.path(save.dir, 'snps_gen.RData'))

