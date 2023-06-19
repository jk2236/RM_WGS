#  Generate n.ped number of random pedigrees

library(vcfR)
library(plyr)

source('./R/helper_func.R')


write.ped.data <- function(meta.file, vcf.header, gt.data,
                           ped.lab, phased, save.name.pre) {
    # This function writes txt and vcf file for simulated ped data
    # meta.file: txt file containing meta data line of original vcf file
    # vcf.header: header of vcf data (one line starting with #)
    # gt.data: genotype data (including first 9 fixed columns)
    # ped.lab: 'child.1', 'par.1', etc
    # phased: if gt.data is phased or not
    # save.name.pre: prefix for save name, identifier for locus

    colnames(gt.data) <- vcf.header

    if (phased == 'p') {
        save.name.post <- 'p'; save.sep='_'
    } else if (phased == 'u') {
        save.name.post <- 'u'; save.sep='_'
    } else if (phased == 'n') {
        save.name.post <- ''; save.sep=''
    } else {
        stop("unsupported phased option.")
    }

    save.name <- paste(save.name.pre, ped.lab, save.name.post, sep=save.sep)

    write.table(gt.data, file=paste(save.name, 'txt', sep='.'),
                sep='\t', row.names=F, col.names=T, quote=F)

    system(paste('cat', meta.file, paste(save.name, 'txt', sep='.'),
                 '>', paste(save.name, 'vcf', sep='.')))

}


gen.vec.to.hap <- function(gen.vec) {
    # This function converts phased genotype vector into haplotypes
    # Each column of hap is a haplotype
    hap <- matrix(unlist(strsplit(gen.vec, split='|', fixed=T), use.names=F),
                  ncol=2, byrow=T)
    return(hap)
}


unphase.geno <- function(phased.gen.mat) {
    # This function takes phased genotype MATRIX and unphase it.
    # To avoid possible bias in BEAGLE phasing & imputation step,
    # the alleles are randomly permutated at each locus during
    # unphasing process.
    # ex) 1|0 1|1 0|0 -> 0/1 1/1 0/0
    unphased.gen.mat <- matrix(NA, nrow=dim(phased.gen.mat)[1],
                               ncol=dim(phased.gen.mat)[2])

    for (i in 1:dim(phased.gen.mat)[2]) {
        tmp.gen.mat <- matrix(unlist(strsplit(phased.gen.mat[,i], split='|', fixed=T),
                                     use.names=F), ncol=2, byrow=T)
        tmp.gen.mat.rand <- apply(tmp.gen.mat, 1, sample, size=2)
        unphased.gen.mat[,i] <- paste(tmp.gen.mat.rand[1,],
                                      tmp.gen.mat.rand[2,], sep='/')
    }

    stopifnot(!any(is.na(unphased.gen.mat)))
    return(unphased.gen.mat)
}


gen.child.gt <- function(n.child=2, par.gen.vec.1, par.gen.vec.2) {
    # This function takes parents' genotype and create a child's genotype.
    # Parents' genotype should be a vector of genotypes (character string).
    # Resulting child's genotype will be in the same form as its parents' format.
    # Input:
    #   n.child: number of children to produce
    #   par.gen.vec.1, par.gen.vec.2
    #                : parent 1 & parent 2's genotype in character vector
    #                  Can be either phased (|) or unphased (/).
    #                  ex) "0/0" "1/0" "0/1" "0/0" "1/0"
    #   phased: whether parents' input genotype is phased or not.
    #
    # Output:
    #   child.gt:
    #       n.loci x n.child matrix.
    #       rows: different loci, cols: different children
    #       i.e., each column is a genotype vector of each child of the parents.
    #       the genotype is in the same format as the input parents' genotype vector.
    stopifnot(length(par.gen.vec.1) == length(par.gen.vec.2))

    par.missing.ind <- union(which(par.gen.vec.1 == '.|.'),
                             which(par.gen.vec.2 == '.|.'))

    par.gen.1 <- gen.vec.to.hap(par.gen.vec.1)
    par.gen.2 <- gen.vec.to.hap(par.gen.vec.2)

    par.1.pick <- sample.int(2, size=n.child, replace=T)
    par.2.pick <- sample.int(2, size=n.child, replace=T)

    child.gt <- matrix(NA, nrow=length(par.gen.vec.1), ncol=n.child)

    for (i in 1:n.child) {
        child.gt[,i] <- paste(par.gen.1[ , par.1.pick[i]],
                              par.gen.2[ , par.2.pick[i]], sep='|')
    }

    child.gt[par.missing.ind, ] <- c('.|.', '.|.')

    return(child.gt)

}



combine.pair <- function(gt.1, gt.2, lab.1, lab.2 ){
    # This function combines a pair's data set. For example, if we're combining
    # par.1 and par.2, the order of the data will be 
    # par.1.1, par.2.1, par.1.2, par.2.2, par.1.3, par.2.3, par.1.4, par.2.4, ... etc
    
    n.ind <- length(lab.1)
    n.locus <- dim(gt.1)[1]
    stopifnot((n.ind == length(lab.2)) & (n.ind == dim(gt.1)[2]))
    stopifnot(all.equal(dim(gt.1), dim(gt.2)))
    
    lab <- c(rbind(lab.1, lab.2))
    gt <- matrix(NA, nrow=n.locus, ncol=n.ind*2)
    for (i in 1:n.ind) {
        gt[,(2*i-1):(2*i)] <- cbind(gt.1[,i], gt.2[,i])
    }
    
    return(list('label'=lab, 'gt'=gt))
}


partition.vcf <- function(id.split, vcf.file, train.save.dir,
                          val.save.dir, save.post='') {
    # Training set: the ones labeled with "1" in id.split and will be 
    #               saved at save.dir/train/ folder
    # Validation set: the ones labeled with "3" in id.split and will be 
    #                 saved at save.dir/valid/ folder
    # Input:
    #   id.split: id indicating which individuals are in training/validation set
    #   vcf.file: vcf file to split
    #   train.save.dir: folder to save the partitioned training vcf file
    #   val.save.dir: folder to save the partitioned validation vcf file
    #   save.post: postfix for save name, identifying individual's role
    #              in the ped, e.g.: '_par', '_child', '_sib_1', etc.
    # Output:
    #   Write a vcf file for training and validation set
    
    save.name.pre <- strsplit(basename(vcf.file), split='\\.')[[1]][1]
    save.name.pre <- paste(save.name.pre, save.post, sep='')
    curr.dir <- getwd()
    vcf.data <- split.vcf.file(vcf.file)
    
    meta.str <- vcf.data$meta
    meta.file <- paste(curr.dir, 'tmp.meta.str.txt', sep='/')
    con <- file(meta.file)
    writeLines(meta.str, con)
    close(con)
    
    header.str <- strsplit(vcf.data$header, split='\t')[[1]]
    label.str <- header.str[-c(1:9)]
    header.str.fixed <- header.str[1:9]
    data.fixed <- vcf.data$data[, 1:9] 
    data.gt <- vcf.data$data[, -c(1:9)]
    
    n.loc <- dim(data.gt)[1]
    n.ind <- dim(data.gt)[2]
    n.iter <- dim(id.split)[2]
    
    for (i in 1:n.iter) {
        curr.save.name <- paste(save.name.pre, '_', i, sep='')
        curr.split <- id.split[,i]
        train.ind <- curr.split == 1
        val.ind <- curr.split == 3
        
        setwd(train.save.dir)
        write.ped.data(meta.file=meta.file,
                       vcf.header=c(header.str.fixed, label.str[train.ind]),
                       gt.data=cbind(data.fixed, data.gt[ ,train.ind]),
                       ped.lab='', phased='n',
                       curr.save.name)
        system('rm *.txt')
        
        setwd(val.save.dir)
        write.ped.data(meta.file=meta.file, 
                       vcf.header=c(header.str.fixed, label.str[val.ind]), 
                       gt.data=cbind(data.fixed, data.gt[ ,val.ind]), 
                       ped.lab='', phased='n', 
                       curr.save.name)
        # system('rm *.txt')
    }
    
    setwd(curr.dir)
    
}



partition.vcf.no.train <- function(id.split, vcf.file, 
                                   val.save.dir, save.post='') {
    # Training set: the ones labeled with "1" in id.split and will be 
    #               saved at save.dir/train/ folder
    # Validation set: the ones labeled with "3" in id.split and will be 
    #                 saved at save.dir/valid/ folder
    # Input:
    #   id.split: id indicating which individuals are in training/validation set
    #   vcf.file: vcf file to split
    #   train.save.dir: folder to save the partitioned training vcf file
    #   val.save.dir: folder to save the partitioned validation vcf file
    #   save.post: postfix for save name, identifying individual's role
    #              in the ped, e.g.: '_par', '_child', '_sib_1', etc.
    # Output:
    #   Write a vcf file for training and validation set
    
    save.name.pre <- strsplit(basename(vcf.file), split='\\.')[[1]][1]
    save.name.pre <- paste(save.name.pre, save.post, sep='')
    curr.dir <- getwd()
    vcf.data <- split.vcf.file(vcf.file)
    
    meta.str <- vcf.data$meta
    meta.file <- paste(curr.dir, 'tmp.meta.str.txt', sep='/')
    con <- file(meta.file)
    writeLines(meta.str, con)
    close(con)
    
    header.str <- strsplit(vcf.data$header, split='\t')[[1]]
    label.str <- header.str[-c(1:9)]
    header.str.fixed <- header.str[1:9]
    data.fixed <- vcf.data$data[, 1:9] 
    data.gt <- vcf.data$data[, -c(1:9)]
    
    n.loc <- dim(data.gt)[1]
    n.ind <- dim(data.gt)[2]
    n.iter <- dim(id.split)[2]
    
    for (i in 1:n.iter) {
        curr.save.name <- paste(save.name.pre, '_', i, sep='')
        curr.split <- id.split[,i]
        train.ind <- curr.split == 1
        val.ind <- curr.split == 3
        
        # setwd(train.save.dir)
        # write.ped.data(meta.file=meta.file, 
        #                vcf.header=c(header.str.fixed, label.str[train.ind]), 
        #                gt.data=cbind(data.fixed, data.gt[ ,train.ind]), 
        #                ped.lab='', phased='n', 
        #                curr.save.name)
        # system('rm *.txt')
        
        setwd(val.save.dir)
        write.ped.data(meta.file=meta.file, 
                       vcf.header=c(header.str.fixed, label.str[val.ind]), 
                       gt.data=cbind(data.fixed, data.gt[ ,val.ind]), 
                       ped.lab='', phased='n', 
                       curr.save.name)
        # system('rm *.txt')
    }
    
    setwd(curr.dir)
    
}



partition.vcf.combine.train <- function(id.split, ped.file.1, ped.file.2, 
                                        train.save.dir) {
    # Training set: the ones labeled with "1" in id.split and will be 
    # saved at save.dir/train/ folder
    # Validation set: the ones labeled with "3" in id.split and will be 
    # saved at save.dir/valid/ folder
    # Input:
    # id.split: id indicating which individuals are in training/validation set
    # train.save.dir: folder to save the partitioned training vcf file
    # val.save.dir: folder to save the partitioned validation vcf file
    # Output:
    # Write a vcf file for training and validation set
    
    stopifnot(strsplit(basename(ped.file.1), split='\\.')[[1]][1] == 
                  strsplit(basename(ped.file.2), split='\\.')[[1]][1])
    save.name.pre <- strsplit(basename(ped.file.1), split='\\.')[[1]][1]
    curr.dir <- getwd()
    data.1 <- split.vcf.file(ped.file.1)
    data.2 <- split.vcf.file(ped.file.2)
    
    meta.str <- data.1$meta
    con <- file(paste(curr.dir, 'tmp.meta.str.txt', sep='/'))
    writeLines(meta.str, con)
    close(con)
    meta.file <- paste(curr.dir, 'tmp.meta.str.txt', sep='/')
    
    header.str.1 <- strsplit(data.1$header, split='\t')[[1]]
    lab.1 <- header.str.1[-c(1:9)]
    header.str.fixed <- header.str.1[1:9]
    data.fixed <- data.1$data[, 1:9] 
    gt.1 <- data.1$data[, -c(1:9)]
    
    header.str.2 <- strsplit(data.2$header, split='\t')[[1]]
    lab.2 <- header.str.2[-c(1:9)]
    data.fixed.2 <- data.2$data[, 1:9] 
    gt.2 <- data.2$data[, -c(1:9)]
    
    stopifnot(all.equal(data.fixed, data.fixed.2))
    
    n.loc <- dim(gt.1)[1]
    n.ind <- dim(gt.1)[2]
    n.iter <- dim(id.split)[2]
    
    for (i in 1:n.iter) {
        curr.save.name <- paste(save.name.pre, '_', i, sep='')
        curr.split <- id.split[,i]
        train.ind <- curr.split == 1
        val.ind <- curr.split == 3
        
        train <- combine.pair(gt.1[ ,train.ind], gt.2[ ,train.ind],
                              lab.1[train.ind], lab.2[train.ind])
        # val <- combine.pair(gt.1[ ,val.ind], gt.2[ ,val.ind],
        #                     lab.1[val.ind], lab.2[val.ind])
        
        setwd(train.save.dir)
        write.ped.data(meta.file, vcf.header=c(header.str.fixed, train$label), 
                       gt.data=cbind(data.fixed, train$gt), 
                       ped.lab='', phased='n', 
                       curr.save.name)
        system('rm *.txt')
        
        # setwd(val.save.dir)
        # write.ped.data(meta.file, vcf.header=c(header.str.fixed, val$label), 
        #                gt.data=cbind(data.fixed, val$gt), 
        #                ped.lab='', phased='n', 
        #                curr.save.name)
        # system('rm *.txt')
    }
    
    setwd(curr.dir)
    
}

