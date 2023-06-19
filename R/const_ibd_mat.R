# This file contains helper functions for the relatedness record matching
# In particular, this file contains functions relevant for computing 
# P(R_Al | R_Bl, M): ibd.term (n.gen x n.gen matrix)
# for a given STR locus l. 
# =======================================================================
# Total match score is determined as: 
#   sum(l=1 to L) [ln(P(R_Al | S_Bl, M)) - ln(P(R_Al)]
# where L is the number of total STR loci in the data set
#
# Denote each term as follow:
#   ln(P(R_Al | S_Bl, M)): cond.ll.tot[[l]] 
#   ln(P(R_Al)): ll.tot[[l]]
# where cond.ll.tot and ll.tot are lists with l-th element is the log likelihood
# of l-th STR locus.
#
# =======================================================================
# P(R_Al | S_Bl, M) = sum(all possible STR genotype of B at locus l)
#                           P(R_Al | R_Bl, M) P(R_Bl | S_Bl)
# Denote each term as follow:
#   n.allele: number of possible distinct alleles of B at locus l 
#   n.gen: number of possible distinct genotypes of B at locus l
#   P(R_Bl | S_Bl): beagle.term (length n.gen vector)
#   P(R_Al | R_Bl, M): ibd.term (n.gen x n.gen matrix)
#                       Rows: B, Columns: A
#
# =======================================================================
# P(R_Al | R_Bl, M) = sum(k=1,9)P(R_Al | C_k, R_Bl)P(C_k | R_Bl, M)
# The sum is over all possible nine identity states C_k, defined in Lange.
# Denote each term as follow:
#   P(R_Al | C_k, R_Bl): prob.ibd.vec (length 9 vector)
#   P(C_k | R_Bl, M): prob.state.mat  (9 x n.allele matrix for homozygote)
#                                     (9 x 1 matrix for heterozygote)

const.prob.state.vec <- function(gen.type, al.freq.vec, del.vec) {
    # This function computes P(C_k | R_Bl, M): prob.state.vec (length 9 vector)
    # Input:
    #   gen.type: either "homo" for homozygote or "hetero" for heterzygote
    #   al.freq.vec: allele frequency vector corresponding to all n.al alleles
    #   del.vec: vector of probabilities of nine condensed identity states
    # Output:
    #   prob.state.mat: P(C_k | R_Bl, M) 9 x n.allele matrix for homozygote
    #                                    9 x 1 matrix for heterozygote
    
    n.al <- length(al.freq.vec)
    fb <- sum(del.vec[1:4]) # imbreeding coefficient of B
    
    if (gen.type == 'homo') {
        prob.state.mat <- matrix(NA, nrow=9, ncol=n.al)
        colnames(prob.state.mat) <- sapply(1:n.al-1, FUN=function(i) 
            paste("a.", i, 'a.', i, sep=""))
        rownames(prob.state.mat) <- sapply(1:9, 
                                           FUN=function(i) paste("D.", i, sep=""))
        prob.state.mat[1:4, ] <- sapply(1:n.al, function(r) {del.vec[1:4] / 
                (al.freq.vec[r] + (1-al.freq.vec[r]) * fb)})
        prob.state.mat[5:9, ] <- sapply(1:n.al, 
                                        function(r) {al.freq.vec[r] * del.vec[5:9] / 
                                                (al.freq.vec[r] + (1-al.freq.vec[r]) * fb)})
    } else if (gen.type == 'hetero') {
        prob.state.mat <- matrix(NA, nrow=9, ncol=1)
        colnames(prob.state.mat) <- 'a.ra.t'
        rownames(prob.state.mat) <- sapply(1:9, 
                                           FUN=function(i) paste("D.", i, sep=""))
        prob.state.mat[1:4, ] <- rep(0, 4)
        prob.state.mat[5:9, ] <- del.vec[5:9] / (1 - fb)
        
    } else {
        stop('unsupported gen.type option.')
    }
    return(prob.state.mat)
}


eval.Hrr.mat <- function(group.id, b.gen.ind, a.gen.ind,
                         al.freq.vec, del.vec, gt.df) {
    # This function evaluates P(R_Al | R_Bl=a_ra_r, Del) = H_rr
    # i.e., the ibd.term with R_Al = group.id given R_Bl is homozygote
    # and a given set of condensed identity state probabilities.
    # Input:
    #   group.id: rr, mm, rm, or mn (R_Al's genotype)
    #   b.gen.ind: index of B's genotype in gt.df (B is assumed to be homozygote with a_r a_r)
    #   a.gen.ind: index of A's genotype in gt.df
    #   al.freq.vec: allele frequency vector for alleles a_0 ... a_n.allele
    #   del.vec: vector of probabilities of nine condensed identity states
    #   gt.df: all genotypes of B in data.frame
    # Output:
    #   (b.ind.rr, a.gen.ind) element of Hrr matrix
    
    b.r.ind <- gt.df[b.gen.ind, 1] + 1
    p.r <- al.freq.vec[b.r.ind]
    
    # P(C_k | R_Bl, M) 9 x n.allele matrix 
    prob.state.mat.homo <- const.prob.state.vec('homo', al.freq.vec, del.vec) 
    
    switch(group.id,
           rr = {#  A = (a_r, a_r)
               prob.ibd.vec.homo <- c(1, p.r, p.r, p.r^2, 
                                      1, p.r, 1, p.r, p.r^2)
           },
           mm = {# A = (a_m, a_m)
               m.ind <- gt.df[a.gen.ind, 1] + 1
               p.m <- al.freq.vec[m.ind]
               prob.ibd.vec.homo <- c(0, p.m, 0, p.m^2, 
                                      0, p.m, 0, 0, p.m^2)
           },
           rm = {# A = (a_r, a_m)
               if (a.gen.ind < b.gen.ind) {m.ind <- gt.df[a.gen.ind, 2] + 1}
               if (a.gen.ind > b.gen.ind) {m.ind <- gt.df[a.gen.ind, 1] + 1}
               p.m <- al.freq.vec[m.ind]
               prob.ibd.vec.homo <- c(0, 0, p.m, 2*p.r*p.m, 
                                      0, 0, 0, p.m, 2*p.r*p.m)
           },
           mn = {# A = (a_m, a_n)
               m.ind <- gt.df[a.gen.ind, 1] + 1
               n.ind <- gt.df[a.gen.ind, 2] + 1
               p.m <- al.freq.vec[m.ind]
               p.n <- al.freq.vec[n.ind]
               prob.ibd.vec.homo <- c(0, 0, 0, 2*p.m*p.n, 
                                      0, 0, 0, 0, 2*p.m*p.n)
           },
           stop('group.id has to be one of "rr, mm, rm, mn"')
    )
    return(sum(prob.ibd.vec.homo * prob.state.mat.homo[, b.r.ind]))
    # this is (b.gen.ind, a.gen.ind) element of Hrr matrix.
}



eval.Trt.mat <- function(group.id, b.gen.ind, a.gen.ind,
                         al.freq.vec, del.vec, gt.df) {
    # This function evaluates P(R_Al | R_Bl=a_ra_t, Del) = T_rt
    # i.e., the ibd.term with R_Al = group.id given R_Bl is heterozygote
    # and a given set of condensed identity state probabilities.
    # Input:
    #   group.id: rr, mm, rm, or mn (R_Al's genotype)
    #   b.r.ind: index of B's allele index (B is assumed to be heterozygote with a_r a_t)
    #   a.gen.ind: index of A's genotype in gt.df
    #   al.freq.vec: allele frequency vector for alleles a_0 ... a_n.allele
    #   del.vec: vector of probabilities of nine condensed identity states
    #   gt.df: all genotypes of B in data.frame
    # Output:
    #   (b.gen.ind, a.gen.ind) element of Hrr matrix
    
    b.r.al <- gt.df[b.gen.ind, 1] # allele 1 of B's genotype
    b.t.al <- gt.df[b.gen.ind, 2] # allele 2 of B's genotype
    p.r <- al.freq.vec[b.r.al + 1]
    p.t <- al.freq.vec[b.t.al + 1]
    
    # P(C_k | R_Bl, M) 9 x n.allele matrix 
    prob.state.mat.hetero <- const.prob.state.vec('hetero', al.freq.vec, del.vec)
    
    switch(group.id,
           rr = {#  A = (a_r, a_r)
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0.5, p.r, 0, 0.5*p.r, p.r^2)
           },
           tt = {# A = (a_t, a_t)
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0.5, p.t, 0, 0.5*p.t, p.t^2)
           },
           mm = {# A = (a_m, a_m)
               m.ind <- gt.df[a.gen.ind, 1] + 1
               p.m <- al.freq.vec[m.ind]
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0, p.m, 0, 0, p.m^2)
           },
           rt = {# A = (a_r, a_t)
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0, 0, 1, 0.5*(p.r+p.t), 2*p.r*p.t)
           },
           rm = {# A = (a_r, a_m)
               m.ind <- gt.df[a.gen.ind, 
                              (which(gt.df[a.gen.ind, ] != b.r.al))] + 1
               if (length(m.ind) != 1) stop('check your input')
               p.m <- al.freq.vec[m.ind]
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0, 0, 0, 0.5*p.m, 2*p.r*p.m)
           },
           tm = {# A = (a_t, a_m)
               m.ind <- gt.df[a.gen.ind, 
                              (which(gt.df[a.gen.ind, ] != b.t.al))] + 1
               if (length(m.ind) != 1) stop('check your input')
               p.m <- al.freq.vec[m.ind]
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0, 0, 0, 0.5*p.m, 2*p.t*p.m)
           },
           mn = {# A = (a_m, a_n)
               m.ind <- gt.df[a.gen.ind, 1] + 1
               n.ind <- gt.df[a.gen.ind, 2] + 1
               p.m <- al.freq.vec[m.ind]
               p.n <- al.freq.vec[n.ind]
               prob.ibd.vec.hetero <- c(0, 0, 0, 0, 
                                        0, 0, 0, 0, 2*p.m*p.n)
           },
           stop('group.id has to be one of "rr, tt, mm, rt, rm, tm, mn"')
    )
    return(sum(prob.ibd.vec.hetero * prob.state.mat.hetero))
    # this is (b.gen.ind, a.gen.ind) element of Hrr matrix.
}



const.ibd.mat <- function(al.freq.vec, del.vec) {
    # This function constructs n.gen x n.gen matrix whose (i,j)-entry is
    # P(R_Al = j | R_Bl = i, M) with Rows: B, Columns: A
    # For each STR locus, this will be computed once and stored to be used later. 
    #
    # Input:
    #   al.freq.vec: allele frequency vector corresponding to all n.al alleles
    #   del.vec: M, vector of probabilities of nine condensed identity states
    # Output:
    #   ibd.mat: n.gen x n.gen matrix
    # Note: 
    #   The alleles from BEAGLE/VCF are numbered from "0". So, if there're
    #   n.al number of alleles, the alleles are: "0", "1", "2", ..., "n.al-1".
    #   The genotype is indexed from 1 to n.gen (= n.al (n.al + 1) / 2) 
    #   in the following order: 
    #   (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...
    #       , (n.al-1,0), (n.al-1,1), ..., (n.al-1, n.al-1)
    #   P(R_Al | R_Bl, M) = sum(k=1,9)P(R_Al | C_k, R_Bl)P(C_k | R_Bl, M)
    #   The sum is over all possible nine identity states C_k, defined in Lange.
    
    # Input format check
    stopifnot(length(del.vec) == 9)
    stopifnot(round(sum(al.freq.vec), digits=12) == 1)
    
    # Initialize variables
    n.al <- length(al.freq.vec)
    n.gen <- n.al * (n.al + 1) / 2
    al.vec <- c(1:n.al) - 1 # allele name vector, starts from "0". 
    gt.df <- ind.to.gt(1:n.gen) # all genotypes of B in data.frame
    fb <- sum(del.vec[1:4]) # imbreeding coefficient of B
    ibd.mat <- matrix(NA, nrow=n.gen, ncol=n.gen)
    names(al.freq.vec) <- sapply(0:(n.al-1), function(r) paste('a.', r, sep=''))
    
    # Initialize Hrr matrix: P(R_Al | R_Bl, Del), n.gen x n.gen matrix
    ibd.mat <- matrix(NA, nrow=n.gen, ncol=n.gen)
    rownames(ibd.mat) <- paste('B', gt.df$a1, gt.df$a2, sep='.')
    colnames(ibd.mat) <- paste('A', gt.df$a1, gt.df$a2, sep='.')
    
    # index of homozygote and heterozygote of B
    b.homo.ind <- sapply(al.vec, function(i) 0.5*(i+1)*(i+2)) 
    b.hetero.ind <- c(1:n.gen)[-b.homo.ind]
    
    # index of homozygote and heterozygote of A
    a.homo.ind <- sapply(al.vec, function(i) 0.5*(i+1)*(i+2)) 
    a.hetero.ind <- c(1:n.gen)[-a.homo.ind]
    
    # P(C_k | R_Bl, M) 9 x n.allele matrix 
    # P(R_Al | C_k, R_Bl): prob.ibd.vec (length 9 vector)
    
    for (al.r in al.vec) {
        r.ind <- al.r + 1 # index of allele a_r in the vector
        # p.r <- al.freq.vec[r.ind] # allele frequency of allele a_r
        
        ## ==== Case 1: B is homozygous with genotype (a_r, a_r) ====
        b.ind.rr <- b.homo.ind[r.ind] # B's genotype index for (a_r, a_r)
        # print(paste("(r,r) =", al.r, al.r))
        
        # group 1: A = (a_r, a_r)
        a.ind.rr <- b.ind.rr # index of a_ra_r genotype
        ibd.mat[b.ind.rr, a.ind.rr] <- eval.Hrr.mat('rr', b.ind.rr, a.ind.rr,
                                                    al.freq.vec, del.vec, gt.df)
        # print(paste('a.ind.rr', a.ind.rr))
        
        # group 2: A = (a_m, a_m)
        a.ind.mm <- a.homo.ind[-r.ind]
        ibd.mat[b.ind.rr, a.ind.mm] <- laply(a.ind.mm, .fun=eval.Hrr.mat, 
                                             group.id='mm', b.gen.ind=b.ind.rr, 
                                             al.freq.vec=al.freq.vec, 
                                             del.vec=del.vec, gt.df=gt.df)
        # print(paste('a.ind.mm', a.ind.mm))
        
        # group 3: A = (a_r, a_m)
        a.ind.rm <- which((gt.df$a1 == al.r | gt.df$a2 == al.r)
                          & !(gt.df$a1 == al.r & gt.df$a2 == al.r))
        ibd.mat[b.ind.rr, a.ind.rm] <- laply(a.ind.rm, .fun=eval.Hrr.mat, 
                                             group.id='rm', b.gen.ind=b.ind.rr, 
                                             al.freq.vec=al.freq.vec, 
                                             del.vec=del.vec, gt.df=gt.df)
        # print(paste('a.ind.rm', a.ind.rm))
        
        # group 4: A = (a_m, a_n)
        a.ind.mn <- which((gt.df$a1 != al.r & gt.df$a2 != al.r)
                          & (gt.df$a1 != gt.df$a2))
        ibd.mat[b.ind.rr, a.ind.mn] <- laply(a.ind.mn, .fun=eval.Hrr.mat, 
                                             group.id='mn', b.gen.ind=b.ind.rr, 
                                             al.freq.vec=al.freq.vec, 
                                             del.vec=del.vec, gt.df=gt.df)
        # print(paste('a.ind.mn', a.ind.mn))
        
        if (length(c(a.ind.rr, a.ind.mm, a.ind.rm, a.ind.mn)) != n.gen) 
            stop('error with matrix indices partition for (a_r,a_r) case')
        
        ## ==== Case 2: B is heterozygous with genotype (a_r, a_t) ====
        if (al.r > 0) {
            b.gen.ind.rt.vec <- (b.homo.ind[r.ind-1] + 1):(b.ind.rr - 1)
            # t.ind.vec <- gt.df[b.ind.rt.vec, 2]
            
            for (b.ind.rt in b.gen.ind.rt.vec) {
                al.t <- gt.df[b.ind.rt, 2] # allele t of B
                t.ind <- al.t + 1
                # print(paste("(r,t) =", al.r, al.t))
                
                # group 1: A = (a_r, a_r)
                a.ind.rr <- b.ind.rr # index of a_ra_r genotype
                ibd.mat[b.ind.rt, a.ind.rr] <- eval.Trt.mat('rr', b.ind.rt, a.ind.rr,
                                                            al.freq.vec, del.vec, gt.df)
                # print(paste('a.ind.rr', a.ind.rr))
                
                # group 2: A = (a_t, a_t)
                a.ind.tt <- which(gt.df$a1 == al.t & gt.df$a2 == al.t) # index of a_ta_t genotype
                ibd.mat[b.ind.rt, a.ind.tt] <- eval.Trt.mat('tt', b.ind.rt, a.ind.tt,
                                                            al.freq.vec, del.vec, gt.df)
                # print(paste('a.ind.tt', a.ind.tt))
                
                # group 3: A = (a_m, a_m)
                a.ind.mm <- b.homo.ind[-c(r.ind, t.ind)] # indices of a_ma_m genotype
                ibd.mat[b.ind.rt, a.ind.mm] <- plyr::laply(a.ind.mm, .fun=eval.Trt.mat, 
                                                           group.id='mm', b.gen.ind=b.ind.rt, 
                                                           al.freq.vec=al.freq.vec, 
                                                           del.vec=del.vec, gt.df=gt.df)
                # print(paste('a.ind.mm', a.ind.mm))
                
                # group 4: A = (a_r, a_t)
                a.ind.rt <- b.ind.rt # index of a_ra_r genotype
                ibd.mat[b.ind.rt, a.ind.rt] <- eval.Trt.mat('rt', b.ind.rt, a.ind.rt,
                                                            al.freq.vec, del.vec, gt.df)
                # print(paste('a.ind.rt', a.ind.rt))
                
                # group 5: A = (a_r, a_m)
                a.ind.rm <- which((gt.df$a1 == al.r | gt.df$a2 == al.r)
                                  & (gt.df$a2 != al.t & gt.df$a1 != gt.df$a2))
                ibd.mat[b.ind.rt, a.ind.rm] <- plyr::laply(a.ind.rm, .fun=eval.Trt.mat, 
                                                           group.id='rm', b.gen.ind=b.ind.rt, 
                                                           al.freq.vec=al.freq.vec, 
                                                           del.vec=del.vec, gt.df=gt.df)
                # print(paste('a.ind.rm', a.ind.rm))
                
                # group 6: A = (a_t, a_m)
                a.ind.tm <- which((gt.df$a1 == al.t | gt.df$a2 == al.t)
                                  & (gt.df$a1 != al.r & gt.df$a1 != gt.df$a2))
                ibd.mat[b.ind.rt, a.ind.tm] <- plyr::laply(a.ind.tm, .fun=eval.Trt.mat, 
                                                           group.id='tm', b.gen.ind=b.ind.rt, 
                                                           al.freq.vec=al.freq.vec, 
                                                           del.vec=del.vec, gt.df=gt.df)  
                # print(paste('a.ind.tm', a.ind.tm))
                
                # group 7: A = (a_m, a_n)
                a.ind.mn <- which((gt.df$a1 != al.r & gt.df$a2 != al.t)
                                  & (gt.df$a1 != al.t & gt.df$a2 != al.r)
                                  & (gt.df$a1 != gt.df$a2))
                if (length(a.ind.mn != 0 )) 
                    ibd.mat[b.ind.rt, a.ind.mn] <- plyr::laply(a.ind.mn, .fun=eval.Trt.mat, 
                                                               group.id='mn', b.gen.ind=b.ind.rt, 
                                                               al.freq.vec=al.freq.vec, 
                                                               del.vec=del.vec, gt.df=gt.df)
                # print(paste('a.ind.mn', a.ind.mn))
                
                
                if (length(c(a.ind.rr, a.ind.tt, a.ind.mm, a.ind.rt, 
                             a.ind.rm, a.ind.tm, a.ind.mn)) != n.gen) 
                    stop('error with matrix indices partition for (a_r,a_t) case')
            }
        }
    }
    return(ibd.mat)
}




