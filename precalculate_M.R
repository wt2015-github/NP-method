#!/usr/bin/env Rscript
# Ting, 2012.12.18
# precalculate M matrix for random walk with restart
# RWR algorithm: P = [I-(1-r)A]^{-1}*r %*% P0  (P0 is initial prob.)
# M matrix named p.single, M = [I-(1-r)A]^{-1}*r
# users can define papameters of adjacent matrix and r by themselves

load('A_adjmat_htri.Rdata');  # adjacent matrix for HTRI network
r <- 0.5;  # RWR restart probability

library('Matrix');

p.single <- RWR_premat(A_adjmat_htri, r);  # precalculate M matrix

save(p.single, file='M_htri_rwr0.5.Rdata'); # presave M matrix


RWR_premat <- function(w, r){

###########################################################################
# inputs:    w: adjmat (weighted matrix)                                  #
#            r: restart probability ((0~1])                               #
# outputs:   p.single: [I-(1-r)A]^{-1}*r                                  #
#            gene: (gene id, then the P0 should have the same order!!)    #
# usage:     output$pre.mat %*% P0 (P0 should have the same gene order!)  #
###########################################################################

  gene <- rownames(w);
  
  #normalize each column for adjacent matrix
  cat('normalizing each column of adjacent matrix.\n');
  col.sum <- colSums(w);
  w.norm <- w;
  for(i in 1:nrow(w)){
	w.norm[i,] <- w[i,]/col.sum;
  }
  w.norm.sp <- as(w.norm,'sparseMatrix');
  
  a <- Diagonal(nrow(w)) - (1-r) * w.norm.sp;
  
  #calculate inverse matrix
  cat('calculating inverse matrix.\n');
  b <- solve(a);
  
  p.single <- b*r;
  
  p.single <- as.matrix(p.single);
  
  rownames(p.single) <- gene;
  colnames(p.single) <- gene;
  
  p.single
}
