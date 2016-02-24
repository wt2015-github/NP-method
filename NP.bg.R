#!/usr/bin/env Rscript
# Ting Wang, 2013.12.10
# NP-method, generate background NPES distributions using gene set permutation analysis

########################
# inputs:
#   genelist.net.file: output of NP.fg.R, 2 columns, first is gene name, second is score
#   p.single.matrix: a N*N matrix precalculated by RWR from single seed gene, each column is prob. and colnames are restarting genes
#   generank.file: output of NP.fg.R, ranked genes in network according to single gene based RWR
#   geneset.net.file: output of NP.fg.R, each column is a geneset, 0/1 represents geneset containing the gene or not
#   perm: number of permutation, could be executed parallelly to save time but attention to the file directions
# outputs:
#   NPES.bg.txt: background NPES
########################

NP.bg <- function(genelist.net.file='genelist.net.txt', p.single.matrix='M_htri_rwr0.5.Rdata', generank.file='gene_p.single.txt', geneset.net.file='geneset.net.txt', perm=1000){
  library('Matrix');
  
  # load single gene based rwr p matrix, p.single
  cat('load p.single matrix\n');
  load(p.single.matrix);
  # load genes ranked by single gene based rwr
  cat('load ranked genes\n');
  generank <- as.character(as.matrix(read.table(generank.file, sep='\t',header=F)));
  #p.single <- p.single[generank, generank];
  
  # load gene list file
  cat('load genelist.net file\n');
  genelist.net <- as.matrix(read.table(genelist.net.file, sep='\t', header=F, row.names=1));
  genelist.net <- genelist.net[,1];
  #genelist.net <- genelist.net[rownames(p.single)];
  
  p.single <- p.single[names(genelist.net), generank];
  
  # load genesets
  cat('load geneset.net file\n');
  geneset.net <- as.matrix(read.table(geneset.net.file, sep='\t', header=T, row.names=1));
  
  ###########################
  # calculate backgound NPES by permuting geneset
  NPES.bg <- c();
  for(i in 1:perm){
    #cat('permutation',i,'\n');
    geneset.net.perm <- geneset.net;
    rownames(geneset.net.perm) <- colnames(p.single)[sample(ncol(p.single), nrow(geneset.net))];
    NPES <- c();
    for(j in 1:ncol(geneset.net.perm)){
      cat('perm',i,'geneset',j,'\n');
      rankgene <- intersect(generank, rownames(geneset.net.perm)[geneset.net.perm[,j]==1]); ## intersect(a,b) order attention!!!
      tmp.matrix <- p.single[,rankgene];
      tmp <- apply(tmp.matrix, 1, cumsum);
      tmp <- tmp/(1:nrow(tmp)); #cumsum based mutiple gene rwr.p
      NPES.multiple <- apply(tmp, 1, function(x) cor(x, genelist.net)); #PCC score
      NPES <- c(NPES, max(NPES.multiple));
    }
    NPES.bg <- cbind(NPES.bg, NPES);
  }
  write.table(NPES.bg, 'NPES.bg.txt',sep='\t',row.names=F,col.names=F);
}

NP.bg(genelist.net.file='genelist.net.txt', p.single.matrix='M_htri_rwr0.5.Rdata', generank.file='gene_p.single.txt', geneset.net.file='geneset.net.txt', perm=1)
