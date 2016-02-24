#!/usr/bin/env Rscript
# Ting Wang, 2013.12.10
# NP-method, get foreground NPES

########################
# inputs:
#   genelist.file: 2 columns, first is gene name, second is score
#   p.single.matrix: a N*N matrix precalculated by RWR, the M matrix
#   geneset.file: each column is a geneset, 0/1 represents geneset containing the gene or not
# outputs:
#   NPES.txt: foreground NPES and size information
#   LE.txt: leading-edge (LE) targets
#   genelist.net.txt: intermediate file, expressions of genes in networks
#   geneset.net.txt: intermediate file, genesets containing the genes in networks
#   gene_p.single.txt: intermediate file, sorted network genes according to single seed based RWR in the farward searching procedure
#   score_rank.txt: intermediate file, sorted targets, fold change and NPES scores of each miRNA target set, every 3 lines are these information of one miRNA target set
########################

NP.fg <- function(genelist.file='DE_test.GSE4107.txt', p.single.matrix='M_htri_rwr0.5.Rdata', geneset.file='TS6.2_hsa.miRFam.mat.txt'){
  library('Matrix');
  
  ####################
  # load input files
  # load single gene based rwr p matrix, p.single
  cat('load p.single matrix\n');
  load(p.single.matrix);
  
  # read gene list file
  cat('load gene list file\n');
  genelist <- as.matrix(read.table(genelist.file, sep='\t', header=T, row.names=1))[,1];
  # get user defined score
  cat('change gene list absolute score\n');
  genelist2 <- abs(genelist);
  
  # map genes on network
  cat('map genes on network\n');
  gene.net <- intersect(rownames(p.single), names(genelist2));
  genelist.net <- genelist2[gene.net];
  
  write.table(genelist.net, 'genelist.net.txt', sep='\t', row.names=T, col.names=F);
  
  p.single <- p.single[gene.net,]; # different from version 1 !!!!
  
  # read foreground genesets
  cat('read foreground genesets\n');
  geneset.tmp <- as.matrix(read.table(geneset.file, sep='\t', header=T, row.names=1));
  geneset <- geneset.tmp[intersect(colnames(p.single), rownames(geneset.tmp)),];
  tmp <- colSums(geneset);
  geneset <- geneset[,tmp>=10]; # gene set size must be not smaller than 10
  
  write.table(geneset, 'geneset.net.txt', sep='\t');
  
  #####################################
  # NP, RWR plus forward searching
  cat('single gene based NP for all genes\n');
  NPES.single <- apply(p.single, 2, function(x) cor(x, genelist.net))
  names(NPES.single) <- colnames(p.single);
  NPES.single <- NPES.single[order(NPES.single, decreasing=T)];
  
  write.table(names(NPES.single), 'gene_p.single.txt', sep='\t',row.names=F, col.names=F); # rank all genes for accelerating following calculations
  
  cat('multiple gene based NP for genesets\n');
  n.set <- ncol(geneset);
  NPES.rank <- c();
  for(i in 1:n.set){
    gene.tmp <- rownames(geneset)[geneset[,i]==1];
    rankgene <- intersect(names(NPES.single), gene.tmp); ## intersect(a,b) order attention!!!
    
    # NPES by increasing genes
    cat('multiple gene based NP for geneset',i,'\n');
    tmp.matrix <- p.single[,rankgene];
    tmp <- apply(tmp.matrix, 1, cumsum);
    tmp <- tmp/(1:nrow(tmp)); #cumsum based mutiple gene rwr.p
    NPES.multiple <- apply(tmp, 1, function(x) cor(x, genelist.net)); #PCC score
    
    NPES.rank <- c(NPES.rank, paste(rankgene, sep='', collapse='\t'));
    NPES.rank <- c(NPES.rank, paste(NPES.single[rankgene], sep='', collapse='\t'));
    NPES.rank <- c(NPES.rank, paste(NPES.multiple, sep='', collapse='\t'));
  }
  
  writeLines(NPES.rank, 'score_rank.txt');
  
  #########################
  # get maximum NPES
  cat('calculate foreground NPES.\n');
  mir <- colnames(as.matrix(read.table('geneset.net.txt',sep='\t',header=T,row.names=1)));
  tmp <- strsplit(readLines('score_rank.txt'), '\t');
  
  NPES <- c();
  for(i in 1:length(mir)){
    NPES.multiple <- as.numeric(tmp[[i*3]]);
    NPES.max <- max(NPES.multiple);
    NPES <- rbind(NPES, c(length(NPES.multiple), which(NPES.multiple==NPES.max), NPES.max));
  }
  rownames(NPES) <- mir;
  colnames(NPES) <- c('size.set','size.maxNPES','maxNPES');
  
  write.table(NPES, 'NPES.txt', sep='\t');
  
  ##########################
  # extract LE targets for each miRs
  cat('extract LE targets for each miRs.\n');
  mir <- rownames(as.matrix(read.table('NPES.txt',sep='\t',header=T, row.names=1)));
  tmp <- strsplit(readLines('score_rank.txt'), '\t');
  
  LE <- c();
  for(i in 1:length(mir)){
    target <- as.numeric(tmp[[(i-1)*3+1]])
    NPES.single <- as.numeric(tmp[[(i-1)*3+2]]);
    NPES.multiple <- as.numeric(tmp[[i*3]]);
    
    LE.place <- which(NPES.multiple == max(NPES.multiple));
    LE.tmp <- target[1:LE.place];
    
    output.tmp <- c(mir[i], paste('LE_',LE.place,sep=''), LE.tmp);
    LE <- c(LE, paste(output.tmp, sep='\t', collapse='\t'));
  }
  
  writeLines(LE, 'LE.txt');
  
}

NP.fg(genelist.file='DE_test.GSE4107.txt', p.single.matrix='M_htri_rwr0.5.Rdata', geneset.file='TS6.2_hsa.miRFam.mat.txt')
