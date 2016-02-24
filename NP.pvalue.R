#!/usr/bin/env Rscript
# Ting Wang, 2013.12.10
# NP-method, normalize NPES and calculate p-value

########################
# inputs:
#   NPES.txt: output of NP.fg.R
#   NPES.bg.txt: output of NP.bg.R, 
# outputs:
#   NPES_result.txt: results of the NPES score, normalized score, p-value and FDR for each miRNA
########################

NPES.fg <- as.matrix(read.table('NPES.txt',sep='\t',header=T,row.names=1));
NPES.fg <- NPES.fg[,3];

NPES.bg <- as.matrix(read.table('NPES.bg.txt',sep='\t',header=F));
rownames(NPES.bg) <- names(NPES.fg);

NPES.zscore <- NPES.fg;
NPES.pvalue <- NPES.fg;
for(i in 1:length(NPES.fg)){
  NPES.zscore[i] <- (NPES.fg[i]-mean(NPES.bg[i,]))/sd(NPES.bg[i,]);
  NPES.pvalue[i] <- sum(NPES.bg[i,] >= NPES.fg[i])/ncol(NPES.bg);
}

NPES.fdr <- p.adjust(NPES.pvalue, method='fdr');

result <- cbind(NPES.fg, NPES.zscore, NPES.pvalue, NPES.fdr);

write.table(result, 'NPES_result.txt',sep='\t');
