# NP-method
A network propagation based method for the inference of perturbed microRNA regulatory networks

R codes (Version 2013.12.13)

[Homepage @ Github](http://wt2015-github.github.io/NP-method/) | [Homepage @ Jin Gu's Lab](http://bioinfo.au.tsinghua.edu.cn/member/jgu/np/) | [Source Code](https://github.com/wt2015-github/NP-method)

## Introduction
This network propagation based method (NP-method), takes advantage of the global network differential information, and carries out a integrated random walk plus forward searching algorithm to infer the perturbed miRNAs as well as their leading-edge target genes using gene differential expression data.

## Details and Usage
### R scripts (should be executed one bye one):
1 . **precalculate_M.R** to precalculate and prestore M matrix for NP-method, skip if using the default M matrix.

2 . **NP.fg.R** to get foreground NPES.

3 . **NP.bg.R** to get background NPES distributions using gene set permutation analysis, to save time users can execute parallelly, but please pay attention to the file directions.

4 . **NP.pvalue.R** to normalize NPES and get p-value using the outputs of NP.fg.R and NP.bg.R.

### Input files:
1 . **A_adjmat_htri.Rdata**: A adjacent matrix for the HTRI network.

2 . **M_htri_rwr0.5.Rdata**: A N*N matrix, named p.single, which is precalculated as the M matrix of NP-method using precalculate.R (with parameters HTRI network and r=0.5).

3 . **TS6.2_hsa.miRFam.mat.txt**: A matrix representing the TargetScan v6.2 annotated conserved targets of miRNA families, each column is a geneset, 0/1 represents geneset containing the gene or not;

4 . **DE_test.GSE4107.txt**: A test gene differential expression file, 2 columns, first is gene name, second is score.

### Output files:
1 . Output of **precalculate_M.R**

* **M matrix file**: containing a matrix variable named p.single.

2 . Outputs of **NP.fg.R**

* **NPES.txt**: foreground NPES and size information

* **LE.txt**: leading-edge (LE) targets

* **genelist.net.txt**: intermediate file, expressions of genes in networks

* **geneset.net.txt**: intermediate file, genesets containing the genes in networks

* **gene_p.single.txt**: intermediate file, sorted network genes according to single seed based RWR in the farward searching procedure

* **score_rank.txt**: intermediate file, sorted targets, fold change and NPES scores of each miRNA target set, every 3 lines are these information of one miRNA target set

3 . Outputs of **NP.bg.R**

* **NPES.bg.txt**: background NPES

4 . Output of **NP.pvalue.R**

* **NPES_result.txt**: results of the NPES score, normalized score, p-value and FDR for each miRNA

**Note**: some NA values in some outputs are reasonable (including the warnings during runing R), because not all genes in network have expression measurements.

## Publication
Ting Wang, Jin Gu, Yanda Li. "Inferring the perturbed microRNA regulatory networks from gene expression data using a network propagation based method." *BMC bioinformatics* 15.1 (2014): 255. PMID: 25069957.

## Contact
[Ting Wang](http://wt2015-github.github.io/) (wang9ting@gmail.com), [Jin Gu](http://bioinfo.au.tsinghua.edu.cn/member/jgu/) (wellgoo@gmail.com).
