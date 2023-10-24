##LCL gene express
library(Rsubread)
library(DESeq)
###
samplenames <- c("NA18486","NA18498","NA18499","NA18501","NA18502","NA18504","NA18505","NA18507","NA18508",
                 "NA18510","NA18511","NA18516","NA18517","NA18519","NA18522","NA18523","NA18852","NA18855",
                 "NA18856","NA18858","NA18861","NA18862","NA18870","NA18907","NA18909","NA18912","NA18913",
                 "NA18916","NA19092","NA19093","NA19098","NA19099","NA19101","NA19102","NA19108","NA19114",
                 "NA19116","NA19119","NA19127","NA19128","NA19130","NA19131","NA19137","NA19138","NA19140",
                 "NA19143","NA19144","NA19147","NA19152","NA19153","NA19159","NA19160","NA19192","NA19193",
                 "NA19200","NA19204","NA19209","NA19222","NA19239","NA19257")
Input_data <- vector()
for(i in 1:length(samplenames)){
  Input_data[i] <- paste0("./Lymphoblastoid_cell_line/",samplenames[i],"/",samplenames[i],"Aligned.sortedByCoord.out.bam")
}

gtf <- "/home/disk3/zhangteng/hg19_GTF/genes.gtf"
gene_count <- featureCounts(Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_name", annot.ext = gtf, isPairedEnd=F, nthreads=20)

###
counts_data <- gene_count$counts
colnames(counts_data) <-samplenames
countMatrix <- sapply(as.matrix(counts_data), as.numeric)
countmatrixgene <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
colnames(countmatrixgene) <- colnames(counts_data)
rownames(countmatrixgene) <- rownames(counts_data)
conds <- factor(colnames(counts_data))
cds <- newCountDataSet(countmatrixgene, conds )
size_factor <-  sizeFactors(estimateSizeFactors( cds )) 
gene_expr <- list(countmatrixgene,size_factor)
names(gene_expr) <- c("gene_reads","size_factor")
save(gene_expr,file = "./Lymphoblastoid_cell_line/gene_expr.Rdata")
