# m6Aexpress-enet
## Introduction
As the most abundant mRNA modification, m6A controls and influences many aspects of mRNA metabolism including the mRNA stability and degradation. However, the role of specific m6A sites in regulating gene expression still remains unclear. In additional, the multicollinearity problem caused by the correlation of methylation level of multiple m6A sites in each gene could influence the prediction performance. To address the above challenges, we propose an enet-regularization negative binomial regression model (called m6Aexpress-enet) to predict which m6A site could potentially regulate its gene expression. 
# The step-by-step implementation of m6Aexpress-enet
## Installed the Required R packages
Before using m6Aexpress-enet, we should firstly install some required R packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsamtools","GenomicAlignments","GenomicRanges",
                       "GenomicFeatures","rtracklayer","exomePeak2", "Rsubread","Biostrings","BSgenome.Hsapiens.UCSC.hg19"))
install.packages("https://www.bioconductor.org/packages/3.8/bioc/src/contrib/DESeq_1.34.1.tar.gz", repos = NULL, type="source")
install.packages("randomForest")
```
## Dowloand SRAMP for predicting m6A sites in single base resolution
1. Go to the webserver: https://www.cuilab.cn/sramp/ and download SRAMP tool

2. Require the Perl (version>5.8)

## Get gene expression
```r
######Import data
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
##### Get annotation file
gtf <- ".hg19_GTF/genes.gtf"
#### Get gene expression
library(Rsubread)
output_dir <- "./Lymphoblastoid_cell_line/gene_express/"
get_gene_expre <- get_gene_express(Input_data,isGTFAnnotationFile=TRUE,GTF.featureType="exon",GTF.attrType="gene_name", annot.ext = gtf, isPairedEnd=F, nthreads=20,output=output_dir)
```

## Peak calling by exomePeak2
```r
######Import data
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

IP_data <- vector()
for(i in 1:length(samplenames)){
  IP_data[i] <- paste0("./Lymphoblastoid_cell_line/IP_samples/",samplenames[i],"_IPAligned.sortedByCoord.out.bam")
}
##### Get annotation file
GENE_ANNO_GTF  <- ".hg19_GTF/genes.gtf"
#### peak calling
library(exomePeak2)
IP_BAM <- c(IP_data)
INPUT_BAM <- c(Input_data)
result =exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff = GENE_ANNO_GTF,
           genome = "hg19",
           paired_end = FALSE,
           parallel = 20,
           save_dir="./Lymphoblastoid_cell_line/exomePeak2_output")
```
## Get the sequence of m6A peak
```r
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
import_file <-  "./Lymphoblastoid_cell_line/exomePeak2_output/Mod.csv"
output_dir <- "./Lymphoblastoid_cell_line/sequences/"
## get the sequence of m6A peak
get_peakseq <- get_peak_seq(f1=import_file,output=output_dir)
```
## Predict m6A sites in single base by SRAMP
```r
#### Obtain m6A sites in single base by SRAMP
cd ./SRAMP/
nohup perl runsramp.pl ./Lymphoblastoid_cell_line/sequences/new_motif_peak_seq.fa ./Lymphoblastoid_cell_line/exomePeak2_output/exomePeak2_result/singlebase_m6Asites.txt mature &

#### Mapping m6A sites to longest transcriptom
f1 <- "./Lymphoblastoid_cell_line/exomePeak2_output/exomePeak2_result/singlebase_m6Asites.txt"
f2 <- "./exomePeak2/exomePeak2_output/Mod.bed"
outputfile <- "./Lymphoblastoid_cell_line/exomePeak2_output/exomePeak2_result/singlebase_m6Asitesmaps.txt"
mapped_LTX_m6A <- mapm6A_LTX(f1=f1,f2=f2,output_file=outputfile) 
```
## Quantifly the methylation level of single-base m6A sites
```r
####Firstly,we should obtain the reads count of each m6A sites in single base
library(GenomicRanges)
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

IP_data <- vector()
for(i in 1:length(samplenames)){
  IP_data[i] <- paste0("./Lymphoblastoid_cell_line/IP_samples/",samplenames[i],"_IPAligned.sortedByCoord.out.bam")
}
ip_bams <- c(IP_data)
input_bams <- c(Input_data)

f1 <- "./exomePeak2_result/singlebase_m6Asitesmaps.txt"
output_file <-  "./Lymphoblastoid_cell_line/exomePeak2_output/exomePeak2_result/IP_Input_SBreadsinfor.Rdata"
get_SBm6A_reads <- getSB_m6A_readscount(ip_bams=ip_bams,input_bams=input_bams,m6A_site_file=f1,output_file=output_file)
```

