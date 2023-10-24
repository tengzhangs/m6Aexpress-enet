library(reshape2)
library(ggplot2)
f1 <- "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/exomePeakoutA/ADDInfo/ADDInfo_ReadsCount.csv"
SB_reads<- read.csv(f1)
f1 <- "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/exomePeakoutB/ADDInfo/ADDInfo_ReadsCount.csv"
SB_reads<- read.csv(f1)
colnames(SB_reads)[1] <- "peak_num"
samplenames <- c("NA18486","NA18498","NA18499","NA18501","NA18502","NA18504","NA18505","NA18507","NA18508",
                 "NA18510","NA18511","NA18516","NA18517","NA18519","NA18522","NA18523","NA18852","NA18855",
                 "NA18856","NA18858","NA18861","NA18862","NA18870","NA18907","NA18909","NA18912","NA18913",
                 "NA18916","NA19092","NA19093","NA19098","NA19099","NA19101","NA19102","NA19108","NA19114",
                 "NA19116","NA19119","NA19127","NA19128","NA19130","NA19131","NA19137","NA19138","NA19140",
                 "NA19143","NA19144","NA19147","NA19152","NA19153","NA19159","NA19160","NA19192","NA19193",
                 "NA19200","NA19204","NA19209","NA19222","NA19239","NA19257")

IP_reads <- SB_reads[,grep("IP",colnames(SB_reads))]
colnames(IP_reads) <- paste0(samplenames,"_IP")

Input_reads <- SB_reads[,-c(1,grep("IP",colnames(SB_reads)))]
colnames(Input_reads) <- paste0(samplenames,"_Input")

f2 <- "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/exomePeakoutA/Mod.csv"
SB_sites <- read.csv(f2)
f2 <- "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/exomePeakoutB/Mod.csv"
SB_sites <- read.csv(f2)
SB_sites_infor <- cbind(SB_sites[,c(1:4,6,13)],IP_reads,Input_reads)
###################
f1 <- "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/singlebase_m6Asitesmaps.txt"
SB_sites <- read.table(f1,header = T)
load("/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/IP_Input_SBreadsinfor.Rdata")
SB_sites_infor <- cbind(SB_sites[,c(1:3,5,8)],IP_Input_SB_readsinfor$IP_Input_SB_readscount[,-1])

###################
Input_reads <- SB_sites_infor[,grep("Input",colnames(SB_sites_infor))]
select_SB <- !apply(Input_reads, 1, function(x) any(x==0))
new_SB_infor <- SB_sites_infor[select_SB,]
IP_reads <- new_SB_infor[,grep("IP",colnames(new_SB_infor))]
Input_reads <- new_SB_infor[,grep("Input",colnames(new_SB_infor))]
load("/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/twokinds_sizefactor.Rdata")
size_factor <- two_sizefactors$size_factor1
size_factor <- two_sizefactors$size_factor2
T1 <- as.numeric(as.character(size_factor))[1:60]
T0 <- as.numeric(as.character(size_factor))[61:120]
all_peak_methylevel <- t( t(IP_reads)/T1 )/ t( t(Input_reads)/T0 )
rownames(all_peak_methylevel) <- NULL
enrichFlag <- apply( all_peak_methylevel,1,function(x){sum(x>1)> 30})
new_SB_sites <- new_SB_infor[enrichFlag,]
OR <- all_peak_methylevel[enrichFlag,]
new_tablegene <- (as.data.frame(table(as.character(new_SB_sites$geneID))))

logOR <- log(OR)
##quantile normalization
logOR.id <- which( rowMeans(logOR) < quantile( rowMeans(logOR), 0.95 ) & rowMeans(logOR) > quantile( rowMeans(logOR), 0.05) )
K_IPe_ij <- apply(logOR[logOR.id,], 2, function(x){
  
  fit <- lm(y~m, data = data.frame(y = x, m=rowMeans(logOR)[logOR.id] ))
  y.est <- predict(fit, newdata =  data.frame(m = rowMeans(logOR)))
  return( y.est - rowMeans(logOR) )
})
new_label <- which(!is.na(rowMeans(K_IPe_ij)))
methy_corrected_SB <- new_SB_sites[new_label,]
LogOR_methy <- as.data.frame(log(OR[new_label,])- K_IPe_ij[new_label,])
for (i in 1:nrow(LogOR_methy)) {
  for (j in 1:ncol(LogOR_methy)) {
    
    if(LogOR_methy[i,j]<=0){
      LogOR_methy[i,j]=0
    }
  }
}

new_select_label <- apply( LogOR_methy,1,function(x){sum(x>0)> 30})
new_select_SB <- methy_corrected_SB[new_select_label,]
new_logOR <- LogOR_methy[new_select_label,]
last_select_SB <-new_select_SB
last_logOR <- new_logOR
##match expression
load("./Lymphoblastoid_cell_line/gene_expr.Rdata")
gene_readscount <- gene_expr$gene_reads
overlap_gene <- intersect(rownames(gene_readscount),last_select_SB$geneID)
match_genereads <- gene_readscount[which(!is.na(match(rownames(gene_readscount),overlap_gene))),]
select_genereads <- match_genereads[apply(match_genereads,1,function(x){sum(x>0)> 30}),]

match_label <- which(!is.na(match(last_select_SB$geneID,rownames(select_genereads))))
match_SB_infor <- last_select_SB[match_label,]
match_logOR <- last_logOR[match_label,]
SB_gene_name <- as.character(match_SB_infor$geneID)
table_gene <- as.data.frame(table(SB_gene_name))
nrow(table_gene[table_gene$Freq>1,])
match_SB_methy <- match_logOR
match_SB_sites <- match_SB_infor
match_genecount <- select_genereads
match_expr_SB_methy_infor <- list(match_expr=select_genereads,SB_sites_infor=match_SB_infor,SB_methylevel=match_logOR)
save(match_expr_SB_methy_infor,file = "./SRAMP/full_RNA_mode/exomePeak2_result/SB_exprmethyinfor_sizefactor.Rdata")


                    
