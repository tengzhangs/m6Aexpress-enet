BED12toGRangesList <- function(filepath, header=FALSE) {
  
  # message
  print("Converting BED12 to GRangesList")
  print("It may take a few minutes")
  
  # read bed file
  a = read.table(filepath,sep="\t",header=header,stringsAsFactors =FALSE)
  # mcols_info =a[,13:length(a[1,])]
  a = a[,1:12]
  
  # get transcripts
  no_tx = length(a[,1])
  tx_id = 1:no_tx;
  tx_name = paste("peak_",1:no_tx,sep="")
  tx_chrom = as.character(a[,1])
  tx_strand = a[,6]
  tx_start = a[,2]+1
  tx_end = a[,3]
  transcripts= data.frame(tx_id,tx_name,tx_chrom,tx_strand,tx_start,tx_end)
  head(transcripts)
  
  
  
  # get genes
  tx_name = tx_name
  gene_id = as.character(a[,4])
  gene_id[is.na(gene_id)]="NA"
  gene=data.frame(tx_name,gene_id)
  
  # 
  splicing <- lapply(1:no_tx, .spliceSingleTrans, a=a, tx_start=tx_start)
  splicing <- .combineListOfSplicing(splicing)
  
  # make txdb
  peaks = suppressWarnings(
    makeTxDb(transcripts=transcripts, 
             splicings=splicing,
             genes=gene))
  
  # generate GRangesList
  tx <- exonsBy(peaks, "tx",use.names=TRUE)
  mcols(tx) <- a
  
  return(tx)
}

.combineListOfSplicing <- function(t){
  
  a <- paste("t[[",1:length(t),"]]", sep="")
  a <- paste(a,collapse =",")
  a <- paste("rbind(",a,")",sep="")
  c <- parse(text=a)
  b <- suppressWarnings(eval(c))
  
  return(b)
}

.spliceSingleTrans <- function(i,a,tx_start) {
  tx = a[i,]
  tx_id = i
  exon_rank=1:as.integer(tx[10])
  
  # get start
  temp = as.integer(strsplit(as.character(tx[12]), ",")[[1]]) + tx_start[i]
  exon_start=temp
  
  # get end
  temp = as.integer(strsplit(as.character(tx[11]), ",")[[1]])
  temp2 = temp + exon_start - 1
  exon_end=temp2
  
  # get CDS
  cds_start = exon_start
  cds_end = exon_end
  
  # get data frame
  splicing_tx = data.frame(tx_id,exon_rank,exon_start,exon_end,cds_start,cds_end)
  return(splicing_tx)
}


f1 <- "./SRAMP/full_RNA_mode/exomePeak2_result/singlebase_m6Asites.txt"
sramp_result <- read.delim(f1,header = T)
##
predict_result <- as.character(sramp_result$Classification)
high_site_label1 <- grep("high",predict_result)
high_site_label2 <- grep("High",predict_result)
high_site_label <- c(high_site_label1,high_site_label2)
high_m6Asites <- sramp_result[high_site_label,]
##conver B12 to GRangesList
fs <- "./exomePeak2/exomePeak2_output/Mod.bed"
allpeak_GRlist <- BED12toGRangesList(filepath = fs,header = F)

fa <- "./exomePeak2/exomePeak2_output/Mod.csv"
all_peak <-  read.csv(fa)
allpeakGR <- GRanges(seqnames = as.character(all_peak$chr),
                     IRanges(start = as.numeric(as.character(all_peak$chromStart)),
                             end = as.numeric(as.character(all_peak$chromEnd))),
                     strand = as.character(all_peak$strand))

mcols(allpeakGR)$gene_name <- as.character(all_peak$geneID)
names(allpeakGR) <- as.character(all_peak$name)

high_m6Asites_peakname <- as.character(high_m6Asites$Seq_ID)
# multi_peak_name <- as.character(high_m6Asites$Seq_ID)[grep("_",as.character(high_m6Asites$Seq_ID))]
select_peakGR <- allpeakGR[which(!is.na(match(names(allpeakGR),unique(high_m6Asites_peakname))))]
select_peakGRList <- allpeak_GRlist[which(!is.na(match(names(allpeak_GRlist),unique(high_m6Asites_peakname))))]
peak_label <- as.character(unique(high_m6Asites_peakname))
pre_m6Asite <- data.frame()
for (i in 1:length(peak_label)) {
  onepeak_GR <- (select_peakGR[which(!is.na(match(names(select_peakGR),peak_label[i])))])
  onepeak_GRList <- unlist(select_peakGRList[which(!is.na(match(names(select_peakGRList),peak_label[i])))])
  select_multi_peak_label <- high_m6Asites_peakname[which(!is.na(match(high_m6Asites_peakname,peak_label[i])))]
  onepeak_m6Asites <- high_m6Asites[which(!is.na(match(high_m6Asites$Seq_ID,select_multi_peak_label))),]
  onepeak_m6Asites <- onepeak_m6Asites[order(onepeak_m6Asites$Position,decreasing = F),]
  one_site_GR <- GRanges(select_multi_peak_label,
                         IRanges(start =as.numeric(as.character(onepeak_m6Asites$Position)),
                                 width = rep(1,nrow(onepeak_m6Asites))),
                         strand =rep(unique(as.character(onepeak_GR@strand)),nrow(onepeak_m6Asites)))
  
  m6Asites_GR <- mapFromTranscripts(one_site_GR,onepeak_GR)
  overlap_ind <- findOverlaps(m6Asites_GR,onepeak_GRList,type = "within")
  if(length(unique(overlap_ind@from))==0){
    next
  }
  if(length(unique(overlap_ind@from))>0){
    overlap_m6AGR <- m6Asites_GR[unique(overlap_ind@from)]
    if(length(overlap_m6AGR)==length(one_site_GR)){
      one_m6Asite <- data.frame(as.data.frame(overlap_m6AGR),peak_num=peak_label[i],onepeak_m6Asites[,-c(1:2)])
      pre_m6Asite <- rbind(pre_m6Asite,one_m6Asite)
    }
    if((length(overlap_m6AGR)<length(one_site_GR))){
      one_m6Asite <- data.frame(as.data.frame(overlap_m6AGR),peak_num=peak_label[i],onepeak_m6Asites[unique(overlap_ind@from),-c(1:2)])
      pre_m6Asite <- rbind(pre_m6Asite,one_m6Asite)
      
    }
    
  }
  
  
}



write.table(pre_m6Asite,file = "./exomePeak2_result/singlebase_m6Asitesmaps.txt",row.names = FALSE)
