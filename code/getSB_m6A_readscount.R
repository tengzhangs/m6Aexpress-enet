##get single base reads count 
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

f1 <- "./SRAMP/full_RNA_mode/exomePeak2_result/singlebase_m6Asitesmaps.txt"
SB_sites <- read.table(f1,header = T)
# f2 <- "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/exomePeakoutA/Mod.csv"
# SB_sites <- read.csv(f2)

SB_sites_GR <- GRanges(seqnames = as.character(SB_sites$seqnames),
                        IRanges(start = as.numeric(as.character(SB_sites$start)),
                                end = as.numeric(as.character(SB_sites$end))),
                        strand = as.character(SB_sites$strand))
names(SB_sites_GR) <- as.character(SB_sites$peak_num)


# SB_sites_GR <- GRanges(seqnames = as.character(SB_sites$chr),
#                       IRanges(start = as.numeric(as.character(SB_sites$chromStart)),
#                               end = as.numeric(as.character(SB_sites$chromEnd))),
#                       strand = as.character(SB_sites$strand))
#bnames(SB_sites_GR) <- as.character(SB_sites$geneID)

bam2counts <- function(bamFile,
                       region,
                       fragLength=100,
                       mapq_value=30){
  
  ### check input
  if (class(region)!="GRanges") stop("region must be a GRanges object")
  if (class(fragLength)!="numeric") stop("fragLength must be numeric")
  
  param <- Rsamtools::ScanBamParam(
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE,
                                  isDuplicate = FALSE),
    what = "mapq")
  aln <- GenomicAlignments::readGAlignments(bamFile,
                                            param = param)
  aln <- aln[GenomicRanges::mcols(aln)$mapq > mapq_value]
  aln <- GenomicRanges::granges(aln)
  # get the chromosome length
  length.chr <- GenomeInfoDb::seqlengths(seqinfo(aln))
  # to make sure the "region + fragLength" is not out-of-bound
  out.bound <- ifelse(
    GenomicRanges::strand(aln)=="+",
    GenomicRanges::start(aln)+ fragLength-1,
    GenomicRanges::end(aln))  >
    length.chr[as.character(GenomicRanges::seqnames(aln))]
  aln <- aln[GenomicRanges::end(aln) > fragLength&!out.bound]
  aln <- GenomicRanges::resize(aln,
                               fragLength,
                               fix = "start")
  total_reads <- length(aln)
  # countOverlaps is strand aware, so remove strand
  GenomicRanges::strand(aln) <- "*"
  counts <- GenomicRanges::countOverlaps(region, aln,type = "within")
  names(counts) <- names(region)
  single_basecount_infor <- list(counts=counts,total_reads=total_reads)
  return(single_basecount_infor)
}

region <- SB_sites_GR
fragLength=100
mapq_value=30
library_size <- vector()
IP_SB_readscounts <- data.frame()
for (i in 1:length(ip_bams)) {
  bamfile <- ip_bams[i]
  one_bamcounts_infor <- bam2counts(bamFile=bamfile,
                                    region=SB_sites_GR,
                                    fragLength=100,
                                    mapq_value=30)
  one_counsts <- one_bamcounts_infor$counts
  library_size[i] <- one_bamcounts_infor$total_reads
  IP_SB_readscounts <- rbind(IP_SB_readscounts,one_counsts)
  
}
IP_SB_readscounts <- t(IP_SB_readscounts)
# rownames(IP_SB_readscounts) <- as.character(SB_sites$peak_num)
colnames(IP_SB_readscounts) <- paste0(samplenames,"_IP")
##Input
library_size_Input <- vector()
Input_SB_readscounts <- data.frame()
for (i in 1:length(input_bams)) {
  bamfile <- input_bams[i]
  one_bamcounts_infor <- bam2counts(bamFile=bamfile,
                                    region=SB_sites_GR,
                                    fragLength=100,
                                    mapq_value=30)
  one_counsts <- one_bamcounts_infor$counts
  library_size_Input[i] <- one_bamcounts_infor$total_reads
  Input_SB_readscounts <- rbind(Input_SB_readscounts,one_counsts)
  
}
Input_SB_readscounts <- t(Input_SB_readscounts)
# rownames(Input_SB_readscounts) <- as.character(SB_sites$peak_num)
colnames(Input_SB_readscounts) <- paste0(samplenames,"_Input")
# IP_SB_reads_infor <- data.frame(peak_label= as.character(SB_sites$peak_num),
#                                 IP_SB_readscounts)
IP_SB_reads_infor <- data.frame(geneID= as.character(SB_sites$geneID),
                                IP_SB_readscounts)
rownames(IP_SB_reads_infor) <- NULL
# Input_SB_reads_infor <- data.frame(peak_label= as.character(SB_sites$peak_num),
#                                    Input_SB_readscounts)
Input_SB_reads_infor <- data.frame(geneID= as.character(SB_sites$geneID),
                                   Input_SB_readscounts)
rownames(Input_SB_reads_infor) <- NULL
####
IP_Input_SB_readscount <- data.frame(IP_SB_reads_infor,Input_SB_reads_infor[,-1])
IP_Input_SB_readsinfor <- list(IP_Input_SB_readscount=IP_Input_SB_readscount,
                               IP_totalreads=library_size,
                               Input_totalreads=library_size_Input)
save(IP_Input_SB_readsinfor,file = "./SRAMP/full_RNA_mode/exomePeak2_result/IP_Input_SBreadsinfor.Rdata")
