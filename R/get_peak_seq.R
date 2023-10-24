library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

###get full RNA mode peak sequence
f1 <- "./exomePeak2_output/Mod.csv"
peak_sites_infor <- read.csv(f1)
peaksites_infor <- peak_sites_infor[,c(1:4,6,13)]
peakGR <- GRanges(seqnames = as.character(peaksites_infor$chr),
                  IRanges(start = as.numeric(as.character(peaksites_infor$chromStart)),
                          end = as.numeric(as.character(peaksites_infor$chromEnd))),
                  strand = as.character(peaksites_infor$strand))
mcols(peakGR)$peak_num <- as.character(peaksites_infor$name)
mcols(peakGR)$gene_name <- as.character(peaksites_infor$geneID)
genome <- BSgenome.Hsapiens.UCSC.hg19
peak_seq <- getSeq(genome, peakGR)
# peak_seq_con <- lapply(peak_seq, unlist)
# peak_seq_con <- DNAStringSet(peak_seq_con)
motif <- c("GGACA", "GGACC", "GGACT", "AGACA", "AGACC", "AGACT", 
           "GAACA", "GAACC", "GAACT", "AAACA", "AAACC", "AAACT",
           "TGACA", "TGACC", "TGACT", "TAACA", "TAACC", "TAACT")
cag_loc <- vmatchPattern(motif[1], peak_seq)
motif_start <- start(cag_loc)

for(i in 2:length(motif)) {
  cag_loc0 <- vmatchPattern(motif[i], peak_seq)
  motif_start0 <- start(cag_loc0)
  motif_start <- mapply(c, motif_start, motif_start0, SIMPLIFY=FALSE)
}

fl <- function(x){
  if(length(x)==0){
    return(0)
  }else{
    return(1)
  }
  
}
motif_select_peak <- which(sapply(motif_start,fl)>0)
motif_peak_line <- as.character(peakGR$peak_num)[motif_select_peak]
##motif peak sequencing select
motif_peak_seq <- peak_seq[motif_select_peak]
motif_m6A_peak <- peakGR[motif_select_peak]
names(motif_peak_seq) <- motif_peak_line

writeXStringSet(motif_peak_seq,"./SRAMP/full_RNA_mode/exomePeak2_result/new_motif_peak_seq.fa")
