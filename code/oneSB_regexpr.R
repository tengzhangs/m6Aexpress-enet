##only one peak gene GLM
library(MASS)
load("./SRAMP/full_RNA_mode/exomePeak2_result/SB_exprmethyinfor_sizefactor.Rdata")
load("./Lymphoblastoid_cell_line/gene_expr.Rdata")
match_SB_methy <- match_expr_SB_methy_infor$SB_methylevel
match_SB_sites <- match_expr_SB_methy_infor$SB_sites_infor
match_genecount <- match_expr_SB_methy_infor$match_expr
match_gene_name <- intersect(unique(as.character(match_SB_sites$geneID)),rownames(match_genecount))  
match_SB_methyinfor <- cbind(match_SB_sites[,1:6],match_SB_methy)
match_gene_expr <- data.frame()
for (i in 1:length(match_gene_name)) {
  one_match_expr <- match_genecount[rownames(match_genecount)==match_gene_name[i],]
  match_gene_expr <- rbind(match_gene_expr,one_match_expr)
}
sample_name  <- c("NA18486","NA18498","NA18499","NA18501","NA18502","NA18504","NA18505","NA18507","NA18508",
                  "NA18510","NA18511","NA18516","NA18517","NA18519","NA18522","NA18523","NA18852","NA18855",
                  "NA18856","NA18858","NA18861","NA18862","NA18870","NA18907","NA18909","NA18912","NA18913",
                  "NA18916","NA19092","NA19093","NA19098","NA19099","NA19101","NA19102","NA19108","NA19114",
                  "NA19116","NA19119","NA19127","NA19128","NA19130","NA19131","NA19137","NA19138","NA19140",
                  "NA19143","NA19144","NA19147","NA19152","NA19153","NA19159","NA19160","NA19192","NA19193",
                  "NA19200","NA19204","NA19209","NA19222","NA19239","NA19257")
colnames(match_gene_expr) <- sample_name
rownames(match_gene_expr) <- match_gene_name
colnames(match_SB_methyinfor)[7:ncol(match_SB_methyinfor)] <- sample_name
##
table_gene <- as.data.frame(table(as.character(match_SB_methyinfor$geneID)))
one_SB_gene <- intersect(as.character(table_gene[table_gene$Freq==1,]$Var1),
                           rownames(match_gene_expr))


one_SB_methy <- data.frame()
one_SB_expr <- data.frame()
for (i in 1:length(one_SB_gene)) {
  gene_SB_methy <- match_SB_methyinfor[match_SB_methyinfor$geneID==one_SB_gene[i],]
  gene_SB_expr <-  match_gene_expr[rownames(match_gene_expr)==one_SB_gene[i],]
  one_SB_methy <- rbind(one_SB_methy,gene_SB_methy)
  one_SB_expr <- rbind(one_SB_expr,gene_SB_expr)
}
##glm NB model 
size_factor <- gene_expr[[2]]
one_SBgene_glm <- data.frame()
one_SB_theta <- vector()
for(i in 1:length(one_SB_gene)){
  one_data <- data.frame(y=t(one_SB_expr[rownames(one_SB_expr)==one_SB_gene[i],]),
                         x=t(one_SB_methy[one_SB_methy$geneID==one_SB_gene[i],7:ncol(one_SB_methy)]))
  colnames(one_data) <- c("expr","methy")
  one_gene_glm <- glm.nb(expr~methy+offset((size_factor)),data = one_data)
  one_glm_summay <- summary(one_gene_glm)
  one_SB_theta[i] <- one_gene_glm$theta
  one_betas_infor <- one_glm_summay[["coefficients"]]
  one_betas_value <- one_betas_infor[,1]
  names(one_betas_value) <- c("beta0","beta1")
  one_betas_SE <- one_betas_infor[,2]
  names(one_betas_SE) <- c("beta0_SE","beta1_SE")
  beta_pvalue <-  one_betas_infor[,4][1:2]
  names(beta_pvalue) <- c("beta0_pvalue",paste0("SB1","_pvalue"))
  
  one_gene_SB_infor <- c("beta0","SB1_beta1")
  gene_name <- rep(one_SB_gene[i],length(one_gene_SB_infor))
  one_gene_result <- data.frame(gene_name=gene_name,
                                SB_infor=one_gene_SB_infor,
                                beta_value=(one_betas_value),
                                beta_pvalue=(beta_pvalue),
                                beta_SE=(one_betas_SE))
  rownames(one_gene_result) <- NULL
  
  one_SBgene_glm <- rbind(one_SBgene_glm,one_gene_result)
}
##sig reg gene
names(one_SB_theta) <- one_SB_gene
oneSB_gene <- one_SBgene_glm[one_SBgene_glm$SB_infor!="beta0",]
onegene_sig <- oneSB_gene[oneSB_gene$beta_pvalue<0.05,]
one_SB_genename <- as.character(onegene_sig$gene_name)
oneSB_sig_regene <- one_SBgene_glm[one_SBgene_glm$gene_name%in%one_SB_genename,]
oneSB_siggene_theta <- one_SB_theta[names(one_SB_theta)%in%one_SB_genename]
oneSB_sig_beta0 <- oneSB_sig_regene[oneSB_sig_regene$SB_infor=="beta0",]
oneSB_sig_beta1 <- oneSB_sig_regene[oneSB_sig_regene$SB_infor=="SB1_beta1",]


save(one_two_SB_reggene,file = "/home/disk3/zhangteng/Lymphoblastoid_cell_line/part_result/single_base_result/SRAMP/full_RNA_mode/exomePeak2_result/one_two_sigSBreg_genesizefactorA.Rdata")

