

r <- getOption("repos");
r["CRAN"] <- "http://cran.rstudio.com/"
options(repos=r)
.libPaths("/home/guanwh/cao00128/Rlibrary/")
# install.packages("bigsnpr")

library("data.table")
library("survival")
library("dplyr")
#library("stringr")
library("survminer")
library("getopt")
library(stringr)
library(genio)
library(bigsnpr)

spec = matrix(c(
  'file.index'      , 'f' , 1, "integer",   "choromosome id, range: 1-22"
), byrow=TRUE, ncol=5)
opt = getopt(spec)
chr <- opt$file.index


filelist <- read.table("/home/guanwh/shared/DeKAF/GWAS/dosage/filelist.txt", header = F)
filelist <- filelist$V1[grepl(paste0("chr", chr, ":"), filelist$V1)]
filelist <- filelist[!filelist == "Minnesota_GoNL_1KG_chr9:45-50.dose"]

geno_1e.5 = c()
geno_1e.4 = c()
geno_1e.3 = c()
geno_1e.2 = c()
geno_5e.2 = c()
bim_1e.5 = c()
bim_1e.4 = c()
bim_1e.3 = c()
bim_1e.2 = c()
bim_5e.2 = c()

gwas <- vector("list",length(filelist))
for(i in 1: length(filelist)){
  print(i)
  gwas[[i]] <- read.csv(paste0("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/", filelist[i], ".csv"),
                   row.names = 1)
  gwas[[i]] <- gwas[[i]][!is.na(gwas[[i]]$beta),]
  gwas[[i]] <- gwas[[i]][gwas[[i]]$ibs_freq>0.05,]
  dosage <- fread(paste0("/home/guanwh/shared/DeKAF/GWAS/dosage/", as.character(filelist[i])))
  dosage <- as.data.frame(dosage)
  
  dosage <- dosage[dosage$rs_id %in% gwas[[i]]$rs_id[which(gwas[[i]]$pvalue<5e-2)],]
  if(nrow(dosage)==0) next
  X_5e.2 = as.matrix(dosage[,14:ncol(dosage)])
  rownames(X_5e.2) = dosage$rs_id
  geno_5e.2 <- rbind(geno_5e.2, X_5e.2)
  bim_5e.2 <- rbind(bim_5e.2, data.frame(chr=as.character(chr),
                                         id=dosage$rs_id,
                                         posg=dosage$position/1e6,
                                         pos=dosage$position,
                                         ref=dosage$A0,
                                         alt=dosage$A1)
  )
  

  
  
  dosage <- dosage[dosage$rs_id %in% gwas[[i]]$rs_id[which(gwas[[i]]$`pvalue`<1e-2)],]
  if(nrow(dosage)==0) next
  X_1e.2 = as.matrix(dosage[,14:ncol(dosage)])
  rownames(X_1e.2) = dosage$rs_id
  geno_1e.2 <- rbind(geno_1e.2, X_1e.2)
  bim_1e.2 <- rbind(bim_1e.2, data.frame(chr=as.character(chr),
                                         id=dosage$rs_id,
                                         posg=dosage$position/1e6,
                                         pos=dosage$position,
                                         ref=dosage$A0,
                                         alt=dosage$A1)
  )
  
  dosage <- dosage[dosage$rs_id %in% gwas[[i]]$rs_id[which(gwas[[i]]$`pvalue`<1e-3)],]
  if(nrow(dosage)==0) next
  X_1e.3 = as.matrix(dosage[,14:ncol(dosage)])
  rownames(X_1e.3) = dosage$rs_id
  geno_1e.3 <- rbind(geno_1e.3, X_1e.3)
  bim_1e.3 <- rbind(bim_1e.3, data.frame(chr=as.character(chr),
                                         id=dosage$rs_id,
                                         posg=dosage$position/1e6,
                                         pos=dosage$position,
                                         ref=dosage$A0,
                                         alt=dosage$A1)
  )
  
  dosage <- dosage[dosage$rs_id %in% gwas[[i]]$rs_id[which(gwas[[i]]$`pvalue`<1e-4)],]
  if(nrow(dosage)==0) next
  X_1e.4 = as.matrix(dosage[,14:ncol(dosage)])
  rownames(X_1e.4) = dosage$rs_id
  geno_1e.4 <- rbind(geno_1e.4, X_1e.4)
  bim_1e.4 <- rbind(bim_1e.4, data.frame(chr=as.character(chr),
                                         id=dosage$rs_id,
                                         posg=dosage$position/1e6,
                                         pos=dosage$position,
                                         ref=dosage$A0,
                                         alt=dosage$A1)
  )
  
  dosage <- dosage[dosage$rs_id %in% gwas[[i]]$rs_id[which(gwas[[i]]$`pvalue`<1e-5)],]
  if(nrow(dosage)==0) next
  X_1e.5 = as.matrix(dosage[,14:ncol(dosage)])
  rownames(X_1e.5) = dosage$rs_id
  geno_1e.5 <- rbind(geno_1e.5, X_1e.5)
  bim_1e.5 <- rbind(bim_1e.5, data.frame(chr=as.character(chr),
                                         id=dosage$rs_id,
                                         posg=dosage$position/1e6,
                                         pos=dosage$position,
                                         ref=dosage$A0,
                                         alt=dosage$A1)
  )
}
gwas_table <- do.call("rbind",gwas)


# test = snp_clumping(
#   G = (matrix(geno_1e.5[order(bim_1e.5$pos),], ncol=dim(dosage)[2])),
#   infos.chr = chr,
#   S = -gwas_table$pvalue[gwas_table$pvalue<1e-5],
#   thr.r2 = 0.2,
#   size = 100/0.2
# )
# write.csv(data.frame(gwas_table[gwas_table$pvalue<1e-5,][test,]), paste0("/scratch.global/cao00128/dekaf_plink/1e.5/clumping_chr_", chr, ".csv"))




if(length(bim_1e.5)>0){
  write_plink(
    file = paste0("/scratch.global/cao00128/dekaf_plink/1e.5/chr_", chr),
    X = matrix(geno_1e.5[order(bim_1e.5$pos),], nrow=nrow(bim_1e.5)),
    bim = bim_1e.5[order(bim_1e.5$pos),],
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
  )
  obj.bed <- bed(paste0("/scratch.global/cao00128/dekaf_plink/1e.5/chr_", chr, ".bed"))
  test = bed_clumping(
    obj.bed,
    # ind.row = gwas_table$pvalue[gwas_table$pvalue<1e-4],
    S = -gwas_table$pvalue[gwas_table$pvalue<1e-5],
    thr.r2 = 0.2,
    size = 100/0.2,
    exclude = NULL,
    ncores = 1
  )
  write.csv(data.frame(gwas_table[gwas_table$pvalue<1e-5,][test,]), paste0("/scratch.global/cao00128/dekaf_plink/1e.5/clumping_chr_", chr, ".csv"))
}
if(length(bim_1e.4)>0){
  write_plink(
    file = paste0("/scratch.global/cao00128/dekaf_plink/1e.4/chr_", chr),
    X = matrix(geno_1e.4[order(bim_1e.4$pos),], nrow=nrow(bim_1e.4)),
    bim = bim_1e.4[order(bim_1e.4$pos),],
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
  )
  obj.bed <- bed(paste0("/scratch.global/cao00128/dekaf_plink/1e.4/chr_", chr, ".bed"))
  test = bed_clumping(
    obj.bed,
    # ind.row = gwas_table$pvalue[gwas_table$pvalue<1e-4],
    S = -gwas_table$pvalue[gwas_table$pvalue<1e-4],
    thr.r2 = 0.2,
    size = 100/0.2,
    exclude = NULL,
    ncores = 1
  )
  write.csv(data.frame(gwas_table[gwas_table$pvalue<1e-4,][test,]), paste0("/scratch.global/cao00128/dekaf_plink/1e.4/clumping_chr_", chr, ".csv"))
}
if(length(bim_1e.3)>0){
  write_plink(
    file = paste0("/scratch.global/cao00128/dekaf_plink/1e.3/chr_", chr),
    X = as.matrix(geno_1e.3[order(bim_1e.3$pos),]),
    bim = bim_1e.3[order(bim_1e.3$pos),],
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
  )
  obj.bed <- bed(paste0("/scratch.global/cao00128/dekaf_plink/1e.3/chr_", chr, ".bed"))
  test = bed_clumping(
    obj.bed,
    # ind.row = gwas_table$pvalue[gwas_table$pvalue<1e-4],
    S = -gwas_table$pvalue[gwas_table$pvalue<1e-3],
    thr.r2 = 0.2,
    size = 100/0.2,
    exclude = NULL,
    ncores = 1
  )
  write.csv(data.frame(gwas_table[gwas_table$pvalue<1e-3,][test,]), paste0("/scratch.global/cao00128/dekaf_plink/1e.3/clumping_chr_", chr, ".csv"))
}
if(length(bim_5e.2)>0){
  write_plink(
    file = paste0("/scratch.global/cao00128/dekaf_plink/5e.2/chr_", chr),
    X = geno_5e.2[order(bim_5e.2$pos),],
    bim = bim_5e.2[order(bim_5e.2$pos),],
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
  )
  obj.bed <- bed(paste0("/scratch.global/cao00128/dekaf_plink/5e.2/chr_", chr, ".bed"))
  test = bed_clumping(
    obj.bed,
    # ind.row = gwas_table$pvalue[gwas_table$pvalue<1e-4],
    S = -gwas_table$pvalue[gwas_table$pvalue<5e-2],
    thr.r2 = 0.2,
    size = 100/0.2,
    exclude = NULL,
    ncores = 1
  )
  write.csv(data.frame(gwas_table[gwas_table$pvalue<5e-2,][test,]), paste0("/scratch.global/cao00128/dekaf_plink/5e.2/clumping_chr_", chr, ".csv"))
}
if(length(bim_1e.2)>0){
  write_plink(
    file = paste0("/scratch.global/cao00128/dekaf_plink/1e.2/chr_", chr),
    X = geno_1e.2[order(bim_1e.2$pos),],
    bim = bim_1e.2[order(bim_1e.2$pos),],
    fam = NULL,
    pheno = NULL,
    verbose = TRUE,
    append = FALSE
  )
  obj.bed <- bed(paste0("/scratch.global/cao00128/dekaf_plink/1e.2/chr_", chr, ".bed"))
  test = bed_clumping(
    obj.bed,
    # ind.row = gwas_table$pvalue[gwas_table$pvalue<1e-4],
    S = -gwas_table$pvalue[gwas_table$pvalue<1e-2],
    thr.r2 = 0.2,
    size = 100/0.2,
    exclude = NULL,
    ncores = 1
  )
  write.csv(data.frame(gwas_table[gwas_table$pvalue<1e-2,][test,]), paste0("/scratch.global/cao00128/dekaf_plink/1e.2/clumping_chr_", chr, ".csv"))
}











### plink2 --pfile all_phase3 --max-alleles 2 --make-bed --out 1000genome
# .libPaths("/home/guanwh/cao00128/Rlibrary/")
# library("data.table")
# genome_bim <- fread("/scratch.global/cao00128/gen03_plink/1000genome.bim")
# 
# for(chr in 1:22){
#   print(chr)
#   chr_bim <- fread(paste0("/scratch.global/cao00128/gen03_plink/chr_", chr, ".bim"))
#   chr_bim$abb <- sub(":.*", "", chr_bim$V2)
#   test <- merge(chr_bim, genome_bim[,2:3], by.y="V2", by.x="abb")
#   
#   write.table(, 
#               paste0(),
#               row.names = F,
#               col.names = F,
#               sep="\t")
# }
# 
# 
# sum=0
# test.bim <- genome_bim[genome_bim$V1==1,]
# for(i in 1:nrow(test.bim)){
#   if(length(grep(test.bim$V2[i], chr_bim$V2))>0){sum=sum+1}
# }


# for(i in 1:22){
#   data.name = paste0("/scratch.global/cao00128/gen03_plink/chr_",
#                      i,
#                      ".bim")
#   data <- read.table(data.name)
#   data$V3 <- data$V4/1e6
#   write.table(data, data.name, row.names = F, col.names = F, sep = "\t")
# }


# data.name = paste0("/scratch.global/cao00128/gen03_plink/prune_chr_",
#                    1:22,
#                    ".prune.in")
# data <- lapply(data.name,function(x){
#   x=read.table(x)
#   x$V1
# })
# snp.list = unlist(data)
# write.table(snp.list, "/home/guanwh/cao00128/ibd/prs/snps_pruned_imputed.txt", row.names = F, col.names = F, sep = "\t")


