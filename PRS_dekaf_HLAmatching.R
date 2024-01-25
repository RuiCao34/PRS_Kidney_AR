.libPaths("/home/guanwh/cao00128/Rlibrary/")

library("data.table")
library("survival")
library("dplyr")
#library("stringr")
library("survminer")
library("getopt")
library(stringr)

spec = matrix(c(
  'file.index'      , 'f' , 1, "integer",   "choromosome id, range: 1-22"
), byrow=TRUE, ncol=5)
opt = getopt(spec)
index <- opt$file.index


gwas.file.list <- list.files(path = "/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/",
                             pattern = "csv",
                             full.names = F)

filelist <- paste0("/home/guanwh/shared/DeKAF/GWAS/dosage/", str_sub(gwas.file.list,1,-5))


gwas.file.list <- paste0("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/",
                         gwas.file.list)




begin=1
end=length(gwas.file.list)
interval=ceiling((end-begin)/100)

clump1e.5 <- list.files("/scratch.global/cao00128/dekaf_plink/1e.5/", "clump",full.names = T)
clump1e.5 <- lapply(clump1e.5, function(x){
  file <- read.csv(x, row.names = 1)
  file$chr <- str_match(x, "_chr_(.*?).csv")[2]
  return(file)
})
clump1e.5 <- do.call("rbind", clump1e.5)
clump1e.4 <- list.files("/scratch.global/cao00128/dekaf_plink/1e.4/", "clump", full.names = T)
clump1e.4 <- lapply(clump1e.4, function(x){
  file <- read.csv(x, row.names = 1)
  file$chr <- str_match(x, "_chr_(.*?).csv")[2]
  return(file)
})
clump1e.4 <- do.call("rbind", clump1e.4)
clump1e.3 <- list.files("/scratch.global/cao00128/dekaf_plink/1e.3/", "clump", full.names = T)
clump1e.3 <- lapply(clump1e.3, function(x){
  file <- read.csv(x, row.names = 1)
  file$chr <- str_match(x, "_chr_(.*?).csv")[2]
  return(file)
})
clump1e.3 <- do.call("rbind", clump1e.3)
clump1e.2 <- list.files("/scratch.global/cao00128/dekaf_plink/1e.2/", "clump", full.names = T)
clump1e.2 <- lapply(clump1e.2, function(x){
  file <- read.csv(x, row.names = 1)
  file$chr <- str_match(x, "_chr_(.*?).csv")[2]
  return(file)
})
clump1e.2 <- do.call("rbind", clump1e.2)
clump5e.2 <- list.files("/scratch.global/cao00128/dekaf_plink/5e.2/", "clump", full.names = T)
clump5e.2 <- lapply(clump5e.2, function(x){
  file <- read.csv(x, row.names = 1)
  file$chr <- str_match(x, "_chr_(.*?).csv")[2]
  return(file)
})
clump5e.2 <- do.call("rbind", clump5e.2)

for(jobid in (begin+(index-1)*interval):min(end,begin+index*interval)){
  gwas.dekaf <- read.csv(gwas.file.list[jobid], row.names = 1)
  gwas.dekaf <- gwas.dekaf[!is.na(gwas.dekaf$beta),]
  
  # gwas.dekaf <- gwas.dekaf[sub(":.*", "", gwas.dekaf$rs_id) %in% prune$V1 | gwas.dekaf$rs_id %in% prune$V1,]
  
  dosage.prs <- fread(filelist[jobid])
  dosage.prs <- as.data.frame(dosage.prs)
  pheno<-read.csv("/home/guanwh/shared/DeKAF/GWAS/phenotype/AR_data_EA_20151208.csv")
  recip_donor_crossmatch_20141119 <- read.csv("/home/guanwh/cao00128/prs_22summer/phenos/recip_donor_crossmatch_20141119.csv")
  pheno <- merge(pheno, recip_donor_crossmatch_20141119, by="pid")
  
  dosage.prs<-dosage.prs[dosage.prs$exp_freq_a1<0.95 & dosage.prs$exp_freq_a1>0.05 & dosage.prs$info>=0.8,]
  sample.geno <- names(dosage.prs)[-c(1:13)]
  sample.geno <- substr(sample.geno, 1, 7)
  
  r_id<-pheno$rgsn[pheno$rgsn %in% sample.geno & pheno$dgsn %in% sample.geno]
  d_id <- pheno$dgsn[pheno$rgsn %in% sample.geno & pheno$dgsn %in% sample.geno]
  snps <- dosage.prs$rs_id
  snps <- snps[snps %in% gwas.dekaf$rs_id]
  geno.info <- dosage.prs[,1:13]
  geno.use <- t(dosage.prs[,-c(1:13)])
  geno.use <- as.data.frame(geno.use)
  rownames(geno.use) <- sample.geno
  colnames(geno.use) <- snps
  geno_r <- geno.use[as.character(r_id),]
  geno_d <- geno.use[as.character(d_id),]
  colnames(geno_d) <- paste0(snps, ".donor")
  rownames(geno_d) <- as.character(r_id)
  data <- merge(pheno, geno_r, by.x="rgsn", by.y="row.names")
  data <- merge(data, geno_d, by.x="rgsn", by.y="row.names")
  
  PRS.score.HLA <- PRS.score <- matrix(0, nrow=nrow(data), ncol=6)
  PRS.score.HLA <- PRS.score <- as.data.frame(PRS.score)
  colnames(PRS.score.HLA) <- colnames(PRS.score) <- c("pid", "IBS1e-5", "IBS1e-4", "IBS1e-3", "IBS1e-2", "IBS5e-2")
  
  
  for(j in 1:length(snps)){
    print(j)
    snpname<-snps[j]
    snpdonor<-paste0(snpname,".donor")
    data$collision2 = abs(data[,snpname]-data[,snpdonor])
    data$collision2[is.na(data$collision2)] = 0
    
    
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-5 & snpname %in% clump1e.5$rs_id){
      PRS.score[,2] = PRS.score[,2] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-4 & snpname %in% clump1e.4$rs_id){
      PRS.score[,3] = PRS.score[,3] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-3 & snpname %in% clump1e.3$rs_id){
      PRS.score[,4] = PRS.score[,4] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-2 & snpname %in% clump1e.2$rs_id){
      PRS.score[,5] = PRS.score[,5] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<5e-2 & snpname %in% clump5e.2$rs_id){
      PRS.score[,6] = PRS.score[,6] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    
    # if(sum(is.na(PRS.score))>0) break
    
    
    if(grepl("chr6:25-30",gwas.file.list[jobid]) | grepl("chr6:30-35",gwas.file.list[jobid])){
      
    }
    
    
  }
  
  PRS.score$pid <- data$rgsn
  write.csv(PRS.score, file = paste0("/home/guanwh/cao00128/ibd/prs/PRS_dekaf_HLAmatching/", as.character(jobid), ".csv"))
  
  
  
  
}