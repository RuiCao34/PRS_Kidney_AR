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

gen03.file.list <- paste0("/home/guanwh/shared/DeKAF/GEN03/imputation/dosage/UMN",
                          str_sub(gwas.file.list,10,-5))


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
  # gwas.dekaf <- gwas.dekaf[gwas.dekaf$ibs_freq>0.05,]
  # gwas.dekaf <- gwas.dekaf[sub(":.*", "", gwas.dekaf$rs_id) %in% prune$V1 | gwas.dekaf$rs_id %in% prune$V1,]
  
  dosage.gen03 <- fread(gen03.file.list[jobid])
  dosage.gen03 <- as.data.frame(dosage.gen03)
  pheno<-read.csv("/home/guanwh/shared/DeKAF/GEN03/pheno/AR_EA_PC_20171012.csv")
  dosage.gen03<-dosage.gen03[dosage.gen03$exp_freq_a1<0.95 & dosage.gen03$exp_freq_a1>0.05 & dosage.gen03$info>=0.8,]
  sample.geno <- names(dosage.gen03)[-c(1:13)]
  r_id<-intersect(as.numeric(sample.geno), as.numeric(sample.geno)-10000)
  r_id<-r_id[!is.na(r_id)]
  r_id<-intersect(r_id, pheno$IID)
  snps <- dosage.gen03$rs_id
  snps <- snps[snps %in% gwas.dekaf$rs_id]
  geno.info <- dosage.gen03[,1:15]
  geno.use <- t(dosage.gen03[,-c(1:15)])
  geno.use <- as.data.frame(geno.use)
  colnames(geno.use) <- snps
  geno_r <- geno.use[as.character(r_id),]
  geno_d <- geno.use[as.character(r_id+10000),]
  colnames(geno_d) <- paste0(snps, ".donor")
  rownames(geno_d) <- as.character(r_id)
  data <- merge(pheno, geno_r, by.x="IID", by.y="row.names")
  data <- merge(data, geno_d, by.x="IID", by.y="row.names")
  
  PRS.score <- matrix(0, nrow=nrow(data), ncol=6)
  PRS.score <- as.data.frame(PRS.score)
  colnames(PRS.score) <- c("IID", "IBS1e-5", "IBS1e-4", "IBS1e-3", "IBS1e-2", "IBS5e-2")
  
  
  for(j in 1:length(snps)){
    # print(j)
    snpname<-snps[j]
    snpdonor<-paste0(snpname,".donor")
    data$collision2 = abs(data[,snpname]-data[,snpdonor])
    data$collision2[is.na(data$collision2)] = 0
    
    
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-5 & snps[j] %in% clump1e.5$rs_id){
      PRS.score[,2] = PRS.score[,2] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-4 & snps[j] %in% clump1e.4$rs_id){
      PRS.score[,3] = PRS.score[,3] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-3 & snps[j] %in% clump1e.3$rs_id){
      PRS.score[,4] = PRS.score[,4] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<1e-2 & snps[j] %in% clump1e.2$rs_id){
      PRS.score[,5] = PRS.score[,5] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    if(gwas.dekaf$pvalue[which(gwas.dekaf$rs_id==snps[j])]<5e-2 & snps[j] %in% clump5e.2$rs_id){
      PRS.score[,6] = PRS.score[,6] + data$collision2 * gwas.dekaf$beta[which(gwas.dekaf$rs_id==snps[j])]
    }
    
    # if(sum(is.na(PRS.score))>0) break
    
    
    if(grepl("chr6:25-30",gwas.file.list[jobid]) | grepl("chr6:30-35",gwas.file.list[jobid])){
      
    }
    
    
  }
  
  PRS.score$IID <- data$IID
  write.csv(PRS.score, file = paste0("/home/guanwh/cao00128/ibd/prs/PRS2GEN03_clump_HLAmatching/", as.character(jobid), ".csv"))
  
  
  
  
}