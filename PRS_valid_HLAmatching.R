# srun -N 1 --ntasks-per-node=1  --mem-per-cpu=32gb -t 10:00:00 -p interactive --pty bash

## whole-genome

.libPaths("/home/guanwh/cao00128/Rlibrary/")
library(data.table)
library("survival")
library(stringr)

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



pheno<-read.csv("/home/guanwh/shared/DeKAF/GEN03/pheno/AR_EA_PC_20171012.csv")
relateness<-read.csv("/home/guanwh/cao00128/ibd/prs/gen03_relateness.csv",row.names = 1)
colnames(relateness)[4]="kinship"
relateness$V1 = as.character(relateness$V1)
# relateness$V1 = paste0(str_sub(relateness$V1,1,1),str_sub(relateness$V1,3,-1))
pheno <- merge(pheno,relateness[,c(1,4)],by.x="IID",by.y="V1")


prsfile <- list.files("/home/guanwh/cao00128/ibd/prs/PRS2GEN03_clump_HLAmatching/", "csv", full.names = T)
PRS_wg <- lapply(prsfile, read.csv, row.names = 1)
prs_sum <- Reduce("+", lapply(PRS_wg, function(x){x[,2:6]}))
PRS_wg <- cbind(PRS_wg[[1]][,1],prs_sum)
colnames(PRS_wg)[1] <- "IID"
gwas.files <- list.files("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/", "csv", full.names = T)
gwas_results <- lapply(gwas.files, function(x){
  file <- read.csv(x, row.names = 1)
  file$chr <- str_match(x, "_chr(.*?):")[2]
  return(file)
})
gwas_results <- do.call("rbind", gwas_results)
gwas_results <- gwas_results[gwas_results$ibs_freq>0.05,]
# prs_snps <- gwas_results[gwas_results$chr==6 & gwas_results$position >= 25759242 & gwas_results$position <= 33534827,]
prs_snps <- gwas_results[sub(":.*", "", gwas_results$rs_id) %in% clump1e.2$rs_id | gwas_results$rs_id %in% clump1e.2$rs_id,]
write.csv(prs_snps,"/home/guanwh/cao00128/ibd/prs/prssnps_HLAmatching_1e2.csv")
prs_snps <- gwas_results[sub(":.*", "", gwas_results$rs_id) %in% clump1e.3$rs_id | gwas_results$rs_id %in% clump1e.3$rs_id,]
write.csv(prs_snps,"/home/guanwh/cao00128/ibd/prs/prssnps_HLAmatching_1e3.csv")
weight <- rep(0, 5)
snp_number <- rep(0, 5)
for(i in 1:length(gwas.files)){
  print(i)
  file <- read.csv(gwas.files[i], row.names = 1)
  file <- file[file$ibs_freq>0.05,]
  # file <- file[sub(":.*", "", file$rs_id) %in% prune$V1 | file$rs_id %in% prune$V1,]
  weight[1] <- weight[1] + sum(file$beta[which(file$pvalue<1e-5 & file$rs_id %in% clump1e.5$rs_id)], na.rm = T)
  weight[2] <- weight[2] + sum(file$beta[which(file$pvalue<1e-4 & file$rs_id %in% clump1e.4$rs_id)], na.rm = T)
  weight[3] <- weight[3] + sum(file$beta[which(file$pvalue<1e-3 & file$rs_id %in% clump1e.3$rs_id)], na.rm = T)
  weight[4] <- weight[4] + sum(file$beta[which(file$pvalue<1e-2 & file$rs_id %in% clump1e.2$rs_id)], na.rm = T)
  weight[5] <- weight[5] + sum(file$beta[which(file$pvalue<5e-2 & file$rs_id %in% clump5e.2$rs_id)], na.rm = T)
  snp_number[1] <- snp_number[1] + sum(file$pvalue<1e-5 & file$rs_id %in% clump1e.5$rs_id, na.rm = T)
  snp_number[2] <- snp_number[2] + sum(file$pvalue<1e-4 & file$rs_id %in% clump1e.4$rs_id, na.rm = T)
  snp_number[3] <- snp_number[3] + sum(file$pvalue<1e-3 & file$rs_id %in% clump1e.3$rs_id, na.rm = T)
  snp_number[4] <- snp_number[4] + sum(file$pvalue<1e-2 & file$rs_id %in% clump1e.2$rs_id, na.rm = T)
  snp_number[5] <- snp_number[5] + sum(file$pvalue<5e-2 & file$rs_id %in% clump5e.2$rs_id, na.rm = T)
}

for(i in 1:5){
  PRS_wg[,i+1] = PRS_wg[,i+1] / weight[i]
}
write.csv(PRS_wg,"/home/guanwh/cao00128/ibd/prs/prs_gen03_scores_HLAmatching_1e3.csv")

for(i in 1:5){
  PRS_wg[,i+1] = scale(PRS_wg[,i+1])
}
HLA_Data <- read.csv("/home/guanwh/cao00128/ibd/prs/dekaf_gen03_demo_20210325.csv")


data.wg <- merge(pheno, PRS_wg, by="IID")
data.wg <- merge(data.wg, HLA_Data[,c(2,47)], by.x="IID", by.y="rgsn")
data.wg$nummismatch <- as.numeric(data.wg$nummismatch)

table2 = matrix(NA, nrow=5, ncol=4)
table2 <- as.data.frame(table2)
colnames(table2) <- c("p cutoff", "2.5% exp coef", "97.5% exp coef", "p value")
table2[,1] = c("1e-5","1e-4","1e-3","1e-2","5e-2")
# covars <- "age_at_tx+genderc+PRA_posneg+nonkidneytx+nummismatch+kinship"
covars <- "age_at_tx+genderc+PRA_posneg+nonkidneytx+nummismatch"
depvar <- "t2ar_ever,ar_ever"
IBS <- c("IBS1e.5","IBS1e.4","IBS1e.3","IBS1e.2","IBS5e.2")

for(i in 1:5){
  res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.wg)  
  res.cox1 <- summary(res.cox1)
  table2[i,2:3] <- res.cox1$conf.int[1,3:4]
  table2[i,4] <- res.cox1$coefficients[1,5]
}
table2 <- cbind(table2, snp_number)
colnames(table2)[5] <- "PRS SNP numbers"
write.csv(table2, "/home/guanwh/cao00128/ibd/prs/dekaf_wg_mismatch_HLAmatching.csv")

i = 3
relateness_gen03<-read.csv("/home/guanwh/cao00128/ibd/prs/gen03_relateness.csv",row.names = 1)
data.wg = merge(data.wg, relateness_gen03, by.x = 'IID', by.y = 'V1')
res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.wg, subset = p > 0.1)  
res.cox1 <- summary(res.cox1)
# table2[i,2:3] <- res.cox1$conf.int[1,3:4]
# table2[i,4] <- res.cox1$coefficients[1,5]

pheno_raw<-read.csv("/home/guanwh/shared/DeKAF/GEN03/pheno/AR_EA_PC_20171012.csv")
donorgender <-pheno_raw$genderc[match(as.character(data.wg$IID+10000),as.character(pheno_raw$IID))] 
sum(donorgender == data.wg$genderc)

pheno_raw2<-read.csv("/home/guanwh/shared/DeKAF/GEN03/pheno/GEN03_eGFR_1year_EA_recip_data_20180124.csv")
donorgender2 <-pheno_raw2$gend_donor[match(as.character(data.wg$IID),as.character(pheno_raw2$IID))] 
table(pheno_raw2$gend_donor[match(as.character(data.wg$IID),as.character(pheno_raw2$IID))],
      pheno_raw2$gender[match(as.character(data.wg$IID),as.character(pheno_raw2$IID))] )

## HLA
PRS_HLA <- read.csv("/home/guanwh/cao00128/ibd/prs/PRS_gen03_HLA_HLAmatching.csv", row.names = 1)
gwas.files <- c("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/Minnesota_GoNL_1KG_chr6:25-30.dose.csv",
                "/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/Minnesota_GoNL_1KG_chr6:30-35.dose.csv")
weight_hla <- rep(0, 5)
snp_number_hla <- rep(0, 5)
for(i in 1:length(gwas.files)){
  print(i)
  file <- read.csv(gwas.files[i], row.names = 1)
  # file <- file[sub(":.*", "", file$rs_id) %in% prune$V1 | file$rs_id %in% prune$V1,]
  file <- file[file$position >= 25759242 & file$position <= 33534827, ]
  file <- file[file$ibs_freq>0.05,]
  weight_hla[1] <- weight_hla[1] + sum(file$beta[which(file$pvalue<1e-5 & file$rs_id %in% clump1e.5$rs_id)], na.rm = T)
  weight_hla[2] <- weight_hla[2] + sum(file$beta[which(file$pvalue<1e-4 & file$rs_id %in% clump1e.4$rs_id)], na.rm = T)
  weight_hla[3] <- weight_hla[3] + sum(file$beta[which(file$pvalue<1e-3 & file$rs_id %in% clump1e.3$rs_id)], na.rm = T)
  weight_hla[4] <- weight_hla[4] + sum(file$beta[which(file$pvalue<1e-2 & file$rs_id %in% clump1e.2$rs_id)], na.rm = T)
  weight_hla[5] <- weight_hla[5] + sum(file$beta[which(file$pvalue<5e-2 & file$rs_id %in% clump5e.2$rs_id)], na.rm = T)
  snp_number_hla[1] <- snp_number_hla[1] + sum(file$pvalue<1e-5 & file$rs_id %in% clump1e.5$rs_id, na.rm = T)
  snp_number_hla[2] <- snp_number_hla[2] + sum(file$pvalue<1e-4 & file$rs_id %in% clump1e.4$rs_id, na.rm = T)
  snp_number_hla[3] <- snp_number_hla[3] + sum(file$pvalue<1e-3 & file$rs_id %in% clump1e.3$rs_id, na.rm = T)
  snp_number_hla[4] <- snp_number_hla[4] + sum(file$pvalue<1e-2 & file$rs_id %in% clump1e.2$rs_id, na.rm = T)
  snp_number_hla[5] <- snp_number_hla[5] + sum(file$pvalue<5e-2 & file$rs_id %in% clump5e.2$rs_id, na.rm = T)
}

for(i in 1:5){
  PRS_HLA[,i+1] = PRS_HLA[,i+1] / weight_hla[i]
}

for(i in 1:5){
  PRS_HLA[,i+1] = scale(PRS_HLA[,i+1])
}
HLA_Data <- read.csv("/home/guanwh/cao00128/ibd/prs/dekaf_gen03_demo_20210325.csv")
data.hla <- merge(pheno, PRS_HLA, by="IID")
data.hla <- merge(data.hla, HLA_Data[,c(2,47)], by.x="IID", by.y="rgsn")
data.hla$nummismatch <- as.numeric(data.hla$nummismatch)


table2 = matrix(NA, nrow=5, ncol=4)
table2 <- as.data.frame(table2)
colnames(table2) <- c("p cutoff", "2.5% exp coef", "97.5% exp coef", "p value")
table2[,1] = c("1e-5","1e-4","1e-3","1e-2","5e-2")
covars <- "age_at_tx+genderc+PRA_posneg+nonkidneytx+nummismatch"
depvar <- "t2ar_ever,ar_ever"
IBS <- c("IBS1e.5","IBS1e.4","IBS1e.3","IBS1e.2","IBS5e.2")

for(i in 3:5){
  res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.hla)  
  res.cox1 <- summary(res.cox1)
  table2[i,2:3] <- res.cox1$conf.int[1,3:4]
  table2[i,4] <- res.cox1$coefficients[1,5]
}
table2 <- cbind(table2, snp_number_hla)
colnames(table2)[5] <- "PRS SNP numbers"
write.csv(table2, "/home/guanwh/cao00128/ibd/prs/dekaf_HLA_mismatch_HLAmatching.csv")

# for(i in 1:5){
#   res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.hla, subset = !(nummismatch==0))
#   print(res.cox1)
#   res.cox1 <- summary(res.cox1)
#   table2[i,2:3] <- res.cox1$conf.int[1,3:4]
#   table2[i,4] <- res.cox1$coefficients[1,5]
# }
# write.csv(table2, "/home/guanwh/cao00128/ibd/prs/dekaf_HLA_mismatch_imputed_nonperfectmatch.csv")


# pdf("/home/guanwh/cao00128/ibd/prs/dekaf_nosplit_pvalue_HLA_ldimputed.pdf")
# par(mfrow=c(1,1))
# main=paste0("SNP p cutoff = ", c("1e-5","1e-4","1e-3","1e-2","5e-2"), " (HLA region)")
# for(i in 1:5){
#   hist(data.hla[,IBS[i]],
#        xlab="PRS",
#        ylab="Count",
#        main=main[i],
#        xlim=c(0,2))
# }
# dev.off()



## non-HLA
# prune <- read.table("/home/guanwh/cao00128/ibd/prs/snps_pruned_imputed.txt")
# pheno<-read.csv("/home/guanwh/shared/DeKAF/GEN03/pheno/AR_EA_PC_20171012.csv")
prsfile <- list.files("/home/guanwh/cao00128/ibd/prs/PRS2GEN03_clump_HLAmatching/", "csv", full.names = T)
PRS_wg <- lapply(prsfile, read.csv, row.names = 1)
prs_sum <- Reduce("+", lapply(PRS_wg, function(x){x[,2:6]}))
PRS_wg <- cbind(PRS_wg[[1]][,1],prs_sum)
colnames(PRS_wg)[1] <- "IID"
gwas.files <- list.files("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/", "csv", full.names = T)
weight <- rep(0, 5)
snp_number <- rep(0, 5)
for(i in 1:length(gwas.files)){
  print(i)
  file <- read.csv(gwas.files[i], row.names = 1)
  file <- file[file$ibs_freq>0.05,]
  # file <- file[sub(":.*", "", file$rs_id) %in% prune$V1 | file$rs_id %in% prune$V1,]
  weight[1] <- weight[1] + sum(file$beta[which(file$pvalue<1e-5 & file$rs_id %in% clump1e.5$rs_id)], na.rm = T)
  weight[2] <- weight[2] + sum(file$beta[which(file$pvalue<1e-4 & file$rs_id %in% clump1e.4$rs_id)], na.rm = T)
  weight[3] <- weight[3] + sum(file$beta[which(file$pvalue<1e-3 & file$rs_id %in% clump1e.3$rs_id)], na.rm = T)
  weight[4] <- weight[4] + sum(file$beta[which(file$pvalue<1e-2 & file$rs_id %in% clump1e.2$rs_id)], na.rm = T)
  weight[5] <- weight[5] + sum(file$beta[which(file$pvalue<5e-2 & file$rs_id %in% clump5e.2$rs_id)], na.rm = T)
  snp_number[1] <- snp_number[1] + sum(file$pvalue<1e-5 & file$rs_id %in% clump1e.5$rs_id, na.rm = T)
  snp_number[2] <- snp_number[2] + sum(file$pvalue<1e-4 & file$rs_id %in% clump1e.4$rs_id, na.rm = T)
  snp_number[3] <- snp_number[3] + sum(file$pvalue<1e-3 & file$rs_id %in% clump1e.3$rs_id, na.rm = T)
  snp_number[4] <- snp_number[4] + sum(file$pvalue<1e-2 & file$rs_id %in% clump1e.2$rs_id, na.rm = T)
  snp_number[5] <- snp_number[5] + sum(file$pvalue<5e-2 & file$rs_id %in% clump5e.2$rs_id, na.rm = T)
}
PRS_HLA <- read.csv("/home/guanwh/cao00128/ibd/prs/PRS_gen03_HLA_HLAmatching.csv", row.names = 1)
gwas.files <- c("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/Minnesota_GoNL_1KG_chr6:25-30.dose.csv",
                "/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/Minnesota_GoNL_1KG_chr6:30-35.dose.csv")
weight_hla <- rep(0, 5)
snp_number_hla <- rep(0, 5)
for(i in 1:length(gwas.files)){
  print(i)
  file <- read.csv(gwas.files[i], row.names = 1)
  file <- file[file$ibs_freq>0.05,]
  # file <- file[sub(":.*", "", file$rs_id) %in% prune$V1 | file$rs_id %in% prune$V1,]
  file <- file[file$position >= 25759242 & file$position <= 33534827, ]
  weight_hla[1] <- weight_hla[1] + sum(file$beta[which(file$pvalue<1e-5 & file$rs_id %in%  clump1e.5$rs_id)], na.rm = T)
  weight_hla[2] <- weight_hla[2] + sum(file$beta[which(file$pvalue<1e-4 & file$rs_id %in%  clump1e.4$rs_id)], na.rm = T)
  weight_hla[3] <- weight_hla[3] + sum(file$beta[which(file$pvalue<1e-3 & file$rs_id %in%  clump1e.3$rs_id)], na.rm = T)
  weight_hla[4] <- weight_hla[4] + sum(file$beta[which(file$pvalue<1e-2 & file$rs_id %in%  clump1e.2$rs_id)], na.rm = T)
  weight_hla[5] <- weight_hla[5] + sum(file$beta[which(file$pvalue<5e-2 & file$rs_id %in%  clump5e.2$rs_id)], na.rm = T)
  snp_number_hla[1] <- snp_number_hla[1] + sum(file$pvalue<1e-5 & file$rs_id %in%  clump1e.5$rs_id, na.rm = T)
  snp_number_hla[2] <- snp_number_hla[2] + sum(file$pvalue<1e-4 & file$rs_id %in%  clump1e.4$rs_id, na.rm = T)
  snp_number_hla[3] <- snp_number_hla[3] + sum(file$pvalue<1e-3 & file$rs_id %in%  clump1e.3$rs_id, na.rm = T)
  snp_number_hla[4] <- snp_number_hla[4] + sum(file$pvalue<1e-2 & file$rs_id %in%  clump1e.2$rs_id, na.rm = T)
  snp_number_hla[5] <- snp_number_hla[5] + sum(file$pvalue<5e-2 & file$rs_id %in%  clump5e.2$rs_id, na.rm = T)
}

PRS_nonHLA = (PRS_wg - PRS_HLA)
PRS_nonHLA[,1]  = PRS_wg[,1] 
for(i in 1:5){
  PRS_nonHLA[,i+1] = (PRS_wg - PRS_HLA)[,i+1] / (weight[i] - weight_hla[i])
}
for(i in 1:5){
  PRS_nonHLA[,i+1] = scale(PRS_nonHLA[,i+1])
}
HLA_Data <- read.csv("/home/guanwh/cao00128/ibd/prs/dekaf_gen03_demo_20210325.csv")
data.nonhla <- merge(pheno, PRS_nonHLA, by="IID")
data.nonhla <- merge(data.nonhla, HLA_Data[,c(2,47)], by.x="IID", by.y="rgsn")
data.nonhla$nummismatch <- as.numeric(data.nonhla$nummismatch)

table2 = matrix(NA, nrow=5, ncol=4)
table2 <- as.data.frame(table2)
colnames(table2) <- c("p cutoff", "2.5% exp coef", "97.5% exp coef", "p value")
table2[,1] = c("1e-5","1e-4","1e-3","1e-2","5e-2")
covars <- "age_at_tx+genderc+PRA_posneg+nonkidneytx+nummismatch"
depvar <- "t2ar_ever,ar_ever"
IBS <- c("IBS1e.5","IBS1e.4","IBS1e.3","IBS1e.2","IBS5e.2")
for(i in 1:5){
  res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.nonhla)  
  res.cox1 <- summary(res.cox1)
  table2[i,2:3] <- res.cox1$conf.int[1,3:4]
  table2[i,4] <- res.cox1$coefficients[1,5]
}
table2 <- cbind(table2, snp_number-snp_number_hla)
colnames(table2)[5] <- "PRS SNP numbers"
write.csv(table2, "/home/guanwh/cao00128/ibd/prs/dekaf_nonHLA_mismatch_HLAmatching.csv")

i=3
relateness_gen03<-read.csv("/home/guanwh/cao00128/ibd/prs/gen03_relateness.csv",row.names = 1)
data.nonhla = merge(data.nonhla, relateness_gen03, by.x = 'IID', by.y = 'V1')
res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.nonhla, subset = p > 0.1)  
summary(res.cox1)



# for(i in 1:5){
#   res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars)), data=data.nonhla, subset = !(nummismatch==0))
#   print(res.cox1)
#   res.cox1 <- summary(res.cox1)
#   table2[i,2:3] <- res.cox1$conf.int[1,3:4]
#   table2[i,4] <- res.cox1$coefficients[1,5]
# }
# sum(data.nonhla$nummismatch==0)
# write.csv(table2, "/home/guanwh/cao00128/ibd/prs/dekaf_nonHLA_mismatch_imputed_nonperfectmatch.csv")


# pdf("/home/guanwh/cao00128/ibd/prs/dekaf_nosplit_pvalue_nonHLA_ldimputed.pdf")
# par(mfrow=c(1,1))
# main=paste0("SNP p cutoff = ", c("1e-5","1e-4","1e-3","1e-2","5e-2"), " (non-HLA region)")
# for(i in 1:5){
#   hist(data.nonhla[,IBS[i]],
#        xlab="PRS",
#        ylab="Count",
#        main=main[i],
#        xlim=c(0,2))
# }
# dev.off()



ibd_gen03 <- read.csv('/home/guanwh/cao00128/ibd/prs/gen03_ibd.csv', row.names = 1)
data <- merge(data.wg,ibd_gen03[,c(1,4)],by.x="IID",by.y="V1")

cor(as.numeric(data$nummismatch), data$IBS1e.5,use="pairwise.complete.obs")
cor(as.numeric(data$nummismatch), data$IBS1e.4,use="pairwise.complete.obs")
cor(as.numeric(data$nummismatch), data$IBS1e.3,use="pairwise.complete.obs")
cor(as.numeric(data$nummismatch), data$IBS1e.2,use="pairwise.complete.obs")
cor(as.numeric(data$nummismatch), data$IBS5e.2,use="pairwise.complete.obs")

covars <- "age_at_tx+genderc+PRA_posneg+nonkidneytx+nummismatch"
depvar <- "t2ar_ever,ar_ever"
data$relateness <- ifelse(data$p>0.4,1,0)
data$nummismatch <- as.numeric(data$nummismatch)
IBS <- c("IBS1e.5","IBS1e.4","IBS1e.3","IBS1e.2","IBS5e.2")
for(i in 1:5){
  res.cox1 <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[i]," + ", covars, " + strata(relateness)")), data=data)  
  print(summary(res.cox1))
  res.cox1 <- summary(res.cox1)
  table2[i,2:3] <- res.cox1$conf.int[1,3:4]
  table2[i,4] <- res.cox1$coefficients[1,5]
}
write.csv(table2, "/home/guanwh/cao00128/ibd/prs/dekaf_nonHLA_mismatch_imputed_stratified.csv")

# cor(data$IBS1e.4, data$prob_ibd)
# data <- merge(data.hla,ibd_gen03,by.x="IID",by.y="V1")
# cor(data$IBS1e.4, data$prob_ibd)
# data <- merge(data.nonhla,ibd_gen03,by.x="IID",by.y="V1")
# cor(data$IBS1e.4, data$prob_ibd)


data <- merge(data.wg,ibd_gen03[,c(1,4)],by.x="IID",by.y="V1")
test <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[2]," + ", covars)), data=merge(data.wg,ibd_gen03[,c(1,4)],by.x="IID",by.y="V1"), subset = p>0.1)  
summary(test)
test <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[2]," + ", covars)), data=merge(data.hla,ibd_gen03[,c(1,4)],by.x="IID",by.y="V1"), subset = p>0.1)  
summary(test)
test <- coxph(as.formula(paste0("Surv(", depvar, ") ~ ", IBS[2]," + ", covars)), data=merge(data.nonhla,ibd_gen03[,c(1,4)],by.x="IID",by.y="V1"), subset = p>0.1)  
summary(test)
#################################### power calculation
library(powerSurvEpi)
X <- rbinom(363, 2, 0.05)
failureFlag <- sample(c(0, 1), 100, prob = c(0.5, 0.5), replace = TRUE)
powerEpi(X, failureFlag = failureFlag, 
         n = 363, theta = 1.5, alpha = 0.05)

set.seed(123456)
X1 <- rnorm(100, mean = 0, sd = 0.3126)
X2 <- sample(c(0, 1), 100, replace = TRUE)
failureFlag <- sample(c(0, 1), 100, prob = c(0.25, 0.75), replace = TRUE)
dat <- data.frame(X1 = X1, X2 = X2, failureFlag = failureFlag)
powerEpiCont(formula = X1 ~ X2,
             dat = dat,
             var.X1 = "X1",
             var.failureFlag = "failureFlag",
             n = 107,
             theta = exp(1),
             alpha = 0.05)


powerEpiCont.default(n = 363,
                     theta = 1.5,
                     sigma2 = 0.05*0.95,
                     psi = 0.168,
                     rho2 =0,
                     alpha = 0.05)