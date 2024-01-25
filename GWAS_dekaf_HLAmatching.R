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

# gwasfile <- list.files(path="/home/guanwh/cao00128/ibd/prs/GWAS_summary/", pattern = "csv")
# gwasfile <- str_sub(gwasfile, 1, -5)
# 
# dosagefile<- filelist <- read.table("/home/guanwh/shared/DeKAF/GWAS/dosage/filelist.txt", header = F)
# dosagefile <- dosagefile$V1
# 
# required_file <- which(dosagefile %in% setdiff(dosagefile, gwasfile))




filelist <- read.table("/home/guanwh/shared/DeKAF/GWAS/dosage/filelist.txt", header = F)
begin=1
end=length(filelist$V1)
interval=ceiling((end-begin)/100)

for(jobid in (begin+(index-1)*interval):min(end,begin+index*interval)){
  print(paste0("jobid:",jobid))
  # if(! jobid %in% required_file) next
  dosage <- fread(paste0("/home/guanwh/shared/DeKAF/GWAS/dosage/", as.character(filelist$V1[jobid])))
  dosage <- as.data.frame(dosage)
  pheno<-read.csv("/home/guanwh/shared/DeKAF/GWAS/phenotype/AR_data_EA_20151208.csv")
  recipient_donor_IDs.tab <- read.delim("/home/guanwh/shared/DeKAF/GWAS/phenotype/recipient_donor_IDs.tab.txt")
  recipient_donor_IDs.tab$recip_ID<-as.character(recipient_donor_IDs.tab$recip_ID)
  recipient_donor_IDs.tab$donor_ID<-as.character(recipient_donor_IDs.tab$donor_ID)
  HLA_Data <- read.csv("/home/guanwh/cao00128/ibd/prs/dekaf_gen03_demo_20210325.csv")
  HLA_Data <- HLA_Data[HLA_Data$nummismatch!='.',]
  HLA_Data$nummismatch = as.numeric(HLA_Data$nummismatch)
  pheno <- merge(pheno, HLA_Data[,c(2,47)], by.x="gseq", by.y="rgsn")
  
  
  dosage<-dosage[dosage$exp_freq_a1<0.95 & dosage$exp_freq_a1>0.05 & dosage$info>=0.8,]
  sample.geno <- names(dosage)[-c(1:13)]
  sample.geno <- substr(sample.geno, 1, 7)
  recipient_donor_IDs.tab<-recipient_donor_IDs.tab[!is.na(match(recipient_donor_IDs.tab$recip_ID, sample.geno)) & !is.na(match(recipient_donor_IDs.tab$donor_ID, sample.geno)), ]
  
  
  snps <- dosage$rs_id
  geno.info <- dosage[,1:13]
  geno.use <- t(dosage[,-c(1:13)])
  geno.use <- as.data.frame(geno.use)
  rownames(geno.use)<-substr(rownames(geno.use), 1, 7)
  colnames(geno.use) <- snps
  
  geno_r <- geno.use[recipient_donor_IDs.tab$recip_ID,]
  geno_d <- geno.use[recipient_donor_IDs.tab$donor_ID,]
  colnames(geno_d) <- paste0(snps, ".donor")
  rownames(geno_d) <- recipient_donor_IDs.tab$recip_ID
  
  data <- merge(pheno, geno_r, by.x="gseq", by.y="row.names")
  data <- merge(data, geno_d, by.x="gseq", by.y="row.names")
  
  data$PRA_posneg<- as.numeric(data$PRA_posneg)
  data$nonkidneytx<-as.numeric(data$nonkidneytx)
  
  
  covars <- "age_at_tx+gender+PRA_posneg+nonkidneytx+nummismatch"
  wholeoutput <- c()
  depvar <- "t2ar_oa,ar_oa"
  
  GWAS.table <- dosage[,1:13]
  GWAS.table <- cbind(GWAS.table, matrix(NA, nrow=nrow(GWAS.table), ncol=6))
  GWAS.table$ibs_freq = NA
  # set.seed(532)
  # test.sample <- sample(1:nrow(data), round(nrow(data)/3))
  # train.sample <- setdiff(1:nrow(data), test.sample)
  # PRS.score <- matrix(0, nrow=261, ncol=7)
  # PRS.score <- as.data.frame(PRS.score)
  # colnames(PRS.score) <- c("gseq", "pid", "IBS1e-5", "IBS1e-4", "IBS1e-3", "IBS1e-2", "IBS5e-2")
  
  
  
  for(j in 1:length(snps)){
    print(j)
    snpname<-snps[j]
    snpdonor<-paste0(snpname,".donor")
    
    ## Method1. Origin
    #genomic collision: GG for recipient and 
    data$collision1 = ifelse((data[,snpname]>1.9 | data[,snpname]<0.1) & data[,snpdonor]>0.1 & data[,snpdonor]<1.9, 1, 0)
    
    ## Method2.IBS Score: absolute difference
    data$collision2 = abs(data[,snpname]-data[,snpdonor])
    
    if(sum(is.na(data$collision2)) / nrow(data) > 0.1) next
    
    # ## Method3.Incomp Score
    # #######################
    # #if diff = 0, score = 0
    # #otherwise score = 1
    # data$collision3 = ifelse(abs(data[,snpname]-data[,snpdonor])<0.2, 0, 1)
    # 
    # ## Method4.AMS
    # #######################
    # #mismatch if D has allele not in R
    # #sum across both alleles in genotype
    # #Score is either 0, 1, or 2
    # data$collision4=0
    # data$collision4[data[,snpdonor]<0.1 & data[,snpname]>1.9]=2
    # data$collision4[data[,snpdonor]>1.9 & data[,snpname]<0.1]=2
    # data$collision4[data[,snpdonor]>0.9 & data[,snpdonor]<1.1 & data[,snpname]>1.9]=1
    # data$collision4[data[,snpdonor]>0.9 & data[,snpdonor]<1.1 & data[,snpname]<0.1]=1
    # data$collision4[is.na(data[,snpdonor]) | is.na(data[,snpname])]=NA
    # 
    # ## Method5.Binary  Mismatch
    # #######################
    # #mismatch if D has allele not in R
    # #Score is either 0 or 1
    # data$collision5=0
    # data$collision5[data[,snpdonor]>0.9 & data[,snpdonor]<1.1 & data[,snpname]>1.9]=1
    # data$collision5[data[,snpdonor]<0.1 & data[,snpname]>1.9]=1
    # data$collision5[data[,snpdonor]>1.9 & data[,snpname]<0.1]=1
    # data$collision5[data[,snpdonor]>0.9 & data[,snpdonor]<1.1 & data[,snpname]<0.1]=1
    # data$collision5[is.na(data[,snpdonor]) | is.na(data[,snpname])]=NA
    
    
    
    
    
    res.tbl <- NULL
    n = dim(data)[1]
    cat(n, " samples in survival\n")
    output<-vector()
    
    # data.train <- data[train.sample,]
    # data.test <- data[test.sample,]
    
    
    for ( i in 2 )   ## 5 collisions
    {
      
      res <- data.frame(name        = geno.info$rs_id[j],
                        #chr         = chr.name,
                        position    = geno.info$position[j], 
                        exp_freq_a1 = geno.info$exp_freq_a1[j], 
                        info        = geno.info$info[j],
                        certainty   = geno.info$certainty[j], 
                        A0          = geno.info$A0[j], 
                        A1          = geno.info$A1[j],
                        #Effective_N = sum(!is.na(snp.geno)),
                        Effective_N = NA,
                        nevent      = NA,
                        beta        = NA, 
                        se          = NA, 
                        z           = NA,
                        pvalue      = NA)
      
      ### cox regression ###
      tryCatch({
        res.cox <- coxph(as.formula(paste0("Surv(", depvar, ") ~ collision", as.character(i), " + ", 
                                           covars)), data=data)  
        # Surv( t2ar_ever,ar_ever ) ~ collision +  age_at_tx+genderc+PC1+PRA_posneg+priorktx+nonkidneytx+nummismatch
        
        res.cox <- summary(res.cox)
        res$Effective_N <- res.cox$n
        res$nevent <- res.cox$nevent
        res$beta   <- round(res.cox$coef[paste0("collision", as.character(i)), "coef"], 4)
        res$se     <- round(res.cox$coef[paste0("collision", as.character(i)), "se(coef)"], 4)
        res$z      <- round(res.cox$coef[paste0("collision", as.character(i)), "z"], 4)
        res$pvalue <- signif(res.cox$coef[paste0("collision", as.character(i)), "Pr(>|z|)"], 4)
        output<-rbind(output,res)
        GWAS.table[j, 14:19] <- output[,8:13]
        
      }, error = function(e) {
        output<-rbind(output,res)
      })
      
      
    } ## for i
    GWAS.table$ibs_freq[j] = mean(data$collision2, na.rm = T)
    
    # data.test$collision2[is.na(data.test$collision2)] <- 0
    # if(is.na(output[1,13])) next
    # if(output[1,13]<1e-5){
    #   PRS.score[,3] = PRS.score[,3] + data.test$collision2 * output[1,10]
    # }
    # if(output[1,13]<1e-4){
    #   PRS.score[,4] = PRS.score[,4] + data.test$collision2 * output[1,10]
    # }
    # if(output[1,13]<1e-3){
    #   PRS.score[,5] = PRS.score[,5] + data.test$collision2 * output[1,10]
    # }
    # if(output[1,13]<1e-2){
    #   PRS.score[,6] = PRS.score[,6] + data.test$collision2 * output[1,10]
    # }
    # if(output[1,13]<2e-2){
    #   PRS.score[,7] = PRS.score[,7] + data.test$collision2 * output[1,10]
    # }
    
    
    
  }## for j
  
  
  colnames(GWAS.table)[14:19] <- c("Effective_N","nevent","beta","se","z","pvalue")
  
  write.csv(GWAS.table, file = paste0("/home/guanwh/cao00128/ibd/prs/GWAS_summary_HLAmatching/", as.character(filelist$V1[jobid]), ".csv"))
  
  # PRS.score$gseq <- data.test$gseq
  # PRS.score$pid <- data.test$pid
  # write.csv(PRS.score, file = paste0("/home/guanwh/cao00128/ibd/prs/PRS.genome/", as.character(filelist$V1[jobid]), ".csv"))
}





# .libPaths("/home/guanwh/cao00128/Rlibrary/")
# 
# library("data.table")
# library("survival")
# library("dplyr")
# library(stringr)
# 
# gwasfile <- list.files(path="/home/guanwh/cao00128/ibd/prs/GWAS_summary/", pattern = "csv")
# gwasfile <- str_sub(gwasfile, 1, -5)
# 
# dosagefile<- filelist <- read.table("/home/guanwh/shared/DeKAF/GWAS/dosage/filelist.txt", header = F)
# dosagefile <- dosagefile$V1
# 
# setdiff(dosagefile, gwasfile)






