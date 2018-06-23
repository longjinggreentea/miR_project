
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
load("miR_clin_Assay_data.RData")

# dataset to be used: clin_miR_g3_pOV_work   miR_ovrlap_clin_g3_pOV_work  miR_clin_g3_OV_patient_ID

# plot for miR data
X11()
par(mfrow=c(1,2))
hist(unlist(miR_ovrlap_clin_g3_pOV_work),xlab="miR TPM", main="miR TPM", freq=FALSE)
hist(log2(unlist(miR_ovrlap_clin_g3_pOV_work)), xlab="Log2(miR TPM)", main="Log2(miR TPM)", freq=FALSE)
# log2 transformation miR data
log.miR.clin.g3.ov <- log2(miR_ovrlap_clin_g3_pOV_work )

# survival object
library(survival)
os.time.g3 <- clin_miR_g3_pOV_work$OS 
names(os.time.g3) <- rownames(clin_miR_g3_pOV_work) 
os.status.g3 <- clin_miR_g3_pOV_work$vital.status.collapsed
names(os.status.g3) <- rownames(clin_miR_g3_pOV_work)
os.g3.event <- (os.status.g3 == "dead")
os.g3.object <- Surv(os.time.g3/30.5, os.g3.event)



# cox function for significant miR screening
my.coxph.fun <- function(x){
  load("miR_g3_OV_survival.RData")
  res.cox <- coxph(os.g3.object ~ x)
  sum.cox <- summary(res.cox)
  return(c(sum.cox$coefficients[c(1,5)], sum.cox$conf.int[c(1, 3, 4)]))
}

sig.miR.g3.ov.screen <- t(apply(log.miR.clin.g3.ov, 1, my.coxph.fun))
sig.miR.g3.ov.screen.df <- data.frame(sig.miR.g3.ov.screen); 
colnames(sig.miR.g3.ov.screen.df)<- c("coef", "p-val","Hazard.Ratio", "lower", "upper")
sig.mature.miR.id <- (rownames(sig.miR.g3.ov.screen.df))[which(sig.miR.g3.ov.screen.df[,"p-val"]<0.05)]

# [1] "hsa-mir-100"   "hsa-mir-101-1" "hsa-mir-124-3" "hsa-mir-1292"  "hsa-mir-140"   "hsa-mir-144"   "hsa-mir-16-1" 
# [8] "hsa-mir-16-2"  "hsa-mir-1975"  "hsa-mir-200a"  "hsa-mir-205"   "hsa-mir-219-2" "hsa-mir-28"    "hsa-mir-410"  
# [15] "hsa-mir-506" 

sig.miR.cox <- sig.miR.g3.ov.screen.df[sig.mature.miR.id,]

sig.miR.cox.HR.lt.1 <- sig.miR.cox[which(sig.miR.cox[,3]< 1), ]
sig.miR.cox.HR.gt.1 <- sig.miR.cox[which(sig.miR.cox[,3]> 1), ]

sig.log.miR <- log.miR.clin.g3.ov[sig.mature.miR.id,]

save(clin_miR_g3_pOV_work, os.g3.object, sig.miR.cox, log.miR.clin.g3.ov, miR_ovrlap_clin_g3_pOV_work, 
     sig.mature.miR.id,
     file="miR_g3_OV_survival.RData")

#------------------------------------------------------------
# Now check the correlation between each pair of sigificant (survival) miR.
cor.test()
cor.test(as.numeric(log.miR.clin.g3.ov["hsa-mir-100",]), as.numeric(log.miR.clin.g3.ov["hsa-mir-101-1",]))
X11()
plot(as.numeric(log.miR.clin.g3.ov["hsa-mir-100",]), as.numeric(log.miR.clin.g3.ov["hsa-mir-101-1",]), 
     xlab="mir-100", ylab="mir-101-1", 
     xlim = c(min(as.numeric(log.miR.clin.g3.ov["hsa-mir-100",])), max(as.numeric(log.miR.clin.g3.ov["hsa-mir-100",]))), 
     ylim = c(min(as.numeric(log.miR.clin.g3.ov["hsa-mir-101-1",])), max(as.numeric(log.miR.clin.g3.ov["hsa-mir-101-1",])))
     )

mat.target <-log.miR.clin.g3.ov[sig.mature.miR.id, ]; n.row <-nrow(mat.target)
# mat.target <-log.miR.clin.g3.ov[rownames(sig.miR.cox.HR.lt.1), ]; n.row <-nrow(mat.target)

miR.1 <- c(); miR.2 <- c(); my.cor <- c(); my.pval<- c()

for (i in (1:(n.row - 1))){
  for (j in ((i+1):n.row)){
    
    cor.reslt <- cor.test(as.numeric(mat.target[i,]), as.numeric(mat.target[j,]))    
    miR.1<- c(miR.1, rownames(mat.target)[i])
    miR.2<- c(miR.2, rownames(mat.target)[j])
    my.cor <- c( my.cor, cor.reslt$estimate)
    my.pval <- c(my.pval, cor.reslt$'p.val')
  }
}
mat.target.cor.df <- data.frame(miR.1, miR.2, my.cor, my.pval )
colnames(mat.target.cor.df) <- c("miR_1", "miR_2", "corr", "pval")

sig.mat.target.cor <- mat.target.cor.df[which((mat.target.cor.df[,4])<0.05), ]








