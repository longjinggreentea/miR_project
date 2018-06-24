###-----------
# work on mRNA RSEM 
setwd("C:/Users/Lixin/Desktop/RPCI_2018Mar/OV_mRNA_seq/OV_mRNA_seq")
OV.TCGA.mRNA.RSEM <- read.table("OV.uncv2.mRNAseq_RSEM_all.txt", header = TRUE, sep ="\t", row.names = 1)
dim(OV.TCGA.mRNA.RSEM)
sample.id.RSEM <- colnames(OV.TCGA.mRNA.RSEM)
sample.id.RSEM <- substr(sample.id.RSEM, 1, 12)
sample.id.RSEM.1 <- gsub("[.]", "-", sample.id.RSEM) # 307

# Update mRNA dataset colnames:
colnames(OV.TCGA.mRNA.RSEM) <- sample.id.RSEM.1

# Recall: qnorm.log.miR.exp.clin.working # 298 382 # ? But where can I find this qnormalized miR dataset? 
####  patient.id.miR <- colnames(qnorm.log.miR.exp.clin.working)  
####  
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
load("miR_clin_Assay_data.RData") 
dim(miR_ovrlap_clin_g3_pOV_work) # 298x382
patient.id.miR <- colnames(miR_ovrlap_clin_g3_pOV_work) 

sample.id.RSEM.miR <- intersect(patient.id.miR, sample.id.RSEM.1) # 244

sub_miR_ovlp_clin_g3_pOV_work<- miR_ovrlap_clin_g3_pOV_work[, sample.id.RSEM.miR] # 298x244

#### Working dataset by combine mRNA and miR (with the common patient ID)
sub.OV.TCGA.mRNA.RSEM <- OV.TCGA.mRNA.RSEM[, sample.id.RSEM.miR] # 20531   245

#### Trim rownames and remove duplicated trimmed rownames
trim.RowName.RSEM <- sapply(rownames(sub.OV.TCGA.mRNA.RSEM), function(x) unlist(strsplit(x, '\\|'))[[1]]) # 20531

dup.RSEM.ind <- which(duplicated(trim.RowName.RSEM)==TRUE)
# length(trim.RowName.RSEM[duplicated(trim.RowName.RSEM)]) # 29

#### Remove duplicated gene/feature: for RSEM expression data

sub.OV.TCGA.mRNA.RSEM.working <- sub.OV.TCGA.mRNA.RSEM[-dup.RSEM.ind, ]
head(rownames(sub.OV.TCGA.mRNA.RSEM.working))
rownames(sub.OV.TCGA.mRNA.RSEM.working) <- sapply(rownames(sub.OV.TCGA.mRNA.RSEM.working),function(x) unlist(strsplit(x, '\\|'))[[1]])
head(rownames(sub.OV.TCGA.mRNA.RSEM.working))
sub.OV.TCGA.mRNA.RSEM.working <- sub.OV.TCGA.mRNA.RSEM.working[-1,]
head(rownames(sub.OV.TCGA.mRNA.RSEM.working))
# length(duplicated(rownames(sub.OV.TCGA.mRNA.RSEM.working)))
dim(sub.OV.TCGA.mRNA.RSEM.working) #  20501   244

#--------------------------------------------------------------------------------------------
#### check missing data in sub.OV.TCGA.mRNA.RSEM.working                                     #|  
# sub.OV.TCGA.mRNA.RSEM.working                                                              #|   
                                                                                             #|   
NAcheck.sub.OV.TCGA.mRNA.RSEM.working <- is.na(sub.OV.TCGA.mRNA.RSEM.working)
mRNA.RSEM.NA.count.byrow <- apply(NAcheck.sub.OV.TCGA.mRNA.RSEM.working, 1, sum)
which(!(mRNA.RSEM.NA.count.byrow==0))

                                                                                             #|
                                                                                             #| 
#### check 0 in sub.OV.TCGA.mRNA.RSEM.working                                                #|
zero.check <- function(x){ifelse((x==0), 1, 0)}
zero.tform.sub.mRNA.RSEM <-zero.check(sub.OV.TCGA.mRNA.RSEM.working)
zero.count.sub.mRNA.RSEM <- apply(zero.tform.sub.mRNA.RSEM, 1, sum)
length(names(zero.count.sub.mRNA.RSEM))
gene.keep.id <- (names(zero.count.sub.mRNA.RSEM))[which(zero.count.sub.mRNA.RSEM==0)] # length(gene.keep.id) 14014 genes
gene.rm.id <- (names(zero.count.sub.mRNA.RSEM))[which(zero.count.sub.mRNA.RSEM!=0)] # length(gene.rm.id) 6487 genes


#### percentage of zero read and replace zero by 0.0001
summary(unlist(sub.OV.TCGA.mRNA.RSEM.working)) 
try5 <- unlist(zero.tform.sub.mRNA.RSEM)
table(try5)
#       0       1 
# 4356302  645942 

nrow.RSEM.sub <- nrow(sub.OV.TCGA.mRNA.RSEM.working)
ncol.RSEM.sub <- ncol(sub.OV.TCGA.mRNA.RSEM.working)
nozero.sub.mRNA.RSEM <- matrix(rep(0, (nrow.RSEM.sub*ncol.RSEM.sub)), nrow.RSEM.sub, ncol.RSEM.sub)

for (i in 1:nrow.RSEM.sub){
      for (j in 1:ncol.RSEM.sub){
           if(zero.tform.sub.mRNA.RSEM[i, j]==0){ 
             nozero.sub.mRNA.RSEM[i,j] <- sub.OV.TCGA.mRNA.RSEM.working[i,j]  
           } else {
             nozero.sub.mRNA.RSEM[i,j] <- 0.0001
           }
         }
     }

dim(nozero.sub.mRNA.RSEM)  
try6 <- unlist(nozero.sub.mRNA.RSEM)
length(which(try6==0.0001))
rownames(nozero.sub.mRNA.RSEM)<-rownames(sub.OV.TCGA.mRNA.RSEM.working)
colnames(nozero.sub.mRNA.RSEM)<-colnames(sub.OV.TCGA.mRNA.RSEM.working)
#--------------------------------------------------------------------------------------------

X11()
par(mfrow=c(2,2))
hist(unlist(nozero.sub.mRNA.RSEM), main="mRNA_TPM", xlab="Log2(mRNA_TPM)")
hist(unlist(log2(nozero.sub.mRNA.RSEM)), main="Log2(mRNA_TPM", xlab="Log2(mRNA_TPM)")
hist(unlist(sub_miR_ovlp_clin_g3_pOV_work), main="miR_TPM", xlab="Log2(miR_TPM)")
hist(unlist(log2(sub_miR_ovlp_clin_g3_pOV_work)), main="Log2(miR_TPM", xlab="Log2(miR_TPM)")

log.nozero.sub.mRNA.RSEM <- log2(nozero.sub.mRNA.RSEM)
log.sub_miR_ovlp_clin_g3_pOV_work <- log2(sub_miR_ovlp_clin_g3_pOV_work)

sc.log.nozero.sub.mRNA.RSEM <- scale(log.nozero.sub.mRNA.RSEM)
sc.log.sub_miR_ovlp_clin_g3_pOV_work <- scale(log.sub_miR_ovlp_clin_g3_pOV_work)

X11()
par(mfrow=c(2,1))
boxplot(sc.log.nozero.sub.mRNA.RSEM, main="mRNA")
boxplot(sc.log.sub_miR_ovlp_clin_g3_pOV_work, main="miR")

save(sc.log.nozero.sub.mRNA.RSEM, sc.log.sub_miR_ovlp_clin_g3_pOV_work, 
     log.nozero.sub.mRNA.RSEM, log.sub_miR_ovlp_clin_g3_pOV_work,
     sub.OV.TCGA.mRNA.RSEM.working, sub_miR_ovlp_clin_g3_pOV_work,
     OV.TCGA.mRNA.RSEM, miR_ovrlap_clin_g3_pOV_work, sample.id.RSEM.miR,
     file="sub_miR_mRNA_RSEM_integration_ASSAY.RDATA")
 # @ C:\Users\Lixin\Desktop\R script downloaded\data\TCGA_miR\OV_clinical_merge_2016_level1\OV_merge_clinical_level_1_2016
 # @ C:\Users\Lixin\Desktop\RPCI_2018Mar\OV_mRNA_seq\OV_mRNA_seq
#--------------------------------------------------------------------------------------------
###  sample.id.RSEM.miR  244
###  sc.log.nozero.sub.mRNA.RSEM  20501x244
###  sc.log.sub_miR_ovlp_clin_g3_pOV_work 298x244

##  variance filtering:
var.RSEM.sub <- apply(sc.log.nozero.sub.mRNA.RSEM, 1, var)
var.RSEM.sub.50.perctl <- quantile(var.RSEM.sub, .50)
var.RSEM.sub.75.perctl <- quantile(var.RSEM.sub, .75)
var.RSEM.sub.80.perctl <- quantile(var.RSEM.sub, .80)

var.select.gene.id <- (names(var.RSEM.sub))[which(var.RSEM.sub > var.RSEM.sub.75.perctl)] # 5125
var.select.sub.mRNA.RSEM <- sc.log.nozero.sub.mRNA.RSEM[var.select.gene.id,] # 5125x244


##  Cox regression:
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
load("miR_g3_OV_survival.RData") # os.g3.object     clin_miR_g3_pOV_work
### subset clin_miR_g3_pOV_work
sub.clin.miR.mRNA.RSEM.g3.pOV <- clin_miR_g3_pOV_work[sample.id.RSEM.miR,] # 244x6
library(survival)
sub.os.time.g3 <- sub.clin.miR.mRNA.RSEM.g3.pOV$OS 
names(sub.os.time.g3) <- rownames(sub.clin.miR.mRNA.RSEM.g3.pOV) 
sub.os.status.g3 <- sub.clin.miR.mRNA.RSEM.g3.pOV$vital.status.collapsed
names(sub.os.status.g3) <- rownames(sub.clin.miR.mRNA.RSEM.g3.pOV)
sub.os.g3.event <- (sub.os.status.g3 == "dead")
sub.os.g3.object <- Surv(sub.os.time.g3/30.5, sub.os.g3.event)
save(sub.os.g3.object, file="sub_os_g3_survival_obj.RData") # @ C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016

my.coxph.fun <- function(x){
  load("sub_os_g3_survival_obj.RData")
  res.cox <- coxph(sub.os.g3.object ~ x)
  sum.cox <- summary(res.cox)
  return(c(sum.cox$coefficients[c(1,5)], sum.cox$conf.int[c(1, 3, 4)]))
}

sig.mRNA.g3.ov.screen <- t(apply(var.select.sub.mRNA.RSEM, 1, my.coxph.fun))
sig.mRNA.g3.ov.screen.df <- data.frame(sig.mRNA.g3.ov.screen); 
colnames(sig.mRNA.g3.ov.screen.df)<- c("coef", "p-val","Hazard.Ratio", "lower", "upper")
sig.sub.mRNA.id <- (rownames(sig.mRNA.g3.ov.screen.df))[which((sig.mRNA.g3.ov.screen.df[,"p-val"]<0.05) & (sig.mRNA.g3.ov.screen.df[,"Hazard.Ratio"]>1))]
# 188
sig.mRNA.cox <- sig.mRNA.g3.ov.screen.df[sig.sub.mRNA.id,]

ordered.sig.mRNA.cox <- sig.mRNA.cox[order(sig.mRNA.cox$'p-val'), ]

save(ordered.sig.mRNA.cox, var.select.sub.mRNA.RSEM, sub.clin.miR.mRNA.RSEM.g3.pOV, file="sig_mRNA_cox.RData")
    ### ordered.sig.mRNA.cox  include 188 mRNA with Hazard Ratio greater than 1 significantly. 
    ### var.select.sub.mRNA.RSEM include 5125 mRNA with variance greater than 75 percentile.
    ### sub.clin.miR.mRNA.RSEM.g3.pOV include clinical data for 244 patients (grade 3 and ovary as primary site)
#--------------------------------------------------------------------------------------------------------------------------
##  correlation: 
###  sc.log.sub_miR_ovlp_clin_g3_pOV_work # 298x244
###  var.select.sub.mRNA.RSEM             # 5125x244

sub.miR.dat <- sc.log.sub_miR_ovlp_clin_g3_pOV_work
sub.mRNA.dat <- var.select.sub.mRNA.RSEM 

identical(colnames(sub.miR.dat), colnames(sub.mRNA.dat)); nmiR<- nrow(sub.miR.dat); nmRNA <- nrow(sub.mRNA.dat)

n.miR <- nrow(sc.log.sub_miR_ovlp_clin_g3_pOV_work)
n.mRNA <- nrow(var.select.sub.mRNA.RSEM)
cor.mat1 <- matrix(rep('0', 2*n.miR*n.mRNA), nrow=n.miR*n.mRNA, ncol=2)
cor.mat2 <- matrix(rep(0, 2*n.miR*n.mRNA), nrow=n.miR*n.mRNA, ncol=2)

# my.pval <- c(); my.cor.est <- c(); miR.name <- c(); mRNA.name <- c()

for (i in (1:nrow(sub.miR.dat))){
  for (j in (1:nrow(sub.mRNA.dat))){
      
         cor.reslt <- cor.test(sub.miR.dat[i,], sub.mRNA.dat[j,])    
         cor.mat1[((i-1)*n.miR)+j, c(1,2)] <- c((rownames(sub.miR.dat))[i], (rownames(sub.mRNA.dat))[j]) 
         cor.mat2[((i-1)*n.miR)+j, c(1,2)] <- c(cor.reslt$estimate, cor.reslt$'p.value')
  }
}

cor.sub.mRNA.miR <- data.frame(cor.mat1[,c(1,2)], cor.mat2[,c(1,2)])
colnames(cor.sub.mRNA.miR) <- c("miR", "mRNA", "Corr", "pval")

save(cor.sub.mRNA.miR, file="miR_mRNA_pear_correlation.RData")
    ### cor.sub.mRNA.miR include the correlation between 298 miR and 5125 mRNA. 


length(which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]< -0.2))) # 168
sig.neg.cor.m <- cor.sub.mRNA.miR[which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]< -0.2)),]



length(which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]< 0) )) # 2947
length(which((cor.sub.mRNA.miR[,4] < 0.001 ) & (cor.sub.mRNA.miR[,3]< -0.2))) # 114
sig.neg.cor.m <- cor.sub.mRNA.miR[which((cor.sub.mRNA.miR[,4] < 0.001 ) & (cor.sub.mRNA.miR[,3]< -0.2)),]

unique(cor.sub.mRNA.miR[which((cor.sub.mRNA.miR[,4] < 0.01 ) & (cor.sub.mRNA.miR[,3]< 0)), 1])

# length(which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]> 0)))
# length(which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]> 0.2))) # 171
# sig.pos.cor.m <- cor.sub.mRNA.miR[which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]> 0.2)),]

cor.sub.mRNA.miR[which((cor.sub.mRNA.miR[,1]=='hsa-mir-100') & ((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]< 0))),]

ordered.cor.est.mRNA.miR<- cor.sub.mRNA.miR[order(cor.sub.mRNA.miR$Corr),]
order100.cor.miR.mRNA <- ordered.cor.est.mRNA.miR[c(1:100),]
           whol.top100.miR <- as.character(order100.cor.miR.mRNA[,1])
           whol.top100.mRNA <- as.character(order100.cor.miR.mRNA[,2])                                                

top100.miR <- unique(as.character(order100.cor.miR.mRNA[,1]))
top100.mRNA <- unique(as.character(order100.cor.miR.mRNA[,2]))

neg.cor.mRNA.miR <- cor.sub.mRNA.miR[(which(cor.sub.mRNA.miR[,3]<0)),] 
ordered.pval.neg.cor.mRNA.miR<- neg.cor.mRNA.miR[order(neg.cor.mRNA.miR$pval),]
order100.pval.neg.cor.miR.mRNA <- ordered.pval.neg.cor.mRNA.miR[c(1:100),]
top100.p.miR <- unique(as.character(order100.pval.neg.cor.miR.mRNA[,1]))
top100.p.mRNA <- unique(as.character(order100.pval.neg.cor.miR.mRNA[,2]))

# ordered.sig.mRNA.cox
intersect(rownames(ordered.sig.mRNA.cox), unique(as.character(order100.pval.neg.cor.miR.mRNA[,2])))

# load("miR_g3_OV_survival.RData")  --> sig.mature.miR.id
intersect(sig.mature.miR.id,unique(as.character(order100.cor.miR.mRNA[,2])))

neg.cor.miR.id.mRNA <-(cor.sub.mRNA.miR[,1])[which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]< 0) )]
   #  2947  miR id (unique 298) significant negatively correlated with mRNA
neg.cor.mRNA.id.miR <-(cor.sub.mRNA.miR[,2])[which((cor.sub.mRNA.miR[,4] < 0.05 ) & (cor.sub.mRNA.miR[,3]< 0) )]
   #  2947  mRNA id (unique 436) significant negatively correlated with miR 


#### negatrive correlation between miR-mRNA 2947 pairs 

intersect(unique(neg.cor.miR.id.mRNA), sig.mature.miR.id)
### 15 common miR 
# [1] "hsa-mir-100"   "hsa-mir-101-1" "hsa-mir-124-3" "hsa-mir-1292"  "hsa-mir-140"   "hsa-mir-144"   "hsa-mir-16-1" 
# [8] "hsa-mir-16-2"  "hsa-mir-1975"  "hsa-mir-200a"  "hsa-mir-205"   "hsa-mir-219-2" "hsa-mir-28"    "hsa-mir-410"  
# [15] "hsa-mir-506"  
intersect(unique(neg.cor.mRNA.id.miR), sig.sub.mRNA.id)
### 16 common mRNA 
# [1] "ADAM21P1" "AGXT"     "ASCL1"    "ATP1A2"   "AQP12A"   "ASXL3"    "ANGPT4"   "AMELX"    "FCAR"     "HS6ST3"   "HSD3B1"  
# [12] "KLK3"     "LMX1B"    "OR8B12"   "PI3"      "SERPINB2"




