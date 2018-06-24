#################################
# STEP 4

# @ C:\Users\Lixin\Desktop\R script downloaded\data\TCGA_miR\OV_clinical_merge_2016_level1\OV_merge_clinical_level_1_2016
# @ C:\Users\Lixin\Desktop\RPCI_2018Mar\OV_mRNA_seq\OV_mRNA_seq


# OV_miRNA_working_h1.R

##### STEP 1: 

##### STEP 2: 

load("miR_g3_OV_survival.RData")
     ### save(clin_miR_g3_pOV_work, os.g3.object, sig.miR.cox, log.miR.clin.g3.ov, miR_ovrlap_clin_g3_pOV_work, 
         ###  sig.mature.miR.id,
         ###  file="miR_g3_OV_survival.RData")

load("correlation_sig_miR.RData")
     ### save(sig.miR.cox.HR.lt.1, sig.miR.cox.HR.gt.1, sig.mat.neg.target.cor, sig.mat.target.cor,
        ###   file = "correlation_sig_miR.RData")

##### STEP 3: 

load("sub_miR_mRNA_RSEM_integration_ASSAY.RDATA")
    ### save(sc.log.nozero.sub.mRNA.RSEM, sc.log.sub_miR_ovlp_clin_g3_pOV_work, 
        ###  log.nozero.sub.mRNA.RSEM, log.sub_miR_ovlp_clin_g3_pOV_work,
        ###  sub.OV.TCGA.mRNA.RSEM.working, sub_miR_ovlp_clin_g3_pOV_work,
        ###  OV.TCGA.mRNA.RSEM, miR_ovrlap_clin_g3_pOV_work, sample.id.RSEM.miR,
        ###  file="sub_miR_mRNA_RSEM_integration_ASSAY.RDATA")

load("miR_mRNA_pear_correlation.RData")
    ###  save(cor.sub.mRNA.miR, file="miR_mRNA_pear_correlation.RData")
         ###  cor.sub.mRNA.miR include the correlation between 298 miR and 5125 mRNA. 

load("sig_mRNA_cox.RData")
     ### save(ordered.sig.mRNA.cox, var.select.sub.mRNA.RSEM, sub.clin.miR.mRNA.RSEM.g3.pOV, file="sig_mRNA_cox.RData")
         ### ordered.sig.mRNA.cox  include 188 mRNA with Hazard Ratio greater than 1 significantly. 
         ### var.select.sub.mRNA.RSEM include 5125 mRNA with variance greater than 75 percentile.
         ### sub.clin.miR.mRNA.RSEM.g3.pOV include clinical data for 244 patients (grade 3 and ovary as primary site)



## Data to be used: 

## Using glmnet package: 
   ###  sc.log.sub_miR_ovlp_clin_g3_pOV_work # 298x244
   ###  var.select.sub.mRNA.RSEM             # 5125x244







## Using glmnet package: 
