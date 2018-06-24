getwd()
  # @  C:\Users\Lixin\Desktop\R script downloaded\data\TCGA_miR\working_code_collection
setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
load("iso_miR_preprocess_OV.RData") # including pre-processed miR data

setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_miRseq_isoform_level3/OV_miRseq_isoform_level3")
load("TCGA_OV_miR_complete.RData") # another version for pre-processed miR data

setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016")
load("OV_clin_merge_raw_plus_working.RData") # clinical data ready for survival analysis
load("miR_clin_Assay_data.RData")

setwd("C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/code_update_from_office")
load("miR_clin_working_data_2018.RData" )

setwd("C:\Users\Lixin\Desktop\R script downloaded\data\TCGA_miR\code_update_from_office")
load("OV_clin_processed_2018.RData") # --> update OV.clin.working --> OV.clin.working$OS.time
### ------------------------------------------------------------------------------------### 
# data to be used: 
#      ready-for-use               miR.collapse.df.working or miR.collapse.df.working.percent  298*461 (YES)
#                                  OV.clin.merge   patient.Id.ov.clin    591                 
#                                  OV.clin   591*5
#      ready-for-use               OV.clin.working  577*5
#                                 OV.clin.working.OS.time 577*6 (YES)
#
#      miR_clin_overlap_patientID  miR_clin_patient_ID    451
#
#     ready-for-use(miR_version)   clin_for_miR   miR_ovrlap_clin  note: be careful about 'TCGA-24-1431'
#                                                             
#   
#     ? ready-foruse(grade=3, primary.site=ovary) miR.exp.clin.working 298*382  clin.miR.working 382*5
#                                                                             clin_for_miR_OS_time 451*6
#
#     ready-for-use(miR )                clin_for_miR_OS_time  # dim 451x6
#                                        miR_ovrlap_clin       # dim 298x451
#     
#     ready-foruse(grade=3, primary.site=ovary, miR)    clin_miR_g3_pOV_work         # dim 382x6
#                                                       miR_ovrlap_clin_g3_pOV_work  # dim 298x382
#                                           in "miR_clin_Assay_data.RData"
#                                           @ "C:/Users/Lixin/Desktop/R script downloaded/data/TCGA_miR/OV_clinical_merge_2016_level1/OV_merge_clinical_level_1_2016"
#
#
#
### ------------------------------------------------------------------------------------### 
setwd("C:/Users/Lixin/Desktop/RPCI_2018Mar/OV_mRNA_seq/OV_mRNA_seq")
# load("OV_TCGA_mRNA_RSEM_miR_intergration.RData", ex3<-new.env())
ls.str(ex3)
load("sub_miR_mRNA_RSEM_integration_ASSAY.RDATA")

