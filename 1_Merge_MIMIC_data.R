library(xlsx)
library(readxl)
library(openxlsx)
library(stringr)
library(dplyr)
library(mice)
library(limma)
#设置工作目录为R脚本所在的路径
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(digits = 15)
source("3_zl_R_Function_v2.r")

#寻找疾病
rc <- read.csv("d_icd_diagnoses.csv")
IBD_ra_ICD <- rc[grep("ulcerative colitis|crohn's disease",tolower(rc$long_title)),]

#导入诊断人群
All_individual <- read.csv("diagnoses_icd.csv")

IBD_indi <- All_individual %>% filter(icd_code %in% IBD_ra_ICD$icd_code)

IBD_all_ICD <- All_individual %>% filter(hadm_id %in% IBD_indi$hadm_id) %>% left_join(rc)
IBD_all_ICD$long_title <- tolower(IBD_all_ICD$long_title)

Diag_count <- table(IBD_all_ICD$long_title) %>% as.data.frame() %>% arrange(desc(Freq))

#诊断疾病
IBD_ra_diag <- diag_disease_fun(data=IBD_all_ICD,
                human_ID="hadm_id",
                diag_col="long_title",
                diag_disease=c("ulcerative colitis|crohn's disease","nicotine"))


#age和sex统计
patient <- read.csv("patients.csv")
IBD_indi_age_sex <- left_join(IBD_indi,patient)

#加载住院信息
admissions <- read.csv("admissions.csv")
IBD_indi_age_sex_ad <- left_join(IBD_indi_age_sex,admissions)

#加载化验
lab_event <- read.csv("labevents.csv")
#save(lab_event,file = "lab_event.RData")
#load("lab_event.RData")
#rm(lab_event)
IBD_lab_event <- lab_event %>% filter(hadm_id %in% IBD_indi$hadm_id)
IBD_lab_event1 <- table(IBD_lab_event$itemid) %>% as.data.frame() %>% arrange(desc(Freq))
IBD_lab_event1 <- IBD_lab_event %>% filter(itemid %in% IBD_lab_event1[1:200,1]) %>% arrange(charttime)

IBD_sin_sub_lab_1 <- extrade_n_lab(data=IBD_lab_event1,
                                   n=1,
                                   human_ID = "subject_id",
                                   lab_col="itemid",
                                   lab_result="value")
#53703
IBD_sin_ham_lab_1 <- extrade_n_lab(data=IBD_lab_event1,
                                   n=1,
                                   human_ID = "hadm_id",
                                   lab_col="itemid",
                                   lab_result="value")

#save(IBD_sin_ham_lab_1,IBD_sin_sub_lab_1,file = "All_sub_ham_lab.RData")
load("All_sub_ham_lab.RData")
#加载用药情况
prescriptions <- read.csv("prescriptions.csv")
IBD_prescriptions <- prescriptions %>% filter(hadm_id %in% IBD_indi$hadm_id|subject_id %in% IBD_indi$subject_id)
IBD_prescriptions <- IBD_prescriptions %>% arrange(starttime)
IBD_prescriptions1 <- table(IBD_prescriptions$drug) %>% as.data.frame() %>% arrange(desc(Freq))
IBD_prescriptions1 <- IBD_prescriptions %>% filter(drug %in% IBD_prescriptions1[1:350,1]) 

IBD_sin_ham_pre_list <- lapply(list(1,2,3),function(x)
                              extrade_n_lab(data=IBD_prescriptions1,
                                            n=x,
                                            human_ID = "hadm_id",
                                            lab_col="drug",
                                            lab_result="form_val_disp"))

IBD_sin_sub_pre_list <- lapply(list(1,2,3),function(x)
                              extrade_n_lab(data=IBD_prescriptions1,
                                            n=x,
                                            human_ID = "subject_id",
                                            lab_col="drug",
                                            lab_result="form_val_disp"))

#微生物情况
microbiologyevents <- read.csv("microbiologyevents.csv")
IBD_microbiologyevents <- microbiologyevents %>% filter(hadm_id %in% IBD_indi$hadm_id|subject_id %in% IBD_indi$subject_id)
IBD_microbiologyevents <- IBD_microbiologyevents %>% arrange(chartdate)

#身高，体重，血压等基础情况
omr <- read.csv("omr.csv")
IBD_omr <- omr %>% filter(subject_id %in% IBD_indi$subject_id) %>% arrange(chartdate)
IBD_omr_sin <- extrade_n_lab(data=IBD_omr,
                             n=1,
                             human_ID = "subject_id",
                             lab_col="result_name",
                             lab_result="result_value")

#将labitem的ID转换成名称
#IBD_sin_sub_lab_1 <- IBD_sin_sub_lab_list[[1]] %>% t() %>% as.data.frame()
#IBD_sin_ham_lab_1 <- IBD_sin_ham_lab_list[[1]] %>% t() %>% as.data.frame()

#读取labitem的ID
d_labitems <- read.csv("d_labitems.csv")

brige_itemID <- colnames(IBD_sin_sub_lab_1) %>% as.data.frame() 
colnames(brige_itemID) <- "itemid"
brige_itemID <- left_join(brige_itemID,d_labitems[,c(1,5)])
brige_itemID[1,2] <- brige_itemID[1,1]
colnames(IBD_sin_sub_lab_1) <- brige_itemID$item_all

brige_itemID <- colnames(IBD_sin_ham_lab_1) %>% as.data.frame() 
colnames(brige_itemID) <- "itemid"
brige_itemID <- left_join(brige_itemID,d_labitems[,c(1,5)])
brige_itemID[1,2] <- brige_itemID[1,1]
colnames(IBD_sin_ham_lab_1) <- brige_itemID$item_all

#融合所有的数据
IBD_indi_age_sex_ad <- IBD_indi_age_sex_ad[,c("subject_id", "hadm_id",  "gender", "anchor_age", "anchor_year","anchor_year_group",
                                              "admittime", "dischtime",  "dod",  "marital_status", "race",  "hospital_expire_flag")]

IBD_indi_age_sex_ad$subject_id <- as.character(IBD_indi_age_sex_ad$subject_id)
IBD_indi_age_sex_ad$hadm_id <- as.character(IBD_indi_age_sex_ad$hadm_id)

#统计随访时间和住院时间
IBD_group_follow <- IBD_indi_age_sex_ad 
IBD_group_follow$Statue <- ifelse(IBD_group_follow$dod=="",0,1)
IBD_group_follow <- IBD_group_follow %>% mutate(Inhospital=as.Date(dischtime)-as.Date(admittime))

tem_data1 <- IBD_group_follow %>% filter(Statue == 1) %>% mutate(svrvi_time=as.Date(dod)-as.Date(admittime))
tem_data2 <- IBD_group_follow %>% filter(Statue != 1) %>% mutate(svrvi_time="no")

IBD_follow <- rbind(tem_data1,tem_data2)

save.image("ALL_LAB_SEX_AGE.RData")
save(IBD_lab_event,IBD_prescriptions,IBD_microbiologyevents,IBD_indi_age_sex_ad,file = "Lab_micro_proce_race.RData")

#分析诊断
diag_list <- list("hadm_id","subject_id")

IBD_ra_diag <- lapply(diag_list,function(x)diag_disease_fun(data=IBD_all_ICD,
                                human_ID=x,
                                diag_col="long_title",
                                diag_disease=c("ulcerative colitis|crohn's disease","ulcerative colitis",
                                               "crohn's disease","nicotine","Diabetes","hypertension",
                                               "hyperlipidemia","rheumatoid arthritis")))


#融合诊断，人口资料，化验
#删除ham_id重复的数据
dupl_ham_id <- duplicated(IBD_follow$hadm_id) 
IBD_follow_uni_ham <- IBD_follow[!dupl_ham_id,] %>% left_join(IBD_omr_sin)

IBD_data_ham <- left_join(IBD_ra_diag[[1]],IBD_follow_uni_ham) %>% left_join(IBD_sin_ham_lab_1)
IBD_data_sub <- left_join(IBD_ra_diag[[2]],IBD_follow_uni_ham) %>% left_join(IBD_sin_sub_lab_1)
IBD_data_ham[IBD_data_ham=="___"] <- NA
IBD_data_sub[IBD_data_sub=="___"] <- NA

#清洗数据
#删除NA值过多的项目
IBD_data_ham_all <- IBD_data_ham %>% filter(!is.na(RDW_Blood_Hematology)) %>% filter(!is.na(Albumin_Blood_Chemistry)) 
IBD_data_ham_all <- na_delet_data(IBD_data_ham_all,nrow(IBD_data_ham_all)*0.2)
colnames(IBD_data_ham_all) %>% dput()

na_count_IBD <- na_count(IBD_data_ham_all,2) %>% as.data.frame() %>% arrange(desc(.))
table(IBD_data_ham_all$`ulcerative colitis|crohn's disease`)
table(IBD_data_ham_all$nicotine)

table(IBD_data_ham_all[,c("ulcerative colitis|crohn's disease","nicotine")])
colnames(IBD_data_ham_all) %>% dput()

index_col <- c("ulcerative colitis|crohn's disease", "ulcerative colitis", 
               "crohn's disease", "nicotine", "Diabetes", "hypertension", "hyperlipidemia", 
               "rheumatoid arthritis", "gender", "anchor_age","race", "Statue", 
               "Inhospital", "svrvi_time", "Weight (Lbs)", "BMI (kg/m2)", "Height (Inches)", 
               "Blood Pressure", "Creatinine_Blood_Chemistry", "Hematocrit_Blood_Hematology", 
               "Hemoglobin_Blood_Hematology", "Platelet Count_Blood_Hematology", 
               "White Blood Cells_Blood_Hematology", "MCH_Blood_Hematology", 
               "MCHC_Blood_Hematology", "MCV_Blood_Hematology", "RDW_Blood_Hematology", 
               "Red Blood Cells_Blood_Hematology", "RDW-SD_Blood_Hematology", 
               "Urea Nitrogen_Blood_Chemistry", "Chloride_Blood_Chemistry", 
               "Bicarbonate_Blood_Chemistry", "Magnesium_Blood_Chemistry", "Phosphate_Blood_Chemistry", 
               "Calcium, Total_Blood_Chemistry", "INR(PT)_Blood_Hematology", 
               "PT_Blood_Hematology", "Bilirubin, Total_Blood_Chemistry", "Asparate Aminotransferase (AST)_Blood_Chemistry", 
               "Alanine Aminotransferase (ALT)_Blood_Chemistry", "Albumin_Blood_Chemistry", 
               "PTT_Blood_Hematology", "Alkaline Phosphatase_Blood_Chemistry")
IBD_data_ham_all <- IBD_data_ham_all[,index_col]
colnames(IBD_data_ham_all) <-  c("IBD", "UC", 
                                 "CD", "nicotine", "Diabetes", "hypertension", "hyperlipidemia", 
                                 "RA", "SEX", "Age","race", "Statue", 
                                 "Inhospital", "svrvi_time", "Weight", "BMI", "Height", 
                                 "BP", "Creatinine", "Hematocrit", 
                                 "Hemoglobin", "Pla", 
                                 "WBC", "MCH", 
                                 "MCHC", "MCV", "RDW", 
                                 "RBC", "RDWSD", 
                                 "BUN", "Chloride", 
                                 "Bicarbonate", "Magnesium", "Phosphate", 
                                 "Calcium", "INR", 
                                 "PT", "Bilirubin", "AST", 
                                 "ALT", "ALB", 
                                 "PTT", "ALP")

BP_split <- strsplit2(IBD_data_ham_all$BP,"/") %>% as.data.frame()
colnames(BP_split) <- c("HBP","LBP")
BP_split$LBP <- ifelse(BP_split$LBP=="",NA,BP_split$LBP)
IBD_data_ham_all <- cbind(IBD_data_ham_all,BP_split)
IBD_data_ham_all <- select(IBD_data_ham_all,-BP)

IBD_data_ham_all[,-1:-14] <- apply(IBD_data_ham_all[,-1:-14],2,as.numeric)
IBD_mice <-mice(IBD_data_ham_all,m=5,seed = 3,method = "rf")
IBD_data_ham_all <- complete(IBD_mice,action = 2)

save(IBD_data_ham_all,file = "IBD_data_ham_all.RData")
save.image("All_mice.Rdata")
