library(tidyverse)
library(dplyr)
library(survival)
library(survminer)

?lung
as_tibble(lung)
lung <- as_tibble(lung)
lung

head(lung)

s <- Surv(lung$time, lung$status)
class(s)
s
head(lung)
?survfit


survfit(s~1)
survfit(Surv(time, status)~1, data=lung)
sfit <- survfit(Surv(time, status)~1, data=lung)
sfit


summary(sfit)
?summary.survfit

# ?summary.survfit
range(lung$time)
seq(0, 1100, 100)

summary(sfit, times=seq(0, 1000, 100))
sfit <- survfit(Surv(time, status)~sex, data=lung)
plot(sfit)


library(survminer)
ggsurvplot(sfit)




ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
           legend.labs=c("Male", "Female"), legend.title="Sex",  
           palette=c("dodgerblue2", "orchid2"), 
           title="Kaplan-Meier Curve for Lung Cancer Survival", 
           risk.table.height=.15)


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("RTCGA.mRNA")

# Load the bioconductor installer. 
# Try http:// if https:// doesn't work.
source("https://bioconductor.org/biocLite.R")

# Install the main RTCGA package
biocLite("RTCGA")

# Install the clinical and mRNA gene expression data packages
biocLite("RTCGA.clinical")
biocLite("RTCGA.mRNA")

library(tidyverse)
library(RTCGA)
infoTCGA()

library(RTCGA.clinical)
?clinical

dim(BRCA.clinical)
names(BRCA.clinical)

# Create the clinical data
clin <- survivalTCGA(BRCA.clinical, OV.clinical, GBM.clinical, 
                     extract.cols="admin.disease_code")
# Show the first few lines
head(clin)

# Tabulate by outcome
xtabs(~admin.disease_code+patient.vital_status, data=clin) %>% addmargins()







pam50_and_surv_only_METABRIC$gse_21653_patient_code[!(pam50_and_surv_only_METABRIC$gse_21653_patient_code %in%  METABRIC_metabolicgroups_survival$METABRIC_patient_code)]



setdiff(METABRIC_metabolicgroups_survival$METABRIC_patient_code, pam50_and_surv_only_METABRIC$gse_21653_patient_code)
METABRIC_metabolicgroups_survival$METABRIC_patient_code[!(METABRIC_metabolicgroups_survival$METABRIC_patient_code %in% pam50_and_surv_only_METABRIC$gse_21653_patient_code)]
METABRIC_metabolicgroups_survival$METABRIC_patient_code[is.na(match(METABRIC_metabolicgroups_survival$METABRIC_patient_code,pam50_and_surv_only_METABRIC$gse_21653_patient_code))]
METABRIC_metabolicgroups_survival$METABRIC_patient_code[!(is.element(METABRIC_metabolicgroups_survival$METABRIC_patient_code,pam50_and_surv_only_METABRIC$gse_21653_patient_code))]












## 1.Select valid cases from the pam50_and_surv_only_gse### based on GSE###_metabolicgroups_survival
METABRIC_surv_info_by_metsubtypes = pam50_and_surv_only_METABRIC[pam50_and_surv_only_METABRIC$gse_21653_patient_code %in% METABRIC_metabolicgroups_survival$METABRIC_patient_code,]
View(METABRIC_surv_info_by_metsubtypes) 









## Sort pam_and_surv_only_gse### as GSE###_metabolicgroups_survival
METABRIC_metabolicgroups_survival_SORTED = METABRIC_metabolicgroup_survival[order(match(METABRIC_metabolicgroup_survival[,1],pam50_and_surv_only_METABRIC[,1])),]
View(METABRIC_metabolicgroups_survival_SORTED)



## Selecting survival info ordered according to metabolic subtypes

### select valid cases from the metabolic subtypes
METABRIC_surv_info_by_metsubtypes = Pam50_and_surv_only_METABRIC[Pam50_and_surv_only_METABRIC$SAMPLE_ID %in% METABRIC_metabolicgoups_survival$patient_code,]
View(METABRIC_surv_info_by_metsubtypes)




METABRIC_readytouse_survinfo = METABRIC_surv_info_by_metsubtypes[order(match(METABRIC_surv_info_by_metsubtypes[,1],METABRIC_metabolicgoups_survival[,1])),]
head(METABRIC_readytouse_survinfo)
View(METABRIC_readytouse_survinfo)

METABRIC_readytouse_survinfo <- cbind(METABRIC_readytouse_survinfo, three_METABRIC_group = METABRIC_metabolicgoups_survival$group)
View(METABRIC_readytouse_survinfo)




# Tabulate by outcome
xtabs(~three_metabric_group+os_months, data=METABRIC_readytouse_survinfo) %>% addmargins()



coxph(Surv(drfs_time_months, drfs_1_event_0_censored)~three_groups, data=global_drfs)
sfit_three_drfs_drfs_global <- survfit(Surv(drfs_time_months, drfs_1_event_0_censored)~three_groups, data=drfs_global)
summary(sfit_os_METABRIC, times=seq(0,365*5,365))

ggsurvplot(sfit_rfs_METABRIC, conf.int=F, pval=TRUE)

ggsurvplot(sfit_three_drfs_drfs_global, conf.int=F, 
           pval=TRUE,
           risk.table = TRUE, 
           risk.table.height=.25,
           palette = "jco", 
           legend.labs = c("grupo 1", "grupo 2","grupo 3"),
           ylab = "Probabilidad de supervivencia",
           xlab = "Tiempo",
           title = "GLOBAL",
           font.title = c(24, "black"))

#####

ggsurvplot(sfit_three_dfs_GSE48390, conf.int=F, 
           pval=TRUE,
           risk.table = TRUE, 
           risk.table.height=.25,
           palette = c("#EFC000FF", "#0073C2FF", "#868686FF"), 
           legend.labs = c("grupo 2", "grupo 1","grupo 3"),
           ylab = "Probabilidad de supervivencia",
           xlab = "Tiempo",
           title = "GSE48390",
           font.title = c(24, "black"))








## A r x c table  Agresti (2002, p. 57) Job Satisfaction
Job <- matrix(c(1,2,1,0, 3,3,6,1, 10,10,14,9, 6,7,12,11), 4, 4,
              dimnames = list(income = c("< 15k", "15-25k", "25-40k", "> 40k"),
                              satisfaction = c("VeryD", "LittleD", "ModerateS", "VeryS")))
fisher.test(Job) # 0.7827
fisher.test(Job, simulate.p.value = TRUE, B = 1e5) # also close to 0.78

##
fish <- matrix(c(7,0,11,3,9, 49, 2,0,2,2, 6, 7,	44,	23,	7), 3, 5,
              dimnames = list(subgrupos = c("SUBGRUPO 1", "SUBGRUPO 2", "SUBGRUPO 3"),
                              pam50 = c("Basal", "HER2", "LUMINAL A", "LUMINAL B", "NORMAL")))
fisher.test(fish) # 0.7827
a = fisher.test(fish, simulate.p.value = TRUE, B = 1e5) # also close to 0.78
a

