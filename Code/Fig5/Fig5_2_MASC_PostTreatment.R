source("~/ws/my.source.R")
options(stringsAsFactors=F)

###################
suppressPackageStartupMessages(library(lme4))



celltype_tmp = "NaiveCD4"
cellstate_tmp = "NaiveCD4-0"

# extract "Flare"-"Treated" sample pairs
sample_meta5=fread_FT("CNA/clinical/Sample_meta_clinical.txt")
Treated_donors=sample_meta5%>%filter(DAgroup=="Treated")%>%pull(donor_id)%>%sort()

# filter metadata
# each celltype
meta_tmp=fread_FT("SNN_clustering/SLE_maindata_4thQC_new_meta_anno.txt")%>%
		     filter(celltype==celltype_tmp)%>%
         filter(donor_id%in%Treated_donors & Status%in% c("Flare","Treated"))%>%
         mutate(pop_cov=ifelse(pop_cov=="European","European","non-European")) # since AA and Hispanic are only in HDA

dim(meta_tmp)%>%print()

# MASC: logistic regression
meta_tmp$object=ifelse(meta_tmp$cellstate==cellstate_tmp,1,0)
table(meta_tmp$object)%>%print()

meta_tmp$Status=factor(meta_tmp$Status,levels=c("Flare","Treated"))
meta_tmp$Sex=factor(meta_tmp$Sex,levels=c("Female","Male"))
meta_tmp$pop_cov=factor(meta_tmp$pop_cov,levels=c("European","non-European"))

# scale continuous variables
meta_tmp$Age_scaled=scale(meta_tmp$Age)

# Covariate: Age, Pop (all female)
T1   = glmer(object~Status+Age_scaled+pop_cov+(1|donor_id)+(1|lib), data=meta_tmp,family=binomial(link="logit"),nAGQ=1,control=glmerControl(optimizer="bobyqa"))
T0   = glmer(object~       Age_scaled+pop_cov+(1|donor_id)+(1|lib), data=meta_tmp,family=binomial(link="logit"),nAGQ=1,control=glmerControl(optimizer="bobyqa"))

wald_tmp=summary(T1)$coefficients%>%as.data.frame()%>%rownames_to_column("Parameter")
colnames(wald_tmp)=c("Parameter","Estimate","se","Z","P")

# use LRT Pvalue
anova_tmp= data.frame(Parameter="anova",Estimate=NA,se=NA,Z=NA,P=anova(T1,T0)[2,8])

res_tmp=rbind(wald_tmp,anova_tmp)%>%mutate(celltype=celltype_tmp)%>%mutate(cellstate=cellstate_tmp)

write.table_FT_2(res_tmp,paste0("MASC/Flare_Treated/cellstate_percelltype_Flare_Treated_",cellstate_tmp,"_agepop_res.txt"))

