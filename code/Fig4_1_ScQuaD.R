source("~/ws/my.source.R")
options(stringsAsFactors=F)

###################
suppressPackageStartupMessages(library(lme4))

## NaiveCD4
celltype_tmp = "NaiveCD4"
job_tmp      = 1

# low-exp filter-passed genes
feature_passed = fread_FT(paste0(celltype_tmp,"_10p_3count_lowexp_filter_gene_wochrYbadRBCPLTgenes.txt"))%>%
                 pull(Gene)
feature_passed = gsub("/","-",feature_passed)

# 200 genes per chunk
start=(200*(job_tmp-1)+1)
if(200*job_tmp>length(feature_passed)){end=length(feature_passed)}else{end=200*job_tmp}

# metadata with target samples
meta_tmp=fread_FT("SLE_maindata_4thQC_new_meta_anno.txt")%>%
          filter(!str_detect(donor_id,pattern="IGTB"))%>%
          filter(celltype==celltype_tmp)%>%
          mutate(DAgroup=ifelse(Status=="Healthy","Healthy",
               ifelse(Status=="Flare","HDA",
               ifelse(Status=="Treated","Others",
               ifelse(Status=="Managed" & is.na(sledaiscore),"Others",
               ifelse(Status=="Managed" & sledaiscore==0,"Inactive",
               ifelse(Status=="Managed" & sledaiscore>=1 & sledaiscore<=4,"LDA",
               ifelse(Status=="Managed" & sledaiscore>=5 & sledaiscore<=9,"MDA",
               ifelse(Status=="Managed" & sledaiscore>=10,"HDA",NA)))))))))%>%
         mutate(pop_cov=ifelse(pop_cov=="European","European","nonEuropean"))%>% # since AA and Hispanic are only in HDA
         filter(DAgroup%in%c("Healthy","Inactive"))

# CNA nbhd correlations
cna_res_tmp=readRDS(paste0(celltype_tmp,"/knn_CLUESsample/InactivevsHC_agesexpopcellnum_res.rds"))
corr_tmp=cna_res_tmp$ncorrs%>%as.data.frame()%>%rownames_to_column("cell_id")
colnames(corr_tmp)[2]="nbhd_corr"


for(iii in start:end){
    gene_tmp =feature_passed[iii]
    paste0(iii,";",gene_tmp)%>%print()
    
    count_tmp=fread_FT(paste0("singlecell_NBmodel/celltype_each_count/",celltype_tmp,"/",gene_tmp,"_count.txt")) # gene count
    data_tmp =inner_join(count_tmp,corr_tmp,by="cell_id")%>%
              inner_join(.,meta_tmp%>%dplyr::select(cell_id,DAgroup,Age,Sex,pop_cov,lib,sample_uuid,nCount_RNA,percent.mt),by="cell_id")

    data_tmp$DAgroup =factor(data_tmp$DAgroup,levels=c("Healthy","Inactive","LDA","MDA","HDA","Others"))
    data_tmp$Sex     =factor(data_tmp$Sex,levels=c("Female","Male"))
    data_tmp$pop_cov =factor(data_tmp$pop_cov,levels=c("European","nonEuropean"))
    
    data_tmp$Age_scaled       =scale(data_tmp$Age)
    data_tmp$lognUMI_scaled   =scale(log(data_tmp$nCount_RNA))
    data_tmp$percent.mt_scaled=scale(data_tmp$percent.mt)
  
    print("Full model start")

    T2   = try(glmer.nb(Count~DAgroup+nbhd_corr+Age_scaled+Sex+pop_cov+lognUMI_scaled+percent.mt_scaled+(1|lib)+(1|sample_uuid), 
                        data=data_tmp,nAGQ=0,control = glmerControl(optimizer = "bobyqa")))
       
    print("Null model 1 start")

    T1   = try(glmer.nb(Count~        nbhd_corr+Age_scaled+Sex+pop_cov+lognUMI_scaled+percent.mt_scaled+(1|lib)+(1|sample_uuid), 
                        data=data_tmp,nAGQ=0,control = glmerControl(optimizer = "bobyqa")))

    print("Null model 2 start")

    T0   = try(glmer.nb(Count~ DAgroup         +Age_scaled+Sex+pop_cov+lognUMI_scaled+percent.mt_scaled+(1|lib)+(1|sample_uuid), 
                        data=data_tmp,nAGQ=0,control = glmerControl(optimizer = "bobyqa")))

    
    if(class(T2)=="try-error"|class(T1)=="try-error"|class(T0)=="try-error"){
      res_tmp = data.frame(parameter=c("(Intercept)","DAgroupInactive","nbhd_corr",
                                       "Age_scaled","SexMale","pop_covnonEuropean",
                                       "lognUMI_scaled","percent.mt_scaled","anova1","anova2"),
                           Full_Estimate=NA,Full_SE=NA,Full_Z=NA,Full_P=NA,
                           gene=gene_tmp)
    }else{

    anova1= anova(T2,T1)[2,8]
    anova2= anova(T2,T0)[2,8]

    res_tmp2=data.frame(summary(T2)$coefficients)%>%rownames_to_column("parameter")
    colnames(res_tmp2)[2:5]=c("Full_Estimate","Full_SE","Full_Z","Full_P")

    anova1_tmp=data.frame(parameter="anova1",Full_Estimate=NA,Full_SE=NA,Full_Z=NA,Full_P=anova1)
    anova2_tmp=data.frame(parameter="anova2",Full_Estimate=NA,Full_SE=NA,Full_Z=NA,Full_P=anova2)

    res_tmp=bind_rows(res_tmp2,anova1_tmp,anova2_tmp)%>%
            mutate(gene=gene_tmp)
    }
    rm(T2);gc()
    rm(T1);gc()
    rm(T0);gc()

    if(iii==start){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}
write.table_FT_2(res_sum,paste0("singlecell_NBmodel/NBmodel_final/",celltype_tmp,"/InactivevsHC/with_percMT/",celltype_tmp,"_chunk_",job_tmp,"_res.txt"))


##############################################################################################################################
##############################################################################################################################

sessionInfo()

R version 4.4.2 (2024-10-31)
Platform: x86_64-conda-linux-gnu
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS/LAPACK: /rshare1/ZETTAI_path_WA_slash_home_KARA/home/imgnaka/miniconda3_new/envs/seurat_v5/lib/libopenblasp-r0.3.28.so;  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=ja_JP.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=ja_JP.UTF-8        LC_COLLATE=ja_JP.UTF-8    
 [5] LC_MONETARY=ja_JP.UTF-8    LC_MESSAGES=ja_JP.UTF-8   
 [7] LC_PAPER=ja_JP.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=ja_JP.UTF-8 LC_IDENTIFICATION=C       

time zone: Asia/Tokyo
tzcode source: system (glibc)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] lme4_1.1-36       Matrix_1.6-5      data.table_1.16.4 lubridate_1.9.4  
 [5] forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2      
 [9] readr_2.1.5       tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.1    
[13] tidyverse_2.0.0  

loaded via a namespace (and not attached):
 [1] gtable_0.3.6     compiler_4.4.2   Rcpp_1.0.14      tidyselect_1.2.1
 [5] splines_4.4.2    scales_1.3.0     boot_1.3-31      lattice_0.22-6  
 [9] R6_2.5.1         generics_0.1.3   rbibutils_2.3    MASS_7.3-64     
[13] nloptr_2.1.1     munsell_0.5.1    minqa_1.2.8      pillar_1.10.1   
[17] tzdb_0.4.0       rlang_1.1.5      stringi_1.8.4    timechange_0.3.0
[21] cli_3.6.3        withr_3.0.2      magrittr_2.0.3   Rdpack_2.6.2    
[25] grid_4.4.2       hms_1.1.3        nlme_3.1-165     lifecycle_1.0.4 
[29] reformulas_0.4.0 vctrs_0.6.5      glue_1.8.0       colorspace_2.1-1
[33] tools_4.4.2      pkgconfig_2.0.3 








