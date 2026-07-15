source("~/ws/my.source.R")
options(stringsAsFactors=F)

###########################################################################
args    = commandArgs(trailingOnly = T)
LIST    = as.character(args[1])
task_id = as.numeric_f(args[2])

## e.g., NaiveCD4, disease-state
celltype_tmp  = "NaiveCD4"
compare_tmp  = "InactivevsHC"
cov_tmp      = "agesexpopcellnum"

## Here I utilized rcna codes at "https://github.com/korsunskylab/rcna/tree/main" with some modification
source("code_R/CNA/CNA_utils/cna_utils_stats.R") # from rcna
source(paste0("code_R/CNA/CNA_utils/cna_reg_",cov_tmp,".R")) # association test with covariates


###################
## 1. Prepare sample-level metadata and NAM-PCs
###################
sample_meta=fread_FT("metadata.txt")%>%
            dplyr::select(sample_uuid,Status_SLE,Status_state,Status_activity,Status_LDA,Status_MDA,European_state,European_activity,Asian_state,Asian_activity,AvsE_SLE,AvsE_HC,Age,Sex,pop_cov,Cellnum_per_sample,DAgroup)%>%
            unique()%>%
            mutate(pop_cov=ifelse(pop_cov=="European","European","non-European")) # since AA and Hispanic are only in HDA
sample_meta$Sex    =factor(sample_meta$Sex,levels=c("Female","Male"))
sample_meta$pop_cov=factor(sample_meta$pop_cov,levels=c("European","non-European"))
sample_meta$DAgroup=factor(sample_meta$DAgroup,levels=c("Healthy","Inactive","LDA","MDA","HDA","Others"))

U = fread_n(paste0("NAM_PCA_",compare_tmp,"_sample_loadings_no_covs.csv"))%>%as.matrix()
V = fread_n(paste0("NAM_PCA_",compare_tmp,"_nbhd_loadings_no_covs.csv"))%>%as.matrix()
sv=fread_FT(paste0("NAM_PCA_",compare_tmp,"_svs_no_covs.csv"))$V2[-1]

# extract sample_id, target phenotype and covariates
if(compare_tmp=="InactivevsHC")   {sample_lim=sample_meta%>%filter(!is.na(Status_state))%>%select(1,3,13:17)}
if(compare_tmp=="HDAvsInactive")  {sample_lim=sample_meta%>%filter(!is.na(Status_activity))%>%select(1,4,13:17)}

colnames(sample_lim)[2]="score"
sample_lim$score=scale(sample_lim$score) # scale phenotype
n=nrow(sample_lim) 

M=Matrix::Diagonal(n = n) # no covs before NAM-PC
r=0                       # no covs before NAM-PC
batches_vec = rep(1, n)   # batch=NULL
Nnull=10000
ks=c(1,2,3,4,5)


###################
## 2. non-null test
###################

# get non-null f-test p-value
.tmp = .minp_stats(sample_lim)
k = .tmp$k; p = .tmp$p; r2 = .tmp$r2
if (k == max(ks)) {
    warning(glue::glue('data supported use of {k} NAM PCs, which is the maximum considered. Consider allowing more PCs by using the "ks" argument.'))        
}

# compute coefficients and r2 with chosen model
.tmp = .reg(sample_lim, k)
yhat = .tmp$qhat; beta = .tmp$beta
r2_perpc = (beta / as.numeric(sqrt(crossprod(sample_lim$score))))**2

# get neighborhood scores with chosen model
ncorrs = V[, 1:k] %*% (sqrt(sv[1:k]) * beta/n)    
rownames(ncorrs) = rownames(V)


###################
## 3. null tests from 10,000 times permutation
###################
# compute final permutation p-value using 10,000 null f-test p-values
# Here I permuted samples maintaining the correlation structure between phenotype and covariates
perm=list()

for(iii in 1:Nnull){
   sample_perm=sample_lim
   set.seed(iii)
   sample_perm$sample_uuid=sample(sample_lim$sample_uuid,replace = FALSE)
   perm[[iii]]=sample_perm
}

.tmp = lapply(perm, .minp_stats)
nullminps = purrr::map_dbl(.tmp, 'p')
nullr2s = purrr::map_dbl(.tmp, 'r2')

pfinal = (sum(nullminps <= p+1e-8) + 1)/(Nnull + 1)
if (sum(nullminps <= p+1e-8) == 0) {
    warning('global association p-value attained minimal possible value. Consider increasing Nnull')
}

# get neighborhood fdrs
# 1000 times permutation
message('computing neighborhood-level FDRs')
set.seed(672)
Nnull2 = sample(1:Nnull,1000,replace=F)%>%sort() # limit to 1000 tests (computational cost)
perm  = perm[Nnull2]
reg_  = lapply(perm,.reg,k=k)

for(iii in 1:1000){
    beta_tmp=reg_[[iii]]$beta
    if(iii==1){gamma_=beta_tmp}else{gamma_=cbind(gamma_,beta_tmp)}
}

nullncorrs = abs(V[, 1:k] %*% (sqrt(sv[1:k])*(gamma_ / n)))

maxcorr = max(abs(ncorrs))
fdr_thresholds = seq(maxcorr/4, maxcorr, maxcorr/4000)
fdr_vals = empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)
fdrs = data.frame(
#         threshold = fdr_thresholds
    threshold = head(fdr_thresholds, -1),
    fdr = fdr_vals, 
    num_detected = purrr::map_dbl(head(fdr_thresholds, -1), function(.t) sum(abs(ncorrs) > .t)) 
)
# find maximal FDR<5% and FDR<10% sets
if (min(fdrs$fdr) > 0.05) {        
    fdr_5p_t = NULL
} else {
    fdr_5p_t = min(subset(fdrs, fdr < 0.05)$threshold)        
}
if (min(fdrs$fdr) > 0.1) {        
    fdr_10p_t = NULL
} else {
    fdr_10p_t = min(subset(fdrs, fdr < 0.1)$threshold)
}


## output
res = list(
    p = pfinal, 
    nullminps=nullminps,
    k=k,
    ncorrs=ncorrs, 
    fdrs=fdrs,
    fdr_5p_t=fdr_5p_t, 
    fdr_10p_t=fdr_10p_t,
    yhat=yhat, 
    ks=ks, 
    beta=beta,
    r2=r2, 
    r2_perpc=r2_perpc,
    nullr2_mean=mean(nullr2s), 
    nullr2_std=sd(nullr2s)
)

saveRDS(res,paste0(compare_tmp,"_",cov_tmp,"_res.rds"))


#########################################################################################################################################
#########################################################################################################################################


########################
## Aggregate and visualize
## aggregate all celltypes
########################
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))

suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(ggrepel))

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
source("~/ws/colors.R")


All_list=fread_FT("tmp_job_list/CNA_celltype_list.txt")

lineage_list=c("CD4","CD8","NK","B","Mono","DC")
celltype_list=All_list$celltype
compare_list=c("InactivevsHC","HDAvsInactive")

list=fread_FT("tmp_job_list/modified_cna_celltype_final.txt")%>%
     left_join(.,All_list,by="celltype") # add lineage,mincell,PC,sigma info

list$celltype=factor(list$celltype,levels=celltype_list)
list$compare=factor(list$compare,levels=compare_list)
list=list[order(list$compare),]
list=list[order(list$celltype),]


## summarize cell state stats, median & FDR 5%

## metadata
meta_full=fread_FT("SNN_clustering/SLE_maindata_4thQC_new_meta_anno.txt")

for(kkk in 1:nrow(list)){
    dir_tmp     =list$dir[kkk]%>%as.character()
    lineage_tmp =list$lineage[kkk]%>%as.character()
    celltype_tmp=list$celltype[kkk]%>%as.character()
    compare_tmp =list$compare[kkk]%>%as.character()
    cov_tmp     =list$cov[kkk]%>%as.character()

    paste0(celltype_tmp,"_",compare_tmp,"_",cov_tmp)%>%print()

    cna_res_tmp=readRDS(paste0("CNA/",compare_tmp,"_",cov_tmp,"_res.rds"))

    fdr_5p_tmp =cna_res_tmp$fdr_5p_t
    if(is.null(fdr_5p_tmp)){fdr_5p_tmp=NA}

    corr_max=max(abs(cna_res_tmp$ncorrs))

    corr_tmp   =cna_res_tmp$ncorrs%>%as.data.frame()%>%rownames_to_column("cell_id")
    colnames(corr_tmp)[2]="corr"

    corr_tmp=corr_tmp%>%inner_join(.,meta_full%>%select(cell_id,celltype,cellstate),by="cell_id")
    dim(corr_tmp)%>%print()

    # cell state median and quantiles
    corr_median2=corr_tmp%>%
                      group_by(cellstate)%>%
                      summarise_at(vars(corr),funs(median=median,quantile1=quantile(.,0.25,type=1),quantile3=quantile(.,0.75,type=1)))%>%
                      mutate(fdr_5p_t=fdr_5p_tmp)%>%mutate(corr_max=corr_max)%>%
                      mutate(lineage=lineage_tmp)%>%mutate(celltype=celltype_tmp)%>%
                      mutate(compare=compare_tmp)%>%mutate(cov=cov_tmp)%>%
                      mutate(type=ifelse(median     > fdr_5p_tmp,"Expand2",
                                  ifelse(quantile3  > fdr_5p_tmp,"Expand1",
                                  ifelse(median     < (-fdr_5p_tmp),"Deplete2",
                                  ifelse(quantile1  < (-fdr_5p_tmp),"Deplete1","Not significant")))))

    if(kkk==1){corr_median2_sum=corr_median2}else{corr_median2_sum=rbind(corr_median2_sum,corr_median2)}
}
write.table_FT_2(corr_median2_sum,paste0("CNA/summary/eachcelltype_localtest_summary_cellstate.txt")) 


##############
# UMAP
##############

# change color scale: max 0.6
myPalette = colorRampPalette(c(RColorBrewer::brewer.pal(9, "RdBu")[9:6], "#DCDCDC", RColorBrewer::brewer.pal(9, "RdBu")[4:1]))
limit=0.6

for(kkk in 1:nrow(list)){
   dir_tmp     =list$dir[kkk]%>%as.character()
   lineage_tmp =list$lineage[kkk]%>%as.character()
   celltype_tmp=list$celltype[kkk]%>%as.character()
   compare_tmp =list$compare[kkk]%>%as.character()
   cov_tmp     =list$cov[kkk]%>%as.character()

   mincell_tmp=list$mincell[kkk]
   PC_tmp     =list$PC[kkk]
   sigma_tmp  =list$sigma[kkk]

   paste0(celltype_tmp," ",compare_tmp)%>%print()

   umap_post=fread_FT(paste0("SNN_clustering/",lineage_tmp,"_",celltype_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_postUMAP.txt"))
   dim(umap_post)%>%print()

   cna_res_tmp=readRDS(paste0("CNA/",compare_tmp,"_",cov_tmp,"_res.rds"))
   fdr_10p_tmp=cna_res_tmp$fdr_10p_t
   fdr_5p_tmp =cna_res_tmp$fdr_5p_t
   corr_tmp   =cna_res_tmp$ncorrs%>%as.data.frame()%>%rownames_to_column("cell_id")
   corr_tmp2=corr_tmp%>%left_join(.,umap_post,by="cell_id")
   colnames(corr_tmp2)[2]="corr"
   dim(corr_tmp2)%>%print()

   if(is.null(fdr_10p_tmp)){fdr_10p_tmp=NA}
   if(is.null(fdr_5p_tmp)){fdr_5p_tmp=NA}

   corr_tmp2=corr_tmp2%>%mutate(corr=ifelse(corr>limit,limit,
                                     ifelse(corr<(-limit),-limit,corr)))%>% # upper/lower limit
                         mutate(corr_sig=ifelse(abs(corr)>fdr_5p_tmp,corr,0))
   corr_tmp2$corr_sig[is.na(corr_tmp2$corr_sig)]=0

   p =  ggplot(corr_tmp2,aes(x=UMAP_1,y=UMAP_2,color=corr))+
        geom_point(size=0.2)+
        theme_void()+
        scale_colour_gradientn(na.value="lightgray",colors = myPalette(100), limits = c(-limit,limit)) + 
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title=element_blank(),
              legend.position="none")

    png(paste0("CNA/umap/",lineage_tmp,"_",celltype_tmp,"_",compare_tmp,"_",cov_tmp,"_localtest_corr_all.png"),height=5, width=5, units = "in",res = 300)
     plot(p)
    dev.off()

   p =  ggplot(corr_tmp2,aes(x=UMAP_1,y=UMAP_2,color=corr_sig))+
        geom_point(size=0.2)+
        theme_void()+
        scale_colour_gradientn(na.value="lightgray",colors = myPalette(100), limits = c(-limit,limit)) + 
        theme(axis.text.x=element_blank(),
              axis.text.y=element_blank(),
              axis.title.x=element_blank(),
              axis.title.y=element_blank(),
              plot.title=element_blank(),
              legend.position="none")

    png(paste0("CNA/umap/",lineage_tmp,"_",celltype_tmp,"_",compare_tmp,"_",cov_tmp,"_localtest_corr_onlyfdr005.png"),height=5, width=5, units = "in",res = 300)
     plot(p)
    dev.off()
}


##############
# Violinplot
##############
corr_median2_sum_all=fread_FT("CNA/summary/eachcelltype_localtest_summary_cellstate.txt")

# anno info
meta_full=fread_FT("SNN_clustering/SLE_maindata_4thQC_new_meta_anno.txt")

compare_lim=c("InactivevsHC","HDAvsInactive")

# change color scale: max 0.6
myPalette = colorRampPalette(c(RColorBrewer::brewer.pal(9, "RdBu")[9:6], "#DCDCDC", RColorBrewer::brewer.pal(9, "RdBu")[4:1]))
limit=0.6

dir.create_p("CNA/violin_stateactivity/")

for(kkk in 1:nrow(list3)){
    dir_tmp     =list3$dir[kkk]%>%as.character()
    celltype_tmp=list3$celltype[kkk]%>%as.character()
    compare_tmp =list3$compare[kkk]%>%as.character()
    cov_tmp     =list3$cov[kkk]%>%as.character()
    paste0(celltype_tmp,";",compare_tmp)%>%print()

    cna_res_tmp=readRDS(paste0("CNA/",compare_tmp,"_",cov_tmp,"_res.rds"))

    corr_median2_tmp=corr_median2_sum_all%>%filter(celltype==celltype_tmp & compare==compare_tmp)

    fdr_5p_tmp =cna_res_tmp$fdr_5p_t
    if(is.null(fdr_5p_tmp)){fdr_5p_tmp=corr_median2_tmp$corr_max[1]} # need to be numeric

    corr_tmp   =cna_res_tmp$ncorrs%>%as.data.frame()%>%rownames_to_column("cell_id")%>%
                mutate(celltype=celltype_tmp)%>%mutate(compare=compare_tmp)
    colnames(corr_tmp)[2]="corr"

    corr_tmp2=corr_tmp%>%inner_join(.,meta_full%>%select(cell_id,celltype,cellstate),by=c("cell_id","celltype"))%>%
                         mutate(cellstate_id=take_factor(cellstate,2,"-")%>%as.numeric())
    dim(corr_tmp2)%>%print()

    #corr_tmp2$cellstate=factor(corr_tmp2$cellstate,levels=rev(sort(unique(corr_tmp2$cellstate))))
    corr_tmp2$cellstate_id=factor(corr_tmp2$cellstate_id,levels=rev(sort(unique(corr_tmp2$cellstate_id))))

    p =  ggplot(corr_tmp2,aes(x=cellstate_id,y=corr))+
         coord_flip()+
         geom_hline(yintercept=fdr_5p_tmp, col="darkgrey",linetype="dashed")+
         geom_hline(yintercept=-fdr_5p_tmp,col="darkgrey",linetype="dashed")+
         geom_hline(yintercept=0,col="darkgrey",linetype="dashed")+
         geom_quasirandom(aes(color=corr),width=0.3,size=0.2) +
         stat_summary(aes(group=cellstate_id),fun=function(x){quantile(x,0.25,type=1)},fun.min=function(x){quantile(x,0.25,type=1)},fun.max=function(x){quantile(x,0.25,type=1)},geom="crossbar",color="darkgrey",width=0.4,lwd=0.3,linetype="dashed")+
         stat_summary(aes(group=cellstate_id),fun=function(x){quantile(x,0.75,type=1)},fun.min=function(x){quantile(x,0.75,type=1)},fun.max=function(x){quantile(x,0.75,type=1)},geom="crossbar",color="darkgrey",width=0.4,lwd=0.3,linetype="dashed")+
         stat_summary(aes(group=cellstate_id),fun=median,fun.min=median,fun.max=median,geom="crossbar",color="grey10",width=0.5,lwd=0.4,linetype="dashed")+
         theme_classic()+
         scale_colour_gradientn(na.value="lightgray",colors=myPalette(100),limits=c(-0.6,0.6)) + 
         theme(axis.text.x=element_blank(),
               axis.text.y=element_text(colour="black",size=20),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(),
               plot.title=element_blank(),
               legend.position="none")+
         scale_y_continuous(limits=c(-0.6,0.6),breaks=seq(-0.4,0.4,by=0.2))

    nrow=length(unique(corr_tmp2$cellstate))

    png(paste0("CNA/violin_stateactivity/",celltype_tmp,"_",compare_tmp,"_",cov_tmp,"_localtest_corr_violin_all_cellstate.png"),height=0.3*nrow, width=4, units = "in",res = 300)
     plot(p)
    dev.off()
}


#########################################################################################################################################
#########################################################################################################################################

sessionInfo()

R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS:   /usr/local/package/r/4.4.0/lib64/R/lib/libRblas.so 
LAPACK: /usr/local/package/r/4.4.0/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

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
[1] grid      stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] ggsci_3.2.0           RColorBrewer_1.1-3    circlize_0.4.16      
 [4] ComplexHeatmap_2.22.0 ggrepel_0.9.6         ggbeeswarm_0.7.2     
 [7] cowplot_1.1.3         gridExtra_2.3         data.table_1.16.4    
[10] lubridate_1.9.4       forcats_1.0.0         stringr_1.5.1        
[13] dplyr_1.1.4           purrr_1.0.2           readr_2.1.5          
[16] tidyr_1.3.1           tibble_3.2.1          ggplot2_3.5.2        
[19] tidyverse_2.0.0      

loaded via a namespace (and not attached):
 [1] generics_0.1.3      shape_1.4.6.1       stringi_1.8.4      
 [4] hms_1.1.3           digest_0.6.37       magrittr_2.0.3     
 [7] timechange_0.3.0    iterators_1.0.14    foreach_1.5.2      
[10] doParallel_1.0.17   GlobalOptions_0.1.2 scales_1.3.0       
[13] codetools_0.2-20    cli_3.6.3           rlang_1.1.5        
[16] crayon_1.5.3        munsell_0.5.1       withr_3.0.2        
[19] tools_4.4.0         parallel_4.4.0      tzdb_0.4.0         
[22] colorspace_2.1-1    BiocGenerics_0.52.0 GetoptLong_1.0.5   
[25] vctrs_0.6.5         R6_2.5.1            png_0.1-8          
[28] stats4_4.4.0        matrixStats_1.5.0   lifecycle_1.0.4    
[31] S4Vectors_0.44.0    IRanges_2.40.1      clue_0.3-66        
[34] vipor_0.4.7         cluster_2.1.6       pkgconfig_2.0.3    
[37] beeswarm_0.4.0      pillar_1.10.1       gtable_0.3.6       
[40] glue_1.8.0          Rcpp_1.0.14         tidyselect_1.2.1   
[43] rjson_0.2.23        compiler_4.4.0 







