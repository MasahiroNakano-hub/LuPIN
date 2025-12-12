source("~/ws/my.source.R")
options(stringsAsFactors=F)

###################
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(presto))
suppressPackageStartupMessages(library(singlecellmethods))
source("~/ws/Immunogenomics.R")


lineage_tmp="CD8"
celltype_tmp="EMCD8"

# celltype-level refmap cell-state annotation
anno_tmp=fread_FT("celltype_refmap_knn_predict_sum.txt")%>%
           filter(celltype==celltype_tmp)

seurat=readRDS(paste0(lineage_tmp,"_4thQC_seurat.rds"))%>%
       subset(cells=anno_tmp$cell_id)

DefaultAssay(seurat)="ADT"
seurat=NormalizeData(seurat,normalization.method="CLR",margin=2)

# limit to cell states with >5 cells
sig_clusters=table(anno_tmp$cell_type_pred_knn)[table(anno_tmp$cell_type_pred_knn)>5]%>%names()%>%sort()

for(iii in 1:length(sig_clusters)){
   cellstate_tmp=sig_clusters[iii]

   cells1=anno_tmp%>%filter(cell_type_pred_knn==cellstate_tmp)%>%pull(cell_id)
   cells2=anno_tmp%>%filter(!cell_type_pred_knn==cellstate_tmp)%>%pull(cell_id)
   
   # report all markers
   # also infer negative markers
   marker_tmp=FindMarkers(seurat$ADT,cells.1=cells1,cells.2=cells2,test.use="wilcox",min.pct=0,logfc.threshold=0)%>%
              mutate(cellnum.1=length(cells1))%>%mutate(cellnum.2=length(cells2))%>%
              mutate(Prot=rownames(.))%>%mutate(cellstate=cellstate_tmp)

   if(iii==1){marker_sum=marker_tmp}else{marker_sum=rbind(marker_sum,marker_tmp)}           
}
write.table_FT_2(marker_sum,paste0(celltype_tmp,"_refmap_SEPs.txt"))



########################
### EMCD8-5 figures
########################

suppressPackageStartupMessages(library(ggrepel))

cellstate_tmp="EMCD8-5"


res_tmp=fread_FT(paste0(celltype_tmp,"_refmap_SEPs_v5.txt"))%>%
        filter(cellstate==cellstate_tmp)

## volcano plot
marker_tmp2=res_tmp%>%
            mutate(key=ifelse(Prot%in%c("HLA-DR-ADT","CD27-ADT","CD38-ADT","CD28-ADT"),"Yes","No"))%>%
            mutate(label=ifelse(abs(avg_log2FC)>0.4 & p_val<1e-50,gsub("-ADT","",Prot),""))
marker_tmp2$key=factor(marker_tmp2$key,levels=c("Yes","No"))

max=max(abs(marker_tmp2$avg_log2FC))
# 1.595288
marker_tmp2$p_val[marker_tmp2$p_val<1e-200]=1e-200 # upper limit

# use avg_log2FC (Seurat)
p =   ggplot(marker_tmp2,aes(x=avg_log2FC,y=-log10(p_val),color=key,label=label))+
      geom_point(size=1)+
      geom_text_repel(size=3,max.overlaps=Inf)+
      theme_classic()+
      scale_color_manual(values=c("red","black")) +
      theme(axis.text.x=element_text(colour="black",size=10),
            axis.text.y=element_text(colour="black",size=10),
            axis.title.x=element_text(colour="black",size=12),
            axis.title.y=element_text(colour="black",size=12),
            plot.title=element_blank(),
            legend.position="none")+
      scale_x_continuous(limits=c(-1.6,1.6))+
      scale_y_continuous(limits=c(0,210))+
      labs(x="logFC",y="-log10(P)")

pdf_3(paste0("summary/",celltype_tmp,"_refmap_",cellstate_tmp,"_SEPs_volcano_logFCSeurat.pdf"),h=2.5,w=4)
     plot(p)
dev.off()


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
 [1] ggrepel_0.9.6           singlecellmethods_0.1.0 presto_1.0.0           
 [4] Rcpp_1.0.14             Seurat_5.1.0            SeuratObject_5.0.2     
 [7] sp_2.1-4                data.table_1.16.4       lubridate_1.9.4        
[10] forcats_1.0.0           stringr_1.5.1           dplyr_1.1.4            
[13] purrr_1.0.2             readr_2.1.5             tidyr_1.3.1            
[16] tibble_3.2.1            ggplot2_3.5.1           tidyverse_2.0.0        

loaded via a namespace (and not attached):
 [1] deldir_2.0-4           pbapply_1.7-2          gridExtra_2.3         
 [4] rlang_1.1.5            magrittr_2.0.3         RcppAnnoy_0.0.22      
 [7] spatstat.geom_3.3-4    matrixStats_1.5.0      ggridges_0.5.6        
[10] compiler_4.4.2         png_0.1-8              vctrs_0.6.5           
[13] reshape2_1.4.4         pkgconfig_2.0.3        fastmap_1.2.0         
[16] promises_1.3.2         tzdb_0.4.0             jsonlite_1.8.9        
[19] goftest_1.2-3          later_1.4.1            spatstat.utils_3.1-2  
[22] irlba_2.3.5.1          parallel_4.4.2         cluster_2.1.8         
[25] R6_2.5.1               ica_1.0-3              stringi_1.8.4         
[28] RColorBrewer_1.1-3     spatstat.data_3.1-4    reticulate_1.40.0     
[31] spatstat.univar_3.1-1  parallelly_1.41.0      lmtest_0.9-40         
[34] scattermore_1.2        tensor_1.5             future.apply_1.11.3   
[37] zoo_1.8-12             sctransform_0.4.1      httpuv_1.6.15         
[40] Matrix_1.6-5           splines_4.4.2          igraph_2.0.3          
[43] timechange_0.3.0       tidyselect_1.2.1       abind_1.4-5           
[46] spatstat.random_3.3-2  codetools_0.2-20       miniUI_0.1.1.1        
[49] spatstat.explore_3.3-4 listenv_0.9.1          lattice_0.22-6        
[52] plyr_1.8.9             shiny_1.10.0           withr_3.0.2           
[55] ROCR_1.0-11            Rtsne_0.17             future_1.34.0         
[58] fastDummies_1.7.5      survival_3.8-3         polyclip_1.10-7       
[61] fitdistrplus_1.2-2     pillar_1.10.1          KernSmooth_2.23-26    
[64] plotly_4.10.4          generics_0.1.3         RcppHNSW_0.6.0        
[67] hms_1.1.3              munsell_0.5.1          scales_1.3.0          
[70] globals_0.16.3         xtable_1.8-4           glue_1.8.0            
[73] lazyeval_0.2.2         tools_4.4.2            RSpectra_0.16-2       
[76] RANN_2.6.2             leiden_0.4.3.1         dotCall64_1.2         
[79] cowplot_1.1.3          grid_4.4.2             colorspace_2.1-1      
[82] nlme_3.1-165           patchwork_1.3.0        cli_3.6.3             
[85] spatstat.sparse_3.1-0  spam_2.11-1            viridisLite_0.4.2     
[88] uwot_0.2.2             gtable_0.3.6           digest_0.6.37         
[91] progressr_0.15.1       htmlwidgets_1.6.4      farver_2.1.2          
[94] htmltools_0.5.8.1      lifecycle_1.0.4        httr_1.4.7            
[97] mime_0.12              MASS_7.3-64 



