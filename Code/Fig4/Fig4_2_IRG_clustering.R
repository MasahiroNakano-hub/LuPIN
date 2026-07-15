################################################################
## 1. IRGs clustering
################################################################

R
source("~/ws/my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(viridis))
source("~/ws/colors.R")

All_list=fread_FT("tmp_job_list/CNA_celltype_list.txt")
celltype_list=All_list$celltype


##############
# 1-2-0. Pascual 100 genes
##############
ISG_geneset_1=fread_FT("geneset_pathways/251026_ISG_pascual_sc2020_v4.txt")%>%pull(genes)
# modified for final figure
# 100 genes 
# exclude CCL3L1 and CCL4L2 from v3

# ScQuaD model sumstats
res_sum =fread_FT("singlecell_NBmodel/NBmodel_final/summary/eachcelltype_eachcompare_with_percMT_res_sum_allFDR_final.txt")
intersect(ISG_geneset_1,unique(res_sum$gene))%>%length()
# 100, all genes included in low-exp filter in at least one celltype


##############
# 8-2-1. disease-state
##############
## clustering by Z score
res_ISG_state=res_ISG%>%filter(compare=="InactivevsHC")%>%
                        dplyr::select(celltype,gene,dys_Z,nbhd_Z)%>%
                        pivot_longer(col=-c(celltype,gene),names_to="sig",values_to="Zscore")%>%
                        mutate(type=paste0(celltype,"_",sig))%>%
                        dplyr::select(type,gene,Zscore)%>%
                        pivot_wider(names_from="type",values_from="Zscore")%>%
                        column_to_rownames("gene")

colnames(res_ISG_state)=gsub("dys_Z","dysregulated",colnames(res_ISG_state))
colnames(res_ISG_state)=gsub("nbhd_Z","abundance",colnames(res_ISG_state))
res_ISG_state[is.na(res_ISG_state)]=0

# scale signatures to make them comparable 
res_ISG_state_scaled=scale(res_ISG_state)

# column order: dysregulated vs abundance
labels_col = colnames(res_ISG_state_scaled)
labels_col = labels_col[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,
                          1,3,5,7,9, 11,13,15,17,19,21,23,25,27,29,31,33,35
                          )]

res_ISG_state_reordered=res_ISG_state_scaled[,labels_col]

# column annotations: no signature annotation
celltype_sig_list=data.frame(celltype_sig=colnames(res_ISG_state_reordered))%>%
                  mutate(celltype=take_factor(celltype_sig,1,"_"))%>%
                  mutate(sig=take_factor(celltype_sig,2,"_"))%>%
                  left_join(.,All_list[,c(1,2)],by="celltype")
celltype_sig_list$lineage=factor(celltype_sig_list$lineage,levels=c("CD4","CD8","NK","B","Mono","DC"))
celltype_sig_list$sig    =factor(celltype_sig_list$sig,levels=c("dysregulated","abundance"))

lineage_col=c("CD4"="#cc0010","CD8"="#8b7355","NK"="#e0d51a","B"="#c4f20b","Mono"="#1d2088","DC"="#a4fbfa")
sig_col=c("dysregulated"="#E64B35FF","abundance"="#4DBBD5FF")
ha = HeatmapAnnotation(lineage=celltype_sig_list$lineage,
                       col=list(lineage=lineage_col),show_annotation_name=F,
                       height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=T,
                       annotation_legend_param=list(title_gp=gpar(fontsize=10),labels_gp=gpar(fontsize=10)))

# colnames: only cell types
res_ISG_state_reordered2=res_ISG_state_reordered
colnames(res_ISG_state_reordered2)=take_factor(colnames(res_ISG_state_reordered),1,"_")

# color: Pascual
col_fun = colorRamp2(c(-1.5,0,1.5), c("cyan","black","yellow"))

p=Heatmap(res_ISG_state_reordered2,
                clustering_distance_rows="euclidean",clustering_method_rows="ward.D2",
                cluster_columns=FALSE,
                col=col_fun,height=unit(18,"cm"),width=unit(10,"cm"),
                row_names_gp=gpar(fontsize=6),row_names_side=c("left"),
                column_names_gp=gpar(fontsize=10),
                bottom_annotation = ha,
                heatmap_legend_param=list(title="Z score (scaled)",at=c(-1.5,0,1.5),
                                          direction="vertical",title_position="topcenter",
                                          grid_height=unit(3,"cm"),legend_width=unit(0.5,"cm"),
                                          title_gp=gpar(fontsize=10),labels_gp=gpar(fontsize=10)))
  
pdf_3("singlecell_NBmodel/ISGs/ISG_pascual_InactivevsHC_Zscore_scaled_Pascual_oriorder_genenames.pdf",w=7,h=11)
 draw(p,heatmap_legend_side="bottom",annotation_legend_side="right")
dev.off()


## reorder genes maintaining cluster structures
## no clustering for signatures
# cluster info based on scaled values
hc_row = hclust(dist(res_ISG_state_scaled), method="ward.D2")

pdf_3(paste0("singlecell_NBmodel/ISGs/ISG_pascual_InactivevsHC_Zscore_scaled_hclust_gene.pdf"),h=5,w=15)
 plot(hc_row)
dev.off()

labels_row = hc_row$labels[hc_row$order]
labels_row = labels_row[c(89:100,87:88,
                          47:64,
                          65:67,68,70:71,69,
                          78:86,72:77,
                          1:46
                          )]

# column order: dysregulated vs abundance
labels_col = colnames(res_ISG_state_scaled)
labels_col = labels_col[c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,
                          1,3,5,7,9, 11,13,15,17,19,21,23,25,27,29,31,33,35
                          )]

res_ISG_state_reordered=res_ISG_state_scaled[labels_row,labels_col]


## No rownames, resize ver
# column annotations: no signature annotation
celltype_sig_list=data.frame(celltype_sig=colnames(res_ISG_state_reordered))%>%
                  mutate(celltype=take_factor(celltype_sig,1,"_"))%>%
                  mutate(sig=take_factor(celltype_sig,2,"_"))%>%
                  left_join(.,All_list[,c(1,2)],by="celltype")
celltype_sig_list$lineage=factor(celltype_sig_list$lineage,levels=c("CD4","CD8","NK","B","Mono","DC"))
celltype_sig_list$sig    =factor(celltype_sig_list$sig,levels=c("dysregulated","abundance"))

lineage_col=c("CD4"="#cc0010","CD8"="#8b7355","NK"="#e0d51a","B"="#c4f20b","Mono"="#1d2088","DC"="#a4fbfa")
sig_col=c("dysregulated"="#E64B35FF","abundance"="#4DBBD5FF")
ha = HeatmapAnnotation(lineage=celltype_sig_list$lineage,
                       col=list(lineage=lineage_col),show_annotation_name=F,
                       height=unit(0.2,"cm"),simple_anno_size_adjust=T,show_legend=T,
                       annotation_legend_param=list(title_gp=gpar(fontsize=10),labels_gp=gpar(fontsize=10)))

# colnames: only cell types
res_ISG_state_reordered2=res_ISG_state_reordered
colnames(res_ISG_state_reordered2)=take_factor(colnames(res_ISG_state_reordered),1,"_")%>%
                                   gsub("CD56brightNK","CD56++NK",.)%>%
                                   gsub("CD16nMono","CD16-Mono",.)%>%
                                   gsub("CD16pMono","CD16+Mono",.)

col_fun = colorRamp2(c(-1.5,0,1.5), c("cyan","black","yellow"))

p=Heatmap(res_ISG_state_reordered2,
                cluster_rows=FALSE,
                cluster_columns=FALSE,
                col=col_fun,height=unit(8,"cm"),width=unit(11,"cm"),
                show_row_names=F,
                column_names_gp=gpar(fontsize=10),
                bottom_annotation = ha,
                heatmap_legend_param=list(title="Z score (scaled)",at=c(-1.5,0,1.5),
                                          direction="vertical",title_position="topcenter",
                                          grid_height=unit(3,"cm"),legend_width=unit(0.5,"cm"),
                                          title_gp=gpar(fontsize=10),labels_gp=gpar(fontsize=10)))
  
pdf_3("singlecell_NBmodel/ISGs/ISG_pascual_InactivevsHC_Zscore_scaled_Pascual_reordered_nogenenames.pdf",h=7,w=7)
 draw(p,heatmap_legend_side="bottom",annotation_legend_side="right")
dev.off()


##############################################################################################################################
##############################################################################################################################

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
 [1] viridis_0.6.5         viridisLite_0.4.2     ggsci_3.2.0          
 [4] RColorBrewer_1.1-3    circlize_0.4.16       ComplexHeatmap_2.22.0
 [7] data.table_1.16.4     lubridate_1.9.4       forcats_1.0.0        
[10] stringr_1.5.1         dplyr_1.1.4           purrr_1.0.2          
[13] readr_2.1.5           tidyr_1.3.1           tibble_3.2.1         
[16] ggplot2_3.5.2         tidyverse_2.0.0      

loaded via a namespace (and not attached):
 [1] generics_0.1.3      shape_1.4.6.1       stringi_1.8.4      
 [4] hms_1.1.3           digest_0.6.37       magrittr_2.0.3     
 [7] timechange_0.3.0    iterators_1.0.14    foreach_1.5.2      
[10] doParallel_1.0.17   GlobalOptions_0.1.2 gridExtra_2.3      
[13] scales_1.3.0        codetools_0.2-20    cli_3.6.3          
[16] rlang_1.1.5         crayon_1.5.3        munsell_0.5.1      
[19] withr_3.0.2         tools_4.4.0         parallel_4.4.0     
[22] tzdb_0.4.0          colorspace_2.1-1    GetoptLong_1.0.5   
[25] BiocGenerics_0.52.0 vctrs_0.6.5         R6_2.5.1           
[28] png_0.1-8           matrixStats_1.5.0   stats4_4.4.0       
[31] lifecycle_1.0.4     S4Vectors_0.44.0    IRanges_2.40.1     
[34] clue_0.3-66         cluster_2.1.6       pkgconfig_2.0.3    
[37] pillar_1.10.1       gtable_0.3.6        glue_1.8.0         
[40] tidyselect_1.2.1    rjson_0.2.23        compiler_4.4.0     

















