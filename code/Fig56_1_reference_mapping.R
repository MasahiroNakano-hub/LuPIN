source("~/ws/my.source.R")
options(stringsAsFactors=F)

###################
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(sctransform))
suppressPackageStartupMessages(library(glmGamPoi))
suppressPackageStartupMessages(library(singlecellmethods))
suppressPackageStartupMessages(library(harmony))
suppressPackageStartupMessages(library(uwot))
source("~/ws/Immunogenomics.R")
suppressPackageStartupMessages(library(symphony))

suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))

suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
source("~/ws/colors.R")


lineage_tmp = "CD4"
celltype_tmp= "NaiveCD4"
mincell_tmp = 10
PC_tmp=20
sigma_tmp=0.25
prune_tmp=0.05


# maindata symphony reference
reference=readRDS(paste0(lineage_tmp,"/",celltype_tmp,"/",lineage_tmp,"_",celltype_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_synphony_reference.rds"))
gene_used=reference$vargenes$symbol

#############################

# 3rd data
# lineage-level refmap cell-type annotation
anno_tmp=fread_FT("SNN_clustering/eachlineage_refmap/lineage_refmap_knn_predict_sum.txt")%>%
         filter(cell_type_pred_knn==celltype_tmp)

# metadata
meta_tmp=fread_FT("SLE_CITEseq_4thQC_final_meta_lineage_anno_full.txt")%>%
         filter(cell_id%in%anno_tmp$cell_id)

# countdata
seurat=readRDS(paste0(lineage_tmp,"_4thQC_seurat.rds"))%>%
       subset(cells=meta_tmp$cell_id)

# Split the dataset into batches
DefaultAssay(seurat) = "RNA"
seurat@meta.data$batch=factor(seurat@meta.data$batch,levels=unique(seurat@meta.data$batch))
seurat.list=SplitObject(seurat, split.by = "batch")

# sctransform for each batch, min_cells=3
# SCT Pearson residuals are not centered nor scaled here
seurat.list=lapply(seurat.list,function(x){SCTransform(x, vst.flavor = "v1", method = "glmGamPoi", do.scale = FALSE, do.center = FALSE, verbose = FALSE, min_cells = 3, return.only.var.genes = FALSE)})

# gene intersect between 2nd and main datasets
gene_kept=intersect(rownames(seurat.list$Batch1$SCT@scale.data),rownames(seurat.list$Batch2$SCT@scale.data))%>%
          intersect(.,gene_used)
gene_kept_df=data.frame(Gene=gene_kept)


# merge each SCT data with reference vargenes
seurat.list2=PrepSCTIntegration(object.list = seurat.list, anchor.features = gene_kept)
seurat_merged2=merge(seurat.list2[[1]], y=seurat.list2[2:length(seurat.list2)],merge.data = TRUE)
dim(seurat_merged2$SCT@scale.data)%>%print()

# metadata: confirm the consistency of cell_id order
meta_tmp=meta_tmp%>%column_to_rownames("cell_id")
meta_tmp=meta_tmp[colnames(seurat_merged2$SCT@scale.data),]
all.equal(colnames(seurat_merged2$SCT@scale.data),rownames(meta_tmp))%>%print() #[1] TRUE

# query mapping
set.seed(42)
query = mapQuery(
    seurat_merged2$SCT@scale.data ,
    meta_tmp, 
    reference,
    vars = c("lib","Donor_ID"), 
    do_normalize = FALSE, # already SCT-normalized
    sigma = sigma_tmp
)
harmony_query=t(query$Z)
write.table_n_2(as.data.frame(harmony_query),"cell_id",paste0(celltype_tmp,"_query_harmony.txt"))


# knn predict
query = knnPredict(query, reference, reference$meta_data$cellstate, k = 5) # label: cellstate
knn_predict_res=query$meta_data%>%dplyr::select(cell_type_pred_knn,cell_type_pred_knn_prob)
write.table_n_2(knn_predict_res,"cell_id",paste0(celltype_tmp,"_knn_predict_res.txt"))


# UMAP projection onto maidata UMAP model
uwot_model_path=paste0(lineage_tmp,"/",celltype_tmp,"/",lineage_tmp,"_",celltype_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_postUMAP_model")
ref_umap_model = uwot::load_uwot(uwot_model_path, verbose = FALSE)
umap_query = uwot::umap_transform(harmony_query, ref_umap_model)
colnames(umap_query) = paste0("UMAP_",1:2)
write.table_n_2(as.data.frame(umap_query),"cell_id",paste0(celltype_tmp,"_query_UMAP.txt"))


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
 [1] ggsci_3.2.0             RColorBrewer_1.1-3      cowplot_1.1.3          
 [4] gridExtra_2.3           symphony_0.1.1          uwot_0.2.2             
 [7] Matrix_1.6-5            harmony_1.2.3           Rcpp_1.0.14            
[10] singlecellmethods_0.1.0 glmGamPoi_1.18.0        sctransform_0.4.1      
[13] Seurat_5.1.0            SeuratObject_5.0.2      sp_2.1-4               
[16] data.table_1.16.4       lubridate_1.9.4         forcats_1.0.0          
[19] stringr_1.5.1           dplyr_1.1.4             purrr_1.0.2            
[22] readr_2.1.5             tidyr_1.3.1             tibble_3.2.1           
[25] ggplot2_3.5.1           tidyverse_2.0.0        

loaded via a namespace (and not attached):
  [1] jsonlite_1.8.9              magrittr_2.0.3             
  [3] spatstat.utils_3.1-2        farver_2.1.2               
  [5] zlibbioc_1.52.0             vctrs_0.6.5                
  [7] ROCR_1.0-11                 spatstat.explore_3.3-4     
  [9] S4Arrays_1.6.0              htmltools_0.5.8.1          
 [11] SparseArray_1.6.0           parallelly_1.41.0          
 [13] KernSmooth_2.23-26          htmlwidgets_1.6.4          
 [15] ica_1.0-3                   plyr_1.8.9                 
 [17] plotly_4.10.4               zoo_1.8-12                 
 [19] igraph_2.0.3                mime_0.12                  
 [21] lifecycle_1.0.4             pkgconfig_2.0.3            
 [23] R6_2.5.1                    fastmap_1.2.0              
 [25] GenomeInfoDbData_1.2.13     MatrixGenerics_1.18.0      
 [27] fitdistrplus_1.2-2          future_1.34.0              
 [29] shiny_1.10.0                digest_0.6.37              
 [31] colorspace_2.1-1            patchwork_1.3.0            
 [33] S4Vectors_0.44.0            tensor_1.5                 
 [35] RSpectra_0.16-2             irlba_2.3.5.1              
 [37] GenomicRanges_1.58.0        progressr_0.15.1           
 [39] spatstat.sparse_3.1-0       timechange_0.3.0           
 [41] httr_1.4.7                  polyclip_1.10-7            
 [43] abind_1.4-5                 compiler_4.4.2             
 [45] withr_3.0.2                 fastDummies_1.7.5          
 [47] MASS_7.3-64                 DelayedArray_0.32.0        
 [49] tools_4.4.2                 lmtest_0.9-40              
 [51] httpuv_1.6.15               future.apply_1.11.3        
 [53] goftest_1.2-3               glue_1.8.0                 
 [55] nlme_3.1-165                promises_1.3.2             
 [57] grid_4.4.2                  Rtsne_0.17                 
 [59] cluster_2.1.8               reshape2_1.4.4             
 [61] generics_0.1.3              gtable_0.3.6               
 [63] spatstat.data_3.1-4         tzdb_0.4.0                 
 [65] hms_1.1.3                   XVector_0.46.0             
 [67] BiocGenerics_0.52.0         spatstat.geom_3.3-4        
 [69] RcppAnnoy_0.0.22            ggrepel_0.9.6              
 [71] RANN_2.6.2                  pillar_1.10.1              
 [73] spam_2.11-1                 RcppHNSW_0.6.0             
 [75] later_1.4.1                 splines_4.4.2              
 [77] lattice_0.22-6              survival_3.8-3             
 [79] deldir_2.0-4                tidyselect_1.2.1           
 [81] miniUI_0.1.1.1              pbapply_1.7-2              
 [83] IRanges_2.40.0              SummarizedExperiment_1.36.0
 [85] scattermore_1.2             stats4_4.4.2               
 [87] Biobase_2.66.0              matrixStats_1.5.0          
 [89] stringi_1.8.4               UCSC.utils_1.2.0           
 [91] lazyeval_0.2.2              codetools_0.2-20           
 [93] cli_3.6.3                   xtable_1.8-4               
 [95] reticulate_1.40.0           munsell_0.5.1              
 [97] GenomeInfoDb_1.42.0         globals_0.16.3             
 [99] spatstat.random_3.3-2       png_0.1-8                  
[101] spatstat.univar_3.1-1       parallel_4.4.2             
[103] dotCall64_1.2               listenv_0.9.1              
[105] viridisLite_0.4.2           scales_1.3.0               
[107] ggridges_0.5.6              crayon_1.5.3               
[109] leiden_0.4.3.1              rlang_1.1.5 













