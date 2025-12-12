source("~/ws/my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(monocle3))

## CD4 monocle3
lineage_tmp="CD4"
mincell_tmp=15;PC_tmp=30;sigma_tmp=0.2;res_tmp=0.15

meta_tmp=fread_FT("SNN_clustering/SLE_maindata_4thQC_new_meta_anno.txt")%>%
         filter(lineage==lineage_tmp)
# 410340

# countdata
seurat=readRDS(paste0("~/ws/2025/250206_SLE_maindata_new_QC/RNA_QC/4thQC/summary/",lineage_tmp,"_4thQC_seurat.rds"))
# 410340
all.equal(meta_tmp$cell_id,seurat@assays[["RNA"]]@counts@Dimnames[[2]])
# TRUE

pca_res=readRDS(paste0("SNN_clustering/eachlineage/",lineage_tmp,"/mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"/",lineage_tmp,"_mincell",mincell_tmp,"_PCA_PC30.rds"))
rownames(pca_res$rotation)=names(pca_res$center)

# umap: use seed0
# flip UMAP1and2 positions for final figure
umap_post=fread_n(paste0("SNN_clustering/eachlineage/",lineage_tmp,"/mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"/",lineage_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_seed0_postUMAP.txt"))
umap_post$UMAP_1=-umap_post$UMAP_1; umap_post$UMAP_2=-umap_post$UMAP_2
umap_post=umap_post[meta_tmp$cell_id,]
all.equal(meta_tmp$cell_id,rownames(umap_post))
# TRUE

### Building the necessary parts for a basic cds

# part one, gene annotations
gene_annotation <- as.data.frame(names(pca_res$center), row.names = names(pca_res$center))
colnames(gene_annotation) <- "gene_short_name"


# part two, cell information
cell_metadata <- as.data.frame(seurat@assays[["RNA"]]@counts@Dimnames[[2]], row.names = seurat@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


# part three, counts sparse matrix
New_matrix <- seurat@assays[["RNA"]]@counts
New_matrix <- New_matrix[names(pca_res$center), ]
expression_matrix <- New_matrix
# 2895 410340


### Construct the basic cds object
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)

### Construct and assign the made up partition
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


### Assign the cluster info
list_cluster <- meta_tmp$cellstate
names(list_cluster) <- meta_tmp$cell_id

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster


### Could be a space-holder, but essentially fills out louvain parameters
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


### Assign UMAP coordinate
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <- umap_post


### Learn graph, this step usually takes a significant period of time for larger samples
print("Learning graph, which can take a while depends on the sample")

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = T)

save_monocle_objects(cds_from_seurat, directory_path='monocle3/allcells/')

# Before we order the cells in pseudotime, we need to identify which is the pricipal node/cells (e.g.,HSC)
# In CD4-lineage, it would be good to set "SOX4+NaiveCD4"

# This is using the GUI
# cds_from_seurat <- order_cells(cds_from_seurat)

root_cell_list <- meta_tmp%>%filter(cellstate=="NaiveCD4-3")%>%pull(cell_id)
cds_from_seurat <- order_cells(cds_from_seurat, root_cells = "AACCGCGTCGACCAGC-1")

## trajectory umap
## color: pseudotime
  p <- plot_cells(cds_from_seurat, 
                color_cells_by = 'pseudotime',
                label_groups_by_cluster=FALSE,
                label_roots =FALSE,
                label_leaves=FALSE,
                label_branch_points=FALSE)

png(paste0("monocle3/",lineage_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_postUMAP_monocle3_PT.png"),height=2.5, width=3.5, units = "in",res = 300)
   plot(p)
dev.off()


##############################################################################################################################
##############################################################################################################################

sessionInfo()

R version 4.4.3 (2025-02-28)
Platform: x86_64-conda-linux-gnu
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS/LAPACK: /rshare1/ZETTAI_path_WA_slash_home_KARA/home/imgnaka/miniconda3_new/envs/monocle3/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

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
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] monocle3_1.4.26             SingleCellExperiment_1.28.0
 [3] SummarizedExperiment_1.36.0 GenomicRanges_1.58.0       
 [5] GenomeInfoDb_1.42.0         IRanges_2.40.0             
 [7] S4Vectors_0.44.0            MatrixGenerics_1.18.0      
 [9] matrixStats_1.5.0           Biobase_2.66.0             
[11] BiocGenerics_0.52.0         Seurat_5.3.1               
[13] SeuratObject_5.2.0          sp_2.2-0                   
[15] data.table_1.17.8           lubridate_1.9.4            
[17] forcats_1.0.1               stringr_1.6.0              
[19] dplyr_1.1.4                 purrr_1.2.0                
[21] readr_2.1.6                 tidyr_1.3.1                
[23] tibble_3.3.0                ggplot2_4.0.1              
[25] tidyverse_2.0.0            

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3      jsonlite_2.0.0          magrittr_2.0.4         
  [4] spatstat.utils_3.2-0    nloptr_2.2.1            farver_2.1.2           
  [7] zlibbioc_1.52.0         vctrs_0.6.5             ROCR_1.0-11            
 [10] minqa_1.2.8             spatstat.explore_3.6-0  htmltools_0.5.8.1      
 [13] S4Arrays_1.6.0          SparseArray_1.6.0       sctransform_0.4.2      
 [16] parallelly_1.45.1       KernSmooth_2.23-26      htmlwidgets_1.6.4      
 [19] ica_1.0-3               plyr_1.8.9              plotly_4.11.0          
 [22] zoo_1.8-14              igraph_2.1.4            mime_0.13              
 [25] lifecycle_1.0.4         pkgconfig_2.0.3         Matrix_1.7-4           
 [28] R6_2.6.1                fastmap_1.2.0           rbibutils_2.4          
 [31] GenomeInfoDbData_1.2.13 fitdistrplus_1.2-4      future_1.68.0          
 [34] shiny_1.12.0            digest_0.6.39           patchwork_1.3.2        
 [37] tensor_1.5.1            RSpectra_0.16-2         irlba_2.3.5.1          
 [40] progressr_0.18.0        spatstat.sparse_3.1-0   timechange_0.3.0       
 [43] httr_1.4.7              polyclip_1.10-7         abind_1.4-8            
 [46] compiler_4.4.3          withr_3.0.2             S7_0.2.1               
 [49] fastDummies_1.7.5       MASS_7.3-65             DelayedArray_0.32.0    
 [52] tools_4.4.3             lmtest_0.9-40           otel_0.2.0             
 [55] httpuv_1.6.16           future.apply_1.20.0     goftest_1.2-3          
 [58] glue_1.8.0              nlme_3.1-168            promises_1.5.0         
 [61] grid_4.4.3              Rtsne_0.17              cluster_2.1.8.1        
 [64] reshape2_1.4.5          generics_0.1.4          gtable_0.3.6           
 [67] spatstat.data_3.1-9     tzdb_0.5.0              hms_1.1.4              
 [70] XVector_0.46.0          spatstat.geom_3.6-1     RcppAnnoy_0.0.22       
 [73] ggrepel_0.9.6           RANN_2.6.2              pillar_1.11.1          
 [76] spam_2.11-1             RcppHNSW_0.6.0          later_1.4.4            
 [79] splines_4.4.3           lattice_0.22-7          survival_3.8-3         
 [82] deldir_2.0-4            tidyselect_1.2.1        miniUI_0.1.2           
 [85] pbapply_1.7-4           reformulas_0.4.2        gridExtra_2.3          
 [88] scattermore_1.2         stringi_1.8.7           UCSC.utils_1.2.0       
 [91] boot_1.3-32             lazyeval_0.2.2          codetools_0.2-20       
 [94] cli_3.6.5               uwot_0.2.4              Rdpack_2.6.4           
 [97] xtable_1.8-4            reticulate_1.44.1       Rcpp_1.1.0             
[100] globals_0.18.0          spatstat.random_3.4-3   png_0.1-8              
[103] spatstat.univar_3.1-5   parallel_4.4.3          dotCall64_1.2          
[106] lme4_1.1-38             listenv_0.10.0          viridisLite_0.4.2      
[109] scales_1.4.0            ggridges_0.5.7          crayon_1.5.3           
[112] rlang_1.1.6             cowplot_1.2.0 







