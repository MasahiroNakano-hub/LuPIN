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

suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(cowplot))


## Ex. NaiveCD4
lineage_tmp = "CD4"
cluster_tmp = "NaiveCD4"
mincell_tmp = 10
PC_tmp   =20
sigma_tmp=0.25
prune_tmp=0.05
res_tmp  =0.35

##############################################################################################################################


## Cell_id list for one cell type
id_tmp=fread_FT(paste0("SNN_clustering/eachcelltype_R2/",lineage_tmp,"/",lineage_tmp,"_anno_R3.txt"))%>%
       filter(anno_R3==cluster_tmp)

# metadata
meta_tmp=fread_FT("~/ws/2025/250206_SLE_maindata_new_QC/RNA_QC/4thQC/summary/SLE_maindata_4thQC_new_meta_clinical.txt")%>%
         filter(cell_id%in%id_tmp$cell_id)

# countdata (SeuratObject)
seurat=readRDS(paste0("~/ws/2025/250206_SLE_maindata_new_QC/RNA_QC/4thQC/summary/",lineage_tmp,"_4thQC_seurat.rds"))%>%
         subset(cells=meta_tmp$cell_id)

# Split the dataset into batches
DefaultAssay(seurat) = "RNA"
seurat@meta.data$batch=factor(seurat@meta.data$batch,levels=unique(seurat@meta.data$batch))
seurat.list=SplitObject(seurat, split.by = "batch")

# sctransform for each batch, opt min_cells
# SCT Pearson residuals are not centered nor scaled here
seurat.list=lapply(seurat.list,function(x){SCTransform(x, vst.flavor = "v1", method = "glmGamPoi", do.scale = FALSE, do.center = FALSE, verbose = FALSE, min_cells = mincell_tmp)})

# pick up hv3000 and exclude bad genes
features = SelectIntegrationFeatures(object.list = seurat.list, nfeatures = 3000)
bad_genes = grep("^MT-|MALAT1|NEAT1|TMSB4X|TMSB10", rownames(seurat$RNA@counts),value=TRUE)
chrY_genes=fread_FF("geneset/hg38_chrY_gene_list.txt")%>%pull(V1)
CC_genes=Seurat::cc.genes
RBC_PLT_genes=fread_FT("geneset/RBCPLT_genes_list_v2.txt")%>%pull(Gene)
IGKL_genes= grep("^IGKC|^IGKJ|^IGKV|^IGLC|^IGLJ|^IGLL|^IGLV", rownames(seurat$RNA@counts),value=TRUE)
features_use=features[!is.element(features,c(bad_genes,chrY_genes,CC_genes$s.genes,CC_genes$g2m.genes,RBC_PLT_genes,IGKL_genes))]
features_df=data.frame(Gene=features_use)


# merge each SCT data
seurat.list2=PrepSCTIntegration(object.list = seurat.list, anchor.features = features_use)
seurat_merged2=merge(seurat.list2[[1]], y=seurat.list2[2:length(seurat.list2)],merge.data = TRUE)
dim(seurat_merged2$SCT@scale.data)%>%print()

# mean sd for symphony
mean=apply(seurat_merged2$SCT@scale.data,1,mean)
sd  =apply(seurat_merged2$SCT@scale.data,1,sd)
meansd=data.frame(symbol=names(mean),mean=mean,stddev=sd)

# scale Pearson rediduals across all batches
# cannot use apply function due to large file size 
exprs_scaled = matrix(0, nrow=ncol(seurat_merged2$SCT@scale.data), ncol=nrow(seurat_merged2$SCT@scale.data) )
rownames(exprs_scaled) = colnames(seurat_merged2$SCT@scale.data)
colnames(exprs_scaled) = rownames(seurat_merged2$SCT@scale.data)

for(iii in 1:nrow(seurat_merged2$SCT@scale.data)){
    print(iii)
    gene_tmp=rownames(seurat_merged2$SCT@scale.data)[iii]
    scale_tmp=scale(seurat_merged2$SCT@scale.data[iii,])
    exprs_scaled[,gene_tmp]=scale_tmp
}

dim(exprs_scaled)%>%print()

# PCA: ~PC30
set.seed(42)
pca_res = irlba::prcomp_irlba(exprs_scaled, 30)
rownames(pca_res$x)=rownames(exprs_scaled)
pca_res2=pca_res$x[,1:PC_tmp]
loadings=pca_res$rotation[,1:PC_tmp]
rownames(loadings)=names(pca_res$center)

# metadata: confirm the consistency of cell_id order
meta_tmp=meta_tmp%>%column_to_rownames("cell_id")
meta_tmp=meta_tmp[rownames(pca_res2),]
all.equal(rownames(pca_res2),rownames(meta_tmp))%>%print() #[1] TRUE

# new harmony with each sigma
# return_object = T (for reference mapping)
set.seed(42)
harmony_obj = RunHarmony(pca_res2, meta_tmp, c("lib","sample_uuid"), theta = c(2,2), lambda = c(1,1),
                               plot_convergence = TRUE, nclust = NULL, sigma=sigma_tmp,
                               early_stop = FALSE,  max_iter = 20, verbose = T, return_object = T)

harmony_res=t(harmony_obj$Z_corr)
rownames(harmony_res)=rownames(pca_res2)
colnames(harmony_res)=paste0("Harmony_",1:PC_tmp)

# build symphony refence
reference = symphony::buildReferenceFromHarmonyObj(
                           harmony_obj,            # output object from HarmonyMatrix()
                           meta_tmp,               # reference cell metadata
                           meansd,                 # gene names, means, and std devs for scaling
                           loadings,               # genes x PCs matrix
                           verbose = TRUE,         # verbose output
                           do_umap = FALSE,         # Set to TRUE only when UMAP model was saved for reference
                           seed=42
                           )

saveRDS(reference,paste0(lineage_tmp,"_",cluster_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_synphony_reference.rds"))


# UMAP post harmony, k=20
# seed=42: recall final UMAP model ret_model = TRUE
# save uwot model
set.seed(42)
umap_post = umap(harmony_res, n_neighbors = 20L, metric = "euclidean", min_dist = .3, ret_model = TRUE)


# clustering
ids_1 = Seurat:::RunModularityClustering(SNN = snn_1, modularity = 1,
        resolution = res_tmp, algorithm = 1, n.start = 20,
        n.iter = 20, random.seed = 0, print.output = FALSE,
        temp.file.location = NULL, edge.file.name = NULL)
# algorithm: 1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm
# Modularity: 1 = standard; 2 = alternative

id_list = data.frame(cell_id=rownames(snn_1),cluster=ids_1)
write.table_FT_2(id_list,paste0(lineage_tmp,"_",cluster_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_res_",res_tmp,"_idlist.txt"))

## See Fig2_2_celltype_clustering_visualization.ipynb for visualization.



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
 [1] cowplot_1.1.3           gridExtra_2.3           uwot_0.2.2             
 [4] Matrix_1.6-5            harmony_1.2.3           Rcpp_1.0.14            
 [7] singlecellmethods_0.1.0 glmGamPoi_1.18.0        sctransform_0.4.1      
[10] Seurat_5.1.0            SeuratObject_5.0.2      sp_2.1-4               
[13] data.table_1.16.4       lubridate_1.9.4         forcats_1.0.0          
[16] stringr_1.5.1           dplyr_1.1.4             purrr_1.0.2            
[19] readr_2.1.5             tidyr_1.3.1             tibble_3.2.1           
[22] ggplot2_3.5.1           tidyverse_2.0.0        

loaded via a namespace (and not attached):
  [1] RColorBrewer_1.1-3          jsonlite_1.8.9             
  [3] magrittr_2.0.3              spatstat.utils_3.1-2       
  [5] farver_2.1.2                zlibbioc_1.52.0            
  [7] vctrs_0.6.5                 ROCR_1.0-11                
  [9] spatstat.explore_3.3-4      S4Arrays_1.6.0             
 [11] htmltools_0.5.8.1           SparseArray_1.6.0          
 [13] parallelly_1.41.0           KernSmooth_2.23-26         
 [15] htmlwidgets_1.6.4           ica_1.0-3                  
 [17] plyr_1.8.9                  plotly_4.10.4              
 [19] zoo_1.8-12                  igraph_2.0.3               
 [21] mime_0.12                   lifecycle_1.0.4            
 [23] pkgconfig_2.0.3             R6_2.5.1                   
 [25] fastmap_1.2.0               GenomeInfoDbData_1.2.13    
 [27] MatrixGenerics_1.18.0       fitdistrplus_1.2-2         
 [29] future_1.34.0               shiny_1.10.0               
 [31] digest_0.6.37               colorspace_2.1-1           
 [33] patchwork_1.3.0             S4Vectors_0.44.0           
 [35] tensor_1.5                  RSpectra_0.16-2            
 [37] irlba_2.3.5.1               GenomicRanges_1.58.0       
 [39] progressr_0.15.1            spatstat.sparse_3.1-0      
 [41] timechange_0.3.0            httr_1.4.7                 
 [43] polyclip_1.10-7             abind_1.4-5                
 [45] compiler_4.4.2              withr_3.0.2                
 [47] fastDummies_1.7.5           MASS_7.3-64                
 [49] DelayedArray_0.32.0         tools_4.4.2                
 [51] lmtest_0.9-40               httpuv_1.6.15              
 [53] future.apply_1.11.3         goftest_1.2-3              
 [55] glue_1.8.0                  nlme_3.1-165               
 [57] promises_1.3.2              grid_4.4.2                 
 [59] Rtsne_0.17                  cluster_2.1.8              
 [61] reshape2_1.4.4              generics_0.1.3             
 [63] gtable_0.3.6                spatstat.data_3.1-4        
 [65] tzdb_0.4.0                  hms_1.1.3                  
 [67] XVector_0.46.0              BiocGenerics_0.52.0        
 [69] spatstat.geom_3.3-4         RcppAnnoy_0.0.22           
 [71] ggrepel_0.9.6               RANN_2.6.2                 
 [73] pillar_1.10.1               spam_2.11-1                
 [75] RcppHNSW_0.6.0              later_1.4.1                
 [77] splines_4.4.2               lattice_0.22-6             
 [79] survival_3.8-3              deldir_2.0-4               
 [81] tidyselect_1.2.1            miniUI_0.1.1.1             
 [83] pbapply_1.7-2               IRanges_2.40.0             
 [85] SummarizedExperiment_1.36.0 scattermore_1.2            
 [87] stats4_4.4.2                Biobase_2.66.0             
 [89] matrixStats_1.5.0           stringi_1.8.4              
 [91] UCSC.utils_1.2.0            lazyeval_0.2.2             
 [93] codetools_0.2-20            cli_3.6.3                  
 [95] xtable_1.8-4                reticulate_1.40.0          
 [97] munsell_0.5.1               GenomeInfoDb_1.42.0        
 [99] globals_0.16.3              spatstat.random_3.3-2      
[101] png_0.1-8                   spatstat.univar_3.1-1      
[103] parallel_4.4.2              dotCall64_1.2              
[105] listenv_0.9.1               viridisLite_0.4.2          
[107] scales_1.3.0                ggridges_0.5.6             
[109] crayon_1.5.3                leiden_0.4.3.1             
[111] rlang_1.1.5                














