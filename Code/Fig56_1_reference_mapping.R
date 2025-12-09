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

















