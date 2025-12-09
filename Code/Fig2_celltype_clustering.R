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
seurat.list=lapply(seurat.list,function(x){SCTransform(x, vst.flavor = "v1", method = "glmGamPoi", vars.to.regress = "percent.mt", do.scale = FALSE, do.center = FALSE, verbose = FALSE, min_cells = mincell_tmp)})

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


###########################################
# Fig2 plot
# NAiveCD4

R4_C0=fread_FT(paste0(lineage_tmp,"_",cluster_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_res_",res_tmp,"_idlist.txt"))%>%
      mutate(anno_final=ifelse(cluster==0,"NaiveCD4-0",
                     ifelse(cluster==1,"NaiveCD4-1",
                     ifelse(cluster==2,"NaiveCD4-2",
                     ifelse(cluster==4,"NaiveCD4-3",
                     ifelse(cluster==3,"NaiveCD4-4",
                     ifelse(cluster==5,"NaiveCD4-5",
                     ifelse(cluster==6,"NaiveCD4-6",NA
                     ))))))))
table(R4_C0$anno_final)
# NaiveCD4-0 NaiveCD4-1 NaiveCD4-2 NaiveCD4-3 NaiveCD4-4 NaiveCD4-5 NaiveCD4-6 
#      81706      58974      40457      11490      14515       5849       3701

umap_C0=left_join(umap_post,R4_C0,by="cell_id")
dim(umap_C0)
# 216692      5
umap_C0$anno_final=factor(umap_C0$anno_final,levels=c(paste0("NaiveCD4-",0:6)))
col_list=brewer.pal(9,"Reds")[3:9]

p = ggplot(umap_C0,aes(x=UMAP_1,y=UMAP_2,color=anno_final))+
    geom_point(size=0.05)+
    theme_void()+
    scale_color_manual(values = col_list) +
    theme(axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title=element_blank(),
          legend.position="none")

png(paste0(lineage_tmp,"_",cluster_tmp,"_mincell",mincell_tmp,"_PC",PC_tmp,"_sigma",sigma_tmp,"_postUMAP_res_",res_tmp,"_cluster.png"),height=5, width=6.5, units = "in",res = 300)
   plot(p)
dev.off()

