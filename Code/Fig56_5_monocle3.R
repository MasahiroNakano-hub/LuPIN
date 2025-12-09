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











