source("~/ws/my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(liana))


#################################
## 1. input data
#################################

show_resources()
#  [1] "Default"          "Consensus"        "Baccin2019"       "CellCall"        
#  [5] "CellChatDB"       "Cellinker"        "CellPhoneDB"      "CellTalkDB"      
#  [9] "connectomeDB2020" "EMBRACE"          "Guide2Pharma"     "HPMR"            
# [13] "ICELLNET"         "iTALK"            "Kirouac2010"      "LRdb"            
# [17] "Ramilowski2015"   "OmniPath"         "MouseConsensus"

db_consensus = select_resource("Consensus")[[1]]
dim(db_consensus)
# [1] 4701   10

## metadata, only CLUES cells
meta_tmp=fread_FT("SLE_maindata_4thQC_new_meta_anno.txt")%>%
         filter(!str_detect(donor_id,pattern="IGTB"))
# 1108721

# countdata (SeuratObject)
seurat_tmp1=readRDS("~/ws/2025/250206_SLE_maindata_new_QC/RNA_QC/2ndQC/SLE_maindata_2ndQC_new_1.rds")%>%subset(cells=meta_tmp$cell_id)
seurat_tmp2=readRDS("~/ws/2025/250206_SLE_maindata_new_QC/RNA_QC/2ndQC/SLE_maindata_2ndQC_new_2.rds")%>%subset(cells=meta_tmp$cell_id)
seurat=merge(seurat_tmp1,seurat_tmp2)
dim(seurat)%>%print()
# 36601 1108721

seurat@meta.data=meta_tmp%>%column_to_rownames("cell_id")

all.equal(rownames(seurat@meta.data),colnames(seurat$RNA@data))
# [1] TRUE

seurat=NormalizeData(seurat)

## Avoid multishread
# future::plan("multisession", workers = 4) # do parallel
# options(future.globals.maxSize = 6000 * 1024^2)


#################################
## 2. Run liana
#################################
liana_test = liana_wrap(seurat,
                        method = c("logfc","sca"),
                        resource = "Consensus",
                        idents_col = "cellstate",
                        min_cells = 10,
                        expr_prop = 0.01, # for ligand
                        base = exp(1))    # Seurat: lognorm with base=e


liana_aggr = liana_test %>%
             liana_aggregate()

liana_aggr2= liana_aggr%>%
             mutate(source_celltype=take_factor(source,1,"-"))%>%
             mutate(target_celltype=take_factor(target,1,"-"))
# 5451435


#################################
## 3. Add low exp filter (0.05) for receptor
#################################
str(liana_test)

dim(liana_test[[1]])
# 5451435      11

min(liana_test[[1]]$ligand.prop) # 0.01000211
min(liana_test[[1]]$receptor.prop) # 0.01000211

lr_genes_raw = union(liana_test[[1]]$ligand.complex, liana_test[[1]]$receptor.complex)
# 1154, including complex receptors
lr_genes_split = str_split(lr_genes_raw, "_") %>% 
  unlist() %>% unique()
# 1087

# low exp filter (0.05) for receptor
# also exclude RBC_PLT_genes
RBC_PLT_genes=fread_FT("geneset/RBCPLT_genes_list_v2.txt")%>%pull(Gene)

liana_NG=    liana_test[[1]]%>%filter(receptor.prop > 0.05)%>%
             mutate(ligand1=take_factor(ligand.complex,1,"_"))%>%
             mutate(ligand2=take_factor(ligand.complex,2,"_"))%>%
             #mutate(ligand3=take_factor(ligand.complex,3,"_"))%>%
             mutate(receptor1=take_factor(receptor.complex,1,"_"))%>%
             mutate(receptor2=take_factor(receptor.complex,2,"_"))%>%
             mutate(receptor3=take_factor(receptor.complex,3,"_"))%>%
             #mutate(receptor4=take_factor(receptor.complex,4,"_"))%>%
             filter(ligand1%in%RBC_PLT_genes | ligand2%in%RBC_PLT_genes | receptor1%in%RBC_PLT_genes | receptor2%in%RBC_PLT_genes | receptor3%in%RBC_PLT_genes)%>%
             dplyr::select(source,target,ligand.complex,receptor.complex,ligand.prop,receptor.prop,ligand.log2FC,receptor.log2FC)
dim(liana_NG)
# 229025      16

liana_NG%>%dplyr::select(ligand.complex,receptor.complex)%>%unique()%>%arrange(ligand.complex)%>%print(n=Inf)
# 88
# Some might include immune-signal wo PLT/RBC, but exclude them.


liana_filter=liana_test[[1]]%>%filter(receptor.prop > 0.05)%>%
             mutate(ligand1=take_factor(ligand.complex,1,"_"))%>%
             mutate(ligand2=take_factor(ligand.complex,2,"_"))%>%
             #mutate(ligand3=take_factor(ligand.complex,3,"_"))%>%
             mutate(receptor1=take_factor(receptor.complex,1,"_"))%>%
             mutate(receptor2=take_factor(receptor.complex,2,"_"))%>%
             mutate(receptor3=take_factor(receptor.complex,3,"_"))%>%
             #mutate(receptor4=take_factor(receptor.complex,4,"_"))%>%
             filter(!ligand1%in%RBC_PLT_genes & !ligand2%in%RBC_PLT_genes & !receptor1%in%RBC_PLT_genes & !receptor2%in%RBC_PLT_genes & !receptor3%in%RBC_PLT_genes)%>%
             dplyr::select(source,target,ligand.complex,receptor.complex,ligand.prop,receptor.prop,ligand.log2FC,receptor.log2FC)
# 2807028
min(liana_filter$receptor.prop) # 0.05000899

liana_aggr3=inner_join(liana_aggr2,liana_filter,by=c("source","target","ligand.complex","receptor.complex"))
dim(liana_aggr3)
# 2807028      16



#################################
## 4. Pick up cell-state-specific interactions
#################################

## sca.LRscore > 0.6 & logfc.logfc_comb > -0.1
liana_aggr3_up=liana_aggr3%>%filter(sca.LRscore>0.6 & logfc.logfc_comb>(-0.1))%>%
               mutate(signal=paste0(ligand.complex,"_",receptor.complex))%>%
               mutate(celltype_signal=paste0(source_celltype,"_",ligand.complex,"_",target_celltype,"_",receptor.complex))
# 1451770


## Pick up the most upregulated cell-state pair as the representative for each unique LR interaction within each of the 27 × 27 cell-type pairs
celltype_list=c("NaiveCD4","CMCD4","EMCD4","CTLCD4","Treg","ProlifCD4",
                "NaiveCD8","CMCD8","EMCD8","gdT","MAIT","dnT","ProlifCD8",
                "CD56dimNK","CD56brightNK","ProlifNK","ILC",
                "NaiveB","MemB","PlasmaB",
                "CD16nMono","CD16pMono","C1QMono",
                "cDC","pDC","ASDC",
                "HSPC")

liana_aggr3_up$source_celltype=factor(liana_aggr3_up$source_celltype,levels=celltype_list)
liana_aggr3_up$target_celltype=factor(liana_aggr3_up$target_celltype,levels=celltype_list)

liana_aggr3_up=liana_aggr3_up%>%arrange(target)%>%arrange(source)%>%
                                arrange(target_celltype)%>%arrange(source_celltype)

signal_unique=unique(liana_aggr3_up$signal)
# 1053
celltype_signal_unique=unique(liana_aggr3_up$celltype_signal)
# 118629

liana_aggr3_up_top=liana_aggr3_up%>%arrange(desc(sca.LRscore))%>%
                   .[!duplicated(.$celltype_signal),]
# 118629

write.table_FT_2(liana_aggr3_up_top,paste0("liana/liana_cellstate_logFC_sca_aggr_filtered_enhanced_cellstate_pairs.txt"))








