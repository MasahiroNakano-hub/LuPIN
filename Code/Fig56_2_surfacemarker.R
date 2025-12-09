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








