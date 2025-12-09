source("~/ws/my.source.R")
options(stringsAsFactors=F)

###################
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(msigdbr))


celltype_tmp = "NaiveCD4"
compare_tmp  = "InactivevsHC"


# Bg genes Z score Top1000 union
res_tmp=fread_FT("singlecell_NBmodel/pathway/ORAtop1000Pval/eachcelltype_eachcompare_dys_nbhd_top1000Pval_gene_union.txt")

# add ENTREZID
entrez_ids=mapIds(org.Hs.eg.db, keys = res_tmp$gene,column = "ENTREZID", keytype = "SYMBOL")%>%
           as.data.frame()%>%rownames_to_column("gene")
colnames(entrez_ids)[2]="ENTREZID"

res_tmp2=left_join(res_tmp,entrez_ids,by="gene")%>%
         .[!is.na(.$ENTREZID),]

paste0("BACK genenum: ",nrow(res_tmp))%>% print()
paste0("BACK analyzed genenum: ",nrow(res_tmp2))%>% print()


##############################
## dysregulated term Top1000 all
##############################
res_tmp_DEG=fread_FT(paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_dys_top1000Pval_genes.txt"))

res_tmp2_DEG=left_join(res_tmp_DEG,entrez_ids,by="gene")%>%
             .[!is.na(.$ENTREZID),]

paste0("DEG genenum: ",nrow(res_tmp_DEG))%>% print()
paste0("DEG analyzed genenum: ",nrow(res_tmp2_DEG))%>% print()


# GSEA (HALLMARK)
HALLMARK = msigdbr(species="Homo sapiens",category="H") %>% 
           dplyr::select(gs_name,entrez_gene)

HALLMARK_res1 = enricher(gene=res_tmp2_DEG$ENTREZID,universe=res_tmp2$ENTREZID,TERM2GENE=HALLMARK,maxGSSize=2000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(HALLMARK_res1,paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_dys_ORAtop1000Pval_HALLMARK.txt"))


# GSEA (KEGG)
KEGG_res1 = enrichKEGG(gene=res_tmp2_DEG$ENTREZID,universe=res_tmp2$ENTREZID,organism='hsa',maxGSSize=2000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(KEGG_res1,paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_dys_ORAtop1000Pval_KEGG.txt"))


##############################
## abundance term Top1000 all
##############################
res_tmp_nbhd=fread_FT(paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_nbhd_top1000Pval_genes.txt"))

res_tmp2_nbhd=left_join(res_tmp_nbhd,entrez_ids,by="gene")%>%
             .[!is.na(.$ENTREZID),]

paste0("nbhd genenum: ",nrow(res_tmp_nbhd))%>% print()
paste0("nbhd analyzed genenum: ",nrow(res_tmp2_nbhd))%>% print()

# GSEA (HALLMARK)
HALLMARK_res2 = enricher(gene=res_tmp2_nbhd$ENTREZID,universe=res_tmp2$ENTREZID,TERM2GENE=HALLMARK,maxGSSize=2000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(HALLMARK_res2,paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_nbhd_ORAtop1000Pval_HALLMARK.txt"))

# GSEA (KEGG)
KEGG_res2 = enrichKEGG(gene=res_tmp2_nbhd$ENTREZID,universe=res_tmp2$ENTREZID,organism='hsa',maxGSSize=2000,pvalueCutoff=1,qvalueCutoff=1)
write.table_FT_2(KEGG_res2,paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_nbhd_ORAtop1000Pval_KEGG.txt"))


##############################
## summarize pathway (no filter)
##############################

HALLMARK_res1_lim=HALLMARK_res1%>%as.data.frame()%>%mutate(sigtype="dys")%>%dplyr::select(sigtype,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue)

HALLMARK_res2_lim=HALLMARK_res2%>%as.data.frame()%>%mutate(sigtype="nbhd")%>%dplyr::select(sigtype,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue)

KEGG_res1_lim=KEGG_res1%>%as.data.frame()%>%mutate(sigtype="dys")%>%dplyr::select(sigtype,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue)

KEGG_res2_lim=KEGG_res2%>%as.data.frame()%>%mutate(sigtype="nbhd")%>%dplyr::select(sigtype,Description,GeneRatio,BgRatio,pvalue,p.adjust,qvalue)

all_sum=bind_rows(HALLMARK_res1_lim,KEGG_res1_lim,HALLMARK_res2_lim,KEGG_res2_lim)
write.table_FT_2(all_sum,paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/",celltype_tmp,"/",celltype_tmp,"_",compare_tmp,"_ORAtop1000Pval_HALLMARK_KEGG_res_sum_all.txt"))



######################################################################################################################################################
######################################################################################################################################################


##############################
## aggregate results
##############################

All_list=fread_FT("tmp_job_list/CNA_celltype_list.txt") # 18 major celltypes
celltype_list=All_list$celltype

list =make_list("singlecell_NBmodel/pathway/ORAtop1000Pval","_ORAtop1000Pval_HALLMARK_KEGG_res_sum_all.txt")
list$celltype=take_factor(list$FILE,1,"_")
list$compare =take_factor(list$FILE,2,"_")
list$celltype=factor(list$celltype,levels=celltype_list)
list$compare =factor(list$compare,levels=c("InactivevsHC","HDAvsInactive"))
list=list[order(list$celltype),]
list=list[order(list$compare),]

for(iii in 1:nrow(list)){
    celltype_tmp=list$celltype[iii]
    compare_tmp =list$compare[iii]
    res_tmp=fread_FT(list$PATH[iii])%>%mutate(celltype=celltype_tmp)%>%mutate(compare=compare_tmp)

   if(iii==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}
res_sum2=res_sum%>%mutate(p.adjust_all=p.adjust(pvalue,method="BH"))%>%
                dplyr::select(8,9,1:5,10)
write.table_FT_2(res_sum2,paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/summary/eachcelltype_eachcompare_ORAtop1000Pval_HALLMARK_KEGG_res_sum_all.txt"))


########################
## representative pathways -logP barplot
########################
suppressPackageStartupMessages(library(ggsci))

All_list=fread_FT("tmp_job_list/CNA_celltype_list.txt")
celltype_list=All_list$celltype

pathway_list=c("HALLMARK_INTERFERON_ALPHA_RESPONSE","HALLMARK_INTERFERON_GAMMA_RESPONSE","HALLMARK_COMPLEMENT",
               "Antigen processing and presentation","Cell adhesion molecules","Cellular senescence","T cell receptor signaling pathway")

res_sum=fread_FT("singlecell_NBmodel/pathway/ORAtop1000Pval/summary/eachcelltype_eachcompare_ORAtop1000Pval_HALLMARK_KEGG_res_sum_all.txt")%>%
        dplyr::select(celltype,compare,sigtype,Description,pvalue,p.adjust_all)

# FDR=0.05 corresponding P
res_sum%>%filter(p.adjust_all<0.05)%>%pull(pvalue)%>%max()
# 0.003329062
res_sum%>%filter(p.adjust_all>=0.05)%>%pull(pvalue)%>%min()
# 0.003335578
FDR_threshold=(0.003329062+0.003335578)/2

res_sum$compare=factor(res_sum$compare,levels=c("InactivevsHC","HDAvsInactive"))
res_sum=res_sum[order(res_sum$compare),]

pattern_all=res_sum[,1:3]%>%unique()

for(iii in 1:length(pathway_list)){
     pathway_tmp=pathway_list[iii]
     print(pathway_tmp)

     res_tmp=res_sum%>%filter(Description==pathway_tmp)%>%
                       full_join(pattern_all,.,by=c("celltype","compare","sigtype")) # show NA
     res_tmp$pvalue[is.na(res_tmp$pvalue)]=1

     res_tmp=res_tmp%>%mutate(logP=-log10(pvalue))%>%
                       mutate(logP_mod=ifelse(logP>lim,lim,logP))

     res_tmp$celltype=factor(res_tmp$celltype,levels=celltype_list)
     res_tmp$compare=factor(res_tmp$compare,levels=c("InactivevsHC","HDAvsInactive"))
     res_tmp$sigtype=factor(res_tmp$sigtype,levels=c("nbhd","dys"))
     levels(res_tmp$sigtype)=c("Abundance signature","Dysregulated signare")

     levels(res_tmp$celltype)=gsub("CD56brightNK","CD56++NK",levels(res_tmp$celltype))%>%
                              gsub("CD16nMono","CD16-Mono",.)%>%
                              gsub("CD16pMono","CD16+Mono",.)   
                              
  p =  ggplot(data=res_tmp,aes(x=celltype,y=logP_mod,fill=sigtype))+
     geom_hline(yintercept=-log10(FDR_threshold),linewidth=0.5,col="darkgrey")+
     geom_bar(stat="identity",position="dodge")+
     theme_classic()+
     scale_fill_npg() +
     facet_wrap(~compare,ncol=2)+
     theme(axis.text.x=element_text(colour="black",size=14,angle=45,hjust=1,vjust=1),
           axis.text.y=element_text(colour="black",size=12),
           axis.title.x=element_blank(),
           axis.title.y=element_text(colour="black",size=14),
           plot.title=element_blank(),
           strip.background=element_blank(),
           strip.text=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=14))+
     scale_y_continuous(limits=c(0,lim))+
     labs(y="-log10(P)")

  pdf_3(paste0("singlecell_NBmodel/pathway/ORAtop1000Pval/summary/eachcelltype_eachcompare_ORAtop1000Pval_",pathway_tmp,"_bar.pdf"),h=4,w=12)
     plot(p)
  dev.off()
}











