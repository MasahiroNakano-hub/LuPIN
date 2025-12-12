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



##############################################################################################################################
##############################################################################################################################

sessionInfo()

R version 4.4.3 (2025-02-28)
Platform: x86_64-conda-linux-gnu
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS/LAPACK: /rshare1/ZETTAI_path_WA_slash_home_KARA/home/imgnaka/miniconda3_new/envs/clusterprofiler/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0

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
 [1] msigdbr_25.1.1         org.Hs.eg.db_3.20.0    AnnotationDbi_1.68.0  
 [4] IRanges_2.40.0         S4Vectors_0.44.0       Biobase_2.66.0        
 [7] BiocGenerics_0.52.0    clusterProfiler_4.14.0 data.table_1.17.8     
[10] lubridate_1.9.4        forcats_1.0.0          stringr_1.5.2         
[13] dplyr_1.1.4            purrr_1.1.0            readr_2.1.5           
[16] tidyr_1.3.1            tibble_3.3.0           ggplot2_4.0.0         
[19] tidyverse_2.0.0       

loaded via a namespace (and not attached):
 [1] tidyselect_1.2.1        farver_2.1.2            blob_1.2.4             
 [4] R.utils_2.13.0          Biostrings_2.74.0       S7_0.2.0               
 [7] lazyeval_0.2.2          fastmap_1.2.0           digest_0.6.37          
[10] timechange_0.3.0        lifecycle_1.0.4         KEGGREST_1.46.0        
[13] tidytree_0.4.6          RSQLite_2.4.3           magrittr_2.0.4         
[16] compiler_4.4.3          rlang_1.1.6             tools_4.4.3            
[19] igraph_2.1.4            ggtangle_0.0.7          curl_7.0.0             
[22] bit_4.6.0               gson_0.1.0              plyr_1.8.9             
[25] RColorBrewer_1.1-3      aplot_0.2.9             BiocParallel_1.40.0    
[28] babelgene_22.9          withr_3.0.2             R.oo_1.27.1            
[31] grid_4.4.3              GOSemSim_2.32.0         enrichplot_1.26.1      
[34] colorspace_2.1-1        GO.db_3.20.0            scales_1.4.0           
[37] cli_3.6.5               crayon_1.5.3            treeio_1.30.0          
[40] generics_0.1.4          ggtree_3.14.0           httr_1.4.7             
[43] reshape2_1.4.4          tzdb_0.5.0              ape_5.8-1              
[46] DBI_1.2.3               qvalue_2.38.0           cachem_1.1.0           
[49] DOSE_4.0.0              zlibbioc_1.52.0         splines_4.4.3          
[52] assertthat_0.2.1        parallel_4.4.3          ggplotify_0.1.2        
[55] XVector_0.46.0          vctrs_0.6.5             yulab.utils_0.2.1      
[58] Matrix_1.7-4            jsonlite_2.0.0          patchwork_1.3.2        
[61] gridGraphics_0.5-1      hms_1.1.3               ggrepel_0.9.6          
[64] bit64_4.6.0-1           glue_1.8.0              codetools_0.2-20       
[67] cowplot_1.2.0           stringi_1.8.7           gtable_0.3.6           
[70] GenomeInfoDb_1.42.0     UCSC.utils_1.2.0        pillar_1.11.1          
[73] rappdirs_0.3.3          fgsea_1.32.2            GenomeInfoDbData_1.2.13
[76] R6_2.6.1                lattice_0.22-7          R.methodsS3_1.8.2      
[79] png_0.1-8               memoise_2.0.1           ggfun_0.2.0            
[82] Rcpp_1.1.0              fastmatch_1.1-6         nlme_3.1-168           
[85] fs_1.6.6                pkgconfig_2.0.3 





