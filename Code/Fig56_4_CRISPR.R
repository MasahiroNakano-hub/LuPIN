
source("~/ws/my.source.R")
options(stringsAsFactors=F)
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))

# CRISPR sample metadata
metadata=fread_FT("Info/sampleID_corresp_v2_all.txt")%>%
         mutate(Donor_ID=paste0("Donor",take_factor(sample_ID,2,"_")))
# 225
KO_Gene_list=unique(metadata$KO_Gene)[1:22]

# each KO cond GLMMres
list       = make_list("RNA/GLMM_DEGs/GLMM_DEGs_2donors_ctrl_v1_TMM","_GLMMres.txt") 
list$KO    = take_factor(list$FILE,1,"_")
list$chunk = take_factor(list$FILE,3,"_")%>%as.numeric()
list$KO    = factor(list$KO,levels=KO_Gene_list)
list=list[order(list$chunk),]
list=list[order(list$KO),]
     
for(kkk in 1:nrow(list)){
     res_tmp=fread_FT(list$PATH[kkk])%>%filter(str_detect(parameter,pattern="typeKO|anova1"))
     if(kkk==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}
# 821698

gene_corresp=fread_FT("RNA/htseq/summary/cellranger_refdata-gex-GRCh38-2020-A_gene_corresp.txt")

# reorder columns
res_sum2=res_sum%>%pivot_wider(names_from="parameter",values_from=c("Full_Estimate","Full_SE","Full_Z","Full_P"))%>%
         inner_join(gene_corresp,.,by=c("gene_id"="gene"))%>%
         mutate(FDR =p.adjust(lrtP,method="BH"))%>%
         dplyr::select(1:4,6,8,10:11,12)
colnames(res_sum2)=c("gene_id","gene_symbol","KO","Estimate","se","Z","waldP","lrtP","FDR")

write.table_FT_2(res_sum2,paste0("RNA/GLMM_DEGs/summary/GLMM_DEGs_2donors_ctrl_v1_TMM_sum.txt"))
# 410849=17863*23


###########################
## Fisher test
###########################
suppressPackageStartupMessages(library(exact2x2))

res_sum2=fread_FT("RNA/GLMM_DEGs/summary/GLMM_DEGs_2donors_ctrl_v1_TMM_sum.txt")


## use NaiveCD4_5 signature genes
NaiveCD4_5_SEG_res=fread_FT("~/ws/2025/250317_SLE_maindata_new_clustering_CNA/SNN_clustering/SEGs_postmerge_v4/summary/cellstates_SEGs_v4_min.pct0.05_logFC0_SLEHC_full_v2.txt")%>%
              filter(cellstate=="NaiveCD4-5")%>%
              mutate(Signif=ifelse(abs(avg_log2FC)>0.1 & p.adjust_all<0.05,"Sig","None"))%>%
              dplyr::select(Gene,avg_log2FC_wo_psuedocount,p_val,Signif)
# 4960

res_compare=inner_join(NaiveCD4_5_SEG_res,res_sum2%>%dplyr::select(gene_symbol,KO,Estimate,se,Z,waldP,lrtP,FDR),by=c("Gene"="gene_symbol"))

# For visualization 
res_compare2=res_compare
res_compare2$avg_log2FC_wo_psuedocount[res_compare2$avg_log2FC_wo_psuedocount<(-2)]=-2
res_compare2$avg_log2FC_wo_psuedocount[res_compare2$avg_log2FC_wo_psuedocount>2]=2
res_compare2$Estimate[res_compare2$Estimate<(-2)]=-2
res_compare2$Estimate[res_compare2$Estimate>2]=2


# visualization 
res_compare3=res_compare2%>%
             mutate(Group=ifelse(Signif=="Sig" & avg_log2FC_wo_psuedocount>0 & Estimate<0 & FDR<0.05,"G1",
                          ifelse(Signif=="Sig" & avg_log2FC_wo_psuedocount>0,"G2",
                          ifelse(Estimate<0 & FDR<0.05,"G3",
                          "G4"))))

res_compare3$Group=factor(res_compare3$Group,levels=c("G1","G2","G3","G4"))
res_compare3=res_compare3[order(res_compare3$Group,decreasing=T),]


# Fisher res
# all TFs
res_compare3_lim=res_compare3%>%filter(Signif=="Sig")
KO_Gene_list_lim=intersect(KO_Gene_list,res_compare3_lim$Gene)

for(iii in 1:length(KO_Gene_list_lim)){
    KO_tmp=KO_Gene_list_lim[iii]

    res_compare3_tmp=res_compare3%>%filter(KO==KO_tmp)

    table_tmp=table(res_compare3_tmp$Group)

    mx1=matrix(c(table_tmp[1],table_tmp[2],table_tmp[3],table_tmp[4]),nrow=2,byrow=T)
    fisher1=fisher.exact(mx1,alternative="greater")
    
    res1=data.frame(KO=KO_tmp,OR=fisher1$estimate,pvalue=fisher1$p.value)

    if(iii==1){res1_sum=res1}else{res1_sum=rbind(res1_sum,res1)}
}
res1_sum=res1_sum%>%mutate(Bonf_P=pvalue*23) # 22(NaiveCD4-5)+1(NaiveCD4-4)

res1_sum_reorder=res1_sum[order(res1_sum$pvalue),]
res1_sum_reorder$KO=factor(res1_sum_reorder$KO,levels=rev(res1_sum_reorder$KO))

p = ggplot(data=res1_sum_reorder, aes(x=KO,y=-log10(pvalue)))+
    geom_bar(stat="identity")+
    geom_hline(yintercept=-log10(0.05/23),col="black",linewidth=0.5,linetype="dashed")+
    coord_flip()+
    theme_classic()+
    #scale_fill_npg()+
    theme(axis.text.x=element_text(colour="black",size=12),
          axis.text.y=element_text(colour="black",size=12),
          axis.title.x=element_text(colour="black",size=14),
          axis.title.y=element_text(colour="black",size=14),
          plot.title=element_blank())+
     labs(x="CRISPR-KO",y="-log10(enrichment P)\n(NaiveCD4-5 signature genes)")

pdf_3(paste0("RNA/GLMM_DEGs/fisher/GLMM_DEGs_2donors_ctrl_v1_TMM_NaiveCD4-5_fisher_bar.pdf"),h=4.5,w=4)
 plot(p)
dev.off()












