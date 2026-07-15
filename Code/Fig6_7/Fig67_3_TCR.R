
################################################################
### V gene level analysis
################################################################
R
source("~/ws/my.source.R")
options(stringsAsFactors=F)

TCR_sum3_QC_beta_clonotype=fread_FT("TCR/cellranger_vdj_vdj_t_filtered_contig_annotations_beta_new_clonotype_meta.txt")
TCR_sum3_QC_alpha_clonotype=fread_FT("TCR/cellranger_vdj_vdj_t_filtered_contig_annotations_alpha_new_clonotype_meta.txt")

# limit to CD4 and CD8 celltypes
use_celltypes=c("NaiveCD4","CMCD4","EMCD4","CTLCD4","Treg",
                "NaiveCD8","CMCD8","EMCD8")

TCR_sum3_QC_beta_clonotype_lim=TCR_sum3_QC_beta_clonotype%>%
                               filter(celltype%in%use_celltypes)
# 74727

TCR_sum3_QC_alpha_clonotype_lim=TCR_sum3_QC_alpha_clonotype%>%
                               filter(celltype%in%use_celltypes)
# 65768


########################
## V gene usage fisher test
## cell state level (across donors)
## within celltype
## use both beta and alpha info
########################
suppressPackageStartupMessages(library(exact2x2))

# beta
for(iii in 1:length(use_celltypes)){
    celltype_tmp=use_celltypes[iii]

    TCR_sum3_QC_beta_clonotype_tmp=TCR_sum3_QC_beta_clonotype_lim%>%
                                   filter(celltype==celltype_tmp)

    cellstate_list=unique(TCR_sum3_QC_beta_clonotype_tmp$cellstate)%>%sort()
    #v_gene_list=unique(TCR_sum3_QC_beta_clonotype_tmp$v_gene)%>%sort()
    v_gene_list2=table(TCR_sum3_QC_beta_clonotype_tmp$v_gene)[table(TCR_sum3_QC_beta_clonotype_tmp$v_gene)>=10]%>%
                 names()%>%sort() # low count Vgene filter

 for(kkk in 1:length(cellstate_list)){
    cellstate_tmp=cellstate_list[kkk]

  for(mmm in 1:length(v_gene_list2)){
    v_gene_tmp=v_gene_list2[mmm]

    TCR_sum3_QC_beta_clonotype_tmp2=TCR_sum3_QC_beta_clonotype_tmp%>%
                                    mutate(cellstate2=ifelse(cellstate==cellstate_tmp,"cellstate_target","cellstate_no"))%>%
                                    mutate(v_gene2=ifelse(v_gene==v_gene_tmp,"vgene_target","vgene_no"))

    TCR_sum3_QC_beta_clonotype_tmp2$cellstate2=factor(TCR_sum3_QC_beta_clonotype_tmp2$cellstate2,levels=c("cellstate_target","cellstate_no"))
    TCR_sum3_QC_beta_clonotype_tmp2$v_gene2   =factor(TCR_sum3_QC_beta_clonotype_tmp2$v_gene2,levels=c("vgene_target","vgene_no"))

    table_tmp=table(TCR_sum3_QC_beta_clonotype_tmp2$v_gene2,TCR_sum3_QC_beta_clonotype_tmp2$cellstate2)

    fisher_res=fisher.exact(table_tmp)

    fisher_tmp=data.frame(celltype=celltype_tmp,cellstate=cellstate_tmp,v_gene=v_gene_tmp,
                          OR=fisher_res$estimate,pvalue=fisher_res$p.value)

    if(mmm==1){fisher_sum=fisher_tmp}else{fisher_sum=rbind(fisher_sum,fisher_tmp)}  
   }
   if(kkk==1){fisher_sum2=fisher_sum}else{fisher_sum2=rbind(fisher_sum2,fisher_sum)}  
  }
 if(iii==1){fisher_sum3_beta=fisher_sum2}else{fisher_sum3_beta=rbind(fisher_sum3_beta,fisher_sum2)}  
}


# alpha
for(iii in 1:length(use_celltypes)){
    celltype_tmp=use_celltypes[iii]

    TCR_sum3_QC_alpha_clonotype_tmp=TCR_sum3_QC_alpha_clonotype_lim%>%
                                   filter(celltype==celltype_tmp)

    cellstate_list=unique(TCR_sum3_QC_alpha_clonotype_tmp$cellstate)%>%sort()
    #v_gene_list=unique(TCR_sum3_QC_alpha_clonotype_tmp$v_gene)%>%sort()
    v_gene_list2=table(TCR_sum3_QC_alpha_clonotype_tmp$v_gene)[table(TCR_sum3_QC_alpha_clonotype_tmp$v_gene)>=10]%>%
                 names()%>%sort()

 for(kkk in 1:length(cellstate_list)){
    cellstate_tmp=cellstate_list[kkk]

  for(mmm in 1:length(v_gene_list2)){
    v_gene_tmp=v_gene_list2[mmm]

    TCR_sum3_QC_alpha_clonotype_tmp2=TCR_sum3_QC_alpha_clonotype_tmp%>%
                                    mutate(cellstate2=ifelse(cellstate==cellstate_tmp,"cellstate_target","cellstate_no"))%>%
                                    mutate(v_gene2=ifelse(v_gene==v_gene_tmp,"vgene_target","vgene_no"))

    TCR_sum3_QC_alpha_clonotype_tmp2$cellstate2=factor(TCR_sum3_QC_alpha_clonotype_tmp2$cellstate2,levels=c("cellstate_target","cellstate_no"))
    TCR_sum3_QC_alpha_clonotype_tmp2$v_gene2   =factor(TCR_sum3_QC_alpha_clonotype_tmp2$v_gene2,levels=c("vgene_target","vgene_no"))

    table_tmp=table(TCR_sum3_QC_alpha_clonotype_tmp2$v_gene2,TCR_sum3_QC_alpha_clonotype_tmp2$cellstate2)

    fisher_res=fisher.exact(table_tmp)

    fisher_tmp=data.frame(celltype=celltype_tmp,cellstate=cellstate_tmp,v_gene=v_gene_tmp,
                          OR=fisher_res$estimate,pvalue=fisher_res$p.value)

    if(mmm==1){fisher_sum=fisher_tmp}else{fisher_sum=rbind(fisher_sum,fisher_tmp)}  
   }
   if(kkk==1){fisher_sum2=fisher_sum}else{fisher_sum2=rbind(fisher_sum2,fisher_sum)}  
  }
 if(iii==1){fisher_sum3_alpha=fisher_sum2}else{fisher_sum3_alpha=rbind(fisher_sum3_alpha,fisher_sum2)}  
}


fisher_sum3=rbind(fisher_sum3_beta,fisher_sum3_alpha)%>%
            mutate(FDR=p.adjust(pvalue,method="BH"))
# 4143

fisher_sum3%>%filter(FDR<0.05)%>%nrow()
# 188

fisher_sum3$cellstate%>%unique()%>%length()
# 50 cell states
# 188/50=3.76

fisher_sum3_sig=fisher_sum3%>%filter(FDR<0.05)
write.table_FT_2(fisher_sum3_sig,paste0("TCR/Vgene/CD4_CD8_Vgene_fisher_res_sum_eachcellstate_10countfilter_FDR0.05.txt"))
# 188


################
## Representative clones Barplot
################

fisher_sum3_sig=fread_FT("TCR/Vgene/CD4_CD8_Vgene_fisher_res_sum_eachcellstate_10countfilter_FDR0.05.txt")

## "TRBV7-9" and "TRBV29-1" in NaiveCD4
celltype_tmp="NaiveCD4"

TCR_sum3_QC_beta_clonotype_tmp=TCR_sum3_QC_beta_clonotype_lim%>%
                               filter(celltype==celltype_tmp)

TCR_sum3_QC_beta_clonotype_tmp2=TCR_sum3_QC_beta_clonotype_tmp%>%
                                mutate(v_gene2=ifelse(v_gene%in%c("TRBV7-9","TRBV29-1"),v_gene,"others"))

TCR_sum3_QC_beta_clonotype_tmp2$v_gene2   =factor(TCR_sum3_QC_beta_clonotype_tmp2$v_gene2,levels=c("TRBV7-9","TRBV29-1","others"))

table_tmp=table(TCR_sum3_QC_beta_clonotype_tmp2$v_gene2,TCR_sum3_QC_beta_clonotype_tmp2$cellstate)

prop_tmp=apply(table_tmp,2,function(x){x/sum(x)*100})%>%as.data.frame()%>%
         rownames_to_column("v_gene2")%>%
         pivot_longer(cols=-v_gene2,names_to="cellstate",values_to="prop")

prop_tmp$v_gene2=factor(prop_tmp$v_gene2,levels=c("TRBV7-9","TRBV29-1","others"))
prop_tmp$cellstate=take_factor(prop_tmp$cellstate,2,"-")

prop_tmp2=prop_tmp%>%filter(!v_gene2=="others")

p = ggplot(data=prop_tmp2,aes(x=cellstate,y=prop))+
    geom_bar(aes(fill=v_gene2),stat="identity",position="stack")+
    scale_fill_manual(values = c(pal_nejm("default")(2),"#cccccc")) +
    theme_classic()+
    theme(axis.text.x=element_text(colour="black",size=14),
          axis.text.y=element_text(colour="black",size=12),
          axis.title.x=element_text(colour="black",size=14),
          axis.title.y=element_text(colour="black",size=14),
          plot.title=element_blank(),
          legend.position="right",
          legend.title=element_blank(),
          legend.text=element_text(colour="black",size=14))+
    labs(x="NaiveCD4 cell state",y="Percentage of cells")

pdf_3("TCR/Vgene/NaiveCD4_TRBV7-9_TRBV29-1_prop_v2.pdf",h=2.5,w=4)
 plot(p)
dev.off()



################################################################
### Clonotype level analysis
## beta chain: CDR3 (amino acid) + V gene
################################################################
R
source("~/ws/my.source.R")
options(stringsAsFactors=F)

# limit to CD4 and CD8 celltypes
use_celltypes=c("NaiveCD4","CMCD4","EMCD4","CTLCD4","Treg",
                "NaiveCD8","CMCD8","EMCD8")

TCR_sum3_QC_beta_clonotype_lim=TCR_sum3_QC_beta_clonotype%>%
                               filter(celltype%in%use_celltypes)
# 74727


########################
## count clonotypes in EMCD8-5
########################
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(viridis))
source("~/ws/colors.R")

cellstate_tmp="EMCD8-5"

res_tmp=TCR_sum3_QC_beta_clonotype%>%filter(cellstate==cellstate_tmp)

donor_list=res_tmp$Donor_ID%>%unique()%>%sort()

## 3-2-1. for each donor count
for(kkk in 1:length(donor_list)){
    donor_tmp=donor_list[kkk]

    res_tmp2=res_tmp%>%filter(Donor_ID==donor_tmp)

    table_tmp=table(res_tmp2$clonotype_new)%>%as.data.frame()
    colnames(table_tmp)=c("clonotype_new","clonenum")

    res_tmp3=left_join(res_tmp2,table_tmp,by="clonotype_new")%>%
             mutate(clonesize=ifelse(clonenum>10,"Verylarge",
                              ifelse(clonenum>5,"Large",
                              ifelse(clonenum>2,"Medium",
                              ifelse(clonenum==2,"Small",
                              ifelse(clonenum==1,"Single",
                              NA))))))
    
    res_tmp3$clonesize=factor(res_tmp3$clonesize,levels=c("Single","Small","Medium","Large","Verylarge"))
    table_tmp2=table_freq(res_tmp3$clonesize)%>%
               mutate(Donor_ID=donor_tmp)%>%
               dplyr::select(Donor_ID,factor,Freq)
    colnames(table_tmp2)[2:3]=c("clonesize","Cellnum")

    if(kkk==1){table_sum=table_tmp2}else{table_sum=rbind(table_sum,table_tmp2)}  
    if(kkk==1){res_sum=res_tmp3}else{res_sum=rbind(res_sum,res_tmp3)}  
}
write.table_FT_2(res_sum,paste0("TCR/clonesize/EMCD8-5_eachclone_colonesize_meta.txt"))


## Prop bar plot
# limit to donors with > 10 clones
# reordered by "single"
donors_lim=table(res_sum$Donor_ID)[table(res_sum$Donor_ID)>10]%>%names()
table_lim=table_sum%>%filter(Donor_ID%in%donors_lim)
prop_lim=table_lim%>%group_by(Donor_ID)%>%reframe(prop=Cellnum/sum(Cellnum)*100)
table_lim2=cbind(table_lim,prop_lim[,2,drop=F])

single_lim2=table_lim2%>%filter(clonesize=="Single")
donor_order=single_lim2[order(single_lim2$prop,decreasing=T),]%>%pull(Donor_ID)

table_lim2$Donor_ID =factor(table_lim2$Donor_ID,levels=donor_order)
levels(table_lim2$Donor_ID)=gsub("Donor_","SLE",levels(table_lim2$Donor_ID))
table_lim2$clonesize=factor(table_lim2$clonesize,levels=rev(c("Single","Small","Medium","Large","Verylarge")))
levels(table_lim2$clonesize)=rev(c("Single (1)","Small (2)","Medium (3-5)","Large (6-10)","Very large (>10)"))

p = ggplot(data=table_lim2,aes(x=Donor_ID,y=prop))+
    geom_bar(aes(fill=clonesize),stat="identity",position="stack")+
    theme_classic()+
    scale_fill_viridis(option="plasma",discrete=T,direction=-1) +
    theme(axis.text.x=element_text(colour="black",,angle=45,hjust=1,vjust=1,size=12),
          axis.text.y=element_text(colour="black",size=12),
          axis.title.x=element_blank(),
          axis.title.y=element_text(colour="black",size=14),
          plot.title=element_blank(),
          legend.position="right",
          legend.title=element_text(colour="black",size=14),
          legend.text=element_text(colour="black",size=14))+
    labs(y="Proportion of cells",fill="Clone size")

pdf_3("TCR/clonesize/EMCD8-5_colonesize_bydonor_>10clones_prop.pdf",h=2.5,w=7)
 plot(p)
dev.off()



########################
## check expanded clones specific to EMCD8-5
########################
suppressPackageStartupMessages(library(exact2x2))

# limit to "Medium","Large","Verylarge"
expanded_clones=res_sum%>%filter(clonesize%in%c("Medium","Large","Verylarge"))%>%dplyr::select(Donor_ID,clonotype_new,clonenum,clonesize)%>%unique()
# 37

## Check alpha chain
for(iii in 1:nrow(expanded_clones)){
    donor_tmp=expanded_clones$Donor_ID[iii]
    clone_tmp=expanded_clones$clonotype_new[iii]

    cell_tmp=TCR_sum3_QC_beta_clonotype_lim%>%
              filter(Donor_ID==donor_tmp & clonotype_new==clone_tmp & cellstate=="EMCD8-5")%>%
              pull(barcode)

    alpha_tmp=TCR_sum3_QC_alpha_clonotype_lim%>%
              filter(barcode%in%cell_tmp)

    alpha_clone=alpha_tmp$clonotype_new%>%unique()%>%paste(.,collapse=", ")

    corresp_tmp=data.frame(Donor_ID=donor_tmp,clonotype=clone_tmp,clonotype_alpha=alpha_clone)    

    if(iii==1){corresp_sum=corresp_tmp}else{corresp_sum=rbind(corresp_sum,corresp_tmp)}  
}


## Fisher test
## two-sided, consistent with V gene analysis
for(iii in 1:nrow(expanded_clones)){
    donor_tmp=expanded_clones$Donor_ID[iii]
    clone_tmp=expanded_clones$clonotype_new[iii]
    num_tmp  =expanded_clones$clonenum[iii]
    size_tmp =expanded_clones$clonesize[iii]

    paste0(donor_tmp,": ",clone_tmp)%>%print()

    check_tmp=TCR_sum3_QC_beta_clonotype_lim%>%
              filter(celltype=="EMCD8")%>%
              filter(Donor_ID==donor_tmp)%>%
              mutate(cellstate_check=ifelse(cellstate=="EMCD8-5","EMCD8-5","Others"))%>%
              mutate(clonotype_check=ifelse(clonotype_new==clone_tmp,clone_tmp,"Others"))

    table_check=table(check_tmp$cellstate_check,check_tmp$clonotype_check)
    fisher_res=fisher.exact(table_check)

    fisher_tmp=data.frame(Donor_ID=donor_tmp,clonotype=clone_tmp,clonenum=num_tmp,clonesize=size_tmp,
                          OR=fisher_res$estimate,pvalue=fisher_res$p.value)

    if(iii==1){fisher_sum=fisher_tmp}else{fisher_sum=rbind(fisher_sum,fisher_tmp)}  
}
fisher_sum$FDR=p.adjust(fisher_sum$pvalue,method="BH")

# add alpha chain info
fisher_sum2=left_join(corresp_sum,fisher_sum,by=c("Donor_ID","clonotype"))
write.table_FT_2(fisher_sum2,paste0("TCR/clonesize/EMCD8-5_expandedclone_fisher_twosided_res.txt"))



##############################################################################################################################
##############################################################################################################################

sessionInfo()

R version 4.4.0 (2024-04-24)
Platform: x86_64-pc-linux-gnu
Running under: Red Hat Enterprise Linux 8.8 (Ootpa)

Matrix products: default
BLAS:   /usr/local/package/r/4.4.0/lib64/R/lib/libRblas.so 
LAPACK: /usr/local/package/r/4.4.0/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0

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
 [1] viridis_0.6.5      viridisLite_0.4.2  ggsci_3.2.0        RColorBrewer_1.1-3
 [5] exact2x2_1.6.9     exactci_1.4-4      testthat_3.2.3     ssanv_1.1         
 [9] data.table_1.16.4  lubridate_1.9.4    forcats_1.0.0      stringr_1.5.1     
[13] dplyr_1.1.4        purrr_1.0.2        readr_2.1.5        tidyr_1.3.1       
[17] tibble_3.2.1       ggplot2_3.5.2      tidyverse_2.0.0   

loaded via a namespace (and not attached):
 [1] gtable_0.3.6     compiler_4.4.0   brio_1.1.5       tidyselect_1.2.1
 [5] gridExtra_2.3    scales_1.3.0     R6_2.5.1         generics_0.1.3  
 [9] munsell_0.5.1    pillar_1.10.1    tzdb_0.4.0       rlang_1.1.5     
[13] stringi_1.8.4    timechange_0.3.0 cli_3.6.3        withr_3.0.2     
[17] magrittr_2.0.3   grid_4.4.0       hms_1.1.3        lifecycle_1.0.4 
[21] vctrs_0.6.5      glue_1.8.0       colorspace_2.1-1 tools_4.4.0     
[25] pkgconfig_2.0.3 













