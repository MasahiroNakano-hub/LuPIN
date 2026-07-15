source("~/ws/my.source.R")
options(stringsAsFactors=F)

# 22 TF list
metadata=fread_FT("~/ws/2025/250725_SLE_CRISPR2/N003/Info/sampleID_corresp_v2_all.txt")%>%
         mutate(Donor_ID=paste0("Donor",take_factor(sample_ID,2,"_")))

KO_Gene_list=unique(metadata$KO_Gene)[1:22]%>%sort()
#  [1] "ARID5B" "ASH1L"  "BACH2"  "ETS1"   "FOXN3"  "FOXO1"  "FOXP1"  "HIVEP2"
#  [9] "KLF12"  "LEF1"   "NCOA2"  "NFKB1"  "NR3C2"  "RFX3"   "RORA"   "RUNX1" 
# [17] "SMYD3"  "STAT4"  "TCF12"  "TSHZ2"  "ZBTB20" "ZEB1"

KO_Gene_df=data.frame(TF=KO_Gene_list)



############################################
### 1. Misg C3 + KEGG
############################################

## C3 GTRD only
res_C3=fread_FT("summary/cellstates_SEGs_v4_min.pct0.05_logFC0.1_SLEHC_AllFDR0.05_ORA_C3GTRD_sum_all.txt")%>%
	   mutate(TF=gsub("_TARGET_GENES","",Description))%>%
	   filter(cellstate%in%c("NaiveCD4-4","NaiveCD4-5","CMCD4-4"))%>%
	   filter(TF%in%KO_Gene_list)%>%
	   mutate(Database="Msig C3GTRD")%>%
	   dplyr::select(TF,cellstate,Database,Description,pvalue)
# ARID5B, ASH1L, FOXN3, RORA, NCOA2

## KEGG
pathway_add=c("FoxO signaling pathway","NF-kappa B signaling pathway","JAK-STAT signaling pathway")

res_KEGG=fread_FT("summary/cellstates_SEGs_v4_min.pct0.05_logFC0.1_SLEHC_AllFDR0.05_ORA_HALLMARK_KEGG_sum_all.txt")%>%
	   	filter(cellstate%in%c("NaiveCD4-4","NaiveCD4-5","CMCD4-4"))%>%
	   	filter(Description%in%pathway_add)%>%
	   	mutate(TF=ifelse(Description=="FoxO signaling pathway","FOXO1",
	   				 ifelse(Description=="NF-kappa B signaling pathway","NFKB1",
	   				 ifelse(Description=="JAK-STAT signaling pathway","STAT4",NA
	   					))))%>%
	   	mutate(Database="KEGG")%>%
	   	dplyr::select(TF,cellstate,Database,Description,pvalue)

res_pathway=rbind(res_C3,res_KEGG)%>%
				mutate(log10P=-log10(pvalue))%>%
				mutate(Bonf_sig=ifelse(pvalue<0.05/16,"Sig","None"))


## Heatmap
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(viridis))
source("~/ws/colors.R")

col_fun = colorRamp2(c(0,5),c("#FEE0D2","#67000D"))

## NaiveCD4-5, CMCD4-4
res_pathway1=res_pathway%>%filter(cellstate%in%c("NaiveCD4-5","CMCD4-4"))%>%
				 dplyr::select(TF,cellstate,log10P)%>%
				 pivot_wider(names_from="cellstate",values_from="log10P")%>%
				 left_join(KO_Gene_df,.,by="TF")%>%
				 column_to_rownames("TF")

# no column annotation
p=Heatmap(res_pathway1,cluster_rows=FALSE,cluster_columns=FALSE,
                col=col_fun,na_col="#DCDCDC",
                height=unit(11,"cm"),width=unit(1,"cm"),rect_gp=gpar(col="gray"),
                row_names_gp=gpar(fontsize=16,fontface="italic"),row_names_side=c("left"),
                column_names_gp=gpar(fontsize=16),
        			 heatmap_legend_param=list(title="-log10(P)",at=c(0,5),
        			 						 title_gp=gpar(fontsize=0),labels_gp=gpar(fontsize=16))
          )

  
pdf_3(paste0("summary/NaiveCD4_5_CMCD4_4_22TF_ORA_C3_KEGG_logP.pdf"),h=6,w=4)
 draw(p,heatmap_legend_side="right")
dev.off()


############################################
### 2. decoupleR
############################################

suppressPackageStartupMessages(library(decoupleR))

##  cell-state SEG full stats (passed low exp filter)
res_sum=fread_FT("summary/cellstates_SEGs_v4_min.pct0.05_logFC0_SLEHC_full_v2.txt")
# 674687

res_sum_NaiveCD4_4=res_sum%>%filter(cellstate=="NaiveCD4-4")%>%
					mutate(SEG=ifelse(avg_log2FC>0.1 & p.adjust_all<0.05, 1, 0))%>%
					dplyr::select(Gene,avg_log2FC,p_val,SEG)%>%
					mutate(logP=-log10(p_val))%>%
					mutate(logP=ifelse(is.infinite(logP),310,logP))%>%
					mutate(multiply=avg_log2FC*logP)%>%
					column_to_rownames("Gene")%>%as.matrix()
# 4582

res_sum_NaiveCD4_5=res_sum%>%filter(cellstate=="NaiveCD4-5")%>%
					mutate(SEG=ifelse(avg_log2FC>0.1 & p.adjust_all<0.05, 1, 0))%>%
					dplyr::select(Gene,avg_log2FC,p_val,SEG)%>%
					mutate(logP=-log10(p_val))%>%
					mutate(logP=ifelse(is.infinite(logP),310,logP))%>%
					mutate(multiply=avg_log2FC*logP)%>%
					column_to_rownames("Gene")%>%as.matrix()
# 4960

res_sum_CMCD4_4=res_sum%>%filter(cellstate=="CMCD4-4")%>%
					mutate(SEG=ifelse(avg_log2FC>0.1 & p.adjust_all<0.05, 1, 0))%>%
					dplyr::select(Gene,avg_log2FC,p_val,SEG)%>%
					mutate(logP=-log10(p_val))%>%
					mutate(logP=ifelse(is.infinite(logP),310,logP))%>%
					mutate(multiply=avg_log2FC*logP)%>%
					column_to_rownames("Gene")%>%as.matrix()
# 5072


## decoupleR
net_dorothea = decoupleR::get_dorothea(organism = 'human', levels = c('A', 'B', 'C'))

## NaiveCD4-5
NaiveCD4_5_dorothea_multiply = decoupleR::run_ulm(mat = res_sum_NaiveCD4_5[,"multiply"], 
                                   				net = net_dorothea, 
                                   				minsize = 5)


## CMCD4-4
CMCD4_4_dorothea_multiply = decoupleR::run_ulm(mat = res_sum_CMCD4_4[,"multiply"], 
                                   				net = net_dorothea, 
                                   				minsize = 5)


## Heatmap
col_fun = colorRamp2(c(0,5),c("#DEEBF7","#08306B"))

## NaiveCD4-5, CMCD4-4
NaiveCD4_5_dorothea_lim=NaiveCD4_5_dorothea_multiply%>%filter(source%in%KO_Gene_list)%>%
								mutate(log10P=-log10(p_value))%>%
								mutate(cellstate="NaiveCD4-5")

CMCD4_4_dorothea_lim   =CMCD4_4_dorothea_multiply%>%filter(source%in%KO_Gene_list)%>%
								mutate(log10P=-log10(p_value))%>%
							   mutate(cellstate="CMCD4-4")

dorothea_lim_1=rbind(NaiveCD4_5_dorothea_lim,CMCD4_4_dorothea_lim)%>%
					dplyr::select(source,cellstate,log10P)%>%
					pivot_wider(names_from="cellstate",values_from="log10P")%>%
				 	left_join(KO_Gene_df,.,by=c("TF"="source"))%>%
				 	column_to_rownames("TF")

# no column annotation
p=Heatmap(dorothea_lim_1,cluster_rows=FALSE,cluster_columns=FALSE,
                col=col_fun,na_col="#DCDCDC",
                height=unit(11,"cm"),width=unit(1,"cm"),rect_gp=gpar(col="gray"),
                row_names_gp=gpar(fontsize=16,fontface="italic"),row_names_side=c("left"),
                column_names_gp=gpar(fontsize=16),
        			 heatmap_legend_param=list(title="-log10(P)",at=c(0,5),
        			 						 title_gp=gpar(fontsize=0),labels_gp=gpar(fontsize=16))
          )

  
pdf_3(paste0("summary/NaiveCD4-5_CMCD4-4_multiply_res_logP.pdf"),h=6,w=4)
 draw(p,heatmap_legend_side="right")
dev.off()



############################################
### 3. Blood-specific ChIP data
############################################

########################
## 3-0. Manual download
## from ChIP Atlas
########################
ll ~/data/ChIP_Atlas

# 13 TFs
# No Blood ChIP data for FOXN3, HIVEP2, KLF12, NR3C2, RFX3, RORA, SMYD3, TSHZ2, ZBTB20
# Few peaks: SRX128782~SRX128784 (RUNX1, Jurkat), RX041288 (STAT4, Tcell)


########################
## 3-1. Make target gene list from each TF
########################

## 3-1-1. TSS to promotor
gtf=~/tool/ref/refdata-gex-GRCh38-2020-A/genes/genes.gtf

mkdir -p ~/tool/ref/refdata-gex-GRCh38-2020-A/bed

# exon tag basic only
cat $gtf |awk 'BEGIN{FS="\t";OFS="\t"}{
    if($3=="transcript" && $9 ~ /tag "basic"/){
      if (match($9, /gene_name "[^"]+"/)) {
       gene = substr($9, RSTART + 11, RLENGTH - 12);
        if ($7 == "+") {
            tss = $4;
         } else {
            tss = $5;
         }
         print $1, tss, $7, gene
      }
    }
}' > tss_temp.txt # 82823

# pick up most upstream TSS for each gene
awk '{if($3=="+") print $0 }' tss_temp.txt | sort -k4,4 -k2,2n | awk '!seen[$4]++' > tss_top.txt
awk '{if($3=="-") print $0 }' tss_temp.txt | sort -k4,4 -k2,2nr | awk '!seen[$4]++' >> tss_top.txt # 36593

# define promotor as TSS -5kb +2kb
cat tss_top.txt | awk 'BEGIN{FS="\t";OFS="\t"}{
    p_start = $2 - 5000;
    p_end = $2 + 2000;

    if(p_start < 0){p_start = 0}
    
	 if($1 ~ /^chr[1-9][0-9]?$/ || $1 == "chrX" || $1 == "chrY" || $1 == "chrM"){
        print $1, p_start, p_end, $4
    }
}' | sort -k1,1 -k2,2n > ~/tool/ref/refdata-gex-GRCh38-2020-A/bed/basic_unique_promoters_minus5kb_plus2kb.bed
# 36564, only auto chr or X,Y,M

rm tss_temp.txt
rm tss_top.txt


## 3-1-2. BED intersect
module use /usr/local/package/modulefiles
module load bedtools/2.31.1

PROMOTER_BED=~/tool/ref/refdata-gex-GRCh38-2020-A/bed/basic_unique_promoters_minus5kb_plus2kb.bed

mkdir -p SNN_clustering/SEGs_postmerge_v4_ChIP_Atlas/Target_gene

find ~/data/ChIP_Atlas -name "*.bed" | while read -r CHIP_FILE; do
    
    TF_NAME=$(basename "$(dirname "$CHIP_FILE")")
    FILENAME=$(basename "$CHIP_FILE" .bed)
    
    echo "Processing ${TF_NAME} (${FILENAME})..."

 bedtools intersect -a "$PROMOTER_BED" -b "$CHIP_FILE" -u | \
    awk '{print $4}' | sort | uniq > SNN_clustering/SEGs_postmerge_v4_ChIP_Atlas/Target_gene/${TF_NAME}_${FILENAME}_genes.txt

    wc -l SNN_clustering/SEGs_postmerge_v4_ChIP_Atlas/Target_gene/${TF_NAME}_${FILENAME}_genes.txt

done


########################
## 3-2. dataset QC
## Peaks overlapped with Promotor (100-2000) 
########################

module load R/4.4.0 
R
source("~/ws/my.source.R")
options(stringsAsFactors=F)

list=make_list("SNN_clustering/SEGs_postmerge_v4_ChIP_Atlas/Target_gene","_genes.txt")
list$TF　=take_factor(list$FILE,1,"_")
list$cond=take_factor(list$FILE,2,"_")
list$ID  =take_factor(list$cond,1,"\\.")
list$Qval=take_factor(list$cond,2,"\\.")
# 65

unique(list$ID)%>%length()
# 43

list$genenum=NA

for(iii in 1:nrow(list)){
	gene_tmp=fread_FF(list$PATH[iii])
	gene_num=nrow(gene_tmp)
	list$genenum[iii]=gene_num
}

list2=list%>%filter(genenum>100 & genenum <2000)%>%.[!duplicated(.$ID),]
# 27
unique(list2$TF)%>%length()
# 13
# OK, at least one dataset retained for each TF

# Excluded IDs
setdiff(unique(list$ID),list2$ID)
#  [1] "SRX015825"   "SRX2900588"  "SRX2900589"  "SRX12209237" "SRX12209238"
#  [6] "SRX12209239" "SRX12209240" "SRX12209241" "SRX12209242" "SRX015826"  
# [11] "SRX1041800"  "SRX1492212"  "SRX2014504"  "SRX151947"   "SRX151948"  
# [16] "SRX151949"
# Too many: SRX015825 (ETS1, Jurkat), SRX1041800,SRX1492212,SRX2014504 (RUNX1, Jurkat)
# Too few:  SRX2900588~89 (NCOA2, THP-1), SRX12209237~42 (NFKB1, CD4T), SRX015826 (RUNX1, Jurkat), SRX151947-49 (STAT4, PBMC)

# write.table_FT_2(list2,paste0("summary/ChIP_Atlas_QCpass_dataset_list.txt"))


########################
## 3-3. Fisher
## NaiveCD4-4, NaiveCD4-5, CMCD4-4 SEGs
########################

suppressPackageStartupMessages(library(exact2x2))

# SEG full stats (passed low exp filter)
res_sum=fread_FT("summary/cellstates_SEGs_v4_min.pct0.05_logFC0_SLEHC_full_v2.txt")
# 674687

res_sum_NaiveCD4_4=res_sum%>%filter(cellstate=="NaiveCD4-4")
# 4582
res_sum_NaiveCD4_4_sig =res_sum_NaiveCD4_4%>%filter(avg_log2FC>0.1 & p.adjust_all<0.05)%>%pull(Gene)
# 106
res_sum_NaiveCD4_4_none=res_sum_NaiveCD4_4%>%filter(!avg_log2FC>0.1 | !p.adjust_all<0.05)%>%pull(Gene)
# 4476

res_sum_NaiveCD4_5=res_sum%>%filter(cellstate=="NaiveCD4-5")
# 4960
res_sum_NaiveCD4_5_sig =res_sum_NaiveCD4_5%>%filter(avg_log2FC>0.1 & p.adjust_all<0.05)%>%pull(Gene)
# 1367
res_sum_NaiveCD4_5_none=res_sum_NaiveCD4_5%>%filter(!avg_log2FC>0.1 | !p.adjust_all<0.05)%>%pull(Gene)
# 3593

res_sum_CMCD4_4=res_sum%>%filter(cellstate=="CMCD4-4")
# 5072
res_sum_CMCD4_4_sig =res_sum_CMCD4_4%>%filter(avg_log2FC>0.1 & p.adjust_all<0.05)%>%pull(Gene)
# 1278
res_sum_CMCD4_4_none=res_sum_CMCD4_4%>%filter(!avg_log2FC>0.1 | !p.adjust_all<0.05)%>%pull(Gene)
# 3794


# Fisher res
# NaiveCD4_5
for(iii in 1:nrow(list2)){
    gene_tmp=fread_FF(list2$PATH[iii])%>%pull(V1)
    TF_tmp  =list2$TF[iii]
    ID_tmp  =list2$ID[iii]
    Qval_tmp=list2$Qval[iii]
    gene_num=list2$genenum[iii]

    G1=intersect(res_sum_NaiveCD4_5_sig,gene_tmp)%>%length()
    G2=intersect(res_sum_NaiveCD4_5_none,gene_tmp)%>%length()
    G3=setdiff(res_sum_NaiveCD4_5_sig,gene_tmp)%>%length()
    G4=setdiff(res_sum_NaiveCD4_5_none,gene_tmp)%>%length()

    mx2=matrix(c(G1,G2,G3,G4),nrow=2,byrow=T)
    fisher2=fisher.exact(mx2,alternative="greater")
    
    res2=data.frame(TF=TF_tmp,dataset_ID=ID_tmp,Qval_cutoff=Qval_tmp,genenum=gene_num,OR=fisher2$estimate,pvalue=fisher2$p.value)

    if(iii==1){res2_sum=res2}else{res2_sum=rbind(res2_sum,res2)}
}

# prioritize one dataset for each TF (manual)
res2_lim=res2_sum%>%filter(genenum>200)%>%.[c(1,2,3,7,8,10,12,14,15,18,20,21,23),]%>%
			mutate(log10P=-log10(pvalue))%>%
			mutate(Bonf_sig=ifelse(pvalue<0.05/13,"Sig","None"))%>%
			mutate(cellstate="NaiveCD4-5")


# CMCD4_4
for(iii in 1:nrow(list2)){
    gene_tmp=fread_FF(list2$PATH[iii])%>%pull(V1)
    TF_tmp  =list2$TF[iii]
    ID_tmp  =list2$ID[iii]
    Qval_tmp=list2$Qval[iii]
    gene_num=list2$genenum[iii]

    G1=intersect(res_sum_CMCD4_4_sig,gene_tmp)%>%length()
    G2=intersect(res_sum_CMCD4_4_none,gene_tmp)%>%length()
    G3=setdiff(res_sum_CMCD4_4_sig,gene_tmp)%>%length()
    G4=setdiff(res_sum_CMCD4_4_none,gene_tmp)%>%length()

    mx3=matrix(c(G1,G2,G3,G4),nrow=2,byrow=T)
    fisher3=fisher.exact(mx3,alternative="greater")
    
    res3=data.frame(TF=TF_tmp,dataset_ID=ID_tmp,Qval_cutoff=Qval_tmp,genenum=gene_num,OR=fisher3$estimate,pvalue=fisher3$p.value)

    if(iii==1){res3_sum=res3}else{res3_sum=rbind(res3_sum,res3)}
}

# consistent dataset with NaiveCD4_5
res3_lim=res3_sum%>%filter(dataset_ID%in%res2_lim$dataset_ID)%>%
			mutate(log10P=-log10(pvalue))%>%
			mutate(Bonf_sig=ifelse(pvalue<0.05/13,"Sig","None"))%>%
			mutate(cellstate="CMCD4-4")


########################
## 3-4. Figure
########################

col_fun = colorRamp2(c(0,5),c("#E5F5E0", "#00441B"))

## NaiveCD4-5, CMCD4-4
res_chip1=rbind(res2_lim,res3_lim)%>%
				 dplyr::select(TF,cellstate,log10P)%>%
				 pivot_wider(names_from="cellstate",values_from="log10P")%>%
				 left_join(KO_Gene_df,.,by="TF")%>%
				 column_to_rownames("TF")

# no column annotation
p=Heatmap(res_chip1,cluster_rows=FALSE,cluster_columns=FALSE,
                col=col_fun,na_col="#DCDCDC",
                height=unit(11,"cm"),width=unit(1,"cm"),rect_gp=gpar(col="gray"),
                row_names_gp=gpar(fontsize=16,fontface="italic"),row_names_side=c("left"),
                column_names_gp=gpar(fontsize=16),
        			 heatmap_legend_param=list(title="-log10(P)",at=c(0,5),
        			 						 title_gp=gpar(fontsize=0),labels_gp=gpar(fontsize=16))
        )

pdf_3(paste0("summary/NaiveCD4_5_CMCD4_4_22TF_Fisher_res_logP.pdf"),h=6,w=4)
 draw(p,heatmap_legend_side="right")
dev.off()


