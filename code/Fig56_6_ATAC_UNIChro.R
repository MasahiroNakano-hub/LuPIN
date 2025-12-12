source("~/ws/my.source.R")
options(stringsAsFactors=F)


###################
# 1. ATAC GLMM
###################
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(lme4))


chunk_tmp=1



metadata=fread_FT("../Info/FACS_metadata.txt")

# lowexp filter TMMnorm
dge_n   = readRDS("ATAC/ATAC_count/summary/count_sum_lowexp0.3_TMM.rds")

# GLMM
offset_tmp = dge_n$samples %>% mutate(norm_libsize=lib.size*norm.factors) %>%
             rownames_to_column("sample_ID") %>% .[,c(1,5)]

# 100 genes per chunk
start=(100*(chunk_tmp-1)+1)
if(100*chunk_tmp>nrow(dge_n$counts)){end=nrow(dge_n$counts)}else{end=100*chunk_tmp}


for(iii in start:end){
    gene_tmp =rownames(dge_n$counts)[iii]
    paste0(iii,";",gene_tmp)%>%print()

    count_tmp=dge_n$counts[gene_tmp,]%>%as.data.frame()%>%rownames_to_column("sample_ID")
    colnames(count_tmp)[2] ="exp"

    data_tmp=count_tmp%>%
             mutate(Library_Sample_ID=take_factor(sample_ID,3,"-"))%>%
             mutate(subset=take_factor(sample_ID,4,"-"))%>%
             left_join(.,metadata,by="Library_Sample_ID")%>%
             left_join(.,offset_tmp,by="sample_ID")

    data_tmp$subset=factor(data_tmp$subset,levels=c("C","T"))


   T1 = try(glmer.nb(exp~subset+(1|Library_Sample_ID)+Age+Sex,offset=log(norm_libsize),data_tmp,control = glmerControl(optimizer = "bobyqa")))
   T0 = try(glmer.nb(exp~       (1|Library_Sample_ID)+Age+Sex,offset=log(norm_libsize),data_tmp,control = glmerControl(optimizer = "bobyqa")))

    if(class(T1)=="try-error"|class(T0)=="try-error"){
      res_tmp = data.frame(parameter=c("(Intercept)","subsetT",
                                       "Age","SexM","anova1"),
                           Full_Estimate=NA,Full_SE=NA,Full_Z=NA,Full_P=NA,
                           gene=gene_tmp)
    }else{

    anova1= anova(T1,T0)[2,8]

    res_tmp1=data.frame(summary(T1)$coefficients)%>%rownames_to_column("parameter")
    colnames(res_tmp1)[2:5]=c("Full_Estimate","Full_SE","Full_Z","Full_P")

    anova1_tmp=data.frame(parameter="anova1",Full_Estimate=NA,Full_SE=NA,Full_Z=NA,Full_P=anova1)

    res_tmp=bind_rows(res_tmp1,anova1_tmp)%>%
            mutate(gene=gene_tmp)
    }
    rm(T1);gc()
    rm(T0);gc()

    if(iii==start){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}
write.table_FT_2(res_sum,paste0("ATAC/GLMM_DARs_lowexp0.3/chunk_",chunk_tmp,"_GLMMres.txt"))



###################
## Aggregate results
###################
suppressPackageStartupMessages(library(ggsci))
suppressPackageStartupMessages(library(ggrepel))

list       = make_list("ATAC/GLMM_DARs_lowexp0.3","_GLMMres.txt") 
list$chunk = take_factor(list$FILE,2,"_")%>%as.numeric()
list=list[order(list$chunk),]
     
for(kkk in 1:nrow(list)){
     res_tmp=fread_FT(list$PATH[kkk])%>%filter(str_detect(parameter,pattern="subsetT|anova"))
     if(kkk==1){res_sum=res_tmp}else{res_sum=rbind(res_sum,res_tmp)}
}
# 28818=14409*2

# reorder columns
res_sum2=res_sum%>%pivot_wider(names_from="parameter",values_from=c("Full_Estimate","Full_SE","Full_Z","Full_P"))%>%
         dplyr::select(1:2,4,6,8:9)
colnames(res_sum2)=c("peak","Estimate","se","Z","waldP","lrtP")

# calculate FDR
# keep NA stats
res_sum2=res_sum2%>%
         mutate(FDR =p.adjust(lrtP,method="BH"))
nrow(res_sum2)
# 14409
write.table_FT_2(res_sum2,paste0("ATAC/GLMM_DARs_lowexp0.3/summary/All_GLMMres_sum.txt"))



## volcano plot
## color: FDR<0.1

marker_sum2=res_sum2%>%
            mutate(key=ifelse(FDR<0.1,"Yes","No"))
marker_sum2$key=factor(marker_sum2$key,levels=c("Yes","No"))

marker_sum2$Estimate[marker_sum2$Estimate<(-2)]=-2
marker_sum2$Estimate[marker_sum2$Estimate>2]=2

p =   ggplot(marker_sum2,aes(x=Estimate,y=-log10(lrtP),color=key))+
      geom_point(size=0.25)+
      #geom_text_repel(size=3.5,max.overlaps=Inf)+
      theme_classic()+
      scale_color_manual(values=c("red","black")) +
      theme(axis.text.x=element_text(colour="black",size=12),
            axis.text.y=element_text(colour="black",size=12),
            axis.title.x=element_text(colour="black",size=14),
            axis.title.y=element_text(colour="black",size=14),
            plot.title=element_blank(),
            legend.position="right",
            legend.title=element_blank(),
            legend.text=element_blank())+
      guides(colour=guide_legend(override.aes=list(size = 3)))+
      labs(x="logFC",y="-log10(P)")

pdf_3(paste0("ATAC/GLMM_DARs_lowexp0.3/summary/All_GLMMres_volcano.pdf"),h=2.5,w=4)
     plot(p)
dev.off()


###################
## HOMER
###################

## Prepare bed for HOMER
# peaks bed (FDR<0.1, Up)
res_sum2_sig2up_bed=res_sum2%>%filter(FDR<0.1 & Estimate>0)%>%
                 mutate(chr=take_factor(peak,1,":"))%>%
                 mutate(pos=take_factor(peak,2,":"))%>%
                 mutate(pos1=take_factor(pos,1,"-"))%>%
                 mutate(pos2=take_factor(pos,2,"-"))%>%
                 dplyr::select(chr,pos1,pos2)
write.table_FF_2(res_sum2_sig2up_bed,paste0("ATAC/GLMM_DARs_lowexp0.3/homer/All_GLMMres_FDR0.1up.bed"))
# 44

q()


cd ATAC/GLMM_DARs_lowexp0.3/homer

genome=fasta/genome.fa
gtf=gtf/genes.gtf

peak=All_GLMMres_FDR0.1up.bed

findMotifsGenome.pl ${peak} ${genome} All_GLMMres_FDR0.1up_motif -size given -p 4


###################
## UNIChro-Seq
###################
suppressPackageStartupMessages(library(lme4))

data_full=fread_FT("/home/ha7477/share/to_Nakano/qATAC_20251029/result/qATAC_data.txt")
# 128=(6+2:samples)*(2:T/C)*(6+2:region)

# limit to rep1
data_full=data_full%>%filter(replicate==1)
# 96=(6:samples)*(2:T/C)*(6+2:region)


ref_SNP <- "chr5_139561796_C_G"
target_SNPs <- setdiff(unique(data_full$SNP),c("chr5_139561796_C_G"))
results_df <- data.frame()

for (target_SNP in target_SNPs){  
    data_prep <- data_full %>% 
      filter(SNP %in% c(target_SNP, ref_SNP)) %>% 
      group_by(Donor, subset, replicate) %>%
      summarise(
        ref_snp = sum(count[SNP == ref_SNP]),  
        target_snp = sum(count[SNP == target_SNP]), 
        .groups = 'drop'
      ) %>%
      mutate(
        condition_binary = ifelse(subset == "Target", 1, 0)
      ) %>% 
      as.data.frame()
    

long_data <- data_prep %>%
      pivot_longer(
        cols = c(ref_snp, target_snp),
        names_to = "snp_type",
        values_to = "count"
      ) %>%
      mutate(
        snp_binary = ifelse(snp_type == "target_snp", 1, 0),
        Donor = factor(Donor, levels = unique(Donor))
      ) %>% rename(condition=subset) %>%
      uncount(count) %>% 
      as.data.frame()
    
    glmer_model <- glmer(condition_binary ~ snp_binary + 
                          (1 + snp_binary|Donor),
                        family = binomial,
                        data = long_data)

    model_summary <- summary(glmer_model)
    coef_table <- model_summary$coefficients
    
    results_df <- rbind(results_df, data.frame(
      target_SNP = target_SNP,
      estimate = coef_table["snp_binary", "Estimate"],
      std_error = coef_table["snp_binary", "Std. Error"],
      z_value = coef_table["snp_binary", "z value"],
      p_value = coef_table["snp_binary", "Pr(>|z|)"],
      AIC = AIC(glmer_model)
    ))
}

results_df$adj_pvalue<-p.adjust(results_df$p_value,method="bonferroni")
# chr12_68181269_A_T: significant


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
 [1] ggrepel_0.9.6     ggsci_3.2.0       lme4_1.1-37       Matrix_1.7-0     
 [5] edgeR_4.4.1       limma_3.62.2      data.table_1.16.4 lubridate_1.9.4  
 [9] forcats_1.0.0     stringr_1.5.1     dplyr_1.1.4       purrr_1.0.2      
[13] readr_2.1.5       tidyr_1.3.1       tibble_3.2.1      ggplot2_3.5.2    
[17] tidyverse_2.0.0  

loaded via a namespace (and not attached):
 [1] generics_0.1.3   stringi_1.8.4    lattice_0.22-6   hms_1.1.3       
 [5] magrittr_2.0.3   grid_4.4.0       timechange_0.3.0 scales_1.3.0    
 [9] reformulas_0.4.1 Rdpack_2.6.4     cli_3.6.3        rlang_1.1.5     
[13] rbibutils_2.3    munsell_0.5.1    splines_4.4.0    withr_3.0.2     
[17] tools_4.4.0      tzdb_0.4.0       nloptr_2.2.1     minqa_1.2.8     
[21] colorspace_2.1-1 locfit_1.5-9.10  boot_1.3-30      vctrs_0.6.5     
[25] R6_2.5.1         lifecycle_1.0.4  MASS_7.3-60.2    pkgconfig_2.0.3 
[29] pillar_1.10.1    gtable_0.3.6     glue_1.8.0       Rcpp_1.0.14     
[33] statmod_1.5.0    tidyselect_1.2.1 nlme_3.1-164     compiler_4.4.0  








