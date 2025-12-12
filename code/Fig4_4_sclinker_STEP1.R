source("~/ws/my.source.R")
options(stringsAsFactors=F)

###################
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(optparse))





## 1. gene score

# ScQuaD sumstats
res_sum=fread_FT("singlecell_NBmodel/NBmodel_final/summary/eachcelltype_eachcompare_with_percMT_res_sum_allFDR_final.txt")
# 416134

list=fread_FT("tmp_job_list/eachcelltype_eachcompare_list.txt")

# Positive statistics
for(iii in 1:nrow(list)){
    celltype_tmp=list$celltype[iii]
    compare_tmp =list$compare[iii]

    res_tmp=res_sum%>%filter(celltype==celltype_tmp & compare==compare_tmp)

    # negative Z score to 0 (P value = 1), consistent with original method
    # min-max normalization
    # https://github.com/kkdey/GSSG/blob/master/code/calc_gene_scores/process_sclinker_output.R
    res_tmp2=res_tmp%>%mutate(dys_Z2 =ifelse(is.na(dys_Z),0,
                                      ifelse(dys_Z<0,0,dys_Z)))%>%
                       mutate(nbhd_Z2 =ifelse(is.na(nbhd_Z),0,
                                       ifelse(nbhd_Z<0,0,nbhd_Z)))%>%
                       mutate(dys_P2 =2*pnorm(abs(dys_Z2),lower.tail=F))%>%
                       mutate(nbhd_P2=2*pnorm(abs(nbhd_Z2),lower.tail=F))%>%
                       mutate(dys_log=-2*log(dys_P2+1e-08))%>%
                       mutate(nbhd_log=-2*log(nbhd_P2+1e-08))%>%
                       mutate(dys_PR =(dys_log-min(dys_log))/(max(dys_log)-min(dys_log)))%>%
                       mutate(nbhd_PR=(nbhd_log-min(nbhd_log))/(max(nbhd_log)-min(nbhd_log)))

    res_tmp2_dys=res_tmp2%>%dplyr::select(gene,dys_PR)
    res_tmp2_dys=res_tmp2_dys[order(res_tmp2_dys$dys_PR,decreasing=T),]

    res_tmp2_nbhd=res_tmp2%>%dplyr::select(gene,nbhd_PR)
    res_tmp2_nbhd=res_tmp2_nbhd[order(res_tmp2_nbhd$nbhd_PR,decreasing=T),]

    write.table_FF_2(res_tmp2_dys, paste0("sclinker/genescore/disease/",celltype_tmp,"_",compare_tmp,"_DEGs_pos_genescore.txt"))
    write.table_FF_2(res_tmp2_nbhd,paste0("sclinker/genescore/disease/",celltype_tmp,"_",compare_tmp,"_nbhd_pos_genescore.txt"))
}






## 2. geneset to probabilistic bed
## e.g., NaiveCD4
list     = fread_FT(LIST)
path_tmp = list[task_id,1] # gene score path
type_tmp = "NaiveCD4"
bed_dir = list[task_id,4]  # bed output directory


enhancer_tissue = "BLD"

dir.create_p(paste0(bed_dir,"/",type_tmp))

source("/home/imgnaka/tool/ref/sclinker/code/GeneSet_toS2G/all_bedgraph_methods_mod.R")

gene_scores = fread_FT(path_tmp) 

scores = gene_scores[,2]
names(scores) = gene_scores[,1]

out = ABC_Road_bedgraph_calc(scores,
                             output_cell = paste0(bed_dir,"/",type_tmp),
                             tissuename = enhancer_tissue,
                             output_bed = paste0("ABC_Road_", enhancer_tissue, ".bed"))

df = fread_FT("/home/imgnaka/tool/ref/sclinker/processed_data/Gene_100kb.txt")
df[which(df[,2] < 0), 2] = 0 # start >=0 
score = gene_scores[match(df$gene, gene_scores[,1]), 2]
score[is.na(score)] = 0

temp = cbind.data.frame(df[,1:3], score)
temp2 = temp[which(temp$score != 0),]
write.table(temp2, file = paste0(bed_dir, "/", type_tmp, "/", "100kb.bed"),
            quote=F, sep = "\t", row.names=F, col.names=F)



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
 [1] optparse_1.7.5    R.utils_2.12.3    R.oo_1.27.0       R.methodsS3_1.8.2
 [5] data.table_1.16.4 lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1    
 [9] dplyr_1.1.4       purrr_1.0.2       readr_2.1.5       tidyr_1.3.1      
[13] tibble_3.2.1      ggplot2_3.5.2     tidyverse_2.0.0  

loaded via a namespace (and not attached):
 [1] gtable_0.3.6     compiler_4.4.0   tidyselect_1.2.1 scales_1.3.0    
 [5] R6_2.5.1         generics_0.1.3   munsell_0.5.1    pillar_1.10.1   
 [9] tzdb_0.4.0       rlang_1.1.5      getopt_1.20.4    stringi_1.8.4   
[13] timechange_0.3.0 cli_3.6.3        withr_3.0.2      magrittr_2.0.3  
[17] grid_4.4.0       hms_1.1.3        lifecycle_1.0.4  vctrs_0.6.5     
[21] glue_1.8.0       colorspace_2.1-1 tools_4.4.0      pkgconfig_2.0.3 
















