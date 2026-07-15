source("~/ws/my.source.R")
options(stringsAsFactors=F)

suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggsci))
source("~/ws/colors.R")


#################################
## 1. hclust based on cell-state prop
#################################

## Sample-level metadata
sample_meta5=fread_FT("CNA/clinical/Sample_meta_clinical.txt")

## Cell-state prop over celltype
## Exclude very rare 2 cell states 
prop_sum=fread_FT("CNA/cellstate_prop_percelltype/cellstate_prop_percelltype.txt")%>%
         filter(celltype%in%celltype_list & !cellstate%in%c("cDC-5","gdT-5"))%>%
         dplyr::select(cellstate,sample_uuid,Prop)%>%
         pivot_wider(names_from="sample_uuid",values_from="Prop")%>%
         column_to_rownames("cellstate")
dim(prop_sum)
# 103 269

prop_scaled=apply(prop_sum,1,scale)%>%t()
colnames(prop_scaled)=colnames(prop_sum)
prop_scaled[prop_scaled>3]=3        # clipping
prop_scaled[prop_scaled<(-3)]=-3    # clipping

cellstates_ISG  =c("NaiveCD4-6","CMCD4-5","EMCD4-8","Treg-4","NaiveCD8-5","CMCD8-6","EMCD8-6","gdT-5","CD56dimNK-4","CD56brightNK-4","NaiveB-5","MemB-7","CD16nMono-5","CD16pMono-4","cDC-6","pDC-3")
cellstates_CD69 =c("NaiveCD4-2","CMCD4-2","Treg-1","NaiveCD8-2","CMCD8-2","MAIT-2","CD56dimNK-2","CD56brightNK-3","NaiveB-3","MemB-4")
cellstates_FOXO1=c("NaiveCD4-5","CMCD4-4","EMCD4-4","NaiveCD8-4","EMCD8-3")
cellstates_focus=c("NaiveCD4-4","CMCD4-3","EMCD4-7","Treg-2","CMCD8-4","EMCD8-5","gdT-0","gdT-1","MemB-6","CD16nMono-3") # expanded


# column annotation
colanno_tmp=data.frame(sample_uuid=colnames(prop_scaled))%>%inner_join(.,sample_meta5,by="sample_uuid")%>%column_to_rownames("sample_uuid")%>%
            mutate('Disease status'=ifelse(Status=="Healthy","Healthy",
                           ifelse(Status=="Flare","HDA",
                           ifelse(Status=="Treated","Treated",
                           ifelse(Status=="Managed" & is.na(sledaiscore),"NA",
                           ifelse(Status=="Managed" & sledaiscore==0,"Inactive",
                           ifelse(Status=="Managed" & sledaiscore>=1 & sledaiscore<=4,"LDA",
                           ifelse(Status=="Managed" & sledaiscore>=5 & sledaiscore<=9,"MDA",
                           ifelse(Status=="Managed" & sledaiscore>=10,"HDA",NA)))))))))%>%
            dplyr::select('Disease status')

colanno_tmp$'Disease status'=factor(colanno_tmp$'Disease status',levels=c("Healthy","Inactive","LDA","MDA","HDA","Treated","NA"))

col3=c("Healthy"="#6BAED6","Inactive"="#FFEDA0","LDA"="#FEB24C","MDA"="#FC4E2A","HDA"="#BD0026","Treated"="#9E9AC8","NA"="white")

ha = HeatmapAnnotation(df=colanno_tmp,
                        col=list('Disease status'=col3),
                        height=unit(0.3,"cm"),simple_anno_size_adjust=T,annotation_name_gp=gpar(fontsize=0),
                        annotation_legend_param=list('Disease status'=list(title_gp=gpar(fontsize=0),labels_gp=gpar(fontsize=16)))
                        )

# row annotation
rowanno_tmp=data.frame(cellstate=rownames(prop_scaled))%>%
            mutate(celltype=take_factor(cellstate,1,"-"))%>%
            mutate(type=ifelse(cellstate%in%cellstates_ISG,"ISG15+MX1+",
                        ifelse(cellstate%in%cellstates_CD69,"JUN+CD69",
                        ifelse(cellstate%in%cellstates_FOXO1,"FOXO1+ARHGAP15+",
                        ifelse(cellstate%in%cellstates_focus,"Others (expand)","NA")))))%>%
            column_to_rownames("cellstate")

rowanno_tmp$celltype=factor(rowanno_tmp$celltype,levels=celltype_list)
levels(rowanno_tmp$celltype)=gsub("CD56brightNK","CD56++NK",levels(rowanno_tmp$celltype))%>%
                             gsub("CD16nMono","CD16-Mono",.)%>%
                             gsub("CD16pMono","CD16+Mono",.)
rowanno_tmp$type    =factor(rowanno_tmp$type,levels=c("ISG15+MX1+","JUN+CD69","FOXO1+ARHGAP15+","Others (expand)","NA"))

celltype_col=c("NaiveCD4"="#921b1b","CMCD4"="#cc0010","EMCD4"="#e66557","CTLCD4"="#ff9994","Treg"="#e6a88f",
               "NaiveCD8"="#995522","CMCD8"="#cc9955","EMCD8"="#eebb77","gdT"="#f77308","MAIT"="#feaf53",
               "CD56dimNK"="#e0d51a","CD56++NK"="#f8e58c",
               "NaiveB"="#00561f","MemB"="#8fc31f",
               "CD16-Mono"="#1d2088","CD16+Mono"="#1e90ff",
               "cDC"="#0bdaf2","pDC"="#a4fbfa"
               )
type_col=c("ISG15+MX1+"="#BC3C29FF","JUN+CD69"="#0072B5FF","FOXO1+ARHGAP15+"="#E18727FF","Others (expand)"="#20854EFF","NA"="white")

ra = rowAnnotation(df=rowanno_tmp,col=list(celltype=celltype_col,type=type_col),
                   width=unit(0.8,"cm"),simple_anno_size_adjust=T,show_annotation_name=F,
                   annotation_legend_param=list(celltype=list(title_gp=gpar(fontsize=0),labels_gp=gpar(fontsize=16)),
                                                type=list(title_gp=gpar(fontsize=0),labels_gp=gpar(fontsize=16))))

# color scale
col_fun = colorRamp2(seq(-1.5,1.5,length=9),viridis(9))


## Fig. 5A
p=Heatmap(prop_scaled,clustering_distance_rows="euclidean",clustering_distance_columns="euclidean",
        clustering_method_rows="ward.D2",clustering_method_columns="ward.D2",
        col=col_fun,width=unit(24,"cm"),height=unit(6,"cm"),
        show_column_names=F,show_row_names=F,
        top_annotation = ha,
        right_annotation = ra,
        heatmap_legend_param=list(title="Frequency (scaled)",at=c(-1.5,0,1.5),title_gp=gpar(fontsize=0,fontface="bold"),labels_gp=gpar(fontsize=16),
                                  direction="horizontal",title_position="topcenter",grid_height=unit(0.6,"cm"),legend_width=unit(4,"cm"))
        )

pdf_3(paste0("CNA/clinical/hclust/allsamples_hclust_final.pdf"),h=8,w=12)
 draw(p,heatmap_legend_side="bottom",annotation_legend_side="right")
dev.off()



#################################
## 2. Association with clinical traits
#################################

# save PGs
hc_prop = hclust(dist(t(prop_scaled)), method="ward.D2")

labels_col = hc_prop$labels[hc_prop$order]

PG1_1=labels_col[213:236]
PG1_2=labels_col[237:269]
PG1_3=labels_col[193:212]
PG1_4=labels_col[161:192]
PG2  =labels_col[1:160]

sample_meta6=sample_meta5%>%
             mutate(PG=ifelse(sample_uuid%in%PG1_1,"PG1_1",
                       ifelse(sample_uuid%in%PG1_2,"PG1_2",
                       ifelse(sample_uuid%in%PG1_3,"PG1_3",
                       ifelse(sample_uuid%in%PG1_4,"PG1_4",
                       ifelse(sample_uuid%in%PG2,"PG2",NA))))))%>%
            mutate(DAgroup2=ifelse(Status=="Healthy","Healthy",
                            ifelse(Status=="Flare","HDA",
                            ifelse(Status=="Treated",NA,
                            ifelse(Status=="Managed" & is.na(sledaiscore),NA,
                            ifelse(Status=="Managed" & sledaiscore==0,"Inactive",
                            ifelse(Status=="Managed" & sledaiscore>=1 & sledaiscore<=4,"LDA",
                            ifelse(Status=="Managed" & sledaiscore>=5 & sledaiscore<=9,"MDA",
                            ifelse(Status=="Managed" & sledaiscore>=10,"HDA",NA)))))))))

table(sample_meta6$PG)
# PG1_1 PG1_2 PG1_3 PG1_4   PG2 
#    24    33    20    32   160


## save CGs
hc_prop2 = hclust(dist(prop_scaled), method="ward.D2")

labels_row = hc_prop2$labels[hc_prop2$order]

CG1=data.frame(CG="CG1",cellstate=labels_row[1:17])
CG2=data.frame(CG="CG2",cellstate=labels_row[33:48])
CG3=data.frame(CG="CG3",cellstate=labels_row[18:32])
CG4=data.frame(CG="CG4",cellstate=labels_row[62:81])
CG5=data.frame(CG="CG5",cellstate=labels_row[82:103])
CG6=data.frame(CG="CG6",cellstate=labels_row[49:61])



## DA group frequency
## Fisher test
table1=table(sample_meta6$PG,sample_meta6$DAgroup2)
table1=table1[,c("Inactive","LDA","MDA","HDA")]

set.seed(0)
fisher.test(table1,simulate.p.value=TRUE, B=10000)
# p-value = 4e-04

prop1=apply(table1,1,function(x){x/sum(x)*100})%>%
      as.data.frame()%>%
      rownames_to_column("DAgroup")%>%
      pivot_longer(col=-DAgroup,names_to="PG",values_to="Prop")

prop1$DAgroup=factor(prop1$DAgroup,levels=c("Inactive","LDA","MDA","HDA"))
prop1$PG     =factor(prop1$PG,levels=c("PG1_1","PG1_2","PG1_3","PG1_4","PG2"))
levels(prop1$PG)=c("1-1","1-2","1-3","1-4","2")

p =  ggplot(data=prop1,aes(x=PG,y=Prop,fill=DAgroup))+
     geom_bar(stat="identity",position="stack")+
     #geom_text(aes(label=Num),size=4,vjust=0)+
     theme_classic()+
     scale_fill_manual(values=c("#FFEDA0","#FEB24C","#FC4E2A","#BD0026"))+
     theme(axis.text.x=element_text(colour="black",size=14,angle=45,hjust=1,vjust=1),
           axis.text.y=element_text(colour="black",size=12),
           axis.title.x=element_text(colour="black",size=14),
           axis.title.y=element_text(colour="black",size=14),
           plot.title=element_blank(),
           legend.position="right",
           legend.title=element_blank(),
           legend.text=element_text(colour="black",size=14))+
     scale_y_continuous(breaks=c(0,50,100))+
     labs(x="SG",y="Proportion (%)")
# "PG" to "SG" since healthy is not patients

pdf_3(paste0("CNA/clinical/hclust/allsamples_hclust_PG_DAgroup_bar.pdf"),h=2.5,w=3.5)
     plot(p)
dev.off()













