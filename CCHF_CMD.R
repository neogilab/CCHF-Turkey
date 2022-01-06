
library(NOISeq)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Network/Prot_Coding_CCHF.txt",header = T,row.names = 1)
head(data)
dataX=as.matrix(data)
norm=tmm(dataX)
head(norm)

write.table(norm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Network/TMM_Normed.txt",sep="\t",col.names = NA,quote = FALSE)

library(edgeR)
library(Biobase)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Coding_RawCount_AT_BT.txt",row.names = 1, check.names = FALSE)
countX=as.matrix(count)
head(countX)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/MetaData.txt",row.names = 1)
expSet=ExpressionSet(countX, phenoData=AnnotatedDataFrame(data=meta))
expSet <- expSet[rowSums(exprs(expSet)) != 0, ] # to remove genes having zeros in all samples
log2cpm <- cpm(exprs(expSet), log = TRUE)

write.table(log2cpm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Coding_lCPM_AT_BT.txt",sep="\t",quote = FALSE,col.names = NA)

library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Coding_lCPM_AT_BT.txt",header=T,row.names = 1)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/lcpm_Umap.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
library(ggrepel)
head(dat)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/lcpm_Umap.txt",check.names = FALSE)
head(dat)
dat$Label <- factor(dat$Label, levels=c("a_BT_S1", "b_AT_S1", "a_BT_S2","b_AT_S2"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/UMAP_BT_AT.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Label)) + geom_point(size=6.5,shape=21,aes(fill=Label))+theme_grey()+
  stat_ellipse(aes(x=V1, y=V2,group=Group),color="#7da7ca",linetype=2,level=0.95)+
  scale_color_manual(labels=c("Severity group 1(Acute)","Severity group 1(Recovered)","Severity group 2(Acute)","Severity group 2(Recovered)"),
                     values=c(a_BT_S1="#996300",b_AT_S1="#cc8400",a_BT_S2="#990000",b_AT_S2="#cc0000"))+
  scale_fill_manual(labels=c("Severity group 1(Acute)", "Severity group 1(Recovered)","Severity group 2(Acute)","Severity group 2(Recovered)"),
                    values=c(a_BT_S1="#b27300",b_AT_S1="#ffa500",a_BT_S2="#b20000",b_AT_S2="#ff0000"))+
  theme(axis.title = element_text(size=15,hjust = 0.5),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),axis.text = element_text(colour = "black"),
        legend.position = c(0.8, 0.109),legend.text = element_text(size = 11),legend.title = element_blank())+labs(x = "UMAP 1",y="UMAP 2",fill="Group",color="Group")
dev.off()

?stat_ellipse

######################## DESEQ2 After -vs- Before
library(DESeq2)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Coding_RawCount_AT_BT.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MetaData.txt")
head(Design)
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~pair+Group)
ds=DESeq(ds)
res=results(ds,independentFiltering = FALSE)
?results
res=results(ds,contrast = c("Group","BT","AT"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Result.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Sev1_AT_BT/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Sev1_AT_BT/design.txt")
head(Design)
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~pair+Group)
ds=DESeq(ds)
res=results(ds,independentFiltering = FALSE)
?results
res=results(ds,contrast = c("Group","BT","AT"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Sev1_AT_BT/Result.txt",sep="\t", quote=FALSE,col.names = NA)


count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Sev2_AT_BT/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Sev2_AT_BT/design.txt")
head(Design)
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~pair+Group)
ds=DESeq(ds)
res=results(ds,independentFiltering = FALSE)
?results
res=results(ds,contrast = c("Group","BT","AT"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Sev2_AT_BT/Result.txt",sep="\t", quote=FALSE,col.names = NA)



################################## HEATMAP
library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Union_LCPM.txt",header = T,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Zscore.txt",header = T, row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/MetaData.txt",row.names = 1,header = T)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/GeneLFC2.txt",header = T,check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 20),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                          title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(SeverityScore = c("S1"="#2e6c92","S2"="#f29393"),Group=c("At_Infection"="#6d5210","12M_After"="#c4941c")))
col_fun2 = colorRamp2(c(3, 1, 0, -1, -3), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun2,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),
           top_annotation  =ha,row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",10),rep("S2",14)))

col_fun_lfc = colorRamp2(c(-5, -3,-2,-1, 0,1,2,3,5), c("#0000cc","#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200","#7f7f00"))

H3=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(5, "cm"),
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(40, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

tt=H1 + H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/HeatMap.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()

Cls=row_order(H1)
row.names(Zscore)

??annotation_legend_param
clu = t(t(row.names(Zscore[row_order(H1)[[2]],])))
write.table(as.data.frame(clu),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/CL2.txt",sep="\t",quote = FALSE,col.names = NA)



################################################

library(NOISeq)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Network/Global/Coding_RawCount_AT_BT.txt",header = T,row.names = 1)
head(data)
dataX=as.matrix(data)
norm=tmm(dataX)
head(norm)

write.table(norm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Network/Global/TMM.txt",sep="\t",col.names = NA,quote = FALSE)


####################################


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Cluste1/Input.txt",header = T)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Cluste1/Pathway.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(high="#fce79a",low="#b3962d",limits=c(3E-11,0.05),breaks=c(3E-10,3E-05,0.000001, 0.0001,0.01))+
  scale_y_discrete(position = "right")+theme_bw()+scale_size(range = c(5,7))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),legend.position = c(9, -0.12),
       legend.box="vertical",legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),
        panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(1,4.8,3.5,5, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow = 1,override.aes = list(size = 5)))
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Cluste2/Input.txt",header = T)
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/Cluste2/Pathway.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(high="#ebf9f6",low="#0d5446",limits=c(3E-15,0.05),breaks=c(3E-10,3E-05,0.00001, 0.0001,0.01))+
  scale_y_discrete(position = "right")+theme_bw()+scale_size(range = c(5,7))+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),legend.position = c(7, -0.12),
        legend.box="vertical",legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),
        panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(0.5,3.7,3.5,5, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow = 1,override.aes = list(size = 5)))
dev.off()

######################## I S G
library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/ISG_LCPM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)


Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/MetaData.txt",row.names = 1)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/ISG_LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(SeverityScore = c("S1"="#ffa500","S2"="#cc0000"),Group=c("At_Infection"="#590059","12M_After"="#cc99cc")))
#col_fun2 = colorRamp2(c(2, 1, 0, -1, -2), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#5d5dff","#4c4cff","#0000cc"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),
           top_annotation  =ha,row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(37, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",10),rep("S2",14)))

col_fun_lfc = colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#006600","#198c19","#4ca64c","white","#ca8383","#c26e6e" ,"#b04545"))
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(1, "cm"),
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 11),height  = unit(37, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

GO=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Untitled Folder/GO_ISG_LFC2.txt",row.names = 1,check.names = FALSE)
col_go = colorRamp2(c(0,1), c("#d8d8d8","red"))
head(LFC)
ha2 = rowAnnotation(foo = anno_mark(at = c(2,4,5,7,8,9,10,11,12,13,17,18,21,22,25,27,28,30,33,36,37,38,39,40,41,42,43,44,45,46,50,51,52,53,55,56,57,58,59,60,61,62,63,65,67,68,69,70,72,73,74,
                                             75,78,80,81,82,85,86,87,88,94,101,103,111,112,116,117,118,120,122,123,124,125,126,128,133,134,135),
                                      labels_gp = gpar(fontsize=15),lines_gp = gpar(col="black"),link_height = unit(7, "mm"),link_width=unit(12, "mm"),labels = c("ADAR","B2M","BST2","CAMK2G","CCL1","CCL2","CCL25","CCL4","CCL5","CCL8","CD44","CIITA","DAPK1","DAPK3","EPRS","FCGR1A","FCGR1B","GAPDH","GBP1","HCK","HFE","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DQB2","HLA-DRA","HLA-DRB1","HLA-E","HLA-F","ICAM1","IFI27","IFI35","IFI6","IFIT1","IFIT2","IFIT3","IFITM1","IFITM2","IFITM3","IFNAR1","IFNAR2","IFNGR2",
                                                                                                                               "IL12RB1","IRF2","IRF4","IRF7","IRF9","ISG15","ISG20","JAK1","JAK2","MT2A","MX1","MX2","NCAM1","OAS1","OAS2","OAS3","OASL","PML","PRKCD","PSMB8","RPL13A","RSAD2","SP100","STAT1","STAT2","SYNCRIP","TRIM21","TRIM22","TRIM25","TRIM26","TRIM5","TRIM8","XAF1","XCL1","XCL2")))

h3= columnAnnotation(foo = anno_mark(at = c(1,2,3,4,5,6,7,8,9,10,11,12,13),
                                      labels_gp = gpar(fontsize=18),lines_gp = gpar(col=NA),link_width=unit(0.8, "mm"),
                                     labels = c("5","5","6","10","12","13","14","14","16","23","29","35","57")))

H3=Heatmap(as.matrix((GO)),col=col_go,cluster_rows=FALSE,cluster_columns = FALSE,name="GO",width  = unit(10, "cm"),show_heatmap_legend = FALSE,
           right_annotation = ha2,show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           bottom_annotation = h3,row_names_gp =gpar(fontsize = 11),height  = unit(37, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

tt=H1 + H2 + H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/HeatMap.pdf",height = 22,width =20)
draw(tt, merge_legend = TRUE)
dev.off()



######


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MAplot_Result.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MAplot.pdf")
ggplot(data, aes(x=log2(baseMean), y=log2FoldChange)) + 
  geom_point(data=subset(data, padj>=0.05),aes(log2(baseMean),y=log2FoldChange),color="#ffffb2",size=1.2)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0),aes(x=log2(baseMean),y=log2FoldChange),color="#007300",size=1.5)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0),aes(x=log2(baseMean),y=log2FoldChange),color="#b20000",size=1.5)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 9),aes(x=log2(baseMean),y=log2FoldChange),color="#b20000",size=3)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,
                   label.size = NA,segment.alpha=0.75,box.padding=0.5,nudge_x = 0.55,nudge_y = 2)+
  scale_fill_manual(values=c(HC_EC="#71abab",Middle="#b69c67",NS="#99aab5"))+guides(fill = guide_legend(title = "Regulation",override.aes = aes(label = "")))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=25),legend.position = "none",
        axis.title.x=element_text(size=25),axis.text=element_text(size=20,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 mean expression",y="Log2 fold change")
dev.off()

head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MAplot.pdf")
ggplot(data, aes(x=log2(baseMean), y=log2FoldChange)) + 
  geom_point(data=subset(data, padj>=0.05),aes(log2(baseMean),y=log2FoldChange),color="#bfbfbf",size=1.2)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0),aes(x=log2(baseMean),y=log2FoldChange,color=abs(log2FoldChange)),size=2.2)+
  scale_color_gradient(low = "#d1ff19", high = "#198c19")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0),shape=21,aes(x=log2(baseMean),y=log2FoldChange,fill=abs(log2FoldChange)),size=2.2,color="transparent")+
  scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 9),aes(x=log2(baseMean),y=log2FoldChange),color="#b20000",size=3.5)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,
                   label.size = NA,segment.alpha=0.75,box.padding=0.5,nudge_x = 0.55,nudge_y = 2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=25),legend.position = "none",
        axis.title.x=element_text(size=25),axis.text=element_text(size=20,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 mean expression",y="Log2 fold change")
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MAplot_2.pdf")
ggplot(data, aes(x=log2(baseMean), y=log2FoldChange)) + 
  geom_point(data=subset(data, padj>=0.05),aes(log2(baseMean),y=log2FoldChange),color="#bfbfbf",size=1.2)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0),aes(x=log2(baseMean),y=log2FoldChange,color=abs(log2FoldChange)),size=2.2)+
  scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0),shape=21,aes(x=log2(baseMean),y=log2FoldChange,fill=abs(log2FoldChange)),size=2.2,color="transparent")+
  scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 9),aes(x=log2(baseMean),y=log2FoldChange),color="#b20000",size=3.5)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,
                   label.size = NA,segment.alpha=0.75,box.padding=0.5,nudge_x = 0.55,nudge_y = 2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=25),legend.position = "none",
        axis.title.x=element_text(size=25),axis.text=element_text(size=20,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 mean expression",y="Log2 fold change")
dev.off()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MAplot_Seveity_BT.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/MAplot_Seveity_BT.pdf")
ggplot(data, aes(x=log2(baseMean), y=log2FoldChange)) + 
  geom_point(data=subset(data, padj>=0.05),aes(log2(baseMean),y=log2FoldChange),color="#999999",size=1.2)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange < 0),aes(x=log2(baseMean),y=log2FoldChange),color="#007300",size=1.5)+
  geom_point(data=subset(data, padj<0.05 & log2FoldChange > 0),aes(x=log2(baseMean),y=log2FoldChange),color="#b20000",size=1.5)+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4,label.size = NA,segment.alpha=0.75,box.padding=0,nudge_x = 0.55,nudge_y = 0.25)+
  scale_fill_manual(values=c(HC_EC="#71abab",Middle="#b69c67",NS="#99aab5"))+guides(fill = guide_legend(title = "Regulation",override.aes = aes(label = "")))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=25),legend.position = "none",
        axis.title.x=element_text(size=25),axis.text=element_text(size=20,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 mean expression",y="Log2 fold change")
dev.off()

########################  PROTEOMICS

library(NormalyzerDE)
?normalyzer
normalyzer(jobName="CCHF",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/metadata.txt",
             dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Raw.txt",
             outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/NORM"
             )


ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Quantile.txt",row.names = 1,check.names = FALSE)
head(ip)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/metadata.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
(groupH48_Treated - groupH48_Untreated) - (groupH24_Treated - groupH24_Untreated)
cont.matrix <-makeContrasts(H48_vs_H24 = ((groupH48_Treated - groupH48_Untreated) - (groupH24_Treated - groupH24_Untreated)),
                            levels = design)

cont.matrix <-makeContrasts(H48_vs_H24 = ((groupH48_Treated - groupH24_Treated) - (groupH48_Untreated - groupH24_Untreated)),
                            levels = design)

fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)

write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_TIME.txt",sep="\t",col.names = NA,quote = FALSE)




ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Quantile_24h.txt",row.names = 1,check.names = FALSE)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/des_24h.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(H24_Untreated_vs_H24_Treated = (groupH24_Treated - groupH24_Untreated),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)

write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Results_24h.txt",sep="\t",col.names = NA,quote = FALSE)


ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Quantile_48h.txt",row.names = 1,check.names = FALSE)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/des_48h.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(H48_Untreated_vs_H48_Treated = (groupH48_Treated - groupH48_Untreated),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)

write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Results_48h.txt",sep="\t",col.names = NA,quote = FALSE)

#################


library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Quantile.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/PCA.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Untreated_24h","Treated_24h","Untreated_48h","Treated_48h"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/PCA_2.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(Untreated_24h="#38725e",Treated_24h="#cc8400",Untreated_48h="#4c8ec1",Treated_48h="#cc5151"))+
  scale_fill_manual(values=c(Untreated_24h="#478f76",Treated_24h="#ffa500",Untreated_48h="#5fb2f2",Treated_48h="#ff6666"))+
  labs(x="PC1, 50.8%",y="PC2, 28.27%")+theme(axis.title = element_text(size=13),axis.text=element_text(size = 12),legend.position = c(0.84,0.1),plot.margin = margin(1,1,1,1, "cm"),
                                              legend.title=element_blank(),legend.text=element_text(size=13),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 5))) 
dev.off()


################ HeatMap DEGs


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/DEG_Quantile.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/HeatmapDEG.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/metadata.txt",row.names = 1,check.names = FALSE)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(`Time point` = c("24h"="#ffa500","48h"="#cc0000"),Condition=c("Untreated"="#cc99cc","Treated"="#590059")))


col_fun2 = colorRamp2(c(3, 1, 0, -1, -3), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
ncol(Zscore)
H1=Heatmap(as.matrix((Zscore)),col=col_fun2,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(4, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),
           top_annotation  =ha,row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",6),rep("S2",6)))

col_fun_lfc = colorRamp2(c(-5, -3,-2,-1, 0,1,2,3,5), c("#0000cc","#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200","#7f7f00"))

H3=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(2, "cm"),
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(40, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

tt=H1+H3
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/HeatMap.pdf",height = 20,width =15)
draw(tt, merge_legend = TRUE)
dev.off()

Cls=row_order(H1)
row.names(Zscore)

??annotation_legend_param
clu = t(t(row.names(Zscore[row_order(H1)[[2]],])))
write.table(as.data.frame(clu),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/CL2.txt",sep="\t",quote = FALSE,col.names = NA)






################## ELISA Violin Plot


Bean=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Elisa/Test/BoxInput.txt",check.names = FALSE)
head(Bean)
P1=ggplot(Bean,aes(x=Group,y=`GM-CSF`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "GM-CSF")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 6, by = 1))+
  scale_color_manual(labels = c("Severity1","Severity2"),values=c(Severity1="#cc8400",Severity2="#cc0000"))+
  scale_fill_manual(labels = c("Severity1","Severity2"),values=c(Severity1="#ffa500",Severity2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=15),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P2=ggplot(Bean,aes(x=Group,y=`TNF-alpha`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
    geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
    labs(y="log2(pg/ml)",title = "TNF-alpha")+
    theme_bw()+scale_y_continuous(breaks = seq(5, 9, by = 1))+
    scale_color_manual(labels = c("Severity1","Severity2"),values=c(Severity1="#cc8400",Severity2="#cc0000"))+
    scale_fill_manual(labels = c("Severity1","Severity2"),values=c(Severity1="#ffa500",Severity2="#ff0000"))+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
          axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

  
P3=ggplot(Bean,aes(x=Group,y=`IL-8`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
    geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
    labs(y="log2(pg/ml)",title = "IL-8")+
    theme_bw()+scale_y_continuous(breaks = seq(3, 8, by = 1.5))+
    scale_color_manual(labels = c("Severity1","Severity2"),values=c(Severity1="#cc8400",Severity2="#cc0000"))+
    scale_fill_manual(labels = c("Severity1","Severity2"),values=c(Severity1="#ffa500",Severity2="#ff0000"))+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
          axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))
  
  
P4=ggplot(Bean,aes(x=Group,y=`IL-10`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
    geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
    labs(y="log2(pg/ml)",title = "IL-10")+
    theme_bw()+scale_y_continuous(breaks = seq(5, 12, by = 2))+
    scale_color_manual(labels = c("Severity score 1","Severity score 2"),values=c(Severity1="#cc8400",Severity2="#cc0000"))+
    scale_fill_manual(labels = c("Severity score 1","Severity score 2"),values=c(Severity1="#ffa500",Severity2="#ff0000"))+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(1.7,0.5),panel.border = element_blank(),
          axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
          axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))
library(ggpubr)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Elisa/Test/Box.pdf",height = 7,width =12)
ggarrange(P1,P2,P3,P4,nrow = 1)+theme(plot.margin = margin(4,7,4,0.5, "cm"))
dev.off()


######################## UMAP Before therapay
library(edgeR)
library(Biobase)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_Prot_Coding.txt",row.names = 1, check.names = FALSE)
countX=as.matrix(count)
head(countX)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_metadata.txt",row.names = 1)
expSet=ExpressionSet(countX, phenoData=AnnotatedDataFrame(data=meta))
expSet <- expSet[rowSums(exprs(expSet)) != 0, ] # to remove genes having zeros in all samples
log2cpm <- cpm(exprs(expSet), log = TRUE)

write.table(log2cpm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_LCPM.txt",sep="\t",quote = FALSE,col.names = NA)



library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_LCPM.txt",row.names = 1)
Art.umap = umap(t(data))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_UMAP.txt",sep="\t",col.names = NA,quote = FALSE)
library(ggplot2)
library(ggrepel)
head(dat)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_UMAP.txt",check.names = FALSE)
head(dat)
dat$Severity <- factor(dat$Severity, levels=c("Severity_1", "Severity_2", "Severity_3"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Severity_UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Severity)) + geom_point(size=6.5,shape=21,aes(fill=Severity))+theme_grey()+
  scale_color_manual(labels=c("Severity group 1","Severity group 2","Severity group 3"),
                     values=c(Severity_1="#cc8400",Severity_2="#cc0000",Severity_3="#7f0000"))+
  scale_fill_manual(labels=c("Severity group 1","Severity group 2","Severity group 3"),
                     values=c(Severity_1="#ffa500",Severity_2="#ff0000",Severity_3="#990000"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),axis.text = element_text(colour = "black",size=15),
        legend.position = c(0.19, 0.09),legend.text = element_text(size = 18),legend.title = element_blank())+labs(x = "UMAP 1",y="UMAP 2",fill="Group",color="Group")
dev.off()


library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_LCPM.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_PCA.txt",check.names = FALSE)
head(dat)
dat$Severity <- factor(dat$Severity, levels=c("Severity_1", "Severity_2", "Severity_3"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/BT_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Severity)) + geom_point(size=6.5,shape=21,aes(fill=Severity))+theme_grey()+
  scale_color_manual(labels=c("Severity score 1","Severity score 2","Severity score 3"),
                     values=c(Severity_1="#cc8400",Severity_2="#cc0000",Severity_3="#7f0000"))+
  scale_fill_manual(labels=c("Severity score 1","Severity score 2","Severity score 3"),
                    values=c(Severity_1="#ffa500",Severity_2="#ff0000",Severity_3="#990000"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),axis.text = element_text(colour = "black",size=15),
        legend.position = c(0.19, 0.09),legend.text = element_text(size = 18),legend.title = element_blank())+labs(x = "PC1, 22.03% variation",y="PC2, 15.24% variation",fill="Group",color="Group")
dev.off()


############################# GO Proteomics Heatmap

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/GO_Quantile.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/GO.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/GO_Zscore",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/GO_Zscore", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/metadata.txt",row.names = 1,check.names = FALSE)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/GO_3.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(`Time point` = c("24h"="#ffa500","48h"="#cc0000"),Condition=c("Untreated"="#cc99cc","Treated"="#590059")))


col_fun2 = 
col_fun2=colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200"))
ncol(Zscore)
H1=Heatmap(as.matrix((Zscore)),col=col_fun2,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),
           top_annotation  =ha,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",6),rep("S2",6)))

col_fun_lfc = colorRamp2(c(-2, -1, 0, 1, 2), c("#004c00","#008000" ,"white","#ab6d2e","#66411b"))

H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(2, "cm"),
           show_row_names = TRUE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 15),height  = unit(40, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#8c8c8c") 

GO=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/X_X_GO.txt",check.names = FALSE,row.names = 1,sep = "\t")

col_go = colorRamp2(c(0,1), c("#d8d8d8","#740001"))
H3=Heatmap(as.matrix((GO)),col=col_go,cluster_rows=FALSE,cluster_columns = FALSE,name="GO",width  = unit(2, "cm"),show_heatmap_legend = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 7),height  = unit(40, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 

rna=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/RNASeq.txt",check.names = FALSE,row.names = 1)
H4=Heatmap(as.matrix((rna)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="RNA_LFC",width  = unit(1, "cm"),row_names_side = "right",show_heatmap_legend = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 15),height  = unit(40, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#8c8c8c") 

tt=H1+H3+H2+H4
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/GO.pdf",height = 20,width =25)
draw(tt, merge_legend = TRUE,auto_adjust = FALSE)
dev.off()

#######################################  pathway bubble & bar graph


tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/METABOLIC_Size.txt")
head(tmp)
tmp$value <- factor(tmp$value, levels = tmp$value)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/Metabolic_24H.pdf")
 ggplot(data=tmp, aes(x=factor(Pathway,levels = unique(Pathway)), y=Perc_24h, fill=Reg_24h))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 75, by = 30))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(1,8,4,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y =element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation at 24h"))

dev.off() 
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/Metabolic_48H.pdf")
ggplot(data=tmp, aes(x=factor(Pathway,levels = unique(Pathway)), y=Perc_48h, fill=Reg_48h))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 75, by = 30))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(1,8,4,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y =element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation at 48h"))

dev.off()
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/ip_Metabolic.txt")
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/MetabolicBubble.pdf")
  ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#ffffff",low="#bd782f",breaks=c(2e-16,2e-10,2e-8,0.000001,0.00001,0.0001),limits=c(2e-16,0.001))+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        axis.ticks = element_blank(),legend.position = c(-4.5, -0.2),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=15,colour = "black"),
        legend.title = element_text(size=15),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(1,7.1,5,0.5, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Proteins",nrow = 1),color=guide_legend(title = "Adj.Pvalue",nrow = 2,override.aes = list(size = 5)))
dev.off()



tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/OTHER_Size.txt")
head(tmp)
tmp$value <- factor(tmp$value, levels = tmp$value)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/NonMetabolic_24H.pdf")
ggplot(data=tmp, aes(x=factor(Pathway,levels = unique(Pathway)), y=Perc_24h, fill=Reg_24h))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 50, by = 20))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(1,8,4,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y =element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation at 24h"))

dev.off() 
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/NonMetabolic_48H.pdf")
ggplot(data=tmp, aes(x=factor(Pathway,levels = unique(Pathway)), y=Perc_48h, fill=Reg_48h))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 75, by = 30))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(1,8,4,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y =element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation at 48h"))

dev.off()
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/ip_Other.txt")
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Pathway/KEGG/PLOT/NonMetabolicBubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=Gene,color=pval))+scale_x_discrete(limits=c("1"))+scale_color_gradient(high="#b9d5f0",low="#135ca4",breaks=c(1e-16,1e-10,0.000001,0.0001),limits=c(1.6E-22,0.0007))+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        axis.ticks = element_blank(),legend.position = c(-3.7, -0.2),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=15,colour = "black"),
        legend.title = element_text(size=15),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(1,7.9,5,0.5, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Genes",nrow = 2),color=guide_legend(title = "Adj.Pvalue",nrow = 2,override.aes = list(size = 5)))
dev.off()


######################### ELISA all protein violin

Bean=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Elisa/Test/log2.tab",check.names = FALSE)
head(Bean)
P1=ggplot(Bean,aes(x=Group,y=`Eotaxin`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "Eotaxin")+
  theme_bw()+scale_y_continuous(breaks = seq(7, 11, by = 1))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=15),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P2=ggplot(Bean,aes(x=Group,y=`G-CSF`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "G-CSF")+
  theme_bw()+scale_y_continuous(breaks = seq(5, 9, by = 1))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P3=ggplot(Bean,aes(x=Group,y=`IFN-a2`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IFN-a2")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 9, by = 2))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P4=ggplot(Bean,aes(x=Group,y=`IFN-y`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IFN-y")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 8, by = 1.5))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))

P5=ggplot(Bean,aes(x=Group,y=`IL-12`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-12")+
  theme_bw()+scale_y_continuous(breaks = seq(1, 5, by = 1))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P6=ggplot(Bean,aes(x=Group,y=`IL-15`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-15")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 6, by = 1))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=15),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P7=ggplot(Bean,aes(x=Group,y=`IL-17a`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-17a")+
  theme_bw()+scale_y_continuous(breaks = seq(1, 7, by = 1.5))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P8=ggplot(Bean,aes(x=Group,y=`IL-1a`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-1a")+
  theme_bw()+scale_y_continuous(breaks = seq(5, 10, by = 1.5))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P9=ggplot(Bean,aes(x=Group,y=`IL-9`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-9")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 3, by = 1.5))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P10=ggplot(Bean,aes(x=Group,y=`IL-1b`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-1b")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 3, by = 0.8))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P11=ggplot(Bean,aes(x=Group,y=`IL-2`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-2")+
  theme_bw()+scale_y_continuous(breaks = seq(-1, 3, by = 0.8))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=15),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P12=ggplot(Bean,aes(x=Group,y=`IL-4`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-4")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 9, by = 2))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P13=ggplot(Bean,aes(x=Group,y=`IL-5`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-5")+
  theme_bw()+scale_y_continuous(breaks = seq(-3, 5, by = 2))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P14=ggplot(Bean,aes(x=Group,y=`IL-6`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IL-6")+
  theme_bw()+scale_y_continuous(breaks = seq(1, 8, by = 2))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P15=ggplot(Bean,aes(x=Group,y=`IP-10`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "IP-10")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 15, by = 3.5))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P16=ggplot(Bean,aes(x=Group,y=`MCP-1`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "MCP-1")+
  theme_bw()+scale_y_continuous(breaks = seq(8, 15, by = 2))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_text(size=15),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P17=ggplot(Bean,aes(x=Group,y=`MIP-1a`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "MIP-1a")+
  theme_bw()+scale_y_continuous(breaks = seq(2, 5, by = 1))+
  scale_color_manual(values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = "none",panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))


P18=ggplot(Bean,aes(x=Group,y=`MIP-1b`,fill=Group,color=Group))+geom_violin(alpha=0.5)+
  geom_jitter(shape=16, size=1,color="black",position=position_jitter(0.05))+
  labs(y="log2(pg/ml)",title = "MIP-1b")+
  theme_bw()+scale_y_continuous(breaks = seq(3, 6, by = 1))+
  scale_color_manual(labels = c("Severity score 1","Severity score 2"),values=c(S1="#cc8400",S2="#cc0000"))+
  scale_fill_manual(labels = c("Severity score 1","Severity score 2"),values=c(S1="#ffa500",S2="#ff0000"))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position = c(1.6,0.6),panel.border = element_blank(),
        axis.ticks.x = element_blank(),axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size = 15),
        axis.text.y = element_text(color="black",size=13),plot.title = element_text(hjust = 0.5,size = 15))
library(ggpubr)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Elisa/Test/All.pdf",height = 10,width =12)
ggarrange(P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,nrow = 4,ncol = 5)+theme(plot.margin = margin(2,0.5,2,0.5, "cm"))
dev.off()

############################### Deseq2 Recovered S1-vs-S2

library(DESeq2)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/Input.txt",row.names = 1, check.names = FALSE)
Design=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/Metadata.txt")
head(Design)
ds=DESeqDataSetFromMatrix(countData = count,colData = Design,design = ~Severity)
ds=DESeq(ds)
res=results(ds,contrast = c("Severity","S2","S1"),independentFiltering = FALSE)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/Result.txt",sep="\t", quote=FALSE,col.names = NA)


library(edgeR)
library(Biobase)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/Input.txt",row.names = 1, check.names = FALSE)
countX=as.matrix(count)
head(countX)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/Metadata.txt",row.names = 1)
expSet=ExpressionSet(countX, phenoData=AnnotatedDataFrame(data=meta))
expSet <- expSet[rowSums(exprs(expSet)) != 0, ] # to remove genes having zeros in all samples
log2cpm <- cpm(exprs(expSet), log = TRUE)

write.table(log2cpm,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/LCPM.txt",sep="\t",quote = FALSE,col.names = NA)

library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/LCPM.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/Metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/PCA.txt",check.names = FALSE)
head(dat)
dat$Severity <- factor(dat$Severity, levels=c("Severity_1", "Severity_2"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Recovered_S1_vs_S2/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Severity)) + geom_point(size=6.5,shape=21,aes(fill=Severity))+theme_grey()+
  scale_color_manual(labels=c("Severity score 1","Severity score 2"),
                     values=c(Severity_1="#cc8400",Severity_2="#cc0000"))+
  scale_fill_manual(labels=c("Severity score 1","Severity score 2"),
                    values=c(Severity_1="#ffa500",Severity_2="#ff0000"))+
  theme(axis.title = element_text(size=20,hjust = 0.5),plot.margin = margin(0.7,0.5,0.7,0.5, "cm"),axis.text = element_text(colour = "black",size=15),
        legend.position = c(0.19, 0.09),legend.text = element_text(size = 18),legend.title = element_blank())+labs(x = "PC1, 41.55% variation",y="PC2, 17.36% variation",fill="Group",color="Group")
dev.off()

####### GO bubble plot

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/GO_Input.txt")
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/GO_Bubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=size,color=Ratio))+scale_x_discrete(limits=c("1"))+scale_color_gradient(low="#b3cde0",high="#005b96")+
  scale_y_discrete(position = "left")+theme_bw()+coord_fixed(ratio = 1.7)+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=13,color="black"),
        axis.ticks = element_blank(),legend.position = c(-3.7, -0.3),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=15,colour = "black"),plot.margin = margin(2.5,3,5,0.5, "cm"),
        legend.title = element_text(size=15),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "-log10(FDR)",nrow = 2),color=guide_legend(title = "Gene-ratio(%)",nrow = 2,override.aes = list(size = 5)))
dev.off()


tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/GO_Bar.txt")
head(tmp)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/GO_Bar.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=Reg))+ylab("% Genes")+scale_y_continuous(breaks = seq(0, 35, by = 10))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(1,8,4,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y = element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_blank(),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation at 24h"))

dev.off() 


################ Time series pathway

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/PATHWY.txt")
head(tmp)
tmp$value <- factor(tmp$value, levels = tmp$value)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/PATHWY_bargraph.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=nGene))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 50, by = 20))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(0.2,8,0.5,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y =element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/Input2.txt")
head(data)
library(ggplot2)
data$Term <- factor(data$Term, levels = data$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/Pathway_bubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=Term,size=nGene,color=pval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(high="#b9d5f0",low="#135ca4",breaks=c(0.0000000001,0.0001,0.001,0.05,0.03),limits=c(0.0000000001,0.04))+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),legend.position = c(8, 0.5),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=15,colour = "black"),plot.margin = margin(0.2,7.2,0,1.5, "cm"),
        legend.title = element_text(size=15),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "# Proteins",nrow = 2),color=guide_legend(title = "Adj.Pvalue",nrow = 2,override.aes = list(size = 5)))
dev.off()

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/24H/PATHWAY.txt")
head(tmp)
tmp$value <- factor(tmp$value, levels = tmp$value)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/24H/PATHWAY1.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=nGene))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 50, by = 20))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(0.2,8,0.5,7, "cm"),panel.grid.major = element_blank(),
        axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 



tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/48H/PATHWAY.txt")
head(tmp)
tmp$value <- factor(tmp$value, levels = tmp$value)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/KEGG_NEW/48H/PATHWAY.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=nGene))+ylab("% Proteins")+scale_y_continuous(breaks = seq(0, 50, by = 20))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=15),plot.margin = margin(0.2,8,0.5,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y = element_blank(),axis.text.x = element_text(size=12,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=15),legend.title = element_text(size=15),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 

######################## Volcanoplots

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/24H_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/24H_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/48H_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/48H_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Time_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Time_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "yellow", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()


############################ Heatmap of Hif - targeted genes 

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/HIF_targetted/Genes_LCPM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/HIF_targetted/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/HIF_targetted/Zscore.txt",sep="\t",quote = FALSE,col.names = NA)


Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/HIF_targetted/Zscore.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/MetaData.txt",row.names = 1)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)
?Heatmap
LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/HIF_targetted/LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(SeverityScore = c("S1"="#ffa500","S2"="#cc0000"),Group=c("At_Infection"="#590059","12M_After"="#cc99cc")))
#col_fun2 = colorRamp2(c(2, 1, 0, -1, -2), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#5d5dff","#4c4cff","#0000cc"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),
           top_annotation  =ha,row_split = 2,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(37, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",10),rep("S2",14)))

col_fun_lfc = colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#006600","#198c19","#4ca64c","white","#ca8383","#c26e6e" ,"#b04545"))
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(1, "cm"),
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 11),height  = unit(37, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 


tt=H1 + H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/HIF_targetted/Heatmap.pdf",height = 22,width =20)
draw(tt, merge_legend = TRUE)
dev.off()

#######################################################

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Untitled Folder/24H_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Untitled Folder/24H_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Untitled Folder/48H_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Untitled Folder/48H_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=abs(logFC),size=-log10(adj.P.Val)))+scale_color_gradient(low = "yellow", high = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=abs(logFC),size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Untitled Folder/TIME_ISG.txt",row.names = 1)
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/VolcanoPlots/Untitled Folder/Time_ISG.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_point(data=subset(data, adj.P.Val>=0.05),aes(x=logFC, y=-log10(adj.P.Val),size=-log10(adj.P.Val)),color="#bfbfbf")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC < 0),aes(x=logFC, y=-log10(adj.P.Val),color=logFC,size=-log10(adj.P.Val)))+scale_color_gradient(high = "yellow", low = "#2c766f")+
  geom_point(data=subset(data, adj.P.Val<0.05 & logFC > 0),shape=21,aes(x=logFC, y=-log10(adj.P.Val),fill=logFC,size=-log10(adj.P.Val)),color="transparent")+scale_fill_gradient(low = "yellow", high = "#b20000")+
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=4.5,label.size = NA,segment.alpha=0.75,
                   box.padding=0,nudge_x = 0.1,nudge_y = 0.2)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),axis.title.y=element_text(size=15),legend.position = "none",
        axis.title.x=element_text(size=15),axis.text=element_text(size=10,color="black"),plot.margin = margin(1.5,1.5,1.5,1, "cm"))+
  labs(x="Log2 fold change",y="-log10 (adj.Pvalue)")+guides(size=guide_legend(override.aes=list(colour="grey")))
dev.off()
?comleheatmap
library(reshape2)
library(ggplot2)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/LineT.txt",check.names = FALSE)
head(data)
M=melt(data)
head(M)
write.table(M,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/LineT.txt",sep = "\t",quote = FALSE,col.names = NA)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Line.pdf")
ggplot(data=data, aes(x=variable, y=value, group=Gene)) +
  geom_line(aes(color=Gene))+
  geom_point()+labs(y="Log2 FoldChange")+theme(axis.title.x = element_blank(),plot.margin = margin(3,3,3,2, "cm"),legend.position = c(0.9,0.18))
dev.off()

geom_point(data=subset(M, variable==`0h`),aes(x=variable, y=value),color="#ffffb2",size=1.2)




################### Organoid data heatmap

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_Lung.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_Lung.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_LungZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_LungZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Design.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)


head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Time = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Time=c("1.Non-Infected (Day1)"="#addbad","2.Conc_10^3 (Day1)"="#5cb85c","3.Conc_10^6 (Day1)"="#2e5c2e",
                                         "4.Non-Infected (Day3)"="#addfee","5.Conc_10^3 (Day3)"="#5bc0de","6.Conc_10^6 (Day3)"="#367385",
                                         "7.Non-Infected (Day6)"="#a3a3ff","8.Conc_10^3 (Day6)"="#6666ff","9.Conc_10^6 (Day6)"="#4747b2")))

MET=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Annot.txt",row.names = 1)
head(MET)
ha2 = columnAnnotation(df = MET,show_annotation_name = FALSE,annotation_legend_param = list(Pathway = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                    col = list(Pathway=c("3.Citrate cycle (TCA cycle)"="#800000","2.Fructose and mannose metabolism"="#f07aaa","1.Glycolysis / Gluconeogenesis"="#FFA500")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 10),height  = unit(38, "cm"),width  = unit(10, "cm"),
           column_split =c(rep("Day1",3),rep("Day3",3),rep("Day6",3)),row_split =c(rep("A",63),rep("B",20),rep("C",23)))


LFC=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/lung.txt",row.names = 1)
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(2, "cm"),show_row_names = TRUE,
           top_annotation = ha2,show_heatmap_legend = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_names_gp =gpar(fontsize = 10),height  = unit(38, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6",
           row_split =c(rep("A",63),rep("B",20),rep("C",23))) 
t=H2+H1


pdf("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_Lung.pdf",height = 20,width =25)
draw(t,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


##########

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_SI.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_SI.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_SI_Z.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_SI_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Design.txt",row.names = 1)


library(ComplexHeatmap)
library(circlize)


head(Zscore)
ha = HeatmapAnnotation(df = sampleinfo,show_annotation_name = FALSE,annotation_name_side = "right",annotation_legend_param = list(Time = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Time=c("1.Non-Infected (Day1)"="#addbad","2.Conc_10^3 (Day1)"="#5cb85c","3.Conc_10^6 (Day1)"="#2e5c2e",
                                         "4.Non-Infected (Day3)"="#addfee","5.Conc_10^3 (Day3)"="#5bc0de","6.Conc_10^6 (Day3)"="#367385",
                                         "7.Non-Infected (Day6)"="#a3a3ff","8.Conc_10^3 (Day6)"="#6666ff","9.Conc_10^6 (Day6)"="#4747b2")))

MET=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Annot.txt",row.names = 1)
head(MET)
ha2 = columnAnnotation(df = MET,show_annotation_name = FALSE,annotation_legend_param = list(Pathway = list(direction = "horizontal",grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 18), labels_gp = gpar(fontsize = 18))),
                       col = list(Pathway=c("3.Citrate cycle (TCA cycle)"="#800000","2.Fructose and mannose metabolism"="#f07aaa","1.Glycolysis / Gluconeogenesis"="#FFA500")))


col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#bf7fbf","#993299","#590059"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),
           top_annotation = ha,heatmap_legend_param =list(grid_width = unit(1, "cm"),grid_height = unit(1, "cm"),title_gp = gpar(fontsize = 15), labels_gp = gpar(fontsize = 15)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 10),height  = unit(38, "cm"),width  = unit(10, "cm"),
           column_split =c(rep("Day1",3),rep("Day3",3),rep("Day6",3)),row_split =c(rep("A",63),rep("B",20),rep("C",23)))


LFC=read.delim("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/SI.txt",row.names = 1)
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="Log2FoldChange",width  = unit(2, "cm"),show_row_names = TRUE,
           top_annotation = ha2,show_heatmap_legend = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_names_gp =gpar(fontsize = 10),height  = unit(38, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6",
           row_split =c(rep("A",63),rep("B",20),rep("C",23))) 
t=H2+H1


pdf("/home/anoop/Desktop/Organoid/BEA20P106_Ujjwal_Neogi/HeatMap/ForUjjwal/Pathway_SI.pdf",height = 20,width =25)
draw(t,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()



library(NormalyzerDE)
?normalyzer
normalyzer(jobName="Prot",designPath = "/home/anoop/Desktop/Others/tmp/Info.txt",
           dataPath = "/home/anoop/Desktop/Others/tmp/Data_Filt.txt",
           outputDir = "/home/anoop/Desktop/Others/tmp/"
)

time.grp <- rep(c(24, 48, 72), 3)

assay <- c(rep("Sample..1", 3),
           rep("Sample..2", 3),
           rep("Sample..3", 3))

groups <- as.factor(assay)

design <- model.matrix( ~0 + groups + time.grp)
?makeContrasts

################################ NEW Analysis as per MS revision


library(piano)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/Names_Result.txt",header = TRUE,row.names = 2)
p <- data[7]
lfc<-data[4]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/PianoResults.txt",sep="\t",col.names = NA)

?runGSA


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/Res.txt",row.names = 1,sep = "\t")
head(data)

geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,segment.alpha=1,
                 box.padding=0.95,nudge_x = 1.8,nudge_y = -0.01)
  

?networkPlot2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/piano.pdf")
ggplot(data, aes(x=rank, y=Per)) + theme_bw()+scale_y_continuous(limits = c(-90, 90), breaks = seq(-80, 90, by = 20))+
  geom_point(data=subset(data, Per>0),aes(x=rank, y=Per, size=lpval),color="#90001c",fill="#ba0024")+
  geom_point(data=subset(data, Per<0),aes(x=rank, y=Per, size=lpval),color="#004000",fill="#008000")+
  scale_size(range = c(1, 5))+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        axis.title.y=element_text(size=15),legend.position = "bottom",axis.ticks.x = element_blank(),panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),axis.text.x =element_text(size=10,color="black"),plot.margin = margin(3,0.5,3,0.5, "cm"))+
  guides(size=guide_legend(override.aes=list(color="grey",fill="grey")))
dev.off()


table_sel <- res[, c('p adj (dist.dir.dn)', 'p adj (mix.dir.dn)', 'p adj (non-dir.)',
                              'p adj (mix.dir.up)', 'p adj (dist.dir.up)')]

head(table_sel)

table_sel[is.na(table_sel)] <- 1  # set NA p-values equal to 1
rownames(table_sel) <- res$Name
table_sel <- table_sel[apply(table_sel, 1, min) < 0.05, ]

library(pheatmap)
pheatmap(-log10(table_sel),  # log-transform p-values
         scale='none',
         cluster_cols=F,
         clustering_distance_rows='euclidean',
         clustering_method='ward.D2',
         breaks=seq(0, 2, len=100),
         border_color='black',
         angle_col=90)


############################################## Paino Acute -vs- Recovered  Heatmap
library(ComplexHeatmap)
library(circlize)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/Ip_HeatMap.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,0.5,0.75,1,1.25,1.5,1.75,2), c("#d8d8d8" ,"#b2b266","#989832","#7f7f00","#ac7fac","#8a4c8a","#691969","#590059"))
sf=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/label.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Label = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                  grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8))),
                       col = list(Label=c("a) p adj (mix.dir.dn)"="#a6a6a6",
                                          "b) p adj (mix.dir.up)"="#8c8c8c","c) p adj (non-dir.)"="#808080",
                                          "d) p adj (dist.dir.dn)"="#666666","e) p adj (dist.dir.up)"="#4c4c4c")))

H1=Heatmap(as.matrix((data)),top_annotation = ha,col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 8),height  = unit(16, "cm"),width  = unit(4, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Piano/Piano_HeatMap.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE)
dev.off()


############################################ IFN signalling proteomics heatmap


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_Quantile.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_Quantile.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_QuantileZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_QuantileZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/metadata2.txt",row.names = 1,check.names = FALSE)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(Condition=c("Untreated"="#cc99cc","Treated"="#590059"),`Time point` = c("24h"="#ffa500","48h"="#cc0000")))


col_fun2=colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200"))

pwy=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_2.txt",row.names = 1,header = TRUE)
H1=Heatmap(as.matrix((Zscore)),col=col_fun2,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 15),row_split = 3,row_title_gp = gpar(fontsize = 0),
           top_annotation  =ha,row_gap = unit(0, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(60, "cm"),width  = unit(20, "cm"),
           column_split =sampleinfo$`Time point`)

col_fun_lfc = colorRamp2(c(-2, -1, 0, 1, 2), c("#004c00","#008000" ,"white","#ab6d2e","#66411b"))

H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(2, "cm"),
           show_row_names = TRUE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 15),height  = unit(60, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#8c8c8c") 


rna=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_LFC_RNASeq.txt",check.names = FALSE,row.names = 1)
H4=Heatmap(as.matrix((rna)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="RNA_LFC",width  = unit(1, "cm"),row_names_side = "right",show_heatmap_legend = FALSE,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 15),height  = unit(60, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#8c8c8c") 

tt=H1+H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN.pdf",height = 28,width =20)
draw(tt, merge_legend = TRUE,auto_adjust = FALSE)
dev.off()

###################
virus=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/virus.txt",row.names = 1,check.names = FALSE)
Vi=log10(virus)
write.table(M,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/Sele_Melt.txt",sep = "\t",col.names = NA,quote = FALSE)


sel=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/Selected.txt",check.names = FALSE)
head(sel)
library(reshape2)
M=melt(sel)
head(M,30)

tt=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/Sele_Melt.txt",header = TRUE)
head(tt)
library(ggplot2)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/Sele_Melt.pdf")
ggplot(tt, aes(x=variable2, y=value, group=Group)) +facet_wrap(~ Gene, ncol = 3,nrow = 3)+
  geom_line(aes(color=Group),size=1.5)+scale_color_manual(values=c("Avg"="#C44354", "R1"="#edc6cb", "R2"="#edc6cb","R3"="#edc6cb"))+
  geom_point(aes(color=Group),size=1)+theme(legend.position = "none",axis.title = element_blank())
dev.off()

#######################


library(piano)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Piano/48H.txt",header = TRUE,row.names = 8)
head(data)
p <- data[5]
lfc<-data[2]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Piano/PianoResults48H.txt",sep="\t",col.names = NA)



##################
library(reshape2)
library(ggplot2)
tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Piano/For_figure.txt")
head(tmp)
M=melt(tmp)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Piano/CommonPathway.pdf")
ggplot(tmp) + 
  geom_point(data=tmp,aes(x=Type,y=Term,size=prot_ratio,color=pvalue))+scale_color_gradient(low="#57a1ab",high="#346066")+
  theme(axis.title = element_blank(),axis.text.y = element_text(size=10,color="black"),
        axis.ticks = element_blank(),legend.position ="bottom",legend.box="vertical",plot.margin = margin(2,5,7,2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Feature ratio",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()


################################################ UV proteomics


library(NormalyzerDE)
?normalyzer
normalyzer(jobName="UV",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/metadat.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Raw_UV.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized"
)


library(impute)
?impute.knn
BiocManager::install("impute")
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Quantile-normalized.txt",row.names = 1,check.names = FALSE)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Imputed.txt",sep="\t",col.names = NA,quote = FALSE)

library(sva)
prot=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Imputed.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(prot)
des=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/metadat.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Batch_Corrected.txt", sep="\t",quote=FALSE,col.names = NA)


library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Batch_Corrected.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/metadat.txt",row.names = 1)
p <- PCAtools::pca(count, metadata = meta, removeVar = 0.1)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Batch_Corrected_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)
head(p$variance)

library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Imputed.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/metadat.txt",row.names = 1)
p <- PCAtools::pca(count, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Imputed_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)



dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Batch_Corrected_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Batch_Corrected_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",Treated="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",Treated="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 27.88% variance",y="PC2, 12.92% variance")+
  theme(axis.title = element_text(size=10),legend.position = c(0.17, 0.88),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Normalized/Imputed_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/CCR6/NonFraction/Analysis/Normalized/Batch_Corrected_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=5,aes(fill=Group,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 8, 22))+
geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 23.95% variance",y="PC2, 11.68% variance")+
  theme(axis.title = element_text(size=10),legend.position = c(0.15, 0.12),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 2,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()

#########
rm(list = ls())
#################################

library(limma)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Batch_Corrected.txt",row.names = 1,check.names = FALSE)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/metadat.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(Mock_vs_UV_24H = (groupUV_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Mock_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_24H = (groupTreated_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Mock_vs_Treated_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_24H = (groupUV_24h - groupTreated_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Treated_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)

packageVersion("impute")
#################

cont.matrix <-makeContrasts(Mock_vs_UV_48H = (groupUV_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Mock_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_48H = (groupTreated_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Mock_vs_Treated_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_48H = (groupUV_48h - groupTreated_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Treated_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)
##################

library(piano)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Named_Mock_vs_Treated_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Kegg/Mock_vs_Treated_24H.txt",sep="\t",col.names = NA)

#############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Named_Mock_vs_Treated_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Kegg/Mock_vs_Treated_48H.txt",sep="\t",col.names = NA)

###############
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Named_Mock_vs_UV_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Kegg/Mock_vs_UV_24H.txt",sep="\t",col.names = NA)
############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Named_Mock_vs_UV_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Kegg/Mock_vs_UV_48H.txt",sep="\t",col.names = NA)

################

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Named_Treated_vs_UV_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Kegg/Treated_vs_UV_24H.txt",sep="\t",col.names = NA)
##############
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Named_Treated_vs_UV_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/Limma/Kegg/Treated_vs_UV_48H.txt",sep="\t",col.names = NA)
#################################

dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Norm/UV/log2-normalized.txt",row.names = 1,check.names = FALSE)
count=log2(dd)
head(count)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/metadata.txt",row.names = 1)
p <- PCAtools::pca(dd, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Norm/UV/log2-normalized_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)

library(ggplot2)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Norm/UV/X.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Log_PCA_2.pdf")
ggplot(dat, aes(x=x, y=y,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",CCHFV="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",CCHFV="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=x,y=y,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 27.88% variance",y="PC2, 12.92% variance")+
  theme(axis.title = element_text(size=10),legend.position ="none",plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()

############################ Limma PD Norma ########################

library(limma)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Batch_Corrected.txt",row.names = 1,check.names = FALSE)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/metadat.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(Mock_vs_UV_24H = (groupUV_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Mock_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_24H = (groupTreated_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Mock_vs_Treated_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_24H = (groupUV_24h - groupTreated_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Treated_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


#################

cont.matrix <-makeContrasts(Mock_vs_UV_48H = (groupUV_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Mock_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_48H = (groupTreated_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Mock_vs_Treated_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_48H = (groupUV_48h - groupTreated_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma/Treated_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)
##################


library(NormalyzerDE)
?normalyzer
normalyzer(jobName="UV",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/metadata.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Batch_Corrected.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Norm"
)

dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Norm/UV/log2-normalized.txt",row.names = 1, check.names = FALSE)
d <- dist(t(dd),method = "manhattan")
head(d)
fit <- cmdscale(d,eig=TRUE, k=2)
head(fit)
x <- fit$points[,1]
y <- fit$points[,2]

plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
     main="Metric MDS", type="n")
text(x, y, labels = row.names(t(dd)), cex=.7)



fit <- stats::cmdscale(d, eig=TRUE, k=2)
x <- fit$points[, 1]
y <- fit$points[, 2]
graphics::plot(x, y, type="n", main=methodnames[i], xlab="", ylab="")
graphics::text(
  fit$points[, 1], 
  fit$points[, 2], 
  col=filterED, 
  labels=filterED
)


d <- stats::dist(scale(t(stats::na.omit((dd))), center=TRUE, scale=TRUE))

#################

############################ Limma PD Norma +Log ########################

library(limma)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/log2-normalized.txt",row.names = 1,check.names = FALSE)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/metadat.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(Mock_vs_UV_24H = (groupUV_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Mock_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_24H = (groupTreated_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Mock_vs_Treated_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_24H = (groupUV_24h - groupTreated_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Treated_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


#################

cont.matrix <-makeContrasts(Mock_vs_UV_48H = (groupUV_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Mock_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_48H = (groupTreated_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Mock_vs_Treated_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_48H = (groupUV_48h - groupTreated_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Treated_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)
##################

library(piano)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Named_Mock_vs_Treated_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Kegg/Mock_vs_Treated_24H.txt",sep="\t",col.names = NA)

#############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Named_Mock_vs_Treated_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Kegg/Named_Mock_vs_Treated_48H.txt",sep="\t",col.names = NA)

###############
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Named_Mock_vs_UV_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Kegg/Named_Mock_vs_UV_24H.txt",sep="\t",col.names = NA)
############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Named_Mock_vs_UV_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Kegg/Named_Mock_vs_UV_48H.txt",sep="\t",col.names = NA)

################

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Named_Treated_vs_UV_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Kegg/Named_Treated_vs_UV_24H.txt",sep="\t",col.names = NA)
##############
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Named_Treated_vs_UV_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/Kegg/Named_Treated_vs_UV_48H.txt",sep="\t",col.names = NA)
#################################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/heatmap/Oxpho_Gene.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/heatmap/Oxpho_Gene.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/heatmap/Oxpho_GeneZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/heatmap/Oxpho_GeneZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/heatmap/metaData.txt",row.names = 1,check.names = FALSE)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)



ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(`Timepoint` = c("24H"="#cbcba9","48H"="#656554"),Condition=c("Mock"="#cc99cc","CCHFV"="#590059","UV"="#4682B4")))

pwy=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/IFN_2.txt",row.names = 1,header = TRUE)

ha2 = rowAnnotation(df = pwy,annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(Complex= c("Complex1"="#16605f","Complex2"="#891b69","Complex3"="#e5ba33","Complex4"="#332a6c","Complex5"="#d17294")))

col_fun2=colorRamp2(c(-2,-1.5,-1, 0,1,1.5,2), c("#4c4cff","#5d5dff","#6f6fff" ,"white","#c1c132","#b9b919","#b2b200"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2A4E6C","#4682B4" ,"#FFFCC9","#FFA700", "#FF5A00"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun2,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,left_annotation = ha2,
  
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),row_split = pwy$Complex,row_title_gp = gpar(fontsize = 0),
           top_annotation  =ha,row_gap = unit(3, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 15),height  = unit(40, "cm"),width  = unit(20, "cm"),
             column_split =sampleinfo$Timepoint)

tt=H1+H2+H4
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Limma_Log/heatmap/Oxpho_Gene.pdf",height = 28,width =20)
draw(H1, merge_legend = TRUE,auto_adjust = FALSE)
dev.off()

##############################################################
library(impute)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/PD_Norm.txt",row.names = 1,check.names = FALSE)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/PD_Norm_Imputed.txt",sep="\t",col.names = NA,quote = FALSE)


library(NormalyzerDE)
?normalyzer
normalyzer(jobName="UV",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/metadata.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/PD_Norm_Imputed.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Norm"
)

######################################################################################## PD Norm 2Replicate

library(limma)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/log2-normalized.txt",row.names = 1,check.names = FALSE)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/metadata.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(Mock_vs_UV_24H = (groupUV_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Mock_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_24H = (groupTreated_24h - groupMock_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Mock_vs_Treated_24H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_24H = (groupUV_24h - groupTreated_24h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Treated_vs_UV_24H.txt",sep="\t",quote = FALSE,col.names = NA)


#################

cont.matrix <-makeContrasts(Mock_vs_UV_48H = (groupUV_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Mock_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Mock_vs_Treated_48H = (groupTreated_48h - groupMock_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Mock_vs_Treated_48H.txt",sep="\t",quote = FALSE,col.names = NA)


cont.matrix <-makeContrasts(Treated_vs_UV_48H = (groupUV_48h - groupTreated_48h),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Treated_vs_UV_48H.txt",sep="\t",quote = FALSE,col.names = NA)

##################

library(piano)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Named_Mock_vs_Treated_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Kegg/Mock_vs_Treated_24H.txt",sep="\t",col.names = NA)

#############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Named_Mock_vs_Treated_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Kegg/Mock_vs_Treated_48H.txt",sep="\t",col.names = NA)

###############
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Named_Mock_vs_UV_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Kegg/Mock_vs_UV_24H.txt",sep="\t",col.names = NA)
############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Named_Mock_vs_UV_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Kegg/Mock_vs_UV_48H.txt",sep="\t",col.names = NA)

################

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Named_Treated_vs_UV_24H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Kegg/Treated_vs_UV_24H.txt",sep="\t",col.names = NA)
##############
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Named_Treated_vs_UV_48H.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/2_Replicate/Kegg/Treated_vs_UV_48H.txt",sep="\t",col.names = NA)
###################################################################

library(sva)
prot=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_Raw.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(na.omit(prot))
des=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/metadata.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PCA_NoImp.txt", sep="\t",quote=FALSE,col.names = NA)


dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/BatchCorrected_NoImp.txt",row.names = 1,check.names = FALSE)

meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/metadata.txt",row.names = 1)
p <- PCAtools::pca(dd, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PCA.txt",sep = "\t",col.names = NA,quote = FALSE)

library(ggplot2)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",CCHFV="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",CCHFV="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 84.70% variance",y="PC2, 6.1% variance")+
  theme(axis.title = element_text(size=10),legend.position =c(0.9,0.12),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()

######


dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/Batch_Corrected_Impu.txt",row.names = 1,check.names = FALSE)

meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/metadata.txt",row.names = 1)
p <- PCAtools::pca(dd, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/Batch_Corrected_ImpuPCA.txt",sep = "\t",col.names = NA,quote = FALSE)

library(ggplot2)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/Batch_Corrected_ImpuPCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/Batch_Corrected_ImpuPCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",CCHFV="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",CCHFV="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 84.59% variance",y="PC2, 6.18% variance")+
  theme(axis.title = element_text(size=10),legend.position =c(0.9,0.12),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()

##################################################################

library(impute)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Raw.txt",row.names = 1,check.names = FALSE)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Imputed.txt",sep="\t",col.names = NA,quote = FALSE)


dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Imputed.txt",row.names = 1,check.names = FALSE)

meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/metadata.txt",row.names = 1)
p <- PCAtools::pca(dd, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Imputed_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)

library(ggplot2)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Imputed_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Imputed_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",CCHFV="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",CCHFV="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 61.44% variance",y="PC2, 27.35% variance")+
  theme(axis.title = element_text(size=10),legend.position =c(0.85,0.92),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()

#######################################################################

library(NormalyzerDE)
?normalyzer
normalyzer(jobName="UV",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/metadata.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/PD_Norm_Raw.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Norm"
)


library(impute)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized.txt",row.names = 1,check.names = FALSE)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed.txt",sep="\t",col.names = NA,quote = FALSE)


library(sva)
prot=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(na.omit(prot))
des=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/metadata.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$group)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(rnaseqCombat,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed_BatchCorr.txt", sep="\t",quote=FALSE,col.names = NA)



dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed_BatchCorr.txt",row.names = 1,check.names = FALSE)

meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/metadata.txt",row.names = 1)
p <- PCAtools::pca(dd, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed_BatchCorrPCA.txt",sep = "\t",col.names = NA,quote = FALSE)

library(ggplot2)
library(ggrepel)
dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed_BatchCorrPCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_BatchCorr/Quantile-normalized_Imputed_BatchCorrPCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",CCHFV="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",CCHFV="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 27.89% variance",y="PC2, 12.91% variance")+
  theme(axis.title = element_text(size=10),legend.position =c(0.15,0.92),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()
########################################################################

library(impute)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Batch1.txt",row.names = 1,check.names = FALSE)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Batch1_imp.txt",sep="\t",col.names = NA,quote = FALSE)


dd=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Batch1_imp.txt",row.names = 1,check.names = FALSE)

meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/metadata2.txt",row.names = 1)
p <- PCAtools::pca(dd, metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Batch1_imp_PCA.txt",sep = "\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Batch1_imp_PCA.txt",row.names = 1)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/Untitled Folder/PD_Norm_NoBatch/PD_Norm_Batch1_imp_PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Cohort)) + geom_point(size=5,aes(fill=Cohort,shape=Time))+theme_gray()+
  scale_shape_manual(values=c(21, 22, 8))+
  scale_color_manual(values=c(Mock="#457979",CCHFV="#cc8400",UV="#0000e5",Pool="#666666"))+
  scale_fill_manual(values=c(Mock="#5ca2a2",CCHFV="#ffa500",UV="#1919ff",Pool="#808080"))+geom_text_repel(data=dat,aes(x=PC1,y=PC2,label = rownames(dat)),size=1.9,box.padding=0.3,show.legend=FALSE,colour="black")+
  labs(x="PC1, 87.81% variance",y="PC2, 5.25% variance")+
  theme(axis.title = element_text(size=10),legend.position =c(0.85,0.9),plot.margin = margin(1,1,1,1, "cm"),
        legend.title=element_blank(),legend.text=element_text(size=10),legend.key.size = unit(0.5, "cm"))+
  guides(shape = guide_legend(nrow = 1,override.aes=list(fill="grey",color="grey")),color = guide_legend(nrow = 2),fill = guide_legend(nrow = 2))
dev.off()

#####################


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/IFN_Avg.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/IFN_Data2.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/IFN_AVgZ.txt",sep="\t",quote = FALSE,col.names = NA)


Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/IFN_AVgZ.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Heatmap/IFN/metadata.txt",row.names = 1,check.names = FALSE)
head(sampleinfo)

Zscore2=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/IFN_LFC.txt", row.names = 1,check.names = FALSE)

library(ComplexHeatmap)
library(circlize)


ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(`Time point` = c("24h"="#ffa500","48h"="#cc0000"),Condition=c("Untreated"="#cc99cc","Treated"="#590059")))

col_fun1 = colorRamp2(c(14,16,18,20,22,24,26,28), c("#7cc34c","#44aa00","#4c4cff","#0000ff","#0000e5","#d67c6d","#c8503c","#a8210a"))
col_fun1 = colorRamp2(c(-0.75, -0.25,0, 0.25,0.75), c("#004c00","#008000","white","#e50000","#7f0000"))

pwy=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/Pwy2.txt",row.names = 1,header = TRUE)
H1=Heatmap(as.matrix((Zscore2)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = TRUE,column_names_side = "top",
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 20),row_split = pwy$Pwy,row_title_gp = gpar(fontsize = 25),
          row_gap = unit(3, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                              title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "LFC",show_row_names = TRUE,row_names_gp=gpar(fontsize = 15),height  = unit(50, "cm"),width  = unit(10, "cm")
           )

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/UV_Proteomics/NEW_DATA/ISG/IFN_DataLFC.pdf",height = 28,width =20)
draw(H1, merge_legend = TRUE,auto_adjust = FALSE)
dev.off()

##################################### Proteomics New Piano

packageVersion("piano")
library(piano)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/LimmaResultsALL_NEW.txt",header = TRUE,row.names = 8)
p <- data[5]
head(p)
lfc<-data[2]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
gsares$geneSetStat
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/24h_48H_piano.txt",sep="\t",col.names = NA,quote = FALSE)


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/Results_24h.txt",header = TRUE,row.names = 8)
p <- data[5]
head(p)
lfc<-data[2]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/Piano_24h.txt",sep="\t",col.names = NA,quote = FALSE)


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/Results_48h.txt",header = TRUE,row.names = 8)
p <- data[5]
head(p)
lfc<-data[2]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/Piano_48h.txt",sep="\t",col.names = NA,quote = FALSE)
##########################################
?runGSA

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/LimmaResultsALL_NEW.txt",header = TRUE,row.names = 8)
t <- data[4]
head(t)
geneSets <- loadGSC(file = "/home/anoop/Desktop/CCR6/NEW_Analysis/Fraction/ART/KEGG/KEGG.gmt")
gsares <- runGSA(geneLevelStats=t,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
head(res)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/LimmaResultsALL_piano.txt",sep="\t",col.names = NA)


##############

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/prot_c1.txt")
head(tmp)
tmp$Term <- factor(tmp$Term, levels = tmp$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/prot_c1.pdf")
ggplot(tmp) +
  geom_point(data=tmp,aes(x=1,y=Term,size=ratio,color=log))+scale_color_gradient(low="#7f7fff",high="#0000e5")+
  scale_x_discrete(limits=c("1"))+
  theme(axis.title = element_blank(),axis.text.y = element_text(size=10,color="black"),axis.text.x = element_blank(),
        axis.ticks = element_blank(),legend.position =c(-3,-0.2),legend.box="vertical",plot.margin = margin(3,8.4,6,2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Feature ratio",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 2,override.aes = list(size = 5)))
dev.off()


tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/prot_c6.txt")
head(tmp)
tmp$Term <- factor(tmp$Term, levels = tmp$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/prot_c6.pdf")
ggplot(tmp) +
  geom_point(data=tmp,aes(x=1,y=Term,size=ratio,color=log))+scale_color_gradient(low="#7fbf7f",high="#006600")+
  scale_x_discrete(limits=c("1"))+
  theme(axis.title = element_blank(),axis.text.y = element_text(size=10,color="black"),axis.text.x = element_blank(),
        axis.ticks = element_blank(),legend.position ="bottom",legend.box="vertical",plot.margin = margin(6,8.4,6,4.2, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Feature ratio",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()



tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/rna_c1.txt")
head(tmp)
tmp$Term <- factor(tmp$Term, levels = tmp$Term)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/rna_c1.pdf")
ggplot(tmp) +
  geom_point(data=tmp,aes(x=1,y=Term,size=ratio,color=log))+scale_color_gradient(low="#ff6666",high="#cc0000")+
  scale_x_discrete(limits=c("1"))+
  theme(axis.title = element_blank(),axis.text.y = element_text(size=10,color="black"),axis.text.x = element_blank(),
        axis.ticks = element_blank(),legend.position =c(-3,-0.2),legend.box="vertical",plot.margin = margin(4,7,6,2.5, "cm"),
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=10,colour = "black"),
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Feature ratio",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()


#####################

library(ComplexHeatmap)
library(circlize)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/pwy4.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data)
col_fun1 = colorRamp2(c(0,1.25,1.5,2,3), c("grey","#ac7fac","#8a4c8a","#691969","#590059"))
sf=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Info.txt",header = TRUE,row.names = 1)
ha = HeatmapAnnotation(df = sf,show_annotation_name = FALSE,annotation_name_side = "right",
                       annotation_legend_param = list(Analysis = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                       col = list(Analysis=c("24hpi"="#a6a6a6",
                                          "48hpi"="#8c8c8c","RNASeq"="#808080",
                                          "Time series"="#666666")))

H1=Heatmap(as.matrix((data)),top_annotation = ha,col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 15),height  = unit(20, "cm"),width  = unit(5, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/pwy.pdf",width = 10,height = 15)
draw(H1,heatmap_legend_side = "left", annotation_legend_side = "left",merge_legend = TRUE)
dev.off()
?gsa

########
library(reshape)
DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Signi/Top_R.txt")
head(DAT)
m=melt(DAT)
head(m)
Re <- cast(DAT)
head(Re)

write.table(Re,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Signi/Top_R_casted.txt",sep = "\t",col.names = NA,quote = FALSE)
?cast

########################
library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Untitled Folder/Comm_FPKM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Untitled Folder/Comm_FPKM.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Untitled Folder/Comm_FPKM_Z.txt",sep="\t",quote = FALSE,col.names = NA)


com=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Comm.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                      grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                       col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))

com=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Comm.txt",row.names = 1)
met=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/meta.txt",row.names = 1)
ha1 = HeatmapAnnotation(df = met,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(cohort = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(cohort=c("Acute"="#b20e66","Recovered"="#738248")))

zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Comm_FPKM_Z.txt",row.names = 1)

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2A4E6C","#4682B4" ,"#FFFCC9","#FFA700", "#FF5A00"))

col_funNew = colorRamp2(c(-2,-1, 0,1,2), c("#ffd700","#ffae19","#e59400","#cc0000","#990000"))

H1=Heatmap(as.matrix((zscore)),cluster_rows=FALSE,col=col_funNew,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           column_dend_height = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",left_annotation = ha,top_annotation = ha1,
           column_names_gp =gpar(fontsize = 3),
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),column_split = met$cohort,
           name = "Zscaled FPKM",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(8, "cm"),width  = unit(5, "cm"))

col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                        grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))

corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/NewTopCasted.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",top_annotation = ha3,left_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(7, "cm"),width  = unit(7, "cm"))

tt=H1+H2

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Casted.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


##############################



library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/Comm_Quantile.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/Comm_Quantile.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/Comm_QuantileZ.txt",sep="\t",quote = FALSE,col.names = NA)


com=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/Comm.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Comm = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(Comm=c("c1"="#4900ff","c2"="#00babf","c3"="#00c9b0","c4"="#00b8ff","c5"="#aa99ff","c6"="#bbff99","c7"="#99ffcc","c8"="#ff99ff")))

met=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/meta.txt",row.names = 1)
ha1 = HeatmapAnnotation(df = met,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Group = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                     grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Group=c("Mock 24h"="#588ebb","Mock 48h"="#738248","Infected 24h"="#2a4e6c","Infected 48h"="red")))

zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/Comm_QuantileZ.txt",row.names = 1)

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2A4E6C","#4682B4" ,"#FFFCC9","#FFA700", "#FF5A00"))
H1=Heatmap(as.matrix((zscore)),cluster_rows=FALSE,col=col_fun1,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",left_annotation = ha,top_annotation = ha1,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),column_split = met$cohort,
           name = "Zscore",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(8, "cm"),width  = unit(5, "cm"))

col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Comm = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                        grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Comm=c("c1"="#4900ff","c2"="#00babf","c3"="#00c9b0","c4"="#00b8ff","c5"="#aa99ff","c6"="#bbff99","c7"="#99ffcc","c8"="#ff99ff")))

corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/NewTop_R_casted.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),col=col_fun2,cluster_rows=FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",left_annotation  = ha,top_annotation = ha3,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(7, "cm"),width  = unit(7, "cm"))

tt=H1+H2

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/results/Species/Heatmap/Top_R_casted.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

#######################################################################

data2=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fig/H1.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data2)
col_fun3 = colorRamp2(c(0,0.25,0.5,0.75,1,1.5,3), c("white","#D69999","#CC7F7F","#B74C4C","#A31919","#990000","#890000"))
sf=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fig/H1_meta.txt",header = TRUE,row.names = 1)

nrow(sf)
ncol(data2)
ha = HeatmapAnnotation(df = sf, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param  = list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),
                                                       title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                       col = list(column=c("a) p adj (dist.dir.dn)"="#A6A6A6",
                                           "b) p adj (mix.dir.dn)"="#8C8C8C","c) p adj (non-dir.)"="#808080",
                                           "d) p adj (mix.dir.up)"="#666666","e) p adj (dist.dir.up)"="#4C4C4C"),Group=c("24hpi"="#586E7C","48hpi"="#6D5C78","RNASeq"="#4682B4")))


H1=Heatmap(as.matrix((data2)),top_annotation = ha,col=col_fun3,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_title_gp = gpar(fontsize=0),row_split = 2,row_gap = unit(0, "cm"),
           row_dend_width = unit(1, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",column_split = sf$Group,column_gap = unit(0.05, "cm"),
           heatmap_legend_param =list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "-log10(Padj)",show_row_names = FALSE,row_names_gp=gpar(fontsize = 8),height  = unit(12, "cm"),width  = unit(3, "cm"))


data3=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fig/H2.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data2)
col_fun3 = colorRamp2(c(0,0.25,0.5,0.75,1,1.5,3), c("white","#D69999","#CC7F7F","#B74C4C","#A31919","#990000","#890000"))
sf2=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fig/H2_Meta.txt",header = TRUE,row.names = 1)

nrow(sf)
ncol(data2)
ha2 = HeatmapAnnotation(df = sf2, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param  = list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),
                                                       title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                       col = list(column=c("a) p adj (dist.dir.dn)"="#A6A6A6",
                                           "b) p adj (mix.dir.dn)"="#8C8C8C","c) p adj (non-dir.)"="#808080",
                                           "d) p adj (mix.dir.up)"="#666666","e) p adj (dist.dir.up)"="#4C4C4C"),Group=c("Time series"="#308f5a")))


H2=Heatmap(as.matrix((data3)),top_annotation = ha2,col=col_fun3,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",
           heatmap_legend_param =list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 6),height  = unit(12, "cm"),width  = unit(1, "cm"))


tt=H1+H2

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fig/H1.pdf")
draw(tt,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE,ht_gap=unit(0.1, "cm"))
dev.off()


###########################################################################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Comm_FPKM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Comm_FPKM.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Comm_FPKM_Z.txt",sep="\t",quote = FALSE,col.names = NA)


com=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Comm.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                        grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))


met=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/meta.txt",row.names = 1)
ha1 = HeatmapAnnotation(df = met,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(cohort = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                     grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(cohort=c("Acute"="#b20e66","Recovered"="#738248")))

zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Comm_FPKM_Z.txt",row.names = 1)

col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7F7F00","#B2B200" ,"#E5E500","white","#BF7FBF","#993299","#590059"))
col_fun1 = colorRamp2(c(2, 1, 0, -1, -2), c("#2A4E6C","#4682B4" ,"#FFFCC9","#FFA700", "#FF5A00"))

col_funNew = colorRamp2(c(-2,-1, 0,1,2), c("#ffd700","#ffae19","#e59400","#cc0000","#990000"))

H1=Heatmap(as.matrix((zscore)),cluster_rows=FALSE,col=col_funNew,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           column_dend_height = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",left_annotation = ha,top_annotation = ha1,
           column_names_gp =gpar(fontsize = 3),
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),column_split = met$cohort,
           name = "Zscaled FPKM",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(8, "cm"),width  = unit(5, "cm"))

col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))
corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Top_R_casted.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha3,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(8, "cm"),width  = unit(6, "cm"))

tt=H1+H2

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Top_15/Comm_FPKM.pdf")
draw(tt,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


########################### Correlation signal and metabolism
library(psych)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Metabolism_FPKM.txt",header = TRUE)
signal=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Signaling_FPKM.txt",header = TRUE)

Res=corr.test(as.matrix(meta),as.matrix(signal),use = "pairwise",method="spearman",adjust="BH")
cor1=melt(Res$r)
pval=melt(Res$p)

head(pval)
write.table(pval,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/RNA_pval.txt",sep="\t",col.names = NA,quote = FALSE)
write.table(cor1,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/RNA_corr.txt",sep="\t",col.names = NA,quote = FALSE)



library(reshape)
DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/RNASeq_R.txt")
head(DAT)
m=melt(DAT)
head(m)
Re <- cast(DAT)
head(Re)

write.table(Re,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/RNASeq_R_casted.txt",sep = "\t",col.names = NA,quote = FALSE)

met=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/MetaGenes6.txt",row.names = 1)
ha = rowAnnotation(df = met,show_annotation_name = FALSE,simple_anno_size = unit(0.13, "cm"),na_col="white",
                   annotation_legend_param = list(direction = "horizontal",grid_width = unit(0.2, "cm"),
                                                              grid_height = unit(0.2, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8)),
                   col = list(TCA.cycle=c("TCA cycle"="#d62976"),Oxphos=c("Oxphos"="#962fbf"),
                              Glycolysis=c("Glycolysis"="#6d2606"),Amino.and.nucleotide.sugar=c("Amino and Nucleotide"="#d0bd09"),
                              N.Glycan=c("N-glycan"="#005900"),Cys.and.meth=c("Cys and meth"="#00cc00")))

sig=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/SingalGenes4.txt",row.names = 1)
ha2 = HeatmapAnnotation(df = sig,show_annotation_name = FALSE,simple_anno_size = unit(0.13, "cm"),na_col="white",
                   annotation_legend_param = list(direction = "horizontal",grid_width = unit(0.2, "cm"),
                                                  grid_height = unit(0.2, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8)),
                   col = list(Notch=c("Notch"="#f2973a"),FoxO=c("FoxO"="#3a144b"),
                              ErbB=c("ErbB"="#3c7c85"),Phosphatidylinositol=c("Phosphatidylinositol"="#67ae88")))

col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))
corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/RNASeq_R_casted.txt",row.names = 1)
H1=Heatmap(as.matrix((corr)),cluster_rows=TRUE,col=col_fun2,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(0.5, "cm"),column_dend_height = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",right_annotation = ha,
           bottom_annotation = ha2,
           heatmap_legend_param =list(grid_width = unit(0.2, "cm"),grid_height = unit(0.2, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(5, "cm"),width  = unit(5, "cm"))



pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/RNASeq_R.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

#####################

meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/MetIP.txt",header = TRUE)
signal=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Signal_ip.txt",header = TRUE)

Res=corr.test(as.matrix(meta),as.matrix(signal),use = "pairwise",method="spearman",adjust="BH")
cor1=melt(Res$r)
pval=melt(Res$p)

head(pval)
write.table(pval,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Prot_pval.txt",sep="\t",col.names = NA,quote = FALSE)
write.table(cor1,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Prot_corr.txt",sep="\t",col.names = NA,quote = FALSE)


library(reshape)
DAT=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Prot_R.txt")
head(DAT)
m=melt(DAT)
head(m)
Re <- cast(DAT)
head(Re)

write.table(Re,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Prot_R_casted.txt",sep = "\t",col.names = NA,quote = FALSE)



met=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Meta_Genes6.txt",row.names = 1)
ha = rowAnnotation(df = met,show_annotation_name = FALSE,simple_anno_size = unit(0.13, "cm"),na_col="white",
                   annotation_legend_param = list(direction = "horizontal",grid_width = unit(0.2, "cm"),
                                                  grid_height = unit(0.2, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8)),
                   col = list(TCA.cycle=c("TCA cycle"="#d62976"),Oxphos=c("Oxphos"="#962fbf"),
                              Glycolysis=c("Glycolysis"="#6d2606"),Amino.and.nucleotide.sugar=c("Amino and Nucleotide"="#d0bd09"),
                              N.Glycan=c("N-glycan"="#005900"),Cys.and.meth=c("Cys and meth"="#00cc00")))

sig=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/SignalGenes4.txt",row.names = 1)
ha2 = HeatmapAnnotation(df = sig,show_annotation_name = FALSE,simple_anno_size = unit(0.13, "cm"),na_col="white",
                        annotation_legend_param = list(direction = "horizontal",grid_width = unit(0.2, "cm"),
                                                       grid_height = unit(0.2, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 8)),
                        col = list(Notch=c("Notch"="#f2973a"),FoxO=c("FoxO"="#3a144b"),
                                   ErbB=c("ErbB"="#3c7c85"),Phosphatidylinositol=c("Phosphatidylinositol"="#67ae88")))

col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))
corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Prot_R_casted.txt",row.names = 1)
H1=Heatmap(as.matrix((corr)),cluster_rows=TRUE,col=col_fun2,cluster_columns = TRUE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(0.5, "cm"),column_dend_height = unit(0.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",right_annotation = ha,
           bottom_annotation = ha2,
           heatmap_legend_param =list(grid_width = unit(0.2, "cm"),grid_height = unit(0.2, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(5, "cm"),width  = unit(5, "cm"))



pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Pwy_Corr/Proteomics/Prot_R_casted.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

###############


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fishers/Results3.txt",header = TRUE,check.names = FALSE)
head(data)
data$BH =p.adjust(data$pvalue,method = "BH")
head(data)
write.table(data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Fishers/ResultsAdj.txt",sep="\t",col.names = NA,quote = FALSE)



###########################################################################################



com=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/Heatmap/All/Comm.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Comm = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                              grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(Comm=c("c1"="#4900ff","c2"="#00babf","c3"="#00c9b0","c4"="#00b8ff","c5"="#aa99ff","c6"="#bbff99","c7"="#99ffcc","c8"="#ff99ff")))

col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Comm = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Comm=c("c1"="#4900ff","c2"="#00babf","c3"="#00c9b0","c4"="#00b8ff","c5"="#aa99ff","c6"="#bbff99","c7"="#99ffcc","c8"="#ff99ff")))

corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/Heatmap/All/All_R_casted.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),col=col_fun2,cluster_rows=FALSE,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",left_annotation  = ha,top_annotation = ha3,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(7, "cm"),width  = unit(7, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/TEST/Heatmap/All/All_R_casted.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()


#########################


com=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Signi/Comm.txt",row.names = 1)
head(com)
ha = rowAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                   annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                   grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                   col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))

ha3 = HeatmapAnnotation(df = com,show_annotation_name = FALSE,simple_anno_size = unit(0.2, "cm"),
                        annotation_legend_param = list(Community = list(direction = "horizontal",grid_width = unit(0.5, "cm"),
                                                                        grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 0), labels_gp = gpar(fontsize = 10))),
                        col = list(Community=c("c1"="#FF0000","c2"="#ffa500","c3"="#4b0082","c4"="#0000ff","c5"="#ffff00","c6"="#ee82ee","c7"="#008000")))


col_fun2 = colorRamp2(c(1, 0.5, 0, -0.5, -1), c("#ff0000","#ff4c4c" ,"white","#6666ff", "#0000ff"))


corr=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Signi/Top_R_casted.txt",row.names = 1)
H2=Heatmap(as.matrix((corr)),cluster_rows=FALSE,col=col_fun2,cluster_columns = FALSE,show_column_names = FALSE,row_title_gp = gpar(fontsize=0),
           row_dend_width = unit(1.5, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "white",top_annotation = ha3,left_annotation = ha,
           heatmap_legend_param =list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),title_gp = gpar(fontsize = 10),
                                      labels_gp = gpar(fontsize = 10)),
           name = "Correlation Coeficient",show_row_names = FALSE,row_names_gp=gpar(fontsize = 15),height  = unit(7, "cm"),width  = unit(7, "cm"))



pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/Heatmap/Signi/Top_R_casted.pdf")
draw(H2,heatmap_legend_side = "right", annotation_legend_side = "bottom",merge_legend = TRUE)
dev.off()

#############
library(reshape2)
d=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Add_fig/Cartoon/tmp1.txt")
M=melt(d,id.vars =rownames(d))
head(M)

############################ RNASeq community pathway bubble plot

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/Bubble_Ip.txt")
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/C1.pdf")
ggplot(data) +
  geom_point(data=subset(data,Comm==1),aes(x=1,y=factor(Term,levels = unique(Term)),size=Ratio,color=log))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#CC9999",high="#800000")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(5,6,5,3, "cm"),
        legend.direction = "vertical",
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol = 1),
         color=guide_legend(title = "-log10(padj)",ncol = 1,override.aes = list(size = 3)))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/C2.pdf")
ggplot(data) +
  geom_point(data=subset(data,Comm==2),aes(x=1,y=factor(Term,levels = unique(Term)),size=Ratio,color=log))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#ffdb99",high="#ffa500")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(4,5,4,1.5, "cm"),
        legend.title = element_text(size=10),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/C3.pdf")
ggplot(data) +
  geom_point(data=subset(data,Comm==3),aes(x=1,y=factor(Term,levels = unique(Term)),size=Ratio,color=log))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#7586a0",high="#1a3662")+coord_fixed(ratio = 2)+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(4,4,8,1, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/C5.pdf")
ggplot(data) +
  geom_point(data=subset(data,Comm==5),aes(x=1,y=factor(Term,levels = unique(Term)),size=Ratio,color=log))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#ffea91",high="#ccaa1c")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=10,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(0,5.5,0.1,1, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))

dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/C6.pdf")
ggplot(data) +
  geom_point(data=subset(data,Comm==6),aes(x=1,y=factor(Term,levels = unique(Term)),size=Ratio,color=log))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#84b987",high="#0a740f")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        legend.direction = "vertical",
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(5,3,5,1, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/AT_vs_BT/Network/results/Species/NEW_GSET/C7.pdf")
ggplot(data) +
  geom_point(data=subset(data,Comm==7),aes(x=1,y=factor(Term,levels = unique(Term)),size=Ratio,color=log))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#b3cffb",high="#4287f5")+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),
        legend.direction = "vertical",
        axis.text.y = element_text(size=12,color="black"),plot.title = element_text(hjust = 0.5),
        axis.ticks = element_blank(),legend.text = element_text(colour = "black"),plot.margin = margin(3,5,2,1, "cm"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank())+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Gene ratio(%)",ncol=1),
         color=guide_legend(title = "-log10(padj)",ncol=1,override.aes = list(size = 5)))
dev.off()

rm(list = ls())
###########


library(NormalyzerDE)
?normalyzer
normalyzer(jobName="CCHF",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/raw_meta.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/Raw.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM"
)

normalyzer(jobName="CCHF",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/pd_meta.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_Norm.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM"
)

library(limma)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Quantile-normalized.txt",row.names = 1,check.names = FALSE)
head(ip)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/raw_meta.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(H7_Untreated_vs_H7_Treated = (groupTreated_Huh7 - groupUntreated_Huh7),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Results_huh7.txt",sep="\t",col.names = NA,quote = FALSE)

cont.matrix <-makeContrasts(sw_Untreated_vs_sw_Treated = (groupTreated_SW13 - groupUntreated_SW13),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Results_SW13.txt",sep="\t",col.names = NA,quote = FALSE)






ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Quantile-normalized.txt",row.names = 1,check.names = FALSE)
head(ip)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/pd_meta.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(H7_Untreated_vs_H7_Treated = (groupTreated_Huh7 - groupUntreated_Huh7),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)

write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Results_huh7.txt",sep="\t",col.names = NA,quote = FALSE)

cont.matrix <-makeContrasts(sw_Untreated_vs_sw_Treated = (groupTreated_SW13 - groupUntreated_SW13),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Results_SW13.txt",sep="\t",col.names = NA,quote = FALSE)


library(piano)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Results_huh7.txt",header = TRUE,row.names = 1)
p <- data[4]
head(p)
lfc<-data[1]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Piano_Results_huh7.txt",sep="\t",col.names = NA)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Results_SW13.txt",header = TRUE,row.names = 1)
p <- data[4]
head(p)
lfc<-data[1]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/RAW_NORM/Piano_Results_SW13.txt",sep="\t",col.names = NA)


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Results_huh7.txt",header = TRUE,row.names = 1)
p <- data[4]
head(p)
lfc<-data[1]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Piano_Results_huh7.txt",sep="\t",col.names = NA)

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Results_SW13.txt",header = TRUE,row.names = 1)
p <- data[4]
head(p)
lfc<-data[1]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/SW13/PD_NORM/Piano_Results_SW13.txt",sep="\t",col.names = NA)

########################## SW13 Datsset2
packageVersion("NormalyzerDE")
library(NormalyzerDE)
?normalyzer
normalyzer(jobName="CCHF",designPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/design.txt",
           dataPath = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/PD_Norm.txt",
           outputDir = "/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED"
)


library(limma)
ip=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Quantile.txt",row.names = 1,check.names = FALSE)
head(ip)
des=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/design.txt",row.names = 1)
head(des)
design <- model.matrix( ~0 + group, data = des)
fit <- lmFit(ip,design)
cont.matrix <-makeContrasts(H7_Untreated_vs_H7_Treated = (groupTreated_Huh7 - groupUntreated_Huh7),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Result_Huh7.txt",sep="\t",col.names = NA,quote = FALSE)

cont.matrix <-makeContrasts(sw_Untreated_vs_sw_Treated = (groupTreated_SW13 - groupUntreated_SW13),
                            levels = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
limma.res <- topTable(fit.cont,n="Inf")
head(limma.res)
write.table(limma.res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Results_SW13.txt",sep="\t",col.names = NA,quote = FALSE)

library(piano)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Huh7.txt",header = TRUE,row.names = 1)
p <- data[4]
head(p)
lfc<-data[1]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Huh7_Piano.txt",sep="\t",col.names = NA)


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Sw13.txt",header = TRUE,row.names = 1)
p <- data[4]
head(p)
lfc<-data[1]
head(lfc)
geneSets <- loadGSC(file = "/home/anoop/Desktop/Transcriptome/New_RNAseq/Project2_ART/With New HC/KEGG_Selected/KEGG.gmt")
gsares <- runGSA(geneLevelStats=p,directions=lfc,adjMethod = "BH",
                 gsc = geneSets,
                 nPerm = 500)
res=GSAsummaryTable(gsares)
write.table(res,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/SW13_Piano.txt",sep="\t",col.names = NA)

library(ggplot2)
tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/H7_All/Pwy_Direction_Cnt.txt")

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/H7_All/Huh7_Bar.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=Reg))+ylab("% Proteins detected")+scale_y_continuous(breaks = seq(0, 100, by = 50))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=9),plot.margin = margin(1,8,3.5,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y =element_blank(),axis.text.x = element_text(size=9,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=9),legend.title = element_text(size=9),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/H7_All/Signi.txt")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/H7_All/Huh7_Bubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=lpval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#b9c7fd",high="#415cc8")+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),legend.position = c(-2, -0.15),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=9,colour = "black"),
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(1,7.1,4,3.4, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()

############

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Sw_All/Pwy_Direction_Cnt.txt")
head(tmp)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Sw_All/SW13_Bar.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=Reg))+ylab("% Proteins detected")+scale_y_continuous(breaks = seq(0, 100, by = 50))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=9),plot.margin = margin(0.6,8,2.5,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y = element_blank(),axis.text.x = element_text(size=9,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=9),legend.title = element_text(size=9),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Sw_All/Signi.txt")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Sw_All/SW13_Bubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=lpval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#ffdb99",high="#cc8400")+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=8,color="black"),
        axis.ticks = element_blank(),legend.position = c(-2, -0.11),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(0.5,6.85,3,3.4, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()

############

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/48H/Pwy_Direction_Cnt.txt")
head(tmp)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/48H/48H_Bar.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=Reg))+ylab("% Proteins detected")+scale_y_continuous(breaks = seq(0, 100, by = 50))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=9),plot.margin = margin(0.25,8,1,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y = element_blank(),axis.text.x = element_text(size=9,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=9),legend.title = element_text(size=9),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/48H/Signi.txt")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/48H/48H_Bubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=lpval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#9fb8c9",high="#2c4f67")+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=8,color="black"),
        axis.ticks = element_blank(),legend.position = "left",legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(0.25,4.4,1,2, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 4)))
dev.off()

#####################

tmp=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/24H/Pwy_Direction_Cnt.txt")
head(tmp)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/24H/24H_Bar.pdf")
ggplot(data=tmp, aes(x=factor(Term,levels = unique(Term)), y=Perc, fill=Reg))+ylab("% Proteins detected")+scale_y_continuous(breaks = seq(0, 100, by = 50))+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+coord_flip()+scale_fill_manual(values = c(Down="#006600",Up="#ce0000"))+theme_bw()+
  theme(axis.title.y = element_blank(),axis.title.x = element_text(size=9),plot.margin = margin(0.25,8,1,7, "cm"),panel.grid.major = element_blank(),
        axis.text.y = element_blank(),axis.text.x = element_text(size=9,color="black"),legend.position = c(2, 0.5),legend.key.size = unit(0.4, "cm"),panel.border = element_blank(),
        legend.text = element_text(size=9),legend.title = element_text(size=9),axis.ticks.y = element_blank())+guides(fill=guide_legend(title="Regulation"))

dev.off() 


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/24H/Signi.txt")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/24H/24H_Bubble.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=lpval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#dfdfcb",high="#797965")+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=8,color="black"),
        axis.ticks = element_blank(),legend.position = "left",legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=8,colour = "black"),
        legend.title = element_text(size=8),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(0.25,4.4,1,2, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 4)))
dev.off()

#####################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/SigniProts.txt",header = T,row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/SW13_HeatMap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/SigniProtsZ.txt",sep="\t",quote = FALSE,col.names = NA)

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/SigniProtsZ.txt",header = T, row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/meta.txt",row.names = 1,header = T)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

ha2 = rowAnnotation(foo = anno_mark(at = c(32,33,34,37,38,39,42,43,48,52,81,114,123,125,126,127,128,129,130,273,299,300,302,303,304,305,306,
                                          307,308,309,310,311,312,313,314,320,354,356,357,456,700,701,727,728,731,732,733,734,735,759,802,895,898,1043,1044,1045,1105,1237,1283,1295,
                                          1331,1390,1418,1419,1424,1425,1476,1477,1621,1659,1660,1661,1662,1858,1868,1870,1871,1911,1931,1946,2000,2001,2002,2003,2089,2090,2091,2092,
                                          2093,2094,2095,2096,2097,2098,2103,2104,2105,2106,2107,2108,2109,2110,2111,2191,2239,2266,2320,2323,2324,2355,2356,2381,2382,2383,2384,2386,
                                          2388,2404,2430,2455,2456,2457,2469,2470,2524,2564,2982,3017,3018,3163,3277,3423,3681,3682,3683,3684,3685,3686),
                                   labels_gp = gpar(fontsize=7),lines_gp = gpar(col=NA),link_width=unit(0, "mm"),labels = c("ACAA1","ACACA","ACACB","ACADM","ACADSB","ACAT1","ACLY","ACO1","ACSF3","ACSS1","ADPGK","AKR1A1","ALDH1B1","ALDH3A1",
                                                                                                                            "ALDH3A2","ALDH7A1","ALDH9A1","ALDOA","ALDOC","ASH1L","ATP5F1D","ATP5MC1","ATP5ME","ATP5MF","ATP5MG","ATP5PB","ATP5PD",
                                                                                                                            "ATP5PF","ATP5PO","ATP6V0A1","ATP6V1A","ATP6V1C1","ATP6V1E1","ATP6V1G1","ATP6V1H","AUH","BCAT1","BCKDHA","BCKDHB","CAMKMT",
                                                                                                                            "COLGALT1","COLGALT2","COX11","COX15","COX5A","COX5B","COX6A1","COX6B1","COX7A2","CS","CYC1","DLAT","DLST","ENO2","ENO3",
                                                                                                                            "ENOSF1","EZH2","GALM","GLO1","GMPPB","GPI","HADHB","HIBADH","HIBCH","HK2","HKDC1","IDH2","IDH3A","KHK","KMT2A","KMT2B",
                                                                                                                            "KMT2C","KMT2D","MCCC2","MDH2","ME1","ME3","MINPP1","MMUT","MPI","MT-CO1","MT-CO3","MT-ND1","MT-ND2","NDUFA10","NDUFA11","NDUFA12","NDUFA2","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9",
                                                                                                                            "NDUFB10","NDUFB11","NDUFB5","NDUFB9","NDUFS1","NDUFS3","NDUFS4","NDUFS5","NDUFS8","NSD3","OGDH",
                                                                                                                            "OXCT1","PC","PCCA","PCCB","PDHA1","PDHB","PFKFB3","PFKFB4","PFKL","PFKP","PGAM4","PGK1","PHYKPL",
                                                                                                                            "PKM","PLOD1","PLOD2","PLOD3","PMM1","PMM2","PPA2","PRDM2","SDHB","SETD1B","SETD7","SORD","SUCLG2",
                                                                                                                            "TKFC","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1")))

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.001, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 10),
                       annotation_legend_param  = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),
                                                       title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 10)),
                       col = list(Group=c("b)Treated"="#6d5210","a)Untreated"="#c4941c")))
col_fun2 = colorRamp2(c(3, 1, 0, -1, -3), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
H1=Heatmap(as.matrix((Zscore)),col=col_fun2,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(1, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),
           top_annotation  =ha,row_split = 2,row_gap = unit(1, "mm"),column_gap = unit(1, "mm"),heatmap_legend_param = list(grid_width = unit(0.5, "cm"),grid_height = unit(0.5, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 10), labels_gp = gpar(fontsize = 10)),
           name = "Z-Score",show_row_names = FALSE,row_names_gp=gpar(fontsize = 5),height  = unit(10, "cm"),width  = unit(5, "cm"),
           column_split =sampleinfo$Group)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/SW13_HeatMap.pdf")
draw(H1, merge_legend = TRUE)
dev.off()


#############

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/Result_Huh7.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/Huh7_NEW_Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.5,nudge_x = 0.5)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC <= -1),aes(x=logFC,y=-log10(adj.P.Val )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC >= 1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC > 0 & logFC < 1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC < 0 & logFC > -1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, adj.P.Val >=.05),aes(x=logFC,y=-log10(adj.P.Val)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-5, 5), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()


########

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/Results_24h.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/Huh7_24H_Old_Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.5,nudge_x = 0.5)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC <= -1),aes(x=logFC,y=-log10(adj.P.Val )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC >= 1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC > 0 & logFC < 1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC < 0 & logFC > -1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, adj.P.Val >=.05),aes(x=logFC,y=-log10(adj.P.Val)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-5, 5), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()


data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/Results_SW13.txt",row.names = 1)
head(data)
library(ggrepel)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Figs/SW13_NEW_Volcano.pdf")
ggplot(data, aes(x=logFC, y=-log10(adj.P.Val))) + 
  geom_label_repel(aes(label = Label),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = 0.6,nudge_x = 0.5)+
  geom_label_repel(aes(label = Label2),segment.color = '#cccccc',fill=NA,color="black",size=3,label.size = NA,nudge_y = -0.6,nudge_x = 0.5)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC <= -1),aes(x=logFC,y=-log10(adj.P.Val )),pch=21,color="#003900",fill="#326632",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC >= 1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,fill="#b20000",color="#8e0000",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC > 0 & logFC < 1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, adj.P.Val <.05 & logFC < 0 & logFC > -1),aes(x=logFC,y=-log10(adj.P.Val)),pch=21,color="#a0a000",fill="#b2b200",size=2)+
  geom_point(data=subset(data, adj.P.Val >=.05),aes(x=logFC,y=-log10(adj.P.Val)),color="#90b4d2",size=1.2)+
  geom_vline(xintercept=1, linetype="dashed",size=0.35)+
  geom_vline(xintercept=-1, linetype="dashed",size=0.35)+scale_x_continuous(limits = c(-5, 5), breaks = seq(-6, 6, by = 2))+
  geom_hline(yintercept=1.3010299957, linetype="dashed",size=0.35)+
  theme(legend.title=element_text(size=8),legend.text=element_text(size=6),legend.key.size=unit(0.7,"line"),
        plot.title = element_text(hjust = 0.5,size =9),panel.grid.minor = element_blank(),
        axis.title=element_text(size=12),axis.text.y=element_text(size=12),axis.text.x=element_text(size=12),plot.margin = margin(2.5,2,2.5,2, "cm"))+
  labs(x="Log2 Fold Change",y="-log10 (Adj.P)")
dev.off()

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Common/Ip.txt")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/Common/Ip.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=lpval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#b9c7fd",high="#415cc8")+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),legend.position = c(-2, -0.15),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=9,colour = "black"),
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(1,7.1,4,3.4, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()



data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/huh7_sw13_comn/ip.txt")
head(data)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/Sw_13_NewData/NORMED/huh7_sw13_comn/ip.pdf")
ggplot(data, aes(y=Term)) + 
  geom_point(data=data,aes(x=1,y=factor(Term,levels = unique(Term)),size=Overlap,color=lpval))+scale_x_discrete(limits=c("1"))+
  scale_color_gradient(low="#fdccc6",high="#c8665b")+
  scale_y_discrete(position = "left")+theme_bw()+
  theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.text.y = element_text(size=9,color="black"),
        axis.ticks = element_blank(),legend.position = c(-2, -0.15),legend.box="vertical",
        legend.key.size = unit(0.7,"line"),legend.text = element_text(size=9,colour = "black"),
        legend.title = element_text(size=9),panel.border = element_blank(),panel.grid.major = element_blank(),plot.margin = margin(3,8,5,3.4, "cm"))+
  guides(size=guide_legend(override.aes=list(colour="grey"),title = "Overlap",nrow = 1),color=guide_legend(title = "-log10(Padj)",nrow = 1,override.aes = list(size = 5)))
dev.off()


######################### NEW Heatmap


library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN_LCPM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN_LCPM_Z.txt",sep="\t",quote = FALSE,col.names = NA)


Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN_LCPM_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/MetaData.txt",row.names = 1)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN_LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(SeverityScore = c("S1"="#ffa500","S2"="#cc0000"),Group=c("At_Infection"="#590059","12M_After"="#cc99cc")))
#col_fun2 = colorRamp2(c(2, 1, 0, -1, -2), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#5d5dff","#4c4cff","#0000cc"))

ha2 = rowAnnotation(foo = anno_mark(at = c(104,120,51,83,82,45,94,90,47,55,124,89,92,91,46,93,102,11,79,118,20,111,112,72,58,110,30,37,28,113,60,59,76,7,123,101,1,
                                           116,122,10,52,61,24,3,97,8,57,67,126,71,34,121,4,69,117,16,105,106,65,62,14,84,107,108,70,74,68,63,66,64,78),
                                    labels_gp = gpar(fontsize=15),lines_gp = gpar(col="black"),link_height = unit(7, "mm"),link_width=unit(12, "mm"),labels = c("NCAM1","TRIM2","RPS27A","CD44","CAMK2G","PLCG1","ICAM1","HLA-DRA","POM121C","UBA52","TRIM3","HLA-DPA1","HLA-E","HLA-DRB1","POM121","HLA-F","JAK1","EIF4G2","B2M","SUMO1","KPNB1","PRKCD","PTPN1","PSMB8","UBE2N","PML","NUP205","NUP54","NUP160","PTPN11","ADAR","ABCE1",
                                                                                                                                                                "STAT2","EIF4A3","TRIM26","IRF9","AAAS","SP100","TRIM25","EIF4G1","SEC13","BST2","NUP107","DDX58","IRF2","EIF4E","UBE2L6","IFIT5",
                                                                                                                                                                "TRIM5","MX2","NUP37","TRIM21","EIF2AK2","ISG20","STAT1","KPNA2","OAS1","OAS2","IFIT2","IFI35","HERC5","GBP1","OAS3","OASL","MX1","RSAD2","ISG15","IFI6","IFIT3","IFIT1","USP18")))

IFN_Meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN_Meta.txt",row.names = 1)
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),row_title_gp =gpar(fontsize = 20),
           top_annotation  =ha,row_split = IFN_Meta$Pwy,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                            title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",10),rep("S2",14)))

col_fun_lfc = colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#006600","#198c19","#4ca64c","white","#ca8383","#c26e6e" ,"#b04545"))
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(1, "cm"),right_annotation = ha2,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 12),height  = unit(37, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 


tt=H1 + H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/Heatmap.pdf",height = 22,width =20)
draw(tt, merge_legend = TRUE)
dev.off()

##########

Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/Signi/IFN_LCPM_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/MetaData.txt",row.names = 1)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/Signi/IFN_LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(SeverityScore = c("S1"="#ffa500","S2"="#cc0000"),Group=c("At_Infection"="#590059","12M_After"="#cc99cc")))
#col_fun2 = colorRamp2(c(2, 1, 0, -1, -2), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#5d5dff","#4c4cff","#0000cc"))

IFN_Meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/Signi/IFN_Meta.txt",row.names = 1)
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),row_title_gp =gpar(fontsize = 20),
           top_annotation  =ha,row_split = IFN_Meta$Pwy,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(35, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",10),rep("S2",14)))

col_fun_lfc = colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#006600","#198c19","#4ca64c","white","#ca8383","#c26e6e" ,"#b04545"))
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(1, "cm"),
           show_row_names = TRUE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 12),height  = unit(37, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 


tt=H1 + H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/Signi/Heatmap.pdf",height = 22,width =20)
draw(tt, merge_legend = TRUE)
dev.off()

##################### Merged proteomics PCA
BiocManager::install("PCAtools")
n

library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Merged.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(na.omit(count), metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/PCA.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Untreated_24h","Treated_24h","Untreated_48h","Treated_48h","Untreated_SW13","Treated_SW13"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=4,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(Untreated_24h="#cc8400",Treated_24h="#b27300",Untreated_48h="#3f75a2",Treated_48h="#23415a",Untreated_SW13="#f6546a",Treated_SW13="#93323f"))+
  scale_fill_manual(values=c(Untreated_24h="#ffa500",Treated_24h="#cc8400",Untreated_48h="#4682b4",Treated_48h="#2a4e6c",Untreated_SW13="#f77687",Treated_SW13="#ac3a4a"))+
  labs(x="PC1, 98.35%",y="PC2, 0.61%")+theme(axis.title = element_text(size=11),axis.text=element_text(size = 11),legend.position ="bottom",plot.margin = margin(4,4.5,4,4.5, "cm"),
                                             legend.title=element_blank(),legend.text=element_text(size=11),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 4))) 
dev.off()



library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Merged.txt",header=TRUE,row.names = 1)
Art.umap = umap(t(na.omit(data)))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/UMAP.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/UMAP.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Untreated_24h","Treated_24h","Untreated_48h","Treated_48h","Untreated_SW13","Treated_SW13"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/UMAP.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(Untreated_24h="#cc8400",Treated_24h="#b27300",Untreated_48h="#3f75a2",Treated_48h="#23415a",Untreated_SW13="#f6546a",Treated_SW13="#93323f"))+
  scale_fill_manual(values=c(Untreated_24h="#ffa500",Treated_24h="#cc8400",Untreated_48h="#4682b4",Treated_48h="#2a4e6c",Untreated_SW13="#f77687",Treated_SW13="#ac3a4a"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=11),axis.text=element_text(size = 11),legend.position ="bottom",plot.margin = margin(4,4.5,4,4.5, "cm"),
                                             legend.title=element_blank(),legend.text=element_text(size=11),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 4))) 
dev.off()


############
n

library(impute)
?impute.knn
BiocManager::install("sva")
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Merged.txt",row.names = 1,check.names = FALSE)
X=as.matrix(data)
KNN=impute.knn(data=X ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)
write.table(KNN$data,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Merged_Imputed",sep="\t",col.names = NA,quote = FALSE)

library(sva)
library(limma)
prot=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Merged.txt",row.names = 1,check.names = FALSE)
mat <- as.matrix(prot)
des=read.table("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/metadata.txt",sep="\t",header=TRUE)
designCombat = model.matrix(~ des$Group)
res1=removeBatchEffect(mat,batch = des$batch,design=designCombat)
head(res1)
rnaseqCombat = ComBat(mat, batch = des$batch, mod = designCombat, par.prior = TRUE, prior.plots = TRUE)
write.table(res1,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Batch_Corrected_NoImp.txt", sep="\t",quote=FALSE,col.names = NA)



library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Batch_Corrected_NoImp.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/metadata.txt",row.names = 1,check.names = FALSE)
p <- pca(na.omit(count), metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/PCA_batch_noImp.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/PCA_batch_noImp.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Untreated_24h","Treated_24h","Untreated_48h","Treated_48h","Untreated_SW13","Treated_SW13"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/PCA_batch_noImp.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=4,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(Untreated_24h="#cc8400",Treated_24h="#b27300",Untreated_48h="#3f75a2",Treated_48h="#23415a",Untreated_SW13="#f6546a",Treated_SW13="#93323f"))+
  scale_fill_manual(values=c(Untreated_24h="#ffa500",Treated_24h="#cc8400",Untreated_48h="#4682b4",Treated_48h="#2a4e6c",Untreated_SW13="#f77687",Treated_SW13="#ac3a4a"))+
  labs(x="PC1, 98.35%",y="PC2, 0.61%")+theme(axis.title = element_text(size=11),axis.text=element_text(size = 11),legend.position ="bottom",plot.margin = margin(4,4.5,4,4.5, "cm"),
                                             legend.title=element_blank(),legend.text=element_text(size=11),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 4))) 
dev.off()



library(umap)
data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/Merged_Imputed",header=TRUE,row.names = 1)
Art.umap = umap(t(na.omit(data)))
head(Art.umap$layout)
write.table(Art.umap$layout,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/UMAP_2.txt",sep="\t",col.names = NA,quote = FALSE)

dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/UMAP_2.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Untreated_24h","Treated_24h","Untreated_48h","Treated_48h","Untreated_SW13","Treated_SW13"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/UMAP_2.pdf")
ggplot(dat, aes(x=V1, y=V2,color=Group)) + geom_point(size=4,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(Untreated_24h="#cc8400",Treated_24h="#b27300",Untreated_48h="#3f75a2",Treated_48h="#23415a",Untreated_SW13="#f6546a",Treated_SW13="#93323f"))+
  scale_fill_manual(values=c(Untreated_24h="#ffa500",Treated_24h="#cc8400",Untreated_48h="#4682b4",Treated_48h="#2a4e6c",Untreated_SW13="#f77687",Treated_SW13="#ac3a4a"))+
  labs(x="UMAP1",y="UMAP2")+theme(axis.title = element_text(size=11),axis.text=element_text(size = 11),legend.position ="bottom",plot.margin = margin(4,4.5,4,4.5, "cm"),
                                  legend.title=element_blank(),legend.text=element_text(size=11),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 4))) 
dev.off()

###############

data2=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/NewHeatmP/Input.txt",header = TRUE,check.names = FALSE,row.names = 1)
head(data2)
col_fun3 = colorRamp2(c(0,0.25,0.5,0.75,1,1.5,3), c("white","#D69999","#CC7F7F","#B74C4C","#A31919","#990000","#890000"))
sf=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/NewHeatmP/H1_meta.txt",header = TRUE,row.names = 1)

nrow(sf)
ncol(data2)
ha = HeatmapAnnotation(df = sf, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),simple_anno_size = unit(0.2, "cm"),
                       annotation_legend_param  = list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),
                                                       title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 8)),
                       col = list(column=c("a) p adj (dist.dir.dn)"="#004C00",
                                           "b) p adj (mix.dir.dn)"="#008000","c) p adj (non-dir.)"="#808080",
                                           "d) p adj (mix.dir.up)"="#ff6666","e) p adj (dist.dir.up)"="#e50000"),Group=c("a) 24hpi"="#586E7C","b) 48hpi"="#6D5C78","c) Time"="#4682B4","d) SW13"="#ffc04c")))


H1=Heatmap(as.matrix((-log10(data2))),top_annotation = ha,col=col_fun3,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_title_gp = gpar(fontsize=0),row_split = 2,row_gap = unit(0, "cm"),
           row_dend_width = unit(1, "cm"),column_title_gp =gpar(fontsize = 0),na_col = "grey",column_split = sf$Group,column_gap = unit(0.05, "cm"),
           heatmap_legend_param =list(grid_width = unit(0.3, "cm"),grid_height = unit(0.3, "cm"),title_gp = gpar(fontsize = 8),
                                      labels_gp = gpar(fontsize = 8)),
           name = "-log10(Padj)",show_row_names = TRUE,row_names_gp=gpar(fontsize = 6),height  = unit(15, "cm"),width  = unit(6, "cm"))


pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/New_Piano/NewHeatmP/H1_meta.pdf")
draw(H1,heatmap_legend_side = "right", annotation_legend_side = "right",merge_legend = TRUE,ht_gap=unit(0.1, "cm"))
dev.off()


##################


library(PCAtools)
count=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/SW13/SW13.txt",row.names = 1,check.names = FALSE)
meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/SW13/meta.txt",row.names = 1,check.names = FALSE)
p <- pca(na.omit(count), metadata = meta, removeVar = 0.1)
head(p$variance)
write.table(p$rotated,file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/SW13/PCA.txt",sep="\t",col.names = NA,quote = FALSE)


dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/SW13/PCA.txt",row.names = 1,check.names = FALSE)
dat$Group <- factor(dat$Group, levels=c("Untreated","Treated"))
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/MERGED_PCA/SW13/PCA.pdf")
ggplot(dat, aes(x=PC1, y=PC2,color=Group)) + geom_point(size=4,aes(fill=Group),shape=21)+
  scale_color_manual(values=c(Untreated="#f6546a",Treated="#93323f"))+
  scale_fill_manual(values=c(Untreated="#f77687",Treated="#ac3a4a"))+
  labs(x="PC1",y="PC2")+theme(axis.title = element_text(size=11),axis.text=element_text(size = 11),legend.position ="bottom",plot.margin = margin(5,5.5,5,5.5, "cm"),
                                             legend.title=element_blank(),legend.text=element_text(size=11),legend.key.size = unit(0.5, "cm"))+
  guides(color=guide_legend(override.aes = list(size = 4))) 
dev.off()


###########

data=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/cnt.txt")
M=melt(data)

pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/PROTEOMICS/cnt.pdf")
ggplot(data=data, aes(x=factor(comp,levels = unique(comp)), y=count,fill=dge))+ylab("% Proteins detected")+
  geom_bar(stat="identity", position=position_dodge(),width = 0.7)+theme_bw()+ylab("# Protein regulated (p.adj<0.05)")+
  theme(plot.margin = margin(5,3.5,5,3.5, "cm"),axis.title.x = element_blank(),legend.title = element_blank())+
  scale_fill_manual(values = c(down="#006600",up="#ce0000"))
dev.off() 


####################

library(gplots)
Dat=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/IFN_LCPM.txt",row.names = 1,check.names = FALSE)
my_palette <- colorRampPalette(c("#003300","#004000","#008000","white","#ea1313","#8c0b0b","#700808"))(n=40)
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/Heatmap.pdf",width = 20,height = 7)
X=heatmap.2(as.matrix(Dat),tracecol=NA,col=my_palette,cexRow=0.8,cexCol = 1,keysize = 1,Rowv = FALSE, na.rm = TRUE,scale="row",
            key.title=NA,dendrogram = "none",Colv = FALSE,lhei=c(1.5,7),margins  =c(7,10),lwid=c(1,11),labRow = FALSE)
dev.off()

write.table(t(X$carpet),file="/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/IFN_LCPM_Z.txt",sep="\t",quote = FALSE,col.names = NA)


Zscore=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/IFN_LCPM_Z.txt", row.names = 1,check.names = FALSE)
sampleinfo <- read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/MetaData.txt",row.names = 1)
head(sampleinfo)
col_fun1 = colorRamp2(c(-2, -1,0, 1,2), c("#004c00","#008000","white","#e50000","#7f0000"))

library(ComplexHeatmap)
library(circlize)

LFC=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/IFN_LFC.txt",check.names = FALSE,row.names = 1)

ha = HeatmapAnnotation(df = sampleinfo, height =unit(0.1, "mm"),annotation_name_side = "left",annotation_name_gp = gpar(fontsize = 0),
                       annotation_legend_param  = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
                       col = list(SeverityScore = c("S1"="#ffa500","S2"="#cc0000"),Group=c("At_Infection"="#590059","12M_After"="#cc99cc")))
#col_fun2 = colorRamp2(c(2, 1, 0, -1, -2), c("#23415a","#3f75a2" ,"white","#d67834", "#7e3f12"))
col_fun1 = colorRamp2(c(3,2,1, 0,-1,-2,-3), c("#7f7f00","#b2b200" ,"#e5e500","white","#5d5dff","#4c4cff","#0000cc"))

ha2 = rowAnnotation(foo = anno_mark(at = c(69,94,146,72,70,84,147,86,101,149,135,71,85,103,105,74,16,78,95,119,138,76,62,145,161,
                                           90,104,116,18,148,100,4,1,167,182,144,36,8,7,67,124,179,73,49,133,117,65,121,3,178,102,77,75,32,170,53,
                                           66,54,87,11,26,63,169,47,88,115,107,92,109,160,22,127,39,168,30,141,81,96,122,155,80,166,48,120,140,151,152,12,125,171,142,162,154,57,108,131),
                                    labels_gp = gpar(fontsize=14),lines_gp = gpar(col="black"),link_height = unit(9, "mm"),
                                    link_width=unit(12, "mm"),labels = c("IFI27","USP18","OAS1","IFIT1","IFI35","ISG15","OAS2","MX1","CD44","OASL","IRF4","IFI6","ISG20","FCGR1A","GBP1",
                                                                         "IFIT3","HERC5","IFITM3","XAF1","HLA-DQB1","IRF7","IFITM1","UBE2L6","NCAM1","STAT1","RSAD2","FCGR1B","HLA-DPB1",
                                                                         "KPNA2","OAS3","CAMK2G","EIF2AK2","AAAS","TRIM21","TRIM8","MT2A","NUP37","EIF4E","EIF4A3","BST2","HLA-E","TRIM6",
                                                                         "IFIT2","POM121C","IRF2","HLA-DQA1","ABCE1","HLA-DRA","DDX58","TRIM5","CIITA","IFITM2","IFIT5","NUP205","TRIM26",
                                                                         "RPS27A","ADAR","SEC13","MX2","EIF4G1","NUP107","UBE2N","TRIM25","PLCG1","PSMB8","HLA-DPA1","GBP3","STAT2","GBP5",
                                                                         "SP100","KPNB1","ICAM1","NUP54","TRIM22","NUP160","JAK1","IFNAR2","B2M","HLA-DRB1","PTPN11","IFNAR1","TRIM2","POM121",
                                                                         "HLA-DQB2","IRF9","PML","PRKCD","EIF4G2","HLA-F","TRIM3","JAK2","SUMO1","PTPN1","UBA52","GBP4","IFNGR2")))


IFN_Meta=read.delim("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/IFN_Meta.txt",row.names = 1)
H1=Heatmap(as.matrix((Zscore)),col=col_fun1,cluster_rows=TRUE,cluster_columns = FALSE,show_column_names = FALSE,
           row_dend_width = unit(3, "cm"),column_title_gp =gpar(fontsize = 0),column_names_gp =gpar(fontsize = 5),row_title_gp =gpar(fontsize = 0),
           top_annotation  =ha,row_split = IFN_Meta$Reg,row_gap = unit(2, "mm"),column_gap = unit(2, "mm"),heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                                                                       title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           name = "Z-Score",show_row_names = TRUE,row_names_gp=gpar(fontsize = 5),height  = unit(40, "cm"),width  = unit(20, "cm"),
           column_split =c(rep("S1",10),rep("S2",14)))

col_fun_lfc = colorRamp2(c(-3,-2,-1, 0,1,2,3), c("#006600","#198c19","#4ca64c","white","#ca8383","#c26e6e" ,"#b04545"))
H2=Heatmap(as.matrix((LFC)),col=col_fun_lfc,cluster_rows=FALSE,cluster_columns = FALSE,name="LFC",width  = unit(1, "cm"),right_annotation = ha2,
           show_row_names = FALSE,column_names_side = "top",heatmap_legend_param = list(grid_width = unit(0.9, "cm"),grid_height = unit(0.9, "cm"),
                                                                                        title_gp = gpar(fontsize = 20), labels_gp = gpar(fontsize = 20)),
           row_names_gp =gpar(fontsize = 11),height  = unit(37, "cm"),column_names_gp =gpar(fontsize = 20),na_col = "#e6e6e6") 


tt=H1 + H2
pdf("/home/anoop/Desktop/Transcriptome/New_RNAseq/P13608/CCHF/New & Final Analysis/Figures/Heatmap/ISG/Changed_2_ISG/IFN/Heatmap.pdf",height = 22,width =20)
draw(tt, merge_legend = TRUE)
dev.off()




library("NormalyzerDE")
BiocManager::install("NormalyzerDE")
n
