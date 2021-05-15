library(clusterProfiler)
library(openxlsx)
library(org.Mm.eg.db)
library(pathview)
library(tidyverse)
trans<-read.xlsx("12LI.xlsx")  ###读入基因表达量的数据
feature<-read.csv2("jointpa_matched_features.csv",header=T,sep=",")##读入匹配特征数据
names(feature)<-c("pathway_name","matched_features")
########把x改成pathway 手动
ALL<-read.xlsx("all.xlsx")  ##读入差异代谢物数据
pathway<-read.csv2("MetaboAnalyst_result_pathway.csv",header=T,sep=",") ##读入差异表达分析的结果
pathway_include<-feature$matched_features
pathway_name<-feature$pathway_name
id<-bitr(trans$gene_name,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
if(length(pathway$pathway)!=length(pathway_name))
{
  stop(print("error,length of both terms not the same"))
}else{
  pathway$include=pathway_include
}
ALL<-data.frame(ALL %>% group_by(MS2.name) %>% filter(MS2.score==max(MS2.score))) 
##筛选出小于0.05的通路
pathway_expect<-pathway[which(pathway$Raw.p<0.05),]
list<-list()
for(i in 1:length(pathway_expect$pathway))
{ 
  list1=c(as.vector(unlist(strsplit(pathway_expect$include[i],split="; "))))
  list[[i]]=list(all=list1,
                 cpd=(unlist(lapply(list1[which(substring(list1,1,3)=="cpd")],function(x){x=substring(x,5)}))),
                 gene=unlist(lapply(list1[which(substring(list1,1,3)=="mmu")],function(x){x=substring(x,5)})))
}
####分离代谢物和基因
#####添加entrezid
trans$entrez=0
for(i in 1:length(id$ENTREZID))
{
  if(id$SYMBOL[i] %in% trans$gene_name)
  {
    trans$entrez[i]=id$ENTREZID[i]
  }else{
    trans$entrez[i]<-NA
  }
}
####数据归一化处理
####不想标准化的话注释掉234行的命令
select1<-function(x)
{
  new_data1<-x[,c("entrez","T1","T2","T3","T4","T5","T6","C1","C2","C3","C4","C5","C6","C","T")]  
  mean=mean((as.numeric(as.vector(unlist(new_data1[,-1])))))
  sd=sd(as.numeric(as.vector(unlist(new_data1[,-1]))))
  new_data1=data.frame(entrez=new_data1$entrez,lapply(new_data1[,-1],function(x){x=(x-mean)/sd}))
  new_data1$TC<-rowMeans(new_data1[,c("C","T")])
  return(new_data1)
}
select2<-function(y)
{
  new_data2<-y[,c("KEGG.ID","T1","T2","T3","T4","T5","T6","C1","C2","C3","C4","C5","C6")]
  names(new_data2)<-c("KEGG.ID","T1","T2","T3","T4","T5","T6","C1","C2","C3","C4","C5","C6")
  mean=mean((as.numeric(as.vector(unlist(new_data2[,-1])))))
  sd=sd(as.numeric(as.vector(unlist(new_data2[,-1]))))
  new_data2=data.frame(entrez=new_data2$KEGG.ID,lapply(new_data2[,-1],function(x){x=(x-mean)/sd}))
  new_data2$C<-rowMeans(new_data2[,c("C1","C2","C3","C4","C5","C6")])
  new_data2$T<-rowMeans(new_data2[,c("T1","T2","T3","T4","T5","T6")])
  new_data2$TC<-rowMeans(new_data2[,c("C","T")])
  return(new_data2)
}
####搜索生成数据框
if(!dir.exists("path_view"))
{dir.create("path_view")}
for(i in 1:length(pathway_expect$pathway))
{ 
  name<-paste0("path_view/",pathway_expect$pathway[i]," - Mus musculus (mouse)")
  if(!dir.exists(name))
  {dir.create(name)}
  data1<-trans[trans$entrez %in% list[[i]]$gene,]
  data2<-ALL[ALL$KEGG.ID %in% list[[i]]$cpd,]
  new_data1<-select1(data1)
  new_data2<-select2(data2)
  assign(paste0(pathway_expect$pathway[i],"-ge"),new_data1)
  assign(paste0(pathway_expect$pathway[i],"-cp"),new_data2)
  write.xlsx(new_data1,file=paste0(name,"/",pathway_expect$pathway[i]," - Mus musculus (mouse)-ge.xlsx"))
  write.xlsx(new_data2,file=paste0(name,"/",pathway_expect$pathway[i]," - Mus musculus (mouse)-cp.xlsx"))
}
pathwayinfo<-read.xlsx("pathwayinfo.xlsx",colNames = FALSE)
list_name=list.files("path_view")
matched_pathway<-pathwayinfo[pathwayinfo$X2 %in% list_name,]
pathway_id<-matched_pathway$X1
for(i in 1:length(matched_pathway$X1))
{
  pathway_id<-matched_pathway$X2[i]
  cpg_index<-paste0("path_view/",matched_pathway$X2[i],"/",matched_pathway$X2[i],"-cp.xlsx")
  ge_index<-paste0("path_view/",matched_pathway$X2[i],"/",matched_pathway$X2[i],"-ge.xlsx")
  ge_data<-read.xlsx(ge_index,rowNames = TRUE)
  ge_data<-data.frame(TC=ge_data$TC,row.names = rownames(ge_data))
  cpg_data<-read.xlsx(cpg_index,rowNames = TRUE)
  cpg_data<-data.frame(TC=cpg_data$TC,row.names=rownames(cpg_data))
  if(length(cpg_data$TC)==0){
    pv.out <- pathview(gene.data =ge_data,pathway.id =matched_pathway$X1[i], species = "mmu", out.suffix =paste0(matched_pathway$X1[i],"combine1"),
                       kegg.native=T,multi.state = T)
  }else
    pv.out <- pathview(gene.data =ge_data, cpd.data = cpg_data,
                       pathway.id =matched_pathway$X1[i], species = "mmu", out.suffix =paste0(matched_pathway$X1[i],"combine1"),
                       kegg.native = T,multi.state = T)
}
##########################z只取了均值，这样只有对照组和肿瘤组的双色对比。
#####富集分析
library(ggplot2)
enrich_pathway<-pathway_expect
enrich_pathway$enrichfactor<-enrich_pathway$Hits/enrich_pathway$Total
enrich_pathway$classify=0
for(i in 1:length(pathway_expect$pathway))
{
  if(is.null(list[[i]]$cpd))
  {
    enrich_pathway$classify[i]="none-cpd"
    
  }else
  {
    enrich_pathway$classify[i]="cpd"
  }
}
size=-log10(as.numeric(enrich_pathway$Raw.p))
factorcpd=factor(enrich_pathway$classify)
enrichplot<-ggplot(data=enrich_pathway,aes(enrichfactor,pathway))+geom_point(aes(color=factorcpd,size=size))+xlab('enrichfactor')+ylab('pathway')
if(!dir.exists("enrichplot"))
{
  dir.create("enrichplot")
}
ggsave(filename="enrichplot/enrichplot.png",enrichplot)
####相关性分析
##层次聚类分析（热图）
#gene_expression<-trans[,c("gene_name","T1","T2","T3","T4","T5","T6","C1","C2","C3","C4","C5","C6")]
##基因不多用上面的
##基因太多筛选<0.001的
library(psych)
library(pheatmap)
library(ggcorrplot)
library(reshape2)
gene_expression1<-trans %>% filter(pvalue<0.001) %>% dplyr::select("gene_name","T1","T2","T3","T4","T5","T6","C1","C2","C3","C4","C5","C6")
compound_expression1<-ALL[,c("KEGG.ID","T1","T2","T3","T4","T5","T6","C1","C2","C3","C4","C5","C6")]
compound_expression<-compound_expression1[which(!is.na(compound_expression1[,1])),]
gene_expression<-t(data.frame(gene_expression1,row.names = 1))
compound_expression$KEGG.ID<-make.names(compound_expression$KEGG.ID,unique=TRUE)
compound_expression<-t(data.frame(compound_expression,row.names = 1))
cor<-corr.test(gene_expression,compound_expression,method="pearson",adjust="none")
rmt<-cor$r
pmt<-cor$p
mycol<-colorRampPalette(c("blue","white","tomato"))(800)
if(!dir.exists("heatmap"))
   {dir.create("heatmap")}
pheatmap(rmt,scale = "none",cluster_row = T, cluster_col = T, border=NA,
         
         display_numbers = pmt, fontsize_number = 1, number_color = "white",
         
         cellwidth = 20, cellheight = 20,color=mycol,filename= "heatmap/heatmap.pdf")
##相关性分析
combine_expression<-cbind(gene_expression,compound_expression)
cor_2<-cor(combine_expression)
corplot_expression<-ggcorrplot(cor_2,tl.cex=1)
if(!dir.exists("cor.plot"))
{dir.create("cor.plot")}
ggsave(filename="cor.plot/corplot.png")
index=1
list_cor={}
index=1
list_cor={}
cor3<-data.frame(cor_2)
cor3<-data.frame(v=names(cor3),cor3)
cor_3<-melt(cor3,id.vars = c("v"))
cor_3<-cor_3[cor_3$value>0.5,]
####
cor_3_bar=cor_3 %>% filter(cor_3$v %in% compound_expression1$KEGG.ID)
cor_3_bar=cor_3_bar %>% filter(cor_3_bar$variable %in% gene_expression1$gene_name)
cor_3_bar<-cor_3_bar %>% arrange(desc(value)) %>% group_by(v)
top_10=0
true_false=duplicated(cor_3_bar$v)
for(top_index in 1:length(cor_3_bar$v))
{
  if(true_false[top_index]=="FALSE")
  {top_10=top_10+1}else
  {next}
  if(top_10==10){break}
}
cor_3_bar_10<-cor_3_bar[1:top_index,]
cor_3_bar_10$genefc=0
cor_3_bar_10$cpdfc=0
cor_3_bar_10$label=0
for(i in 1:top_index)
{
  cor_3_bar_10$genefc[i]=trans[which(trans$gene_name==cor_3_bar_10$variable[i]),]$log2FoldChange
  cor_3_bar_10$cpdfc[i]=ALL[which(ALL$KEGG.ID==cor_3_bar_10$v[i]),]$LOG_FOLDCHANGE
  cor_3_bar_10$label[i]=paste0(cor_3_bar_10$v[i],"-",cor_3_bar_10$variable[i])
}
dat1<-c(rep(0,2*top_index))
for(i in 1:(top_index))
{
  dat1[2*i-1]=cor_3_bar_10$genefc[i]
  dat1[2*i]=cor_3_bar_10$cpdfc[i]
}
dat_cor_10<-data.frame(GROUP=c(rep(cor_3_bar_10$label,each=2)),data=dat1,sub=c(rep(c("gene","cpd"),top_index)))
barplot2<-ggplot(data=dat_cor_10,aes(GROUP,data))+geom_bar(aes(fill = sub), stat="identity", position="dodge", width=.5)+theme(axis.text.x = element_text(angle = 70, hjust = 1))+ylab("log2foldchange")+ylim(-1,3)
if(!dir.exists("cpd-trans-bar"))
{
  dir.create("cpd-trans-bar")
}
ggsave("cpd-trans-bar/barplot.png",barplot2)
write.xlsx(dat_cor_10,"cpd-trans-bar/dat1.xlsx")
####fold change barplot
write.xlsx(cor_3,"cor.plot/cor_3.xlsx")
output_list<-cor_3$v[!duplicated(cor_3$v)]
gene_list<-output_list[output_list %in% gene_expression1$gene_name]
cpd_list<-output_list[output_list %in% compound_expression1$KEGG.ID]
write.xlsx(gene_list,"cor.plot/gene.xlsx")
write.xlsx(cpd_list,"cor.plot/cpd.xlsx")
######Venn图和barplot图
####ropls分析
#library(ropls)
###使用all
ALL_venn_file<-read.csv2("jointpa_matched_features_ALL.csv",sep=",",header=T)
list_ALL<-{}
for(i in 1:length(ALL_venn_file$X))
{ 
  list1=c(as.vector(unlist(strsplit(ALL_venn_file$matched_features[i],split="; "))))
  list_ALL[[i]]=list(all=list1,
                 cpd=(unlist(lapply(list1[which(substring(list1,1,3)=="cpd")],function(x){x=substring(x,5)}))),
                 gene=unlist(lapply(list1[which(substring(list1,1,3)=="mmu")],function(x){x=substring(x,5)})))
}
class<-c(rep(0,length(list_ALL)))
for(i in 1:length(list_ALL))
{
  if(is.null(list_ALL[[i]]$cpd))
  {
    class[i]<-"gene"
  }
  if(is.null(list_ALL[[i]]$gene))
  {
    class[i]<-"cpd"
  }
}
library(VennDiagram)
if(!dir.exists("venn_plot"))
{
  dir.create("venn_plot")
}
venn.plot <- venn.diagram(
  list(cpd = 1:(length(class[which(class %in% c("cpd","0"))])), gene=(length(class[which(class %in% c("cpd"))])+1):length(class)), 
  "venn_plot/Venn_3set_simple.jpeg",fill=c("red","blue")
)
##barplot
barplot_file<-read.csv2("jointpa_matched_features.csv",sep=",",header=T)
list_ALL1<-{}
for(i in 1:length(barplot_file$X))
{ 
  list1=c(as.vector(unlist(strsplit(barplot_file$matched_features[i],split="; "))))
  list_ALL1[[i]]=list(all=list1,
                     cpd=(unlist(lapply(list1[which(substring(list1,1,3)=="cpd")],function(x){x=substring(x,5)}))),
                     gene=unlist(lapply(list1[which(substring(list1,1,3)=="mmu")],function(x){x=substring(x,5)})))
}
class1<-c(rep(0,length(list_ALL1)))
for(i in 1:length(list_ALL1))
{
  if(is.null(list_ALL1[[i]]$cpd))
  {
    class1[i]<-"gene"
  }
  if(is.null(list_ALL1[[i]]$gene))
  {
    class1[i]<-"cpd"
  }
}
long=c(rep(0,length(class1[which(class1 %in% c("0"))])))
j=1
index_all=c(rep(0,length(class1[which(class1 %in% c("0"))])))
for(i in 1:length(list_ALL1))
{
  if(class1[i]=="0")
  {
    long[j]=length(list_ALL1[[i]]$all)
    index_all[j]=i
    j=j+1
  }
}
long_all<-setNames(long,index_all)

if(length(long_all)>10){
  long_10<-order(long_all,decreasing=TRUE)[1:10]
}else{
  long_10<-order(long_all,decreasing=TRUE)[1:length(long_10)]}
long_all<-long_all[long_10]

data<-c(rep(0,2*length(long_10)))
for(i in 1:length(long_10))
{
  data[2*i-1]=length(list_ALL1[[as.numeric(names(long_all[i]))]]$cpd)
  data[2*i]=length(list_ALL1[[as.numeric(names(long_all[i]))]]$gene)
}
dat<-data.frame(GROUP=rep(ALL_venn_file$X[as.numeric(names(long_all))],each=2),sub=c(rep(c("cpd","gene"),length(long_10))),data=data)
GROUP_SORT=reorder(dat$GROUP,-dat$data)
barplot<-ggplot(data=dat,aes(GROUP_SORT,data))+geom_bar(aes(fill = sub), stat="identity", position="dodge", width=.5)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+ylab("number of features")
if(!dir.exists("barplot"))
{
  dir.create("barplot")
}
ggsave("barplot/barplot.png",barplot)
ggsave("barplot/barplot.pdf",barplot)
write.xlsx(dat,"barplot/pathway_rank10.xlsx")



