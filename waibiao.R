library(openxlsx)
down=read.xlsx("table1.xlsx") #######输入的数据
names(down)[1]="NONGDU"
downnongdu<-down$NONGDU[1:length(down$NONGDU)/2]
half<-(length(down$NONGDU)/2)
for(i in 2:length(names(down)))
{
  for(j in 1:half)
  {
    if(down[j,i]==0)
    {
      down[j+half,i]<-NA
    }
  }
}
new_data=down[(half+1):(half*2),]
t=data.frame(new_data,row.names = 1)
t1=t(t)
t2=as.data.frame(t1,row.name=F)
yuansu<-names(down)[2:length(names(down))]
t3=as.data.frame(cbind(yuansu,t2))
write.xlsx(t3,"t3.xlsx")
init=read.xlsx("table2.xlsx") ###输入数据
length_init=length(init$V1)
nullinit=init[,2]
nullinit=length(nullinit[which(nullinit=="/")])
for(i in 4:length(names(init)))
{
  if(i==4)
  {
    nongdu=c(rep(0,length_init))
    hanliang=c(rep(0,length_init-nullinit),rep("/",nullinit))
    data=data.frame(init[i],nongdu,hanliang)
  }else
  {
    nongdu=c(rep(0,length_init))
    hanliang=c(rep(0,length_init-nullinit),rep("/",nullinit))
    data=cbind(data,init[i],nongdu,hanliang)
  }
}
data=data.frame(init[1:3],data)
write.xlsx(data,"pre_t3.xlsx")
library(ggplot2)
library(broom)
gz=read.xlsx("t3.xlsx")
number<-length(gz$yuansu)
multi=c(rep(1,number)) ######手动设置
gz2=read.xlsx("t3.xlsx",rowNames = FALSE)
name1=names(gz2)[-1]
name1=strsplit(name1,split="_")
standard<-c()
for(i in 1:length(name1))
{
  standard=c(standard,name1[[i]][2])
}
standard<-as.numeric(standard)
#standard=as.numeric(names(gz)[-1])
QC_value="QC_50"
dilute_value=1000
number=length(gz2$yuansu)
coef1=c(rep(0,number))
coef2=c(rep(0,number))
rsquare=c(rep(0,number))
rsquare_ad=c(rep(0,number))
sz=c(rep(0,number))
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 3),
                        b = format(unname(coef(m)[2]), digits = 3),
                        r2 = format(summary(m)$r.squared, digits = 5)))
  as.character(as.expression(eq));
}
##表格显示函数
choose1=function(m){
  m=na.omit(m)
  length=length(m$x)
  fw=m$x[length]
  for(i in 1:number)
  {
    for(j in 1:length(standard))
    {
      if(fw==standard[j]*multi[i])
      {
        n=fw*1.1
      }
    }
  }
  return(n)
}
choose2=function(m)
{
  m=na.omit(m)
  length=length(m$y)
  fw2=m$y[length]
  return(fw2)
}
##定义一个作图范围
options(scipen = 10)
for(i in 1:number)
{
  y=gz[i,2:(length(standard)+1)]
  y=as.vector(unlist(y))
  x=standard*multi[i]
  datam=data.frame(x,y)
  datam=na.omit(datam)
  pre=lm(y~x,data=datam)
  limit=choose1(datam)
  ylimit=choose2(datam)
  filename=paste(gz[i,1],"-plot.png")
  plot=ggplot(data=datam,aes(x=x,y=y))+geom_point(size=2)+geom_smooth(method="lm")+ geom_text(size=6,x = limit/2, y = ylimit, label = lm_eqn(datam), parse = TRUE)+xlim(c(0,limit))+xlab("浓度")+ylab("峰面积")+ggtitle(gz[i,1])
  filename=paste(gz[i,1],"-plot.png")
  ggsave(filename,plot=plot)
  summary=summary(pre)
  rsquare[i]=summary$r.squared
  rsquare_ad[i]=summary$adj.r.squared
  pre=tidy(pre)
  coef2[i]=pre$estimate[1]
  coef1[i]=pre$estimate[2]
}
#生成图片
coef1=round(coef1,digits=4)
coef2=round(coef2,digits=4)
rsquare=round(rsquare,digits=3)
rsquare_ad=round(rsquare_ad,digits=3)
for(i in 1:number)
{
  sz[i]=paste0("y=",coef2[i],"+",coef1[i],"·x")
}
gz$nihe=sz
gz$slope=coef1
gz$intercept=coef2
gz$rsq=rsquare
gz$rsqad=rsquare_ad
write.xlsx(gz,"predict.xlsx")
ndpredict=read.xlsx("pre_t3.xlsx")
totallen=length(ndpredict$V1)
QC=ndpredict[which(ndpredict$V3==QC_value),]
QClength=length(QC$V1)
prelen=totallen-QClength
chengyang=as.numeric(ndpredict$V2[1:prelen])
for(i in 1:number)
{
  j=3*i+1
  if(i==1)
  {
    qqq=ndpredict[,j]
    data1=as.numeric(qqq)
    nongdu=(data1-coef2[i])/coef1[i]
    nongdu1=nongdu[1:prelen]
    hanliang=(nongdu1*0.6/chengyang)*dilute_value
    hanliang=c(hanliang,c(rep("/",QClength)))
    data=data.frame(data1,nongdu,hanliang)
  }
  else
  {
    qqq=ndpredict[,j]
    data1=as.numeric(qqq)
    nongdu=(data1-coef2[i])/coef1[i]
    nongdu1=nongdu[1:prelen]
    hanliang=(nongdu1*0.6/chengyang)*dilute_value
    hanliang=c(hanliang,c(rep("/",QClength)))
    data=cbind(data,data1,nongdu,hanliang)
  }
}
data=cbind(ndpredict[,1:3],data)
RSDdata=data[which(data$V3==QC_value),]
RSD=c(rep(0,number*3))
for(i in 1:number)
{
  k=3*i+2
  RSD[3*i+1]=round(sd(RSDdata[1:3,k])/mean(RSDdata[1:3,k])*100,digits=2)
}
#计算相对标准偏差
RSD[1]="RSD%"
name=c("峰面积","浓度(ng/mL)","样品含量(pg/g)")
name=c("序号","称样量（g）","上机编号",rep(name,number))
names(data)=name
data=rbind(data,RSD)
name2<-names(down)[-1]
name2<-c(0,0,0,rep(name2,each=3))
data<-rbind(name2,data)
write.xlsx(data,"predic_init.xlsx")

