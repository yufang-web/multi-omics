######内标一点法
####输入文件
library(openxlsx)
library(tidyverse)
library(broom)
library(ggplot2)
inner_table<-read.xlsx("16.xlsx",rowNames=FALSE) ###改名字
colnames(inner_table)[1]="sample"
####设置参数
std_dl<-500
in_dl<-50
inner_num=1
if(inner_num != 1){
num<-length(names(inner_table))-(inner_num+1)   ###改数字
inner<-inner_table[,c("D9_1")] ###改参数
inner_Table<-inner_table[,1:(num+1)]
###看情况调整
#####function
###判断这一列属于哪个一个内标
use_num<-c(rep(0,num))
###3个内标
analyze1<-function(x)
{
  x=as.vector(unlist(x))
  long=length(x)
  use<-c(rep(0,long))
  for(i in 1:length(x))
  {
    for(j in 1:length(inner))
    {
      if(x[i]<=inner[length(inner[,1]),1]){use[i]=1
      }else if((inner[length(inner[,1]),j]<x[i])&(x[i]<=inner[length(inner[,1]),j]))
      {
        use=j+1
      }else{use=length(inner)}
    }
  }
}
use_num<-analyze1(inner_Table[length(inner_Table$sample),-1])
###函数遇到别的情况也要改
dilution_result<-matrix(nrow=(length(inner_Table$sample)-1),ncol=(length(inner_Table)-1))
f_test<-c(rep(0,(length(inner_Table)-1)))
for(i in 1:(length(inner_Table)-1))
{for(j in 2:(length(inner_Table$sample)-1))
  {
  f=(std_dl/inner_Table[1,i+1])/(in_dl/(inner[1,use_num[i]]))
  dilution_result[j-1,i]<-(f*in_dl)/inner[j,use_num[i]]*inner_Table[j,i+1]

}
}
f_test<-(std_dl/inner_Table[1,-1])/(in_dl/(inner[1,use_num]))
dilu_test<-c(rep(0,length(inner_Table$sample)-1))
dilu_test<-as.vector(unlist(inner[1,use_num]))
dilution_result<-data.frame(dilution_result)
names(dilution_result)<-names(inner_Table)[-1]
dilution_result<-cbind(sample=inner_table$sample[2:(length(inner_table$sample))],dilution_result)
dilution_result[length(dilution_result$sample),-1]<-inner_Table[length(inner_Table$sample),-1]
biao_test<-as.vector(unlist(inner_Table[1,][-1]))
F_test<-rbind(f_test,biao_test,dilu_test)
F_test<-data.frame(t(F_test))
F_test<-cbind(sample=names(f_test),F_test)
F_test<-F_test[order(F_test$X3),]
result<-dilution_result$sample
null<-c(rep(0,length(dilution_result$sample)))
for(i in 1:(length(dilution_result)-1))
{
  result<-cbind(result,inner_Table[2:(length(inner_Table$sample)),i+1],dilution_result[,i+1],null)
}
result<-cbind(result,inner[-1,])
name=c("sample",rep(names(inner_Table)[-1],each=3),c(rep("inner",inner_num)))
for(i in 1:(length(names(inner_Table)[-1])-1))
{
  name[3*i]=""
  name[3*i+1]=""
}
RSD=c("RSD",rep(0,3*length(names(inner_Table)[-1])),c(rep(0,inner_num)))
result<-data.frame(result)
names(result)<-name
xuhao<-1:(length(result))
result<-rbind(result,RSD,xuhao)
write.xlsx(result,"result.xlsx")
names(F_test)<-c("物质名称","相对系数（f）","峰面积","内标峰面积")
write.xlsx(F_test,"f.xlsx")
}else{
  num<-length(names(inner_table))-(inner_num+1)   ###改数字
  inner<-inner_table[,c("D9_1")] ###改参数
  inner_Table<-inner_table[,1:(num+1)]
  ###看情况调整
  #####function
  ###判断这一列属于哪个一个内标
  use_num<-c(rep(0,num))
  use_num<-c(rep(1,num))
  ###函数遇到别的情况也要改
  dilution_result<-matrix(nrow=(length(inner_Table$sample)-1),ncol=(length(inner_Table)-1))
  f_test<-c(rep(0,(length(inner_Table)-1)))
  for(i in 1:(length(inner_Table)-1))
  {for(j in 2:(length(inner_Table$sample)-1))
  {
    f=(std_dl/inner_Table[1,i+1])/(in_dl/(inner[use_num[i]]))
    dilution_result[j-1,i]<-(f*in_dl)/inner[j]*inner_Table[j,i+1]
    
  }
  }
  f_test<-(std_dl/inner_Table[1,-1])/(in_dl/(inner[use_num]))
  dilu_test<-c(rep(0,length(inner_Table$sample)-1))
  dilu_test<-as.vector(unlist(inner[use_num]))
  dilution_result<-data.frame(dilution_result)
  names(dilution_result)<-names(inner_Table)[-1]
  dilution_result<-cbind(sample=inner_table$sample[2:(length(inner_table$sample))],dilution_result)
  dilution_result[length(dilution_result$sample),-1]<-inner_Table[length(inner_Table$sample),-1]
  biao_test<-as.vector(unlist(inner_Table[1,][-1]))
  F_test<-rbind(f_test,biao_test,dilu_test)
  F_test<-data.frame(t(F_test))
  F_test<-cbind(sample=names(f_test),F_test)
  F_test<-F_test[order(F_test$X3),]
  result<-dilution_result$sample
  null<-c(rep(0,length(dilution_result$sample)))
  for(i in 1:(length(dilution_result)-1))
  {
    result<-cbind(result,inner_Table[2:(length(inner_Table$sample)),i+1],dilution_result[,i+1],null)
  }
  result<-cbind(result,inner[-1])
  name=c("sample",rep(names(inner_Table)[-1],each=3),"inner")
  for(i in 1:(length(names(inner_Table)[-1])-1))
  {
    name[3*i]=""
    name[3*i+1]=""
  }
  RSD=c("RSD",rep(0,3*length(names(inner_Table)[-1]),0))
  result<-data.frame(result)
  names(result)<-name
  xuhao<-1:(length(result))
  result<-rbind(result,RSD,xuhao)
  write.xlsx(result,"result.xlsx")
  names(F_test)<-c("物质名称","相对系数（f）","峰面积","内标峰面积")
  write.xlsx(F_test,"f.xlsx")
}

