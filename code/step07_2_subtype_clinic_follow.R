data <- read.csv("D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\treatment\\clinical\\data_forR.csv")  
data$ID <- 1:47
library(reshape2)  
data_long<-melt(data=data,id.vars=c('ID','subtype'), #保留不变的变量
                measure.vars=c('time1','time2'), #想要转换的变量
                variable.name='time', #转换后的分类变量名
                value.name='HDRS') #转换后的数值变量名
library(ez) #调用包“ez”
result <- ezANOVA(data=data_long,dv=HDRS,wid=ID,within=.(time),
                  between=subtype,detailed=T) #重复测量方差分析
write.csv(data.frame(result),file="D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\treatment\\clinical\\result_rANOVA.csv",row.names = FALSE)
