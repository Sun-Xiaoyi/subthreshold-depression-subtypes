
library(reshape2) 
library(ez)

for (i in 1:220){
  data <- read.csv(paste0('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\treatment\\brain\\data_forR',i,'.csv'))
  data$ID <- 1:47
  data_long <- melt(data=data,id.vars=c('ID','subtype'), #保留不变的变量
                measure.vars=c('time1','time2'), #想要转换的变量
                variable.name='time', #转换后的分类变量名
                value.name='Deviation') #转换后的数值变量名
  result <- ezANOVA(data=data_long,dv=Deviation,wid=ID,within=.(time),
                    between=subtype,detailed=T) #重复测量方差分析
  write.csv(data.frame(result),file=paste0('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\treatment\\brain\\result_rANOVA',i,'.csv'),row.names = FALSE)
}
