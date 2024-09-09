# Identifying StD subtypes based on individual FCS deviations

# install NbClust
library("NbClust")

#determine N
data <- read.table('D:\\Data_Chen\\With_DIDA_all_HC\\res_norm\\res_norm\\z_base.txt',header = FALSE)
index_all = list("kl","ch","hartigan","cindex","db","silhouette","duda","pseudot2","beale","ratkowsky","ball","ptbiserial","gap","frey","mcclain","gamma","gplus","tau","dunn","sdindex","sdbw")
for (i in 1:21){
  res <- NbClust(data, distance = "euclidean", min.nc = 2, max.nc = 10, method = "kmeans",index=c(index_all[i]))
  write.csv(res$Best.partition,paste0('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\partition_',i,'.csv'))
  write.csv(res$Best.nc,paste0('D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\nc_',i,'.csv'))
}

# result from kmeans
KMClu2<-kmeans(x=data,centers=2,nstart=10)
save(KMClu2,file='D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\res_KMClus2.Rdata')
write.csv(KMClu2$cluster,'D:\\Data_Chen\\With_DIDA_all_HC\\subtype\\KMClus2.csv')