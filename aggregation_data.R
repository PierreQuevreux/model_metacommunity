path="results/"
balance<-NULL
balance_minus<-NULL
balance_plus<-NULL
biomass_community_CV<-NULL
biomass_community_mean<-NULL
biomass_species_CV<-NULL
biomass_species_mean<-NULL
biomass_species_equilibrium<-NULL
biomass_species_toital_CV<-NULL
biomass_species_total_mean<-NULL
correlation<-NULL
covariance<-NULL
metacommunity<-NULL
press<-NULL
time_series<-NULL

file<-NULL

nSlice=10

for (i in 0:(nSlice-1)){
  file_name=paste(path,"balance_",i,".txt",sep="")
  if (file.access(file_name)==0){
    file<-read.table(file_name,sep=';',header=T)
    balance<-rbind(balance,file)
  }
  file_name=paste(path,"biomass_species_equilibrium_",i,".txt",sep="")
  if (file.access(file_name)==0){
    file<-read.table(file_name,sep=';',header=T)
    biomass_species_equilibrium<-rbind(biomass_species_equilibrium,file)
  }
  file_name=paste(path,"correlation_",i,".txt",sep="")
  if (file.access(file_name)==0){
    file<-read.table(file_name,sep=';',header=T)
    correlation<-rbind(correlation,file)
  }
  file_name=paste(path,"covariance_",i,".txt",sep="")
  if (file.access(file_name)==0){
    file<-read.table(file_name,sep=';',header=T)
    covariance<-rbind(covariance,file)
  }
}
write.table(balance,paste(path,"balance.txt",sep=""),sep=';',row.names=F)
write.table(biomass_species_equilibrium,paste(path,"biomass_species_equilibrium.txt",sep=""),sep=';',row.names=F)
write.table(correlation,paste(path,"correlation.txt",sep=""),sep=';',row.names=F)
write.table(covariance,paste(path,"covariance.txt",sep=""),sep=';',row.names=F)