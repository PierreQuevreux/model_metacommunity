path="results/"

IDsimu=0 # number of the slice of simulation
n_slice=1

nReplicates=1
d_min=-5
d_max=5
n_dispersal=(d_max-d_min)*10

nSpecies=4
nCommunity=2

m=c(0.0065,0.065,0.65,6.5,65)
g=1
r=0
D=1
e=0.65
a=c(1/6.5,1/0.65,1/0.065)
FR=1
h=0

# models ---------
params_matrix<-read.table(paste(path,"parameters_matrix.txt",sep=""),sep=';',header=TRUE)
models<-read.table(paste("models.txt",sep=""),sep=';',header=TRUE)
models$model<-paste(models$perturbation,models$dispersal,sep="-")
# find the ID of perturbations matrix
params<-params_matrix[params_matrix$type=="perturbation",-1]
names(params)[1]="perturbation"
models<-merge(models,params,by=c("sigma_exo","sigma_demo","sigma_env","perturbation"))
models$perturbation<-NULL
names(models)[dim(models)[2]]="perturbation"
# find the ID of dispersal matrix
params<-params_matrix[params_matrix$type=="dispersal",-c(1,which(names(params_matrix)%in%c("sigma_exo","sigma_demo","sigma_env")))]
names(params)[1]="dispersal"
models<-merge(models,params,by=c("dispersal"))
models$dispersal<-NULL
names(models)[dim(models)[2]]="dispersal"
# find the ID of perturbations matrix
params<-params_matrix[params_matrix$type=="asymmetry",-c(1,which(names(params_matrix)%in%c("sigma_exo","sigma_demo","sigma_env")))]
names(params)[1]="asymmetry"
models<-merge(models,params,by=c("asymmetry"))
models$asymmetry<-NULL
names(models)[dim(models)[2]]="asymmetry"
rm(params)
models$model_ID<-seq(1:dim(models)[1])
models<-models[,which(names(models)%in%c("model","perturbation","dispersal","asymmetry","sigma_exo","sigma_demo","sigma_env","model_ID"))]
models<-models[,c(5:7,1:4,8)]

# parameter table ---------
simu_ID<-seq(1:nReplicates)
params<-expand.grid(simu_ID=simu_ID,
                    m=m,
                    g=g,
                    r=r,
                    D=D,
                    e=e,
                    a=a,
                    FR=FR,
                    h=h,
                    model_ID=models$model_ID)
params<-merge(params,models,by=c("model_ID"))
params$ma=params$m*params$a
params<-params[params$ma>0.05 & params$ma<=15,]
params$ea=params$e*params$a
params$model_ID<-NULL
params$simu_ID<-c(1:dim(params)[1])

# split the variables into nSlice sub tables
nSimu=dim(params)[1] # number of simulations
slice_seq=seq(1:n_slice)
slice_seq[1:n_slice]=nSimu %/% n_slice
slice_seq[n_slice]=slice_seq[n_slice]+nSimu-slice_seq[n_slice]*n_slice

######
# Save the parameters

IDsimu_0=0
for (i in 1:n_slice){
  
  # Simulation informations
  nParams=dim(params)[2] # number of parameters
  nThreads=10 # number of threads used
  params_data<-c(slice_seq[i], # number of simulations in the slice
                 nParams,
                 nSpecies,
                 nCommunity,
                 nThreads)
  
  # sample of params
  paramsbis<-params[(i-1)*slice_seq[1]+(1:slice_seq[i]),]
  
  IDsimu = IDsimu_0+i-1 # ID of the slice
  write.table(paramsbis,paste(path,"parameters_",IDsimu,".txt",sep=""),sep=";",row.names = FALSE)
  write.table(params_data,paste(path,"parameters_data_",IDsimu,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)
}

