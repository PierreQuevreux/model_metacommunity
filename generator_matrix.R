# This program generate the dispersion matrix that only contain the factors d that are the specific dispersion relative to the
# selfe regulation coefficient

#library(blockmatrix)

path="results/"
nSpecies=4
nCommunity=2

# dispersal rates (log10 scale)
d_min=-5
d_max=5
n_dispersal = (d_max-d_min)*10

###### Functions generating the matrix ###### 

generate_dispersion_matrix <- function(nSpecies,nCommunity,d) {
  P<-diag(d)
  Pi<-P/(nCommunity-1)
  
  P_final=NULL
  for (i in 1:nCommunity){
    P_row=NULL
    for (j in 1:nCommunity){
      if (i==j){
        P_row<-cbind(P_row,-P)
      }
      else{
        P_row<-cbind(P_row,Pi)
      }
    }
    P_final<-rbind(P_final,P_row)
  }
  
  test=0
  for (i in 1:(nSpecies*nCommunity)){
    test=test+sum(P_final[i,])
  }
  if (test==0){
    print("Dispersion matrix correct")
  }
  else {
    print("Mass conservation not respected")
  }
  return(P_final)
}

generate_perturbation_matrix_block <- function(nSpecies,nCommunity,p,sigma) {
  P <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
  speciesID = 0
  for (i in 1:dim(p)[1]){
    speciesID = p[i,1]+(p[i,2]-1)*nSpecies
    P[speciesID,speciesID] = sigma
  }
  return(P)
}

generate_perturbation_matrix <- function(P_exo,P_demo,P_env) {
  Z <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
  # exogenous perturbation
  P <- cbind(P_exo,Z)
  P <- cbind(P,Z)
  # demographic perturbation
  P_row <- cbind(Z,P_demo)
  P_row <- cbind(P_row,Z)
  P <- rbind(P,P_row)
  # environmental perturbation
  P_row <- cbind(Z,Z)
  P_row <- cbind(P_row,P_env)
  P <- rbind(P,P_row)
  return(P)
}

###### Perturbation matrix ###### 

# demo - 1 - TL1 ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.001
p=rbind(c(1,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
P_env <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=1
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# demo - 1 - TL2 ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.001
p=rbind(c(2,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
P_env <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=2
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# demo - 1 - TL3 ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.001
p=rbind(c(3,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
P_env <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=3
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# demo - 1 - TL4 ####
### demographic perturbation
sigma=0.001
p=rbind(c(4,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=4
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# demo - 1 - H (all trophic levels) ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
P_env <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=0
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# demo - 2 - H (all trophic levels in the two patches) ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1),c(1,2),c(2,2),c(3,2),c(4,2)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
P_env <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=5
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# env - 1 - H ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
P_demo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### environmental perturbation
sigma=0.001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_env<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=6
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# env - 1 - H - demo - 2 - H (all species in the two patches) - sigma ratio = 0.01 ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1),c(1,2),c(2,2),c(3,2),c(4,2)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
sigma=0.00001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_env<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=7
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

# env - 1 - H - demo - 2 - H (all species in the two patches) - sigma ratio = 100 ####
### exogenous perturbation
P_exo <- matrix(0, nSpecies*nCommunity, nSpecies*nCommunity)
### demographic perturbation
sigma=0.00001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1),c(1,2),c(2,2),c(3,2),c(4,2)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_demo<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### environmental perturbation
sigma=0.001
p=rbind(c(1,1),c(2,1),c(3,1),c(4,1)) # coordinates of disturbed trophic levels (i=trophic level and j=community ID)
P_env<-generate_perturbation_matrix_block(nSpecies,nCommunity,p,sigma)
### final matrix
perturbation_matrix<-generate_perturbation_matrix(P_exo,P_demo,P_env)
IDmatrix=8
write.table(perturbation_matrix,paste(path,"matrix_perturbation_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

###### Dispersal matrix ###### 

# H (all trophic levels) ####
v=rep(1,nSpecies)
d=0*v
IDmatrix=0
for (i in seq(from=d_min, to=d_max, by=(d_max-d_min)/n_dispersal)){
  d = c(d,10^i*v)
}
d = matrix(d,
           nrow=n_dispersal+2,
           ncol=nSpecies,
           byrow = TRUE)
for (i in 1:(n_dispersal+2)){
  dispersion_matrix<-generate_dispersion_matrix(nSpecies,nCommunity,d[i,])
  IDd=IDmatrix*(n_dispersal+2) + i-1
  write.table(dispersion_matrix,paste(path,"matrix_dispersal_",IDd,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)
}

# TL1 ####
v=c(1,0,0,0)
d=0*v
IDmatrix=1
for (i in seq(from=d_min, to=d_max, by=(d_max-d_min)/n_dispersal)){
  d = c(d,10^i*v)
}
d = matrix(d,
           nrow=n_dispersal+2,
           ncol=nSpecies,
           byrow = TRUE)
for (i in 1:(n_dispersal+2)){
  dispersion_matrix<-generate_dispersion_matrix(nSpecies,nCommunity,d[i,])
  IDd=IDmatrix*(n_dispersal+2) + i-1
  write.table(dispersion_matrix,paste(path,"matrix_dispersal_",IDd,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)
}

# TL4 ####
v=c(0,0,0,1)
d=0*v
IDmatrix=4
for (i in seq(from=d_min, to=d_max, by=(d_max-d_min)/n_dispersal)){
  d = c(d,10^i*v)
}
d = matrix(d,
           nrow=n_dispersal+2,
           ncol=nSpecies,
           byrow = TRUE)
for (i in 1:(n_dispersal+2)){
  dispersion_matrix<-generate_dispersion_matrix(nSpecies,nCommunity,d[i,])
  IDd=IDmatrix*(n_dispersal+2) + i-1
  write.table(dispersion_matrix,paste(path,"matrix_dispersal_",IDd,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)
}

###### Asymmetry matrix ######

n_params=6
# m , g , r , D , e , a
d = matrix(c(1,1,1,1,1,1,
             1,1,1,1,1,1),
           nrow=nCommunity,
           ncol=n_params,
           byrow = TRUE)
IDmatrix=0
write.table(d,paste(path,"matrix_asymmetry_",IDmatrix,".txt",sep=""),sep=";",col.names = FALSE,row.names = FALSE)

###### Parameter matrix ######

perturbation<-c("exo_0_0_demo_1_TL1_env_0_0",
                "exo_0_0_demo_1_TL2_env_0_0",
                "exo_0_0_demo_1_TL3_env_0_0",
                "exo_0_0_demo_1_TL4_env_0_0",
                "exo_0_0_demo_1_H_env_0_0",
                "exo_0_0_demo_2_H_env_0_0",
                "exo_0_0_demo_0_0_env_1_H",
                "exo_0_0_demo_2_H_env_1_H",
                "exo_0_0_demo_2_H_env_1_H")
perturbation_ID<-c(1:4,
                   0,
                   5:8)
dispersal<-c(rep("disp_H",n_dispersal+2),
             rep("disp_TL1",n_dispersal+2),
                 rep("disp_TL4",n_dispersal+2))
dispersal_ID<-c(0*(n_dispersal+2) + seq(0,(n_dispersal+1)),
                1*(n_dispersal+2) + seq(0,(n_dispersal+1)),
                4*(n_dispersal+2) + seq(0,(n_dispersal+1)))
asymmetry<-"symmetric"
asymmetry_ID<-0
sigma_exo<-0
sigma_demo<-c(rep(0.001,8),0.00001)
sigma_env<-c(rep(0,6),0.001,0.00001,0.001)

params<-data.frame(type="perturbation",
                   name=perturbation,
                   matrix_ID=perturbation_ID,
                   sigma_exo=sigma_exo,
                   sigma_demo=sigma_demo,
                   sigma_env=sigma_env)
databis<-data.frame(type="dispersal",
                    name=dispersal,
                    matrix_ID=dispersal_ID,
                    sigma_exo=NaN,
                    sigma_demo=NaN,
                    sigma_env=NaN)
params<-rbind(params,databis)
databis<-data.frame(type="asymmetry",
                    name=asymmetry,
                    matrix_ID=asymmetry_ID,
                    sigma_exo=NaN,
                    sigma_demo=NaN,
                    sigma_env=NaN)
params<-rbind(params,databis)
rm(databis)
write.table(params,paste(path,"parameters_matrix.txt",sep=""),sep=";",col.names = TRUE,row.names = FALSE)
