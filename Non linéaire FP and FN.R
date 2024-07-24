### Etude sur les faux positifs et faux négatifs

rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")   #  répertoire à indiquer
library(nimble)
library(MASS)
library(latex2exp)
library(gridExtra)
library(abind)

save <- TRUE # sauvegarde des données et résultats : TRUE si oui, FALSE sinon

### Fonction du modèle :

f <- nimbleFunction(
  run = function(phi=double(0),psi1=double(0),psi2=double(0),t=double(0)) {
    result <- psi1 / (1+exp(-(t-phi)/psi2))
    return(result)
    returnType(double(0))
  })

### Paramètres 
{
p = 100
gamma = 200
n = 50
nchains = 2        # nombre de chaînes de Markov lancées
niter = 10000      # nombre d'itérations dans le MCMC
nburnin = 5000     # Temps de chauffe du MCMC
seuil = 1          # méthode de sélection (seuil = 0 pour la méthode avec IC)

# paramètres sur les priors
nu0 = 0.05          
nu1 = 10^4
a = 1
b = p
df = 25
r = 1 
s = 0.1

multivariateNodesAsScalars = TRUE
sigma = sqrt(30)
psi1 = 200
psi2 = 300
J = 10
t = seq(150,3000,length=J)
mu = 1200
beta_true = c(100,50,20,rep(0,p-3))

Cas_Corr = 4     # Cas_Corr = 0 correspond à l'indépendance
rho_sigma = 0.99
Matr_rho = matrix(0,nrow=p,ncol=p)
for(j in 1:p) {Matr_rho[,j] <- rho_sigma^(abs(1:p-j))}
A <- matrix(0,nrow=3,ncol=p-3)
A[3,] <- rho_sigma^(abs(3-4:p))
Corr1 <- rbind(cbind(diag(3),matrix(0,nrow=3,ncol=p-3)),cbind(matrix(0,nrow=p-3,ncol=3),Matr_rho[4:p,4:p]))
Corr2 <- rbind(cbind(diag(3),A),cbind(t(A),diag(p-3))) 
Corr3 <- rbind(cbind(Matr_rho[1:3,1:3],matrix(0,nrow=3,ncol=p-3)),cbind(matrix(0,nrow=p-3,ncol=3),diag(p-3)))
Sigma <- (Cas_Corr == 0)*diag(p)+(Cas_Corr == 1)*Corr1+(Cas_Corr == 2)*Corr2+(Cas_Corr == 3)*Corr3+(Cas_Corr == 4)*Matr_rho
}

Gener_Data <- function(V,beta,sigma,gamma){
  eps = rnorm(n*J,mean = 0,sd = sigma)
  Eps <- matrix(eps,nrow=n,ncol=J,byrow=T)
  ksi <- rnorm(n,0,sd = sqrt(gamma))
  phi <- rep(mu,n) + apply(V,1,function(x) {t(beta) %*% x}) + ksi
  x <- matrix(0,nrow=n,ncol=J)
  for(j in 1:J)  {x[,j] <- f(phi,psi1,psi2,t[j])}
  y_data <- x + Eps
  return(y_data)
}
Misclass_Quant <- function(Cmcmc){
  mcmc <- runMCMC(Cmcmc,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  Summ <- mcmc$summary$all.chains
  if(seuil == 0){
    lower <- Summ[1:p,4]
    upper <- Summ[1:p,5]
    negli <- (0 >= lower) & (0 <= upper)
    FP <- sum((beta_true == 0) & !negli)
    FN <- sum((beta_true != 0) & negli)
  } else{
  FP <- sum((beta_true == 0) & (abs(Summ[1:p,1]) > seuil))
  FN <- sum((beta_true != 0) & (abs(Summ[1:p,1]) < seuil))}
  Result <- c(FP,FN)
  return(Result)
}  # retourne le couple (FP,FN)

Inits <- list(beta = rep(1,p))


# Définition et choix du modèle
SpikeSlabNorm <- nimbleCode({
  alpha ~ dbeta(a,b)
  for(l in 1:p) {
    delta[l] ~ dbern(alpha)
    beta[l] ~ dnorm(0,var = (1-delta[l])*nu0 + delta[l]*nu1)
  }
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})
SpikeSlabStudent <- nimbleCode({
  alpha ~ dbeta(a,b)
  for(l in 1:p) {
    delta[l] ~ dbern(alpha)
    beta[l] ~ dt(mu = 0,tau = 1/((1-delta[l])*nu0 + delta[l]*nu1),df = df)}  
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})
Dirac_Norm <- nimbleCode({
  alpha ~ dbeta(a,b)
  for(l in 1:p) {
    delta[l] ~ dbern(alpha)
    slab[l] ~  dnorm(0,var=nu1)
    beta[l] <- delta[l]*slab[l]}   
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})
PriorHorseshoe <- nimbleCode({
  omega_tau ~ dinvgamma(1/2,1/gamma^2)
  tau2 ~ dinvgamma(1/2,1/omega_tau)
  for(l in 1:p){
    omega_lambda[l] ~ dinvgamma(1/2,1)
    lambda2[l] ~ dinvgamma(1/2,1/omega_lambda[l])
    beta[l] ~ dnorm(0, var = lambda2[l]*tau2)}             
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})
PriorHorseshoePlus <- nimbleCode({
  #tau ~ dunif(0,1)
  #tau0 ~ dunif(0,1)
  omega_tau ~ dinvgamma(1/2,1)
  tau2 ~ dinvgamma(1/2,1/omega_tau)
  for(l in 1:p){ 
    omega_eta[l] ~ dinvgamma(1/2,1)
    eta2[l] ~ dinvgamma(1/2,1/omega_eta[l])
    omega_lambda[l] ~ dinvgamma(1/2,1/(eta2[l]))
    #omega_lambda[l] ~ dinvgamma(1/2,1/(eta2[l]*tau0^2))
    lambda2[l] ~ dinvgamma(1/2,1/omega_lambda[l])
    beta[l] ~ dnorm(0,var = lambda2[l]*tau2)}
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})
PriorLaplace <- nimbleCode({
  lambda2 ~ dgamma(shape = r, rate = s)
  for(l in 1:p){BERN[l] ~ dbern(0.5)
    EXP[l] ~ dexp(rate = sqrt(lambda2/gamma))
    beta[l] <- (2*BERN[l]-1)*EXP[l]}
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})

Constnorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,nu0=nu0,nu1=nu1,a=a,b=b,mu=mu,sigma=sigma,gamma=gamma,t=t)
Conststud <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,nu0=nu0,nu1=nu1,a=a,b=b,mu=mu,df=df,sigma=sigma,gamma=gamma,t=t)
ConstHS <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t)
Constlapl <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,s=s,r=r)

Num_Model <- 1  # à changer selon le modèle voulu

List_Code <- c(SpikeSlabNorm,SpikeSlabStudent,Dirac_Norm,PriorHorseshoe,PriorHorseshoePlus,PriorLaplace)
Names_Model <- c('S&S Continuous N', 'S&S Continous T','S&S Dirac N','Horseshoe','Horseshoe+','Laplace')
List_Const <- list(Constnorm,Conststud,Constnorm,ConstHS,ConstHS,Constlapl)
Code <- List_Code[[Num_Model]]
Const <- List_Const[[Num_Model]]
Name <- Names_Model[Num_Model]
c(n,p,nu0,gamma,Name,Cas_Corr)


### Tests en séquentiel (ne pas exécuter si on souhaite faire du parallèle) : -----------------------------------

# Génération du premier dataset
set.seed(1)
{
  V <- mvrnorm(n,rep(0,p),Sigma)   # Vi : ième ligne de V
  #det(t(V)%*%V)
  #eigen(t(V)%*%V)$values
  eps = rnorm(n*J,mean = 0,sd = sigma)
  Eps <- matrix(eps,nrow=n,ncol=J,byrow=T)
  ksi <- rnorm(n,0,sd = sqrt(gamma))
  phi <- rep(mu,n) + apply(V,1,function(x) {t(beta_true) %*% x}) + ksi
  
  x <- matrix(0,nrow=n,ncol=J)
  for(j in 1:J){
    x[,j] <- f(phi,psi1,psi2,t[j])}
  
  Df_V <- as.data.frame(V)
  y_data <- x + Eps
}
Data <- list(y=y_data,V=V)

Model <- nimbleModel(code = Code,const = Const, data = Data,inits = Inits)
Cmodel <- compileNimble(Model)
conf <- configureMCMC(model = Cmodel,monitors = 'beta')
Build <- buildMCMC(conf)
Cmcmc <- compileNimble(Build)

Misclass_Quant(Cmcmc) 

### Calcul en parallèle : -------------------------------

library(doParallel)
ncores <- detectCores() - 1        # nombre de coeurs utilisés
cluster <- makeCluster(ncores,outfile = 'output_FPandFN.txt')   # type = 'FORK' n'est pas disponible sur Windows
registerDoParallel(cluster)

nbDataset <- 100
nechant = nbDataset %/% ncores
rem_echant <- nbDataset %% ncores
seed_clust <- 1:ncores
seed_clust_rem <- 100 * 1:rem_echant

# Pré-génération de toutes les données pour le calcul parallèle  (et non au fur et à mesure)
set.seed(1)
All_V <- array(0,dim = c(n,p,nbDataset))
All_y <- array(0,dim = c(n,J,nbDataset))
for(z in 1:nbDataset){
  All_V[,,z] <- mvrnorm(n,rep(0,p),Sigma) 
  All_y[,,z] <- Gener_Data(All_V[,,z],beta_true,sigma,gamma)}

# sauvegarde des datasets
if(save){
Covar_NamefileRDS <- paste('Covar_FPandFNs_n=',n,'_p=',p,'_nu0=',nu0,'_ds=',nbDataset,'_Corr=',Cas_Corr,'.rds',sep='')
Y_NamefileRDS <- paste('Y_data_FPandFN_n=',n,'_p=',p,'_nu0=',nu0,'_ds=',nbDataset,'_Corr=',Cas_Corr,'.rds',sep='')
saveRDS(All_V,file = Covar_NamefileRDS)
saveRDS(All_y,file = Y_NamefileRDS)}

clusterExport(cluster,ls())
clusterApply(cluster,seed_clust,fun = function(z){set.seed(z)})  # reproductibilité pour le calcul parallèle

for(z in 1:ncores){
  Covar <- array(0,dim = c(n,p,nechant))
  Y <- array(0,dim = c(n,J,nechant))
  for(k in 1:nechant){                                  # assigne une série de datasets à chaque cluster
    Covar[,,k] <- All_V[,,(z-1)*nechant+k]
    Y[,,k] <- All_y[,,(z-1)*nechant+k]  }
  Data <- list(y=Y[,,1],V=Covar[,,1])
  clusterExport(cluster[z],c("Covar","Y","Data"))
}

t_init <- Sys.time()

{
  Out <- clusterEvalQ(cluster,{
    library(MASS)
    library(nimble)

    Model <- nimbleModel(code = Code,const = Const, data = Data,inits = Inits)
    Cmodel <- compileNimble(Model)
    conf <- configureMCMC(model = Cmodel,monitors = 'beta')
    Build <- buildMCMC(conf)
    Cmcmc <- compileNimble(Build)
    Result <- array(dim = c(nechant,2))
    for(z in 1:nechant){
      print(paste0("Itération n°",z))
      V <- Covar[,,z]
      y_data <- Y[,,z]
      Data <- list(y=y_data,V=V)
      Cmodel$setData(Data)
      Result[z,] <- Misclass_Quant(Cmcmc)
    }
    return(Result)
})
  if(rem_echant != 0){
  clusterApply(cluster[1:rem_echant],seed_clust_rem,fun = function(z){set.seed(z)})
  for(z in 1:rem_echant){
    Covar_rem <- array(0,dim = c(n,p,1))
    Y_rem <- array(0,dim = c(n,J,1))
    Covar_rem[,,1] <- All_V[,,nbDataset-rem_echant+z]
    Y_rem[,,1] <- All_y[,,nbDataset-rem_echant+z]
    clusterExport(cluster[z],c('Covar_rem','Y_rem'))
  }
  Out_rem <- clusterEvalQ(cluster[1:rem_echant],{
    
    Result <- array(dim = c(1,2))
    V <- Covar_rem[,,1]
    y_data <- Y_rem[,,1]
    Data <- list(y=y_data,V=V)
    Cmodel$setData(Data)
    Result[1,] <- Misclass_Quant(Cmcmc)
    return(Result)
})}
  
t_final = difftime(Sys.time(),t_init,units='secs')
stopCluster(cluster)
}
t_final

Matr <- NULL
for(z in 1:ncores) {Matr <- rbind(Matr,Out[[z]])}
if(rem_echant != 0){
for(z in 1:rem_echant) {Matr <- rbind(Matr,Out_rem[[z]])}}

c(max(Matr[,1]),max(Matr[,2]),max(Matr[,1]+Matr[,2]))   # Tous les résultats sont dans 'Matr'
Misc <- rep(0,4)  # Exact, FP but not FN, FN but not FP, FP and FN

for(ds in 1:nbDataset){
  if(Matr[ds,1]+Matr[ds,2] == 0){Misc[1] <- Misc[1] + 1}
  if(Matr[ds,1] != 0 & Matr[ds,2] == 0){Misc[2] <- Misc[2] + 1}
  if(Matr[ds,1] == 0 & Matr[ds,2] != 0){Misc[3] <- Misc[3] + 1}
  if(Matr[ds,1] != 0 & Matr[ds,2] != 0){Misc[4] <- Misc[4] + 1}
}
Misc
Freq_Misc <- round((Misc/nbDataset)*100,2)
Freq_Misc


# Sauvegarde des résultats
if(save){
Res_NamefileRDS <- paste('results_FPandFN_',Name,'_n=',n,'_p=',p,'_nu0=',nu0,'_ds=',nbDataset,'_Corr=',Cas_Corr,'.rds',sep='')
saveRDS(Matr, file = Res_NamefileRDS)}

### Spike and Slab Normal  ----

### Indépendance
### variation en nu0 : 

# nu0 = 10^-6 : 7  7 47 39 -> jusqu'à 76 erreurs !
# nu0 = 0.001 : 3 20 13 64 -> jusqu'à 42 erreurs !
# nu0 = 0.01 : 29 44  7 20 -> jusqu'à 46 erreurs !
# nu0 = 0.05 : 63 32  3  2
# nu0 = 0.1 : 73 26  0  1
# nu0 = 0.3 : 77 21  2  0
# nu0 = 0.5 : 77 22  1  0
# nu0 = 0.6 : 76 24  0  0
# nu0 = 0.7 : 83 17  0  0
# nu0 = 0.8 : 72 28  0  0
# nu0 = 1 : 70 30  0  0
# nu0 = 1.5 : 45 55  0  0
# nu0 = 5 : 0 100   0   0
# nu0 = 10 : 0 100   0   0  -> jusqu'à 56 erreurs !

val_nu0 = c(0.01,0.05,0.1,0.3,0.5,0.6,0.7,0.8,1,1.5)
val_exact = c(29,63,73,77,77,76,83,72,70,45)
val_fpnotfn = c(44,32,26,21,22,24,17,28,30,55)
val_fnnotfp = c(7,3,0,2,1,0,0,0,0,0)
val_fpandfn = c(20,2,1,0,0,0,0,0,0,0)

plot(val_nu0,val_exact,type='l',col='blue',lwd=2,ylim = c(0,85),xlab=TeX('$\\nu_0$'),ylab='% dataset',
     main = TeX("Types d'erreurs dans le S&S Continuous Normal en fonction de $\\nu_0$"))
lines(val_nu0,val_fpnotfn,col = 'red',lwd=2)
lines(val_nu0,val_fnnotfp,col = 'green',lwd=2)
lines(val_nu0,val_fpandfn,col = 'orange',lwd=2)
grid()
legend('topright',legend = c('Exact','FP not FN','FN not FP','FP and FN'),col=c('blue','red','green','orange'),lwd=2)

### Corrélation entre toutes covariables : 

# nu0 = 0.01 : 0  1 53 46  (max 21 erreurs)
# nu0 = 0.05 : 2  0 51 47 (max 21 erreurs)
# nu0 = 0.1 : 1  0 54 45  (max 10 erreurs)
# nu0 = 0.3 : 2  2 55 41 (max 9 erreurs)
# nu0 = 0.5 : 1  4 60 35
# nu0 = 0.6 : 2  1 53 44 (max 7 erreurs)
# nu0 = 0.7 : 3  2 54 41 (max 13 erreurs)
# nu0 = 0.8 : 2  4 52 42
# nu0 = 1 : 3  4 56 37
# nu0 = 1.5 : 2  4 47 47 (max 6 erreurs)
# nu0 = 5 : 1  5 34 60 (max 9 erreurs)
# nu0 = 10 : 0 20  2 78 (max 25 erreurs)

val_nu0 = c(0.01,0.05,0.1,0.3,0.5,0.6,0.7,0.8,1,1.5)
val_exact = c(0, 2, 1 ,2,1,2,3,2,3,2)
val_fpnotfn = c(1, 0,0,2,4,1,2,4,4,5)
val_fnnotfp = c(53,51 ,54,55,60,53,54,52,56,47)
val_fpandfn = c(46,47 ,45,41,35,44,41,42,37,47)

plot(val_nu0,val_exact,type='l',col='blue',lwd=2,ylim = c(0,65),xlab=TeX('$\\nu_0$'),ylab='% dataset',
     main = TeX("Types d'erreurs dans le S&S Continuous Normal en fonction de $\\nu_0$"))
lines(val_nu0,val_fpnotfn,col = 'red',lwd=2)
lines(val_nu0,val_fnnotfp,col = 'green',lwd=2)
lines(val_nu0,val_fpandfn,col = 'orange',lwd=2)
grid()
legend(x=c(1.15,1.55),y=c(15,35),legend = c('Exact','FP not FN','FN not FP','FP and FN'),col=c('blue','red','green','orange'),lwd=2)


### Dirac Normal ----

#### Indép ####

#  82 18  0  0  (max 3 erreurs)

#### Corr ####

#   4  5 63 28   (max 10 erreurs)

### HS ----

#### Indép ####

#  3 97  0  0 (max 22 erreurs)

#### Corr ####

#  1 44  9 46  (max 45 erreurs)

### HS+ ----

#### Indép ####

#  4 96  0  0  (max 17 erreurs)

#### Corr ####

#   1 36 10 53  (max 33 erreurs)

### Laplace ----

#### Indép ####

#  0 100   0   0 (max 68 erreurs)

#### Corr ####

# 0 100   0   0  (max 88 erreurs)

######## Indep
df <- matrix(0,nrow = 4,ncol=6)
df <- as.data.frame(df)
colnames(df) <- c('Model','Exact','FP not FN','FN not FP','FP and FN','#Max error')
df[,1] <- c('S&S Dirac Normal','Horseshoe','Horseshoe+','Laplace')
df[1,-1] <- c(82, 18,  0,  0,3)
df[2,-1] <- c(3, 97,  0,  0,22)
df[3,-1] <- c(4, 96,  0,  0, 17)
df[4,-1] <- c(0, 100,   0,   0,68)
df
myt <- ttheme_default(
  core = list(fg_params=list(cex=1.0)),
  colhead = list(fg_params=list(cex=1.0)),
  rowhead = list(fg_params=list(cex=1.0)))
grid.table(df,theme=myt,rows=NULL)


###### Corr

df <- matrix(0,nrow = 4,ncol=6)
df <- as.data.frame(df)
colnames(df) <- c('Model','Exact','FP not FN','FN not FP','FP and FN','#Max error')
df[,1] <- c('S&S Dirac Normal','Horseshoe','Horseshoe+','Laplace')
df[1,-1] <- c( 4,  5, 63, 28,10)
df[2,-1] <- c( 1, 44 , 9 ,46 ,45)
df[3,-1] <- c(1, 36 ,10 ,53 ,33)
df[4,-1] <- c(0, 100,   0,   0,88)
df
myt <- ttheme_default(
  core = list(fg_params=list(cex=1.0)),
  colhead = list(fg_params=list(cex=1.0)),
  rowhead = list(fg_params=list(cex=1.0)))
grid.table(df,theme=myt,rows=NULL)








### Spike and Slab Student ----
### Résultats en grande dimension (n=50, p=100) : ----------------------


########## Indépendance ##########

##### df = 1 (loi de Cauchy) : 

# nu0 = 1 : 0 100   0   0   (max 29 erreurs)
# nu0 = 0.1 : 5 95  0  0   (max 9 erreurs)
# nu0 = 0.01 : 43 57  0  0  (max 4 erreurs)
# nu0 = 0.001 : 79 10  9  2
# nu0 = 0.0001 : 47  1 45  7   (max 4 erreurs)

df <- matrix(0,nrow = 5,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.0001,47 , 1, 45,  7)
df[2,] <- c(0.001,79, 10 , 9,  2)
df[3,] <- c(0.01,43 ,57 , 0,  0 )
df[4,] <- c(0.1,5, 95,  0 , 0 )
df[5,] <- c(1,0, 100 ,  0 ,  0)
df

##### df = 5 : 

# nu0 = 1 : 22 78  0  0 (max 5 erreurs)
# nu0 = 0.1 : 82 17  0  1  (max 2 erreurs)
# nu0 = 0.01 : 38  6 39 17  (max 5 erreurs)
# nu0 = 5 : 0 100   0   0   -> jusqu'à 38 erreurs !

df <- matrix(0,nrow = 4,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.01,38 , 6 ,39, 17)
df[2,] <- c(0.1,82, 17,  0 , 1)
df[3,] <- c(1,22, 78,  0,  0 )
df[4,] <- c(5,0, 100 ,  0 ,  0 )
df

##### df = 10 :

# nu0 = 1 : 54 45  1  0
# nu0 = 0.1 : 68 21  5  6  (max 1 erreur)
# nu0 = 0.01 : 14  3 56 27
# nu0 = 10 : 0 100   0   0  -> jusqu'à 58 erreurs ! 
# nu0 = 5 : 0 100   0   0  -> jusqu'à 37 erreurs !

df <- matrix(0,nrow = 5,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.01,14 , 3, 56, 27)
df[2,] <- c(0.1,68, 21,  5,  6)
df[3,] <- c(1,54 ,45 , 1,  0)
df[4,] <- c(5,0 ,100  , 0  , 0)
df[5,] <- c(10,0, 100 ,  0 ,  0)
df

##### df = 25 : 

# nu0 = 1 : 62 35  2  1  (max 3 erreurs)
# nu0 = 0.1 : 59 18 18  5  (max 4 erreurs)
# nu0 = 0.01 : 5  0 49 46  (max 17 erreurs)
# nu0 = 5 : 0 100   0   0  -> jusqu'à 33 erreurs !

df <- matrix(0,nrow = 4,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.01,5 , 0, 49 ,46)
df[2,] <- c(0.1,59, 18 ,18,  5)
df[3,] <- c(1,62 ,35 , 2,  1 )
df[4,] <- c(5,0, 100 ,  0 ,  0 )
df


########## Corrélation entre toutes les covariables ##########

##### df = 1 (loi de Cauchy) : 

# nu0 = 1 : 0 55  0 45    -> jusqu'à 37 erreurs !
# nu0 = 0.1 : 3 10 10 77  (max 13 erreurs, FN très fréquents)
# nu0 = 0.01 : 1  0 55 44  (max 7 erreurs, FN très fréquents)
# nu0 = 0.001 : 0  0 83 17  (max 5 erreurs, FN très fréquents)

df <- matrix(0,nrow = 4,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.001,0 , 0 ,83, 17 )
df[2,] <- c(0.01,1 , 0 ,55 ,44)
df[3,] <- c(0.1,3, 10, 10, 77  )
df[4,] <- c(1,0 ,55 , 0, 45 )
df

##### df = 5 : 

# nu0 = 1 : 2  2 52 44  (max 9 erreurs, FN très fréquents)
# nu0 = 0.1 :  1  0 72 27  (max 6 erreurs, FN très fréquents)
# nu0 = 0.01 : 0  0 72 28  (max 6 erreurs, FN très fréquents)
# nu0 = 5 : 0 19  5 76   -> jusqu'à 18 erreurs !

df <- matrix(0,nrow = 4,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.01,0,  0 ,72, 28)
df[2,] <- c(0.1,1,  0 ,72, 27)
df[3,] <- c(1,2  ,2 ,52 ,44 )
df[4,] <- c(5,0 ,19,  5 ,76 )
df

##### df = 10 :

# nu0 = 1 : 0  2 52 46  (max 7 erreurs, FN très fréquents)
# nu0 = 0.1 : 0  0 71 29
# nu0 = 10 : 0 37  0 63  -> jusqu'à 37 erreurs !
# nu0 = 5 : 1 10 17 72  (max 12 erreurs, FN très fréquents)

df <- matrix(0,nrow = 4,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.1,0 , 0, 71, 29)
df[2,] <- c(1,0,  2 ,52, 46)
df[3,] <- c(5,1, 10, 17 ,72  )
df[4,] <- c(10,0 ,37,  0, 63 )
df

##### df = 25 : 

# nu0 = 1 : 1  0 55 44   (max 9 erreurs, FN très fréquents)
# nu0 = 0.1 : 0  1 71 28  (max 8 erreurs, FN très fréquents)
# nu0 = 0.01 : 0  0 73 27  (max 6 erreurs, FN très fréquents)
# nu0 = 5 : 1 3 34 62  (max 11 erreurs, FN très fréquents)

df <- matrix(0,nrow = 4,ncol=5)
df <- as.data.frame(df)
colnames(df) <- c('nu_0','Exact','FP not FN','FN not FP','FP and FN')
df[1,] <- c(0.01,0 , 0, 73, 27)
df[2,] <- c(0.1,0,  1 ,71, 28)
df[3,] <- c(1,1 , 0, 55, 44 )
df[4,] <- c(5,1,3, 34, 62 )
df


