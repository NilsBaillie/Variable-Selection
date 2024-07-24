### Variations des hyperparamètres

### étude pour les Spike and Slab continus Normal et Student

rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")  # répertoire à indiquer
library(latex2exp)
library(nimble)
library(MASS)
library(abind)
library(doParallel)

save <- TRUE   # sauvegarde des données et résultats : TRUE si oui, FALSE sinon


### Choix des hyperparamètres : -----------------------
nu0 = 1
df = 10


### Fonctions et paramètres : -------------------------

f <- nimbleFunction(
  run = function(phi=double(0),psi1=double(0),psi2=double(0),t=double(0)) {
    result <- psi1 / (1+exp(-(t-phi)/psi2))
    return(result)
    returnType(double(0))
  })

Misc_quantities <- function(seuil,Summ,beta){
  nbVar <- dim(Summ)[1]
  if(seuil == 0){
    lower <- Summ[1:nbVar,4]
    upper <- Summ[1:nbVar,5]
    negli <- (0 >= lower) & (0 <= upper)
    TP <- sum((beta != 0) & !negli)
    TN <- sum((beta == 0) & negli)
    FP <- sum((beta == 0) & !negli)
    FN <- sum((beta != 0) & negli)
  } else{
    TP <- sum((beta != 0) & (abs(Summ[1:nbVar,1]) > seuil))           # 'vrai' & 'estim'
    TN <- sum((beta == 0) & (abs(Summ[1:nbVar,1]) < seuil))
    FP <- sum((beta == 0) & (abs(Summ[1:nbVar,1]) > seuil))
    FN <- sum((beta != 0) & (abs(Summ[1:nbVar,1]) < seuil))}
  Misc <- (FP+FN)/nbVar
  Sensi <- TP/(TP+FN)
  Speci <- TN/(TN+FP)
  Preci <- TP/(TP+FP)
  return(list(Misc=Misc,Sensi=Sensi,Speci=Speci,Preci=Preci,FP=FP,FN=FN))
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

# Aire sous une courbe
compute_auc <- function(x,y){
  # tri des valeurs
  tri_indices <- order(x)
  x_tri <- x[tri_indices]
  y_tri <- y[tri_indices]
  # règle des trapèzes
  auc <- sum(diff(x_tri) * (y_tri[-1] + y_tri[-length(y_tri)]) / 2)
  return(auc)
}


# Paramètres généraux
{
  n = 50
  p = 100
  nchains = 2        # nombre de chaînes de Markov lancées
  niter = 15000      # nombre d'itérations dans le MCMC
  nburnin = 5000     # Temps de chauffe du MCMC
  nbModels = 5 
  sigma = sqrt(30)
  gamma = 200
  psi1 = 200
  psi2 = 300
  J = 10
  t = seq(150,3000,length=J)
  mu = 1200
  beta_true = c(100,50,20,rep(0,p-3))
  
  # corrélations
  rho_sigma = 0
  Matr_rho = matrix(0,nrow=p,ncol=p)
  for(j in 1:p) {Matr_rho[,j] <- rho_sigma^(abs(1:p-j))}
  Sigma <- (rho_sigma==0)*diag(p)+(rho_sigma!=0)*Matr_rho
  
  # autres hyperparamètres sur les priors
  nu1 = 10^4
  a = 1      # hyperparamètres sur Beta dans les S&S
  b = p
}

Constnorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)
Conststud<- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,nu0=nu0,nu1=nu1,mu=mu,df=df,sigma=sigma,gamma=gamma,t=t)

Inits <- vector('list',nchains)
Inits[[1]] <- list(beta = rep(0,p))
Inits[[2]] <- list(beta = rep(10,p))

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

Num_Model <- 1  # à changer selon le modèle voulu

List_Code <- c(SpikeSlabNorm,SpikeSlabStudent)
Names_Model <- c('S&S N', 'S&S T')
List_Const <- list(Constnorm,Conststud)
Code <- List_Code[[Num_Model]]
Const <- List_Const[[Num_Model]]
Name <- Names_Model[Num_Model]
c(n,p,nu0,df,Name)


### Calcul en parallèle sur plusieurs datasets : ----------------------------------------------

ncores <- detectCores() - 1        # nombre de coeurs utilisés
#ncores <- 5
#cluster <- makeCluster(ncores) 
cluster <- makeCluster(ncores,outfile = 'output_hyperparam.txt')   # si on veut suivre l'avancement des itérations
registerDoParallel(cluster)

nbDataset <- 100                   # nombre de datasets souhaité
nechant = nbDataset %/% ncores
rem_echant <- nbDataset %% ncores
seed_clust <- 1:ncores
seed_clust_rem <- 100 * 1:rem_echant


# Pré-génération de toutes les données pour le calcul parallèle 
set.seed(1)
All_V <- array(0,dim = c(n,p,nbDataset))
All_y <- array(0,dim = c(n,J,nbDataset))
for(z in 1:nbDataset){
  All_V[,,z] <- mvrnorm(n,rep(0,p),Sigma) 
  All_y[,,z] <- Gener_Data(All_V[,,z],beta_true,sigma,gamma)
}

if(save){
  # sauvegarde des datasets
  Covar_NamefileRDS <- paste('Covar_FPandFNs_n=',n,'_p=',p,'_nu0=',nu0,'_ds=',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
  Y_NamefileRDS <- paste('Y_data_FPandFN_n=',n,'_p=',p,'_nu0=',nu0,'_ds=',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
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
    library(nimble)
    library(MASS)
    library(coda)
    library(gridExtra)
    library(abind)
    Model <- nimbleModel(code = Code,const = Const, data = Data,inits = Inits)
    Cmodel <- compileNimble(Model)
    conf <- configureMCMC(model = Cmodel,monitors = 'beta')
    Build <- buildMCMC(conf)
    Cmcmc <- compileNimble(Build)
    Summ <- array(0,dim=c(p,5,nechant))
    for(z in 1:nechant){
      print(paste0("Itération n°",z))
      V <- Covar[,,z]
      y_data <- Y[,,z]
      Data <- list(y=y_data,V=V)
      Cmodel$setData(Data)
      MCMC <- runMCMC(Cmcmc,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
      Summ[,,z] <- MCMC$summary$all.chains
    }
    return(Summ)
  })
  if(rem_echant != 0){
    clusterApply(cluster[1:rem_echant],seed_clust_rem,fun = function(z){set.seed(z)})
    for(z in 1:rem_echant){
      Covar_rem <- array(0,dim = c(n,p,1))
      Y_rem <- array(0,dim = c(n,J,1))
      Covar_rem[,,1] <- All_V[,,nbDataset-rem_echant+z]
      Y_rem[,,1] <- All_y[,,nbDataset-rem_echant+z]
      clusterExport(cluster[z],c('Covar_rem','Y_rem'))}
    
    Out_rem <- clusterEvalQ(cluster[1:rem_echant],{
      V <- Covar_rem[,,1]
      y_data <- Y_rem[,,1]
      Data <- list(y=y_data,V=V)
      Cmodel$setData(Data)
      Summ <- array(0,dim=c(p,5,1))
      MCMC <- runMCMC(Cmcmc,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
      Summ[,,1] <- MCMC$summary$all.chains
      return(Summ)
    })}
  
  t_final = difftime(Sys.time(),t_init,units='secs')
  stopCluster(cluster)
}
t_final

Res <- NULL
for(z in 1:ncores) {Res<- abind(Res,Out[[z]],along=3)}    

if(rem_echant != 0){
  for(z in 1:rem_echant) {Res <- abind(Res,Out_rem[[z]],along=3) }}


if(Num_Model == 1){filename <- paste('Var_Hyperparam_',Name,'_nu0=',nu0,'_Corr',rho_sigma,'.rds',sep='')} else 
  if(Num_Model == 2){filename <- paste('Var_Hyperparam_',Name,'_nu0=',nu0,'_df=',df,'_Corr',rho_sigma,'.rds',sep='')}

# sauvegarde des summary (sous forme de tableaux 3d, 'pile de matrices')
if(save){saveRDS(object = Res,file = filename)}

# récupération des summary si ils ont déjà été calculés
Res <- readRDS(file = filename)

Final <- apply(Res,2,rbind)
seuil = seq(0.01,20,by=0.01)
nb_seuil = length(seuil)
Misc_rate <- array(0,nb_seuil)

for(k in 1:nb_seuil){Misc_rate[k] <- Misc_quantities(seuil[k],Final,beta_true)$Misc}

# plot misc. rate
plot(seuil,Misc_rate,type='l',col='blue',ylab = 'Misc. rate',main="Taux d'erreur en fonction du seuil",
     lwd=2)
grid()

seuil = c(seq(0.001,1,by=0.0001),seq(1,200,by=0.1))
nb_seuil = length(seuil)
Sensi <- array(0,nb_seuil)
Speci <- array(0,nb_seuil)
Preci <- array(0,nb_seuil)

for(k in 1:nb_seuil){
  Sensi[k] <- Misc_quantities(seuil[k],Final,beta_true)$Sensi
  Speci[k] <- Misc_quantities(seuil[k],Final,beta_true)$Speci
  Preci[k] <- Misc_quantities(seuil[k],Final,beta_true)$Preci}

# ROC
plot(1-Speci,Sensi,type='l',col='blue',main='Courbes ROC',ylab = 'Sensibility',xlab = '1-Specificity',
     lwd=2)
lines(seq(0,1,length=nb_seuil),seq(0,1,length=nb_seuil),col='gray',lwd=2)
grid()

# PR
plot(Sensi,Preci,type='l',col='blue',main='Courbes PR',ylab = 'Precision',xlab = 'Recall',
     lwd=2)
lines(seq(0,1,length=nb_seuil),rep(sum(beta_true!=0)/p,length=nb_seuil),col='gray',lwd=2)
grid()


################################################ Plots avec les différentes courbes ###############################################

### S&S normal :  --------------------- 

nu0 = c(0.01,0.1,0.5,1,1.5,5)
len_nu0 <- length(nu0)
filename <- paste('Var_Hyperparam_S&S N_nu0=',nu0,'_Corr',rho_sigma,'.rds',sep = '')

RES <- array(dim=c(p,5,nbDataset,len_nu0))
Final <- array(dim=c(p*nbDataset,5,len_nu0))
seuil = seq(0.01,20,by=0.01)
nb_seuil = length(seuil)
Misc_rate <- array(dim = c(nb_seuil,len_nu0))

for(z in 1:len_nu0){
RES[,,,z] <- readRDS(filename[z])
Final[,,z] <- apply(RES[,,,z],2,rbind)
}

for(z in 1:len_nu0)
{
for(k in 1:nb_seuil){
  Misc_rate[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Misc}}

### NB : légendes adaptées au format pdf

firstpdf <- paste('Var_hyper_Misc.rate_S&S_N_Corr',rho_sigma,'.pdf',sep='')
pdf(firstpdf,width = 5,height = 6)
plot(seuil,Misc_rate[,1],type='l',col='blue',ylab = 'Misc. rate',main="Taux d'erreur en fonction du seuil",
     lwd=2)
lines(seuil,Misc_rate[,2],col='red',lwd=2)
lines(seuil,Misc_rate[,3],col='magenta',lwd=2)
lines(seuil,Misc_rate[,4],col='green',lwd=2)
lines(seuil,Misc_rate[,5],col='orange',lwd=2)
lines(seuil,Misc_rate[,6],col='black',lwd=2)
grid()
if(rho_sigma==0){
legend(x=c(10,18),y=c(0.22,0.35),legend = c(TeX('$\\nu_0 = 0.01$'),TeX('$\\nu_0 = 0.1$'),TeX('$\\nu_0 = 0.5$'),TeX('$\\nu_0 = 1$'),
       TeX('$\\nu_0 = 1.5$'),TeX('$\\nu_0 = 5$')),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)}
if(rho_sigma==0.9){
  legend(x=c(10,18),y=c(0.12,0.2),legend = c(TeX('$\\nu_0 = 0.01$'),TeX('$\\nu_0 = 0.1$'),TeX('$\\nu_0 = 0.5$'),TeX('$\\nu_0 = 1$'),
                                              TeX('$\\nu_0 = 1.5$'),TeX('$\\nu_0 = 5$')),
         col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)}
dev.off()


seuil = c(seq(10^-8,0.1,by=0.00001),seq(0.1,1,by=0.01),seq(1,200,by=0.1))
nb_seuil = length(seuil)
Sensi <- array(dim = c(nb_seuil,len_nu0))
Speci <- array(dim = c(nb_seuil,len_nu0))
Preci <- array(dim = c(nb_seuil,len_nu0))

for(z in 1:len_nu0)
{
  for(k in 1:nb_seuil){
    Speci[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Speci
    Sensi[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Sensi
    Preci[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Preci
}}


secondpdf <- paste('Var_hyperROC_PR_S&S_N_Corr',rho_sigma,'.pdf',sep='')
pdf(secondpdf,width = 6,height = 6)
{plot(1-Speci[,1],Sensi[,1],type='l',col='blue',main='Courbes ROC',ylab = 'Sensibility',xlab = '1-Specificity',
     ylim = c(0,1),
     lwd=2)
lines(1-Speci[,2],Sensi[,2],col='red',lwd=2)
lines(1-Speci[,3],Sensi[,3],col='magenta',lwd=2)
lines(1-Speci[,4],Sensi[,4],col='green',lwd=2)
lines(1-Speci[,5],Sensi[,5],col='orange',lwd=2)
lines(1-Speci[,6],Sensi[,6],col='black',lwd=2)
lines(seq(0,1,length=nb_seuil),seq(0,1,length=nb_seuil),col='gray',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
lines(rep(0,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
grid()
if(rho_sigma==0){
legend(x=c(0.6,1),y=c(0.09,0.55),
        legend = c(TeX('$\\nu_0 = 0.01$'),TeX('$\\nu_0 = 0.1$'),TeX('$\\nu_0 = 0.5$'),TeX('$\\nu_0 = 1$'),
                  TeX('$\\nu_0 = 1.5$'),TeX('$\\nu_0 = 5$'),'Random','Règle parfaite'),
        col = c('blue','red','magenta','green','orange','black','gray','black'),lty=c(1,1,1,1,1,1,1,2),lwd=2)
}

if(rho_sigma==0.9){
legend(x=c(0.6,1),y=c(0.09,0.55),
       legend = c(TeX('$\\nu_0 = 0.01$'),TeX('$\\nu_0 = 0.1$'),TeX('$\\nu_0 = 0.5$'),TeX('$\\nu_0 = 1$'),
                  TeX('$\\nu_0 = 1.5$'),TeX('$\\nu_0 = 5$'),'Random','Règle parfaite'),
       col = c('blue','red','magenta','green','orange','black','gray','black'),lty=c(1,1,1,1,1,1,1,2),lwd=2)
}}


Preci[is.nan(Preci)] <- 1

{plot(Sensi[,1],Preci[,1],type='l',col='blue',main='Courbes PR',ylab = 'Precision',xlab = 'Recall',
     ylim = c(0,1),
     lwd=2)
lines(Sensi[,2],Preci[,2],col='red',lwd=2)
lines(Sensi[,3],Preci[,3],col='magenta',lwd=2)
lines(Sensi[,4],Preci[,4],col='green',lwd=2)
lines(Sensi[,5],Preci[,5],col='orange',lwd=2)
lines(Sensi[,6],Preci[,6],col='black',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(sum(beta_true!=0)/p,length=nb_seuil),col='gray',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
lines(rep(1,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
grid()
if(rho_sigma==0){
  legend(x=c(0,0.45),y=c(0.22,0.7),
        legend = c(TeX('$\\nu_0 = 0.01$'),TeX('$\\nu_0 = 0.1$'),TeX('$\\nu_0 = 0.5$'),TeX('$\\nu_0 = 1$'),
                    TeX('$\\nu_0 = 1.5$'),TeX('$\\nu_0 = 5$'),'Tous significatifs','Règle parfaite'),
         col = c('blue','red','magenta','green','orange','black','gray','black'),lty=c(1,1,1,1,1,1,1,2),lwd=2)
}
if(rho_sigma==0.9){
legend(x=c(0,0.45),y=c(0.22,0.7),
       legend = c(TeX('$\\nu_0 = 0.01$'),TeX('$\\nu_0 = 0.1$'),TeX('$\\nu_0 = 0.5$'),TeX('$\\nu_0 = 1$'),
                  TeX('$\\nu_0 = 1.5$'),TeX('$\\nu_0 = 5$'),'Tous significatifs','Règle parfaite'),
       col = c('blue','red','magenta','green','orange','black','gray','black'),lty=c(1,1,1,1,1,1,1,2),lwd=2)
}}

dev.off()

# ROC-AUC

compute_auc(1-Speci[,1],Sensi[,1])
compute_auc(1-Speci[,2],Sensi[,2])
compute_auc(1-Speci[,3],Sensi[,3])
compute_auc(1-Speci[,4],Sensi[,4])
compute_auc(1-Speci[,5],Sensi[,5])
compute_auc(1-Speci[,6],Sensi[,6])

# PR-AUC

compute_auc(Sensi[,1],Preci[,1])
compute_auc(Sensi[,2],Preci[,2])
compute_auc(Sensi[,3],Preci[,3])
compute_auc(Sensi[,4],Preci[,4])
compute_auc(Sensi[,5],Preci[,5])
compute_auc(Sensi[,6],Preci[,6])



### S&S student :  --------------------- 

df = c(1,5,20,100)
len_df = length(df)
filename <- paste('Var_Hyperparam_S&S T_nu0=1_df=',df,'_Corr',rho_sigma,'.rds',sep = '')

RES <- array(dim=c(p,5,nbDataset,len_df))
Final <- array(dim=c(p*nbDataset,5,len_df))
seuil = seq(0.01,20,by=0.01)
nb_seuil = length(seuil)
Misc_rate <- array(dim = c(nb_seuil,len_df))

for(z in 1:len_df){
  RES[,,,z] <- readRDS(filename[z])
  Final[,,z] <- apply(RES[,,,z],2,rbind)
}
for(z in 1:len_df)
{
  for(k in 1:nb_seuil){
    Misc_rate[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Misc}}

### NB : légendes adaptées au format pdf

firstpdf <- paste('Var_hyper_Misc.rate_S&S_T_Corr',rho_sigma,'.pdf',sep='')
pdf(firstpdf,width = 5,height = 6)
plot(seuil,Misc_rate[,1],type='l',col='blue',ylab = 'Misc. rate',main="Taux d'erreur en fonction du seuil",
     xlim=c(0,5),
     lwd=2)
lines(seuil,Misc_rate[,2],col='red',lwd=2)
lines(seuil,Misc_rate[,3],col='orange',lwd=2)
lines(seuil,Misc_rate[,4],col='green',lwd=2)
grid()
if(rho_sigma==0){
  legend(x=c(3,4.8),y=c(0.17,0.4),legend = c(TeX('$\\nu = 1$'),TeX('$\\nu = 5$'),TeX('$\\nu = 20$'),TeX('$\\nu = 100$')),
         col = c('blue','red','orange','green'),lty=1,lwd=2)}
if(rho_sigma==0.9){
  legend(x=c(3,4.8),y=c(0.17,0.4),legend = c(TeX('$\\nu = 1$'),TeX('$\\nu = 5$'),TeX('$\\nu = 20$'),TeX('$\\nu = 100$')),
         col = c('blue','red','orange','green'),lty=1,lwd=2)}
dev.off()


seuil = c(seq(10^-8,0.1,by=0.00001),seq(0.1,1,by=0.01),seq(1,200,by=0.1))
nb_seuil = length(seuil)
Sensi <- array(dim = c(nb_seuil,len_df))
Speci <- array(dim = c(nb_seuil,len_df))
Preci <- array(dim = c(nb_seuil,len_df))

for(z in 1:len_df)
{
  for(k in 1:nb_seuil){
    Speci[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Speci
    Sensi[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Sensi
    Preci[k,z] <- Misc_quantities(seuil[k],Final[,,z],beta_true)$Preci
  }}


secondpdf <- paste('Var_hyperROC_PR_S&S_T_Corr',rho_sigma,'.pdf',sep='')
pdf(secondpdf,width = 6,height = 6)
{plot(1-Speci[,1],Sensi[,1],type='l',col='blue',main='Courbes ROC',ylab = 'Sensibility',xlab = '1-Specificity',
      ylim = c(0,1),
      lwd=2)
  lines(1-Speci[,2],Sensi[,2],col='red',lwd=2)
  lines(1-Speci[,3],Sensi[,3],col='orange',lwd=2)
  lines(1-Speci[,4],Sensi[,4],col='green',lwd=2)
  lines(seq(0,1,length=nb_seuil),seq(0,1,length=nb_seuil),col='gray',lwd=2)
  lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
  lines(rep(0,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
  grid()
  legend(x=c(0.6,1),y=c(0.2,0.55),
        legend = c(TeX('$\\nu = 1$'),TeX('$\\nu = 5$'),TeX('$\\nu = 20$'),TeX('$\\nu = 100$'),'Random','Règle parfaite'),
        col = c('blue','red','orange','green','gray','black'),lty=c(1,1,1,1,1,2),lwd=2)
}
  

Preci[is.nan(Preci)] <- 1

{plot(Sensi[,1],Preci[,1],type='l',col='blue',main='Courbes PR',ylab = 'Precision',xlab = 'Recall',
      ylim = c(0,1),
      lwd=2)
  lines(Sensi[,2],Preci[,2],col='red',lwd=2)
  lines(Sensi[,3],Preci[,3],col='orange',lwd=2)
  lines(Sensi[,4],Preci[,4],col='green',lwd=2)
  lines(seq(0,1,length=nb_seuil),rep(sum(beta_true!=0)/p,length=nb_seuil),col='gray',lwd=2)
  lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
  lines(rep(1,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
  grid()
  legend(x=c(0,0.45),y=c(0.35,0.7),
           legend = c(TeX('$\\nu = 1$'),TeX('$\\nu = 5$'),TeX('$\\nu = 20$'),TeX('$\\nu = 100$'),'Tous significatifs','Règle parfaite'),
           col = c('blue','red','orange','green','gray','black'),lty=c(1,1,1,1,1,2),lwd=2)}

dev.off()

# ROC-AUC

compute_auc(1-Speci[,1],Sensi[,1])
compute_auc(1-Speci[,2],Sensi[,2])
compute_auc(1-Speci[,3],Sensi[,3])
compute_auc(1-Speci[,4],Sensi[,4])

# PR-AUC

compute_auc(Sensi[,1],Preci[,1])
compute_auc(Sensi[,2],Preci[,2])
compute_auc(Sensi[,3],Preci[,3])
compute_auc(Sensi[,4],Preci[,4])


