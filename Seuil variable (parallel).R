##### Seuil variable, version parallèle

# On veut déterminer comment évolue le taux de mauvais classement pour chacun des priors
# lorsqu'on fait varier la valeur seuil
# On affiche les courbes ROC et PR (Precision-Recall)


rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")  # répertoire à indiquer
library(latex2exp)
library(nimble)
library(MASS)
library(abind)
library(doParallel)

save <- FALSE   # sauvegarde des données et résultats : TRUE si oui, FALSE sinon (évite d'écraser les données si elles sont déjà là)

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
  mcc <- (TP*TN-FP*FN)
  mcc <- mcc /sqrt((TP+FP))
  mcc <- mcc/sqrt((TP+FN))
  mcc <- mcc/sqrt((TN+FP))
  mcc <- mcc/sqrt((TN+TP))
  return(list(Misc=Misc,Sensi=Sensi,Speci=Speci,Preci=Preci,mcc=mcc))
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
  rho_sigma = 0.6
  Matr_rho = matrix(0,nrow=p,ncol=p)
  for(j in 1:p) {Matr_rho[,j] <- rho_sigma^(abs(1:p-j))}
  Sigma <- (rho_sigma==0)*diag(p)+(rho_sigma!=0)*Matr_rho
  
  # hyperparamètres sur les priors
  df = 10                     
  nu0 = 1          
  nu1 = 10^4
  a = 1      # hyperparamètres sur Beta dans les S&S
  b = p
  r = 1      # hyperparamètres sur Gamma dans Laplace
  s = 0.1
}

# Choix de la version utilisée du prior : (par défaut, HSvers = 2 et HSpvers = 2)
HSvers <- 2
HSpvers <- 2

ConstNorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)
Conststud<- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,nu0=nu0,nu1=nu1,mu=mu,df=df,sigma=sigma,gamma=gamma,t=t)
ConstHS <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t)
Constlapl <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,s=s,r=r)

Inits <- vector('list',nchains)
Inits[[1]] <- list(beta = rep(0,p))
Inits[[2]] <- list(beta = rep(10,p))


### Calcul en parallèle sur plusieurs datasets : ----------------------------------------------

ncores <- detectCores() - 1        # nombre de coeurs utilisés
#ncores <- 5
#cluster <- makeCluster(ncores) 
cluster <- makeCluster(ncores,outfile = 'output_seuil_var.txt')   # si on veut suivre l'avancement des itérations
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


# noms des fichiers ----
{filenorm <- paste('Norm_summary_ds',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
  filestud <- paste('Stud_summary_ds',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
  fileDirac <- paste('Dirac_summary_ds',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
  fileHS <- paste('HS_summary_ds',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
  fileHSp <- paste('HSp_summary_ds',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
  filelapl <- paste('Lapl_summary_ds',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')}

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
  SpikeSlab_Norm <- nimbleCode({
    alpha ~ dbeta(a,b)
    for(l in 1:p) {
      delta[l] ~ dbern(alpha)
      beta[l] ~ dnorm(0,var=(1-delta[l])*nu0 + delta[l]*nu1)}    # nu0 : spike, nu1 : slab
    for(i in 1:n){
      phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
      for(j in 1:J){
        y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
      }}
  })
  Modelnorm <- nimbleModel(code = SpikeSlab_Norm,const = ConstNorm, data = Data,inits = Inits)
  Cnorm <- compileNimble(Modelnorm)
  confnorm <- configureMCMC(model = Cnorm, monitors = 'beta')
  buildnorm <- buildMCMC(confnorm)        # conjugate sampler for beta
  Cmcmcnorm <- compileNimble(buildnorm)
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
  ModelDirac <- nimbleModel(code = Dirac_Norm,const = ConstNorm, data = Data,inits = Inits)
  CDirac <- compileNimble(ModelDirac)
  confDirac <- configureMCMC(model = CDirac, monitors = 'beta')
  buildDirac <- buildMCMC(confDirac)
  CmcmcDirac <- compileNimble(buildDirac)
  SpikeSlab_student <- nimbleCode({
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
  Modelstud <- nimbleModel(code = SpikeSlab_student,const = Conststud, data = Data,inits = Inits)
  Cstud <- compileNimble(Modelstud)
  confstud <- configureMCMC(model = Cstud, monitors = 'beta')
  buildstud <- buildMCMC(confstud)        # RW sampler for beta
  Cmcmcstud <- compileNimble(buildstud)
  if(HSvers == 1){
    PriorHorseshoe <- nimbleCode({
      tau ~ T(dt(mu = 0, tau = 1/gamma^2, df = 1),0,10^300)
      for(l in 1:p){
        lambda[l] ~ T(dt(mu = 0, tau = 1, df = 1),0,10^300)
        beta[l] ~ dnorm(0, sd = lambda[l]*tau)}             
      for(i in 1:n){
        phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
        for(j in 1:J){
          y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
        }}
    })
  } else if(HSvers == 2){
     PriorHorseshoe <- nimbleCode({
      omega_tau ~ dinvgamma(1/2,1)
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
  } else if(HSvers == 3){
    PriorHorseshoe <- nimbleCode({
      norm_tau ~ dnorm(0,var = 1)
      IG_tau ~ dinvgamma(1/2,1/2)
      tau <- abs(norm_tau*sqrt(IG_tau))
      for(l in 1:p){
        norm_lambda[l] ~ dnorm(0,var = 1)
        IG_lambda[l] ~ dinvgamma(1/2,1/2)
        lambda[l] <- abs(norm_lambda[l]*sqrt(IG_lambda[l]))
        beta[l] ~ dnorm(0, sd = lambda[l]*tau)}             
      for(i in 1:n){
        phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
        for(j in 1:J){
          y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
        }}
    })
  }
  Horseshoe <- nimbleModel(code = PriorHorseshoe,constants = ConstHS,data = Data,inits = Inits)
  CHS <- compileNimble(Horseshoe)
  confHS <- configureMCMC(model = CHS,monitors = 'beta')  # conjugate sampler for beta
  BuildHS <- buildMCMC(confHS)
  CmcmcHS <- compileNimble(BuildHS)
  if(HSpvers == 1){
    PriorHorseshoePlus <- nimbleCode({
      tau ~ dunif(0,1)
      for(l in 1:p){ eta[l] ~ T(dt(0,tau = 1, df = 1),0,10^300)
        lambda[l] ~ T(dt(0,tau = 1/(eta[l]*tau)^2,df = 1),0,10^300)
        beta[l] ~ dnorm(0,sd = lambda[l])}
      for(i in 1:n){
        phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
        for(j in 1:J){
          y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
        }}
    })
  } else if(HSpvers == 2){
    PriorHorseshoePlus <- nimbleCode({
      #tau ~ dunif(0,1)
      omega_tau ~ dinvgamma(1/2,1)
      tau2 ~ dinvgamma(1/2,1/omega_tau)
      for(l in 1:p){ 
        omega_eta[l] ~ dinvgamma(1/2,1)
        eta2[l] ~ dinvgamma(1/2,1/omega_eta[l])
        omega_lambda[l] ~ dinvgamma(1/2,1/(eta2[l]))
        lambda2[l] ~ dinvgamma(1/2,1/omega_lambda[l])
        beta[l] ~ dnorm(0,var = lambda2[l]*tau2)}
      for(i in 1:n){
        phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
        for(j in 1:J){
          y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
        }}
    })
  } else if(HSpvers == 3){
    PriorHorseshoePlus <- nimbleCode({
      omega_tau ~ dinvgamma(1/2,1)
      tau2 ~ dinvgamma(1/2,1/omega_tau)
      for(l in 1:p){ 
        omega_eta[l] ~ dinvgamma(1/2,1)
        eta2[l] ~ dinvgamma(1/2,1/omega_eta[l])
        omega_lambda[l] ~ dinvgamma(1/2,1/(eta2[l]*tau2))
        lambda2[l] ~ dinvgamma(1/2,1/omega_lambda[l])
        beta[l] ~ dnorm(0,var = lambda2[l])}
      for(i in 1:n){
        phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
        for(j in 1:J){
          y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
        }} 
    })}
  HorseshoePlus <- nimbleModel(code = PriorHorseshoePlus,constants = ConstHS,data = Data,inits = Inits)
  CHSp <- compileNimble(HorseshoePlus)
  confHSp <- configureMCMC(model = CHSp,monitors = 'beta')  # conjugate sampler for beta
  BuildHSp <- buildMCMC(confHSp)
  CmcmcHSp <- compileNimble(BuildHSp) 
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
  Laplace <- nimbleModel(code = PriorLaplace,constants = Constlapl,data = Data,inits = Inits)
  Clapl <- compileNimble(Laplace)
  conflapl <- configureMCMC(model = Clapl, monitors = 'beta')
  Buildlapl <- buildMCMC(conflapl)
  Cmcmclapl <- compileNimble(Buildlapl)
  
  Summnorm  <- array(0,dim=c(p,5,nechant))
  Summstud <- array(0,dim=c(p,5,nechant))
  SummDirac <- array(0,dim=c(p,5,nechant))
  SummHS <- array(0,dim=c(p,5,nechant))
  SummHSp <- array(0,dim=c(p,5,nechant))
  Summlapl <- array(0,dim=c(p,5,nechant))
  
  for(z in 1:nechant){
    print(paste0("Itération n°",z))
    V <- Covar[,,z]
    y_data <- Y[,,z]
    Data <- list(y=y_data,V=V)
    Cnorm$setData(Data)
    Cstud$setData(Data)
    CDirac$setData(Data)
    CHS$setData(Data)
    CHSp$setData(Data)
    Clapl$setData(Data)
    MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    Summnorm[,,z] <- MCMCnorm$summary$all.chains
    MCMCstud <- runMCMC(Cmcmcstud,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    Summstud[,,z] <- MCMCstud$summary$all.chains
    MCMCDirac <- runMCMC(CmcmcDirac,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    SummDirac[,,z] <- MCMCDirac$summary$all.chains
    MCMCHS <- runMCMC(CmcmcHS,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    SummHS[,,z] <- MCMCHS$summary$all.chains
    MCMCHSp <- runMCMC(CmcmcHSp,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    SummHSp[,,z] <- MCMCHSp$summary$all.chains
    MCMClapl <- runMCMC(Cmcmclapl,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    Summlapl[,,z] <- MCMClapl$summary$all.chains
  }
  return(list(Norm=Summnorm,Stud=Summstud,Dirac=SummDirac,HS=SummHS,HSp=SummHSp,Lapl=Summlapl))
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
    Cnorm$setData(Data)
    Cstud$setData(Data)
    CDirac$setData(Data)
    CHS$setData(Data)
    CHSp$setData(Data)
    Clapl$setData(Data)
    Summnorm  <- array(0,dim=c(p,5,1))
    Summstud <- array(0,dim=c(p,5,1))
    SummDirac <- array(0,dim=c(p,5,1))
    SummHS <- array(0,dim=c(p,5,1))
    SummHSp <- array(0,dim=c(p,5,1))
    Summlapl <- array(0,dim=c(p,5,1))
    MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    Summnorm[,,1] <- MCMCnorm$summary$all.chains
    MCMCstud <- runMCMC(Cmcmcstud,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    Summstud[,,1] <- MCMCstud$summary$all.chains
    MCMCDirac <- runMCMC(CmcmcDirac,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    SummDirac[,,1] <- MCMCDirac$summary$all.chains
    MCMCHS <- runMCMC(CmcmcHS,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    SummHS[,,1] <- MCMCHS$summary$all.chains
    MCMCHSp <- runMCMC(CmcmcHSp,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    SummHSp[,,1] <- MCMCHSp$summary$all.chains
    MCMClapl <- runMCMC(Cmcmclapl,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
    Summlapl[,,1] <- MCMClapl$summary$all.chains
    return(list(Norm=Summnorm,Stud=Summstud,Dirac=SummDirac,HS=SummHS,HSp=SummHSp,Lapl=Summlapl))
})}

t_final = difftime(Sys.time(),t_init,units='secs')
stopCluster(cluster)
}
t_final

# mise en forme des résultats
{Resnorm <- NULL
Resstud <- NULL
ResDirac <- NULL
ResHS <- NULL
ResHSp <- NULL
Reslapl <- NULL}

for(z in 1:ncores) {
  Resnorm <- abind(Resnorm,Out[[z]]$Norm,along=3)
  Resstud <- abind(Resstud,Out[[z]]$Stud,along=3)
  ResDirac <- abind(ResDirac,Out[[z]]$Dirac,along=3)
  ResHS <- abind(ResHS,Out[[z]]$HS,along=3)
  ResHSp <- abind(ResHSp,Out[[z]]$HSp,along=3)
  Reslapl <- abind(Reslapl,Out[[z]]$Lapl,along=3)
}    
if(rem_echant != 0){
for(z in 1:rem_echant) {
  Resnorm <- abind(Resnorm,Out_rem[[z]]$Norm,along=3)
  Resstud <- abind(Resstud,Out_rem[[z]]$Stud,along=3)
  ResDirac <- abind(ResDirac,Out_rem[[z]]$Dirac,along=3)
  ResHS <- abind(ResHS,Out_rem[[z]]$HS,along=3)
  ResHSp <- abind(ResHSp,Out_rem[[z]]$HSp,along=3)
  Reslapl <- abind(Reslapl,Out_rem[[z]]$Lapl,along=3)
}}


# sauvegarde des summary (sous forme de tableaux 3d, 'pile de matrices')
if(save){
saveRDS(object = Resnorm,file = filenorm)
saveRDS(object = Resstud,file = filestud)
saveRDS(object = ResDirac,file = fileDirac)
saveRDS(object = ResHS,file = fileHS)
saveRDS(object = ResHSp,file = fileHSp)
saveRDS(object = Reslapl,file = filelapl)}


# récupération des summary si ils ont déjà été calculés
Resnorm <- readRDS(file = filenorm)
Resstud <- readRDS(file = filestud)
ResDirac <- readRDS(file = fileDirac)
ResHS <-  readRDS(file = fileHS)
ResHSp <- readRDS(file = fileHSp)
Reslapl <- readRDS(file = filelapl)

# transformation en tableau 2d
{Finalnorm <- apply(Resnorm,2,rbind)
Finalstud <- apply(Resstud,2,rbind)
FinalDirac <- apply(ResDirac,2,rbind)
FinalHS <- apply(ResHS,2,rbind)
FinalHSp <- apply(ResHSp,2,rbind)
Finallapl <- apply(Reslapl,2,rbind)}


seuil = seq(0.01,20,by=0.01)
nb_seuil = length(seuil)
{
Misc_rate_norm <- array(0,nb_seuil)
Misc_rate_stud <- array(0,nb_seuil)
Misc_rate_Dirac <- array(0,nb_seuil)
Misc_rate_HS <- array(0,nb_seuil)
Misc_rate_HSp <- array(0,nb_seuil)
Misc_rate_lapl <- array(0,nb_seuil)

mcc_norm <- array(0,nb_seuil)
mcc_stud <- array(0,nb_seuil)
mcc_Dirac <- array(0,nb_seuil)
mcc_HS <- array(0,nb_seuil)
mcc_HSp <- array(0,nb_seuil)
mcc_lapl <- array(0,nb_seuil)
}


for(k in 1:nb_seuil){
    Misc_rate_norm[k] <- Misc_quantities(seuil[k],Finalnorm,beta_true)$Misc
    Misc_rate_stud[k] <- Misc_quantities(seuil[k],Finalstud,beta_true)$Misc
    Misc_rate_Dirac[k] <- Misc_quantities(seuil[k],FinalDirac,beta_true)$Misc
    Misc_rate_HS[k] <- Misc_quantities(seuil[k],FinalHS,beta_true)$Misc
    Misc_rate_HSp[k] <- Misc_quantities(seuil[k],FinalHSp,beta_true)$Misc
    Misc_rate_lapl[k] <- Misc_quantities(seuil[k],Finallapl,beta_true)$Misc
    
    mcc_norm[k] <- Misc_quantities(seuil[k],Finalnorm,beta_true)$mcc
    mcc_stud[k] <- Misc_quantities(seuil[k],Finalstud,beta_true)$mcc
    mcc_Dirac[k] <- Misc_quantities(seuil[k],FinalDirac,beta_true)$mcc
    mcc_HS[k] <- Misc_quantities(seuil[k],FinalHS,beta_true)$mcc
    mcc_HSp[k] <- Misc_quantities(seuil[k],FinalHSp,beta_true)$mcc
    mcc_lapl[k] <- Misc_quantities(seuil[k],Finallapl,beta_true)$mcc
}
  
### Plot : -------------------------
firstpdf = paste('Misc.RateCorr',rho_sigma,'.pdf',sep='')

pdf(firstpdf,width = 6,height = 6)

{plot(seuil,Misc_rate_norm,type='l',col='blue',ylab = 'Misc. rate',main="Taux d'erreur en fonction du seuil",
     #xlim = c(1,2),ylim=c(0,0.1),
     #xlim = c(0,2),
     #ylim = c(0,0.05),
     #ylim=c(0,1.2),
     lwd=2)
lines(seuil,Misc_rate_stud,col='red',lwd=2)
lines(seuil,Misc_rate_Dirac,col='magenta',lwd=2)
lines(seuil,Misc_rate_HS,col='green',lwd=2)
lines(seuil,Misc_rate_HSp,col='orange',lwd=2)
lines(seuil,Misc_rate_lapl,col='black',lwd=2)
grid()
legend('topright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)}
#abline(h=0.02)


{plot(seuil,Misc_rate_norm,type='l',col='blue',ylab = 'Misc. rate',main="Taux d'erreur en fonction du seuil",
     xlim = c(0,2),
     lwd=2)
lines(seuil,Misc_rate_stud,col='red',lwd=2)
lines(seuil,Misc_rate_Dirac,col='magenta',lwd=2)
lines(seuil,Misc_rate_HS,col='green',lwd=2)
lines(seuil,Misc_rate_HSp,col='orange',lwd=2)
lines(seuil,Misc_rate_lapl,col='black',lwd=2)
grid()
legend('topright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)}
dev.off()

# Coefficient de corrélation de Matthews
plot(seuil,mcc_norm,type='l',col='blue',ylab = 'MCC',main="MCC en fonction du seuil",
     ylim = c(0,1),
     lwd=2)
lines(seuil,mcc_stud,col='red',lwd=2)
lines(seuil,mcc_Dirac,col='magenta',lwd=2)
lines(seuil,mcc_HS,col='green',lwd=2)
lines(seuil,mcc_HSp,col='orange',lwd=2)
lines(seuil,mcc_lapl,col='black',lwd=2)
grid()
legend('bottomright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)


### Courbes PR et ROC : -----------------------

seuil = c(seq(10^-8,0.1,by=0.00001),seq(0.1,1,by=0.01),seq(1,200,by=0.1))
nb_seuil = length(seuil)

{
Sensi_norm <- array(0,nb_seuil)
Sensi_stud <- array(0,nb_seuil)
Sensi_Dirac <- array(0,nb_seuil)
Sensi_HS <- array(0,nb_seuil)
Sensi_HSp <- array(0,nb_seuil)
Sensi_lapl <- array(0,nb_seuil)
Speci_norm <- array(0,nb_seuil)
Speci_stud <- array(0,nb_seuil)
Speci_Dirac <- array(0,nb_seuil)
Speci_HS <- array(0,nb_seuil)
Speci_HSp <- array(0,nb_seuil)
Speci_lapl <- array(0,nb_seuil)
Preci_norm <- array(0,nb_seuil)
Preci_stud <- array(0,nb_seuil)
Preci_Dirac <- array(0,nb_seuil)
Preci_HS <- array(0,nb_seuil)
Preci_HSp <- array(0,nb_seuil)
Preci_lapl <- array(0,nb_seuil)
}


for(k in 1:nb_seuil){
  Sensi_norm[k] <- Misc_quantities(seuil[k],Finalnorm,beta_true)$Sensi
  Sensi_stud[k] <- Misc_quantities(seuil[k],Finalstud,beta_true)$Sensi
  Sensi_Dirac[k] <- Misc_quantities(seuil[k],FinalDirac,beta_true)$Sensi
  Sensi_HS[k] <- Misc_quantities(seuil[k],FinalHS,beta_true)$Sensi
  Sensi_HSp[k] <- Misc_quantities(seuil[k],FinalHSp,beta_true)$Sensi
  Sensi_lapl[k] <- Misc_quantities(seuil[k],Finallapl,beta_true)$Sensi
  
  Speci_norm[k] <- Misc_quantities(seuil[k],Finalnorm,beta_true)$Speci
  Speci_stud[k] <- Misc_quantities(seuil[k],Finalstud,beta_true)$Speci
  Speci_Dirac[k] <- Misc_quantities(seuil[k],FinalDirac,beta_true)$Speci
  Speci_HS[k] <- Misc_quantities(seuil[k],FinalHS,beta_true)$Speci
  Speci_HSp[k] <- Misc_quantities(seuil[k],FinalHSp,beta_true)$Speci
  Speci_lapl[k] <- Misc_quantities(seuil[k],Finallapl,beta_true)$Speci
  
  Preci_norm[k] <- Misc_quantities(seuil[k],Finalnorm,beta_true)$Preci
  Preci_stud[k] <- Misc_quantities(seuil[k],Finalstud,beta_true)$Preci
  Preci_Dirac[k] <- Misc_quantities(seuil[k],FinalDirac,beta_true)$Preci
  Preci_HS[k] <- Misc_quantities(seuil[k],FinalHS,beta_true)$Preci
  Preci_HSp[k] <- Misc_quantities(seuil[k],FinalHSp,beta_true)$Preci
  Preci_lapl[k] <- Misc_quantities(seuil[k],Finallapl,beta_true)$Preci
}

### ROC : Sensibilité en fonction de 1-Spécificité

secondpdf = paste('ROC_PR_Corr',rho_sigma,'.pdf',sep='')
pdf(secondpdf,width = 6,height = 6)

plot(1-Speci_norm,Sensi_norm,type='l',col='blue',main='Courbes ROC',ylab = 'Sensibility',xlab = '1-Specificity',
     #xlim=c(0,0.01),
     #xlim=c(0.9,1),
     #xlim=c(0,0.01),ylim=c(0.98,1),  # zoom 
     ylim = c(0,1),
     lwd=2)
lines(1-Speci_stud,Sensi_stud,col='red',lwd=2)
lines(1-Speci_Dirac,Sensi_Dirac,col='magenta',lwd=2)
lines(1-Speci_HS,Sensi_HS,col='green',lwd=2)
lines(1-Speci_HSp,Sensi_HSp,col='orange',lwd=2)
lines(1-Speci_lapl,Sensi_lapl,col='black',lwd=2)
lines(seq(0,1,length=nb_seuil),seq(0,1,length=nb_seuil),col='gray',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
lines(rep(0,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
grid()
legend('bottomright',
       legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace','Random','Règle parfaite'),
       col = c('blue','red','magenta','green','orange','black','gray','black'),lty=c(1,1,1,1,1,1,1,2),lwd=2)
#abline(h=0.9)



### PR : Précision en fonction de Sensibilité :

Preci_norm[is.nan(Preci_norm)] <- 1
Preci_stud[is.nan(Preci_stud)] <- 1
Preci_Dirac[is.nan(Preci_Dirac)] <- 1
Preci_HS[is.nan(Preci_HS)] <- 1
Preci_HSp[is.nan(Preci_HSp)] <- 1
Preci_lapl[is.nan(Preci_lapl)] <- 1

plot(Sensi_norm,Preci_norm,type='l',col='blue',main='Courbes PR',ylab = 'Precision',xlab = 'Recall',
     #xlim=c(0.95,1),
     lwd=2)
lines(Sensi_stud,Preci_stud,col='red',lwd=2)
lines(Sensi_Dirac,Preci_Dirac,col='magenta',lwd=2)
lines(Sensi_HS,Preci_HS,col='green',lwd=2)
lines(Sensi_HSp,Preci_HSp,col='orange',lwd=2)
lines(Sensi_lapl,Preci_lapl,col='black',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(sum(beta_true!=0)/p,length=nb_seuil),col='gray',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
lines(rep(1,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
grid()
legend(x=c(0,0.44),y=c(0.12,0.6),
       legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace','Tous significatifs','Règle parfaite'),
       col = c('blue','red','magenta','green','orange','black','gray','black'),lty=c(1,1,1,1,1,1,1,2),lwd=2)
# légende adaptée au pdf ! 
dev.off()


### Aires sous la courbe
compute_auc <- function(x, y) {
  # tri des valeurs
  tri_indices <- order(x)
  x_tri <- x[tri_indices]
  y_tri <- y[tri_indices]
  # règle des trapèzes
  auc <- sum(diff(x_tri) * (y_tri[-1] + y_tri[-length(y_tri)]) / 2)
  return(auc)}

# ROC-AUC

compute_auc(1-Speci_norm,Sensi_norm)
compute_auc(1-Speci_stud,Sensi_stud)
compute_auc(1-Speci_Dirac,Sensi_Dirac)
compute_auc(1-Speci_HS,Sensi_HS)
compute_auc(1-Speci_HSp,Sensi_HSp)
compute_auc(1-Speci_lapl,Sensi_lapl)

# PR-AUC
compute_auc(Sensi_norm,Preci_norm)
compute_auc(Sensi_stud,Preci_stud)
compute_auc(Sensi_Dirac,Preci_Dirac)
compute_auc(Sensi_HS,Preci_HS)
compute_auc(Sensi_HSp,Preci_HSp)
compute_auc(Sensi_lapl,Preci_lapl)


### Sensibilité et Spécificité :

plot(seuil,Sensi_norm,type='l',col='blue',main='Sensibility',
    #xlim = c(0,50),
    ylim = c(0,1),
    lwd=2)
lines(seuil,Sensi_stud,col='red',lwd=2)
lines(seuil,Sensi_Dirac,col='magenta',lwd=2)
lines(seuil,Sensi_HS,col='green',lwd=2)
lines(seuil,Sensi_HSp,col='orange',lwd=2)
lines(seuil,Sensi_lapl,col='black',lwd=2)
grid()
legend('topright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)


plot(seuil,Speci_norm,type='l',col='blue',main='Specificity',
    #xlim = c(0,5),
    ylim = c(0,1),
    lwd=2)
lines(seuil,Speci_stud,col='red',lwd=2)
lines(seuil,Speci_Dirac,col='magenta',lwd=2)
lines(seuil,Speci_HS,col='green',lwd=2)
lines(seuil,Speci_HSp,col='orange',lwd=2)
lines(seuil,Speci_lapl,col='black',lwd=2)
grid()
legend('bottomright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)

plot(seuil,0.5*(Sensi_norm+Speci_norm),type='l',col='blue',main='Average Sensitivity / Specificity',
     #xlim = c(0,5),
     ylim = c(0,1),
     lwd=2)
lines(seuil,0.5*(Sensi_stud+Speci_stud),col='red',lwd=2)
lines(seuil,0.5*(Sensi_Dirac+Speci_Dirac),col='magenta',lwd=2)
lines(seuil,0.5*(Sensi_HS+Speci_HS),col='green',lwd=2)
lines(seuil,0.5*(Sensi_HSp+Speci_HSp),col='orange',lwd=2)
lines(seuil,0.5*(Sensi_lapl+Speci_lapl),col='black',lwd=2)
grid()
legend('bottomright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)


### G-mean = sqrt(Sensi * Speci)  : 

plot(seuil,sqrt(Sensi_norm*Speci_norm),type='l',col='blue',main='G-mean',
     #xlim = c(0,5),
     ylim = c(0,1),
     lwd=2)
lines(seuil,sqrt(Sensi_stud*Speci_stud),col='red',lwd=2)
lines(seuil,sqrt(Sensi_Dirac*Speci_Dirac),col='magenta',lwd=2)
lines(seuil,sqrt(Sensi_HS*Speci_HS),col='green',lwd=2)
lines(seuil,sqrt(Sensi_HSp*Speci_HSp),col='orange',lwd=2)
lines(seuil,sqrt(Sensi_lapl*Speci_lapl),col='black',lwd=2)
grid()
legend('topright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','black'),lty=1,lwd=2)

