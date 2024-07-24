### Simulation par un modèle non linéaire avec tous les priors, version en parallèle


rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")  # à modifier 
library(latex2exp)
library(nimble)
library(MASS)
library(coda)
library(gridExtra)
library(abind)

save <- TRUE  # sauvegarde des données et résultats : TRUE si oui, FALSE sinon

### Fonctions et paramètres : -------------------------

f <- nimbleFunction(
  run = function(phi=double(0),psi1=double(0),psi2=double(0),t=double(0)) {
    result <- psi1 / (1+exp(-(t-phi)/psi2))
    return(result)
    returnType(double(0))
  })

alpha <- 0.05
IS <- function(alpha,lower,upper,true_val)   # fonction auxiliaire pour le calcul du MIS
{return(((upper-lower) + (2/alpha)*(lower-true_val)*(true_val < lower) + (2/alpha)*(true_val-upper)*(true_val > upper)))}

length <- 10^4
CRPS_custom <- function(beta,samples,length){
  P <- ecdf(samples)
  interv <- seq(beta - 10^3, beta + 10^3, length = length)
  CRPS_sing <- sum((P(interv)-(interv >= beta))^2)
  return(CRPS_sing)
}
Misc_quantities <- function(seuil,Summ,beta){ 
  if(seuil == 0){
    lower <- Summ[1:p,4]
    upper <- Summ[1:p,5]
    negli <- (0 >= lower) & (0 <= upper)
    TP <- sum((beta != 0) & !negli)
    TN <- sum((beta == 0) & negli)
    FP <- sum((beta == 0) & !negli)
    FN <- sum((beta != 0) & negli)
  } else{
    TP <- sum((beta != 0) & (abs(Summ[1:p,1]) > seuil))           # 'vrai' & 'estim'
    TN <- sum((beta == 0) & (abs(Summ[1:p,1]) < seuil))
    FP <- sum((beta == 0) & (abs(Summ[1:p,1]) > seuil))
    FN <- sum((beta != 0) & (abs(Summ[1:p,1]) < seuil))}
  Misc <- (FP+FN)/p
  Sensi <- TP/(TP+FN)
  Speci <- TN/(TN+FP)
  return(list(Misc=Misc,Sensi=Sensi,Speci=Speci))
}
ComputeModel <- function(Cmodel,Cmcmc,Data){
  Cmodel$setData(Data)
  start_time = Sys.time()
  MCMC <- runMCMC(Cmcmc,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE,inits = Inits)
  Time <- difftime(Sys.time(),start_time,units='secs')
  Summ <- MCMC$summary$all.chains
  Samples <- NULL
  list_mcmc_objects <- vector('list',nchains)
  for(z in 1:nchains) {
    Samples <- rbind(Samples,MCMC$samples[[z]])
    list_mcmc_objects[[z]] <- mcmc(MCMC$samples[[z]])
  }
  list_mcmc <- mcmc.list(list_mcmc_objects)
  ess <-  mean(effectiveSize(list_mcmc))
  ESS_per_sec <- ess / as.numeric(Time)
  #ind_VIF <- which(abs(Summ[,1]) > seuil)
  #formula <- parse(text = paste("V",ind_VIF[1],sep="","~."))[[1]]
  #modelVIF <- lm(formula = formula, data = Df_V[,ind_VIF])
  #res_VIF <- mean(vif(modelVIF))
  RMSE <- sqrt(mean((Summ[1:p,1] - beta_true )^2))
  Errors <- Misc_quantities(seuil,Summ,beta_true)
  Misc <- Errors$Misc
  Sensi <- Errors$Sensi
  Speci <- Errors$Speci
  #Misc <- 1-mean((beta_true == 0) == (abs(Summ[1:p,1]) < seuil))
  #TP <- sum((beta_true != 0) & (abs(Summ[1:p,1]) > seuil))           # 'vrai' & 'estim'
  #TN <- sum((beta_true == 0) & (abs(Summ[1:p,1]) < seuil))
  #FP <- sum((beta_true == 0) & (abs(Summ[1:p,1]) > seuil))
  #FN <- sum((beta_true != 0) & (abs(Summ[1:p,1]) < seuil))
  #Sensi <- TP/(TP+FN)
  #Speci <- TN/(TN+FP)
  list_crps <- numeric(p)
  for(l in 1:p) { list_crps[l] <- CRPS_custom(beta_true[l],Samples[,l],length)}
  res_crps <- mean(list_crps)
  res_MIS <- mean(IS(alpha,Summ[1:p,4],Summ[1:p,5],beta_true))
  result <- c(RMSE,res_crps,Misc,Sensi,Speci,res_MIS,ESS_per_sec#,res_VIF)
  )
  return(result)
}   # setData et runMCMC compris
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
  seuil = 1          # méthode de sélection (seuil = 0 pour la méthode avec IC)
  nbModels = 6 + 2*(p <= n)
  nbCrit = 7
  multivariateNodesAsScalars = TRUE  # tous les beta sont obtenus par une normale à 1 dim et non multivariée
                                     # pour les modèles utilisant dmnorm
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
  
  # paramètres sur les priors
  df = 10   
  nu0 = 1          
  nu1 = 10^4
  g = n
  a = 1      # hyperparamètres sur Beta dans les S&S
  b = p
  r = 1      # hyperparamètres sur Gamma dans Laplace
  s = 0.1
}

# Pour les priors avec plusieurs versions :

HSvers <- 2
HSpvers <- 2

ConstNorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)
Conststud<- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,nu0=nu0,nu1=nu1,mu=mu,df=df,sigma=sigma,gamma=gamma,t=t)
Constgslab <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t,g=g)
Constgprior <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,g=g)
ConstHS <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t)
Constlapl <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,s=s,r=r)

Inits <- vector('list',nchains)
Inits[[1]] <- list(beta = rep(0,p))
Inits[[2]] <- list(beta = rep(10,p))

### Calcul en parallèle (tous les priors, données pré-générées) : ----------------------------------

library(doParallel)
ncores <- detectCores() - 1        # nombre de coeurs utilisés
#ncores <- 5
#cluster <- makeCluster(ncores)
cluster <- makeCluster(ncores,outfile = 'output_allprior.txt')   # si on veut suivre l'avancement des itérations
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
  All_y[,,z] <- Gener_Data(All_V[,,z],beta_true,sigma,gamma)
}

# sauvegarde des données générées 
if(save){
Covar_NamefileRDS <- paste('Covar_AllPriors_n=',n,'_p=',p,'_ds=',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
Y_NamefileRDS <- paste('Y_data_AllPriors_n=',n,'_p=',p,'_ds=',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
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
    library(car)
    library(verification)
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
    buildDirac <- buildMCMC(confDirac)     # conjugate sampler for beta
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
    if(p <= n){
      SpikeSlab_gslab <- nimbleCode({
        alpha ~ dbeta(a,b)
        for(l in 1:p){delta[l] ~ dbern(alpha)}
        D_half[1:p,1:p] <- diag(sqrt(nu0*(1-delta[1:p]) + nu1*delta[1:p]))     # nu0 : spike, nu1 : slab
        Mean[1:p] <- rep(0,p)
        Rgslab[1:p,1:p] <- g*inverse( t(V[1:n,1:p]) %*% V[1:n,1:p])
        CovMatr[1:p,1:p] <- D_half[1:p,1:p] %*% Rgslab[1:p,1:p] %*% D_half[1:p,1:p]
        beta[1:p] ~ dmnorm(mean = Mean[1:p],cov = CovMatr[1:p,1:p])
        for(i in 1:n){
          phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
          for(j in 1:J){
            y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
          }}
      })
      Modelgslab <- nimbleModel(code = SpikeSlab_gslab,const = Constgslab, data = Data,inits = Inits)
      Cgslab <- compileNimble(Modelgslab)
      confgslab <- configureMCMC(model = Cgslab, monitors = 'beta',multivariateNodesAsScalars = multivariateNodesAsScalars)
      buildgslab <- buildMCMC(confgslab)                # RW sampler for beta
      Cmcmcgslab <- compileNimble(buildgslab)
      Gprior <- nimbleCode({
        Mean[1:p] <- rep(0,p)
        CovMatr[1:p,1:p] <- gamma*g*inverse( t(V[1:n,1:p]) %*% V[1:n,1:p])
        beta[1:p] ~ dmnorm(mean = Mean[1:p],cov = CovMatr[1:p,1:p])
        for(i in 1:n){
          phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
          for(j in 1:J){
            y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
          }}
      })
      Modelgprior <- nimbleModel(code = Gprior,const = Constgprior, data = Data,inits = Inits)
      Cgprior <- compileNimble(Modelgprior)
      confgprior <- configureMCMC(model = Cgprior, monitors = 'beta',multivariateNodesAsScalars = multivariateNodesAsScalars)
      buildgprior <- buildMCMC(confgprior)              # RW sampler for beta
      Cmcmcgprior <- compileNimble(buildgprior)
    }
    
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
      })
    } else if(HSpvers == 4){
      PriorHorseshoePlus <- nimbleCode({
        tau ~ dunif(0,1)
        for(l in 1:p){ 
          norm_eta[l] ~ dnorm(0,var = 1)
          IG_eta[l] ~ dinvgamma(1/2,1/2)
          eta[l] <- abs(norm_eta[l]*sqrt(IG_eta[l]))
          norm_lambda[l] ~ dnorm(0,sd = tau*eta[l])
          IG_lambda[l] ~ dinvgamma(1/2,1/2)
          lambda[l] <- abs(norm_lambda[l]*sqrt(IG_lambda[l]))
          beta[l] ~ dnorm(0,sd = lambda[l])}
        for(i in 1:n){
          phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
          for(j in 1:J){
            y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
          }} 
      }) 
    }
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
    Result <- array(dim = c(nbModels,nbCrit,nechant))
    
    for(z in 1:nechant){
      print(paste0("Itération n°",z))
      V <- Covar[,,z]
      y_data <- Y[,,z]
      Data <- list(y=y_data,V=V)
      res_norm <- ComputeModel(Cnorm,Cmcmcnorm,Data)
      res_stud <- ComputeModel(Cstud,Cmcmcstud,Data)
      res_Dirac <- ComputeModel(CDirac,CmcmcDirac,Data)
      if(p <= n){
        res_gslab <- ComputeModel(Cgslab,Cmcmcgslab,Data)
        res_gprior <- ComputeModel(Cgprior,Cmcmcgprior,Data)}
      res_HS <- ComputeModel(CHS,CmcmcHS,Data)
      res_HSp <- ComputeModel(CHSp,CmcmcHSp,Data)
      res_lapl <- ComputeModel(Clapl,Cmcmclapl,Data)
      if(p <= n){Result[,,z] <- rbind(res_norm,res_stud,res_gslab,res_gprior,res_HS,res_HSp,res_lapl)}
      else {Result[,,z] <- rbind(res_norm,res_stud,res_Dirac,res_HS,res_HSp,res_lapl)}
    }
    return(Result)
    
  }) # end cluster evalQ
  
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
    
    Result <- array(dim = c(nbModels,nbCrit,1))
    V <- Covar_rem[,,1]
    y_data <- Y_rem[,,1]
    Data <- list(y=y_data,V=V)
    res_norm <- ComputeModel(Cnorm,Cmcmcnorm,Data)
    res_stud <- ComputeModel(Cstud,Cmcmcstud,Data)
    res_Dirac <- ComputeModel(CDirac,CmcmcDirac,Data)
    if(p <= n){
      res_gslab <- ComputeModel(Cgslab,Cmcmcgslab,Data)
      res_gprior <- ComputeModel(Cgprior,Cmcmcgprior,Data)}
    res_HS <- ComputeModel(CHS,CmcmcHS,Data)
    res_HSp <- ComputeModel(CHSp,CmcmcHSp,Data)
    res_lapl <- ComputeModel(Clapl,Cmcmclapl,Data)
    if(p <= n){Result[,,1] <- rbind(res_norm,res_stud,res_gslab,res_gprior,res_HS,res_HSp,res_lapl)}
    else {Result[,,1] <- rbind(res_norm,res_stud,res_Dirac,res_HS,res_HSp,res_lapl)}
    return(Result)
  })}
  t_final = difftime(Sys.time(),t_init,units='secs')
  stopCluster(cluster)
}

t_final

Matr3D <- NULL
for(z in 1:ncores) {Matr3D <- abind(Matr3D,Out[[z]],along=3)} 
if(rem_echant != 0){
for(z in 1:rem_echant) {Matr3D <- abind(Matr3D,Out_rem[[z]],along=3)}}


Final <- apply(Matr3D,c(1,2),mean)       # Matrice des moyennes sur les datasets (ie. sur la 3ème dimension)

CritNames <- c('RMSE','CRPS','Misc. rate','Sensitivity','Specificity','MIS','ESS/s')

if(p <= n){
  ModelNames <- c('S&S Continuous Normal','S&S Continuous Student','S&S Dirac Normal','S&S Continuous g-slab','g-prior','Horseshoe','Horseshoe+','Laplace')
} else {ModelNames <- c('S&S Continuous Normal','S&S Continuous Student','S&S Dirac Normal','Horseshoe','Horseshoe+','Laplace')}

df = matrix(0,nrow = nbModels,ncol=nbCrit+1)
df <- data.frame(df)
colnames(df) <- c('Model',CritNames)

df[,1] <- ModelNames
df[,2] <- round(Final[,1],4)
df[,3] <- round(Final[,2],3)
df[,4] <- round(Final[,3],3)
df[,5] <- round(Final[,4],3)
df[,6] <- round(Final[,5],3)
df[,7] <- round(Final[,6],3)
df[,8] <- round(Final[,7],1)
#df[,11] <- round(Final[,10],4)

df
myt <- ttheme_default(
  core = list(fg_params=list(cex=1.0)), 
  colhead = list(fg_params=list(cex=1.0)),
  rowhead = list(fg_params=list(cex=1.0)))
grid.table(df,theme=myt,rows=NULL)


# Sauvegarde des résultats 
if(save){
Res_NamefileRDS <- paste('results_AllPriors_n=',n,'_p=',p,'_ds=',nbDataset,'_Corr=',rho_sigma,'.rds',sep='')
saveRDS(Matr3D, file = Res_NamefileRDS)   
Result <- readRDS(file = Res_NamefileRDS)}                                  


### Histogrammes : 

#pdf('name_pdf.pdf')  # à renommer
if(p <= n){par(mfrow=c(2,4))} else {par(mfrow=c(2,3))}
for(j in 1:nbCrit){
  for(i in 1:nbModels) {
    if(j==3){
      hist(Result[i,j,],main=CritNames[j],xlab=ModelNames[i],breaks=seq(0,1,length=20))}
    else if(j==4){
      hist(Result[i,j,],main=CritNames[j],xlab=ModelNames[i],breaks=seq(0,1,length=20))}
    else if(j == 5){
      hist(Result[i,j,],main=CritNames[j],xlab=ModelNames[i],breaks=seq(0,1,length=20))}
    else {hist(Result[i,j,],main=CritNames[j],xlab=ModelNames[i])}
  } 
  #plot.new()
}
#dev.off()
