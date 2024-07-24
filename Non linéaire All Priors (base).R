### Simulation par un modèle non linéaire avec tous les priors, version de base

rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")  # à modifier 
library(latex2exp)
library(nimble)
library(MASS)
library(coda)
library(gridExtra)

### Fonctions et paramètres : -------------------------

f <- nimbleFunction(
  run = function(phi=double(0),psi1=double(0),psi2=double(0),t=double(0)) {
    result <- psi1 / (1+exp(-(t-phi)/psi2))
    return(result)
    returnType(double(0))
  })

# nimble_ginv <- nimbleFunction( run = function(M=double(2),n=integer(0)){
#     decomp_M <- svd(M)
#     Sigma <- diag(decomp_M$d)
#     U <- decomp_M$u
#     V <- decomp_M$v
#     for(i in 1:n){if(Sigma[i,i] != 0){Sigma[i,i] <- Sigma[i,i]^(-1)}}
#     result <- V %*% t(Sigma) %*% t(U)
#     return(result)
#     returnType(double(2)) })

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
  return(list(Misc=Misc,Sensi=Sensi,Speci=Speci))}
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

Gelman_Rubin <- function(MCMC){
  list_mcmc_objects <- vector('list',nchains)
  for(z in 1:nchains) {list_mcmc_objects[[z]] <- mcmc(MCMC$samples[[z]])}
  list_mcmc <- mcmc.list(list_mcmc_objects)
  GelRub <- gelman.diag(list_mcmc)
  return(round(GelRub$psrf[,1],2))
}  # estimation ponctuelle du critère de Gelman-Rubin pour chaque coefficient


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
  rho_sigma = 0
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

# Génération du premier dataset
set.seed(1)
#set.seed(3)   # pour plot
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
  
  y_data <- x + Eps
}

# par(mfrow=c(1,1))
# plot(t,x[1,],type='l',ylab = "",lwd=2,cex.lab=1.5)
# title(ylab = TeX('$f(\\varphi_i,t)$'),line=2.2,cex.lab=1.5)
# points(t,y_data[1,],pch=3)
# for(i in 2:7){
#  lines(t,x[i,],col=i,lwd=2)
#  points(t,y_data[i,],col=i,pch=3)
#  }
# grid()
# legend(x=c(2200,2900),y=c(40,100),legend=c("Individu 1",".",".",".", "Individu n"),
#        col=c("black","white","white","white", 6), lty=1, cex=0.8,lwd=2)

# Pour les priors avec plusieurs versions :

HSvers <- 2
HSpvers <- 2

ConstNorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)
Conststud<- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,nu0=nu0,nu1=nu1,mu=mu,df=df,sigma=sigma,gamma=gamma,t=t)
Constgslab <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t,g=g)
Constgprior <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,g=g)
ConstHS <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t)
#ConstHSp <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,tau=tau)
Constlapl <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,mu=mu,sigma=sigma,gamma=gamma,t=t,s=s,r=r)

Inits <- vector('list',nchains)
Inits[[1]] <- list(beta = rep(0,p))
Inits[[2]] <- list(beta = rep(10,p))

Data <- list(y=y_data,V=V)
#Data <- list(y=All_y[,,1],V=All_V[,,1])

### Définition des priors : ------------------------------

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
#confnorm <- configureMCMC(model = Cnorm, monitors = 'beta',onlyRW = TRUE)
#confnorm <- configureMCMC(model = Cnorm, monitors = c('beta','delta'))
#confnorm$printSamplers()                  # conjugate_dnorm_dnorm_linear sampler for beta
buildnorm <- buildMCMC(confnorm)        
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
#confDirac$printSamplers()
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
#confstud$printSamplers()                  # RW sampler for beta
buildstud <- buildMCMC(confstud)        
Cmcmcstud <- compileNimble(buildstud)

if(p <= n){
  SpikeSlab_gslab <- nimbleCode({
    alpha ~ dbeta(a,b)
    for(l in 1:p){delta[l] ~ dbern(alpha)}
    D_half[1:p,1:p] <- diag(sqrt(nu0*(1-delta[1:p]) + nu1*delta[1:p]))     # nu0 : spike, nu1 : slab
    Mean[1:p] <- rep(0,p)
    #Rgslab[1:p,1:p] <- g*nimble_ginv( t(V[1:n,1:p]) %*% V[1:n,1:p],p)
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
  #confgslab$printSamplers()                    # RW sampler for beta
  buildgslab <- buildMCMC(confgslab)               
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
  confgprior$printSamplers()                   # RW sampler for beta
  buildgprior <- buildMCMC(confgprior)             
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
confHS <- configureMCMC(model = CHS,monitors = 'beta') 
#confHS <- configureMCMC(model = CHS,monitors = c('beta','tau')) 
#confHS <- configureMCMC(model = CHS,monitors = c('beta','tau2')) 
#confHS$printSamplers()                                # conjugate_dnorm_dnorm_linear sampler for beta
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
confHSp <- configureMCMC(model = CHSp,monitors = 'beta')  
#confHSp <- configureMCMC(model = CHSp,monitors = c('tau','beta'))  
#confHSp <- configureMCMC(model = CHSp,monitors = c('tau2','beta'))
#confHSp <- configureMCMC(model = CHSp,monitors = c('tau2','beta','tau0'))  
#confHSp$printSamplers()                             # conjugate_dnorm_dnorm_linear sampler for beta
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
#conflapl$printSamplers()                          # RW Sampler for EXP
Buildlapl <- buildMCMC(conflapl)
Cmcmclapl <- compileNimble(Buildlapl)


### Summary MCMC  : --------------------------

set.seed(1)

MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
MCMCnorm$summary$all.chains

MCMCDirac <- runMCMC(CmcmcDirac,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
MCMCDirac$summary$all.chains

MCMCstud <- runMCMC(Cmcmcstud,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
MCMCstud$summary$all.chains
if(p <= n){
  MCMCgslab <- runMCMC(Cmcmcgslab,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
  MCMCgslab$summary$all.chains
  MCMCgprior <- runMCMC(Cmcmcgprior,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
  MCMCgprior$summary$all.chains}
MCMCHS <- runMCMC(CmcmcHS,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
MCMCHS$summary$all.chains
MCMCHSp <- runMCMC(CmcmcHSp,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
MCMCHSp$summary$all.chains
MCMClapl <- runMCMC(Cmcmclapl,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
MCMClapl$summary$all.chains


### Chaînes de Markov :  ------------------------------------
Titles <- paste("beta",1:p,sep = "")

par(mfrow=c(1,1))
plot(MCMCnorm$samples$chain1[,1],type='l',main='',ylab='',xlab = 'Itérations',ylim = c(0,120))
lines(MCMCnorm$samples$chain2[,1],col='red')
title(ylab = TeX('$\\beta_1$'),line = 2.2,cex.lab = 1.5)
grid()
legend('bottomright',legend = c('chain 1','chain 2'),col = c("black","red"),lty=1)
plot(MCMCnorm$samples$chain1[,4],type='l',main='',ylab='',xlab = 'Itérations',ylim = c(-10,10))
lines(MCMCnorm$samples$chain2[,4],col='red')
title(ylab = TeX('$\\beta_4$'),line = 2.2,cex.lab = 1.5)
grid()
legend('bottomright',legend = c('chain 1','chain 2'),col = c("black","red"),lty=1)


par(mfrow=c(3,3))
for(l in 1:p){
  plot(MCMCnorm$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Normal continuous S&S')
  lines(MCMCnorm$samples$chain2[,l],col='red')
}
for(l in 1:p){
  plot(MCMCDirac$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Normal Dirac S&S')
  lines(MCMCDirac$samples$chain2[,l],col='red')
}
for(l in 1:p){
  plot(MCMCstud$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Stud S&S')
  lines(MCMCstud$samples$chain2[,l],col='red')
}
if(p <= n){
for(l in 1:p){
  plot(MCMCgslab$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'G-slab S&S')
  lines(MCMCgslab$samples$chain2[,l],col='red')
}
for(l in 1:p){
  plot(MCMCgprior$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'G-prior')
  lines(MCMCgprior$samples$chain2[,l],col='red')
}}
for(l in 1:p){
  plot(MCMCHS$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Horseshoe')
  lines(MCMCHS$samples$chain2[,l],col='red')
}
for(l in 1:p){
  plot(MCMCHSp$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Horseshoe+')
  lines(MCMCHSp$samples$chain2[,l],col='red')
}
for(l in 1:p){
  plot(MCMClapl$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Laplace')
  lines(MCMClapl$samples$chain2[,l],col='red')
}

### Distributions a posteriori : ---------------------------------


par(mfrow=c(3,3))
for(l in 1:p){
  plot(density(MCMCnorm$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'Normal continuous S&S')
  lines(density(MCMCnorm$samples$chain2[,l]),col='red')
  grid()
}
for(l in 1:p){
  plot(density(MCMCDirac$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'Normal Dirac S&S')
  lines(density(MCMCDirac$samples$chain2[,l]),col='red')
  grid()
}
for(l in 1:p){
  plot(density(MCMCstud$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'Stud S&S')
  lines(density(MCMCstud$samples$chain2[,l]),col='red')
  grid()
}
if(p <= n){
  for(l in 1:p){
    plot(density(MCMCgslab$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'G-slab S&S')
    lines(density(MCMCgslab$samples$chain2[,l]),col='red')
    grid()
  }
  for(l in 1:p){
    plot(density(MCMCgprior$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'G-prior')
    lines(density(MCMCgprior$samples$chain2[,l]),col='red')
    grid()
  }}
for(l in 1:p){
  plot(density(MCMCHS$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'Horseshoe')
  lines(density(MCMCHS$samples$chain2[,l]),col='red')
  grid()
}
for(l in 1:p){
  plot(density(MCMCHSp$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'Horseshoe+')
  lines(density(MCMCHSp$samples$chain2[,l]),col='red')
  grid()
}
for(l in 1:p){
  plot(density(MCMClapl$samples$chain1[,l]),type='l',main=Titles[l],ylab="",xlab = 'Laplace')
  lines(density(MCMClapl$samples$chain2[,l]),col='red')
  grid()
}

### Critère de Gelman-Rubin : ------------------------------

max(Gelman_Rubin(MCMCnorm))
sum(Gelman_Rubin(MCMCnorm) > 1.05)
max(Gelman_Rubin(MCMCDirac))        # critère mal adapté pour le Dirac car la convergence est souvent bonne (visuellement)
sum(Gelman_Rubin(MCMCDirac) > 1.05)
max(Gelman_Rubin(MCMCstud))
sum(Gelman_Rubin(MCMCstud) > 1.05)
max(Gelman_Rubin(MCMCHS))
sum(Gelman_Rubin(MCMCHS) > 1.05)
max(Gelman_Rubin(MCMCHSp))
sum(Gelman_Rubin(MCMCHSp) > 1.05)
max(Gelman_Rubin(MCMClapl))
sum(Gelman_Rubin(MCMClapl) > 1.05)


### Comparaisons sur un dataset : ------------------------------

set.seed(1)
res_norm <- ComputeModel(Cnorm,Cmcmcnorm,Data)
res_Dirac <- ComputeModel(CDirac,CmcmcDirac,Data)
res_stud <- ComputeModel(Cstud,Cmcmcstud,Data)
if(p <= n){
  res_gslab <- ComputeModel(Cgslab,Cmcmcgslab,Data)
  res_gprior <- ComputeModel(Cgprior,Cmcmcgprior,Data)}
res_HS <- ComputeModel(CHS,CmcmcHS,Data)
res_HSp <- ComputeModel(CHSp,CmcmcHSp,Data)
res_lapl <- ComputeModel(Clapl,Cmcmclapl,Data)

CritNames <- c('Model','RMSE','CRPS','Misc. rate','Sensitivity','Specificity','MIS','ESS/s')

if(p <= n){
  ModelNames <- c('S&S Continuous Normal','S&S Continuous Student','S&S Dirac Normal','S&S Continuous g-slab','g-prior','Horseshoe','Horseshoe+','Laplace')
  Final <- rbind(res_norm,res_stud,res_Dirac,res_gslab,res_gprior,res_HS,res_HSp,res_lapl)
} else {ModelNames <- c('S&S Continuous Normal','S&S Continuous Student','S&S Dirac Normal','Horseshoe','Horseshoe+','Laplace')
Final <- rbind(res_norm,res_stud,res_Dirac,res_HS,res_HSp,res_lapl)
}

df = matrix(0,nrow = nbModels,ncol=nbCrit+1)
df <- data.frame(df)
colnames(df) <- CritNames

df[,1] <- ModelNames
df[,2] <- round(Final[,1],4)
df[,3] <- round(Final[,2],3)
df[,4] <- round(Final[,3],3)
df[,5] <- round(Final[,4],3)
df[,6] <- round(Final[,5],3)
df[,7] <- round(Final[,6],3)
df[,8] <- round(Final[,7],1)
df
par(mfrow=c(1,1))
plot.new()
myt <- ttheme_default(
  core = list(fg_params=list(cex=1.0)),
  colhead = list(fg_params=list(cex=1.0)),
  rowhead = list(fg_params=list(cex=1.0)))
grid.table(df,theme=myt,rows=NULL)


#pacf(MCMCnorm$samples$chain1[,1])
#pacf(MCMCnorm$samples$chain2[,1])
#pacf(MCMCDirac$samples$chain1[,1])
#pacf(MCMCDirac$samples$chain2[,1])

