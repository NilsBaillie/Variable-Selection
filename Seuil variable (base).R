##### Seuil variable, version de base

# On veut déterminer comment évolue le taux de mauvais classement pour chacun des priors
# lorsqu'on fait varier la valeur seuil
# On affiche les courbes ROC et PR (Precision-Recall)

rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")  # à modifier 
library(latex2exp)
library(nimble)
library(MASS)
library(abind)


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
  return(list(Misc=Misc,Sensi=Sensi,Speci=Speci,Preci=Preci))
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
  
  # paramètres sur les priors
  df = 10                    
  nu0 = 1          
  nu1 = 10^4
  a = 1      # hyperparamètres sur Beta dans les S&S
  b = p
  r = 1      # hyperparamètres sur Gamma dans Laplace
  s = 0.1
}

# Génération du dataset
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
  y_data <- x + Eps
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
Data <- list(y=y_data,V=V)


# Définition des priors : ------------------------

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
}

HorseshoePlus <- nimbleModel(code = PriorHorseshoePlus,constants = ConstHS,data = Data,inits = Inits)
CHSp <- compileNimble(HorseshoePlus)
confHSp <- configureMCMC(model = CHSp,monitors = 'beta')  
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

# Échantillons a posteriori : ---------------------------------
set.seed(1)
MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
Summnorm <- MCMCnorm$summary$all.chains
plot(1:p,Summnorm[1:p,1],xlab = 'indice',ylab = "",
     col=ifelse(abs(Summnorm[1:p,1]) < seuil,"blue","red"),
     #ylim=c(-1,20),
     pch=16)
title(ylab = TeX('$\\beta$'),line=2.2,cex.lab=1.5)
grid()
legend("topright",legend = c("Covariables significatives","Covariables négligeables"),col=c('red','blue'),pch=16)
MCMCDirac <- runMCMC(CmcmcDirac,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
SummDirac <- MCMCDirac$summary$all.chains
MCMCstud <- runMCMC(Cmcmcstud,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
Summstud <- MCMCstud$summary$all.chains
MCMCHS <- runMCMC(CmcmcHS,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
SummHS <- MCMCHS$summary$all.chains
MCMCHSp <- runMCMC(CmcmcHSp,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
SummHSp <- MCMCHSp$summary$all.chains
MCMClapl <- runMCMC(Cmcmclapl,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
Summlapl <- MCMClapl$summary$all.chains


seuil <- seq(0.01,2,by=0.01)
nb_seuil <- length(seuil)

Misc_rate_norm <- array(0,nb_seuil)
Misc_rate_stud <- array(0,nb_seuil)
Misc_rate_Dirac <- array(0,nb_seuil)
Misc_rate_HS <- array(0,nb_seuil)
Misc_rate_HSp <- array(0,nb_seuil)
Misc_rate_lapl <- array(0,nb_seuil)

for(k in 1:nb_seuil){
  Misc_rate_norm[k] <- Misc_quantities(seuil[k],Summnorm,beta_true)$Misc
  Misc_rate_stud[k] <- Misc_quantities(seuil[k],Summstud,beta_true)$Misc
  Misc_rate_Dirac[k] <- Misc_quantities(seuil[k],SummDirac,beta_true)$Misc
  Misc_rate_HS[k] <- Misc_quantities(seuil[k],SummHS,beta_true)$Misc
  Misc_rate_HSp[k] <- Misc_quantities(seuil[k],SummHSp,beta_true)$Misc
  Misc_rate_lapl[k] <- Misc_quantities(seuil[k],Summlapl,beta_true)$Misc
}

### Plot : ------------------------------

plot(seuil,Misc_rate_norm,type='l',col='blue',ylab = 'Misc. rate',main="Taux d'erreur en fonction du seuil",
     #xlim = c(1,2),ylim=c(0,0.05),
     #xlim = c(0,0.5),
     lwd=2)
lines(seuil,Misc_rate_stud,col='red',lwd=2)
lines(seuil,Misc_rate_Dirac,col='magenta',lwd=2)
lines(seuil,Misc_rate_HS,col='green',lwd=2)
lines(seuil,Misc_rate_HSp,col='orange',lwd=2)
lines(seuil,Misc_rate_lapl,col='purple',lwd=2)
grid()
legend('topright',legend = c('S&S Cont. Norm','S&S Cont. Stud','S&S Dirac Norm','HS','HS+','Laplace'),
       col = c('blue','red','magenta','green','orange','purple'),lty=1,lwd=2)



### Courbes ROC : -----------------------

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
  Sensi_norm[k] <- Misc_quantities(seuil[k],Summnorm,beta_true)$Sensi
  Sensi_stud[k] <- Misc_quantities(seuil[k],Summstud,beta_true)$Sensi
  Sensi_Dirac[k] <- Misc_quantities(seuil[k],SummDirac,beta_true)$Sensi
  Sensi_HS[k] <- Misc_quantities(seuil[k],SummHS,beta_true)$Sensi
  Sensi_HSp[k] <- Misc_quantities(seuil[k],SummHSp,beta_true)$Sensi
  Sensi_lapl[k] <- Misc_quantities(seuil[k],Summlapl,beta_true)$Sensi
  Speci_norm[k] <- Misc_quantities(seuil[k],Summnorm,beta_true)$Speci
  Speci_stud[k] <- Misc_quantities(seuil[k],Summstud,beta_true)$Speci
  Speci_Dirac[k] <- Misc_quantities(seuil[k],SummDirac,beta_true)$Speci
  Speci_HS[k] <- Misc_quantities(seuil[k],SummHS,beta_true)$Speci
  Speci_HSp[k] <- Misc_quantities(seuil[k],SummHSp,beta_true)$Speci
  Speci_lapl[k] <- Misc_quantities(seuil[k],Summlapl,beta_true)$Speci
  Preci_norm[k] <- Misc_quantities(seuil[k],Summnorm,beta_true)$Preci
  Preci_stud[k] <- Misc_quantities(seuil[k],Summstud,beta_true)$Preci
  Preci_Dirac[k] <- Misc_quantities(seuil[k],SummDirac,beta_true)$Preci
  Preci_HS[k] <- Misc_quantities(seuil[k],SummHS,beta_true)$Preci
  Preci_HSp[k] <- Misc_quantities(seuil[k],SummHSp,beta_true)$Preci
  Preci_lapl[k] <- Misc_quantities(seuil[k],Summlapl,beta_true)$Preci
}

plot(1-Speci_norm,Sensi_norm,type='l',col='blue',main='Courbes ROC',lwd=2,xlim = c(0,1),ylim=c(0,1))
lines(1-Speci_stud,Sensi_stud,col='red',lwd=2)
lines(1-Speci_Dirac,Sensi_Dirac,col='magenta',lwd=2)
lines(1-Speci_HS,Sensi_HS,col='green',lwd=2)
lines(1-Speci_HSp,Sensi_HSp,col='orange',lwd=2)
lines(1-Speci_lapl,Sensi_lapl,col='black',lwd=2)
lines(seq(0,1,length=nb_seuil),seq(0,1,length=nb_seuil),col='gray',lwd=2)
lines(seq(0,1,length=nb_seuil),rep(1,nb_seuil),col='black',lwd=1,lty=2)
lines(rep(0,nb_seuil),seq(0,1,length=nb_seuil),col='black',lwd=1,lty=2)
grid()

### Courbes PR : -----------------------

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

### Aires sous la courbe
compute_auc <- function(x,y){
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
