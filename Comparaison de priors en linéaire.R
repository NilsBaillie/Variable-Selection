### Comparaison de priors en linéaire

rm(list=objects())
graphics.off()
library(latex2exp)
library(MASS)
library(nimble)
library(coda)
library(car)
library(verification)
library(gridExtra)
library(foreach)
library(doParallel)
library(abind)
library(BoomSpikeSlab)

### Paramètres utilisés : ----------------------------------------------

n = 100            # nombre d'observations
p = 15             # nombre de variables
nbzero = 5         # nombre de paramètres imposés à zéro
q = 15             # nombre de jeux de données générés
nchains = 2        # nombre de chaînes de Markov lancées
niter = 10000      # nombre d'itérations dans le MCMC
nburnin = 5000     # Temps de chauffe du MCMC
seuil = 0.1        # Seuil en-dessous duquel le theta est estimé à zéro
theta = rep(0,nbzero) 
theta = c(theta,seq(0.1,2,length=p-nbzero))


IS <- function(alpha,lower,upper,true_val)   # fonction auxiliaire pour le calcul du MIS
{return(((upper-lower) + (2/alpha)*(lower-true_val)*(true_val < lower) + (2/alpha)*(true_val-upper)*(true_val > upper)))}
alpha <- 0.05
situation <- 'independant'
#situation <- 'few correlations'
#situation <- 'many correlations'
#situation <- 'multicollinearity'

nbModels = 5 + (situation != 'multicollinearity')
nbCrit = 10
mu = 0               # intercept de la régression linéaire (coeff du terme constant)
sigmaX = 1           # facteur de variance sur les données générées dans X
sigma = 1            # écart-type du bruit eps (pas toujours le même)

-------------------------------
  
### Première version : pas de précompilation
### Pas à jour ! Utiliser la 2ème version !---------------------------------------
##### Compute MCMC function : 
ComputeMCMC = function(){
  
theta = rep(0,nbzero) 
theta = c(theta,seq(0.1,2,length=p-nbzero))  
X = matrix(0,nrow=n,ncol=p)
if(situation == 'independant'){
  for(j in 1:p)
  {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # variables hétéroscédastiques, normales et indépendantes
}

if(situation == 'few correlations'){
  for(j in 1:(p-2))
  {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # 2 dernières variables fonctions des 2 premières
  X[,p-1] <- X[,1]^2
  X[,p] <- sin(X[,2])
}

if(situation == 'multicollinearity')
{
  for(j in 1:(p-1))
  {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}    # p_ème variable = somme de toutes les autres
  X[,p] <- sum(X[,1:(p-1)])
}

eps = rnorm(n,0,sd = sigma)
Y_data = c(mu*rep(1,n) + X %*% theta + eps)
X <- as.data.frame(X)

### Spike & Slab : -------------------------------------------------

### Spike continu normal

SpikeNorm <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2) 
  prob ~ dbeta(a,b)       # probabilité d'acceptation a posteriori
  for(j in 1:p)
  {
    delta[j] ~ dbern(prob)
    theta[j] ~ dnorm(mean = 0,sd = sqrt(V*((1-r)*(delta[j] == 1)+r)))  # Spike & Slab
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)}   # rajouter [1:p] à theta pour avoir la liste
})

ConstN <- list(n=n,p=p,mu=mu,a=1,b=1,V=1,r=10^(-4),X=X) 
Data <- list(Y=Y_data)
Inits <- list(sigma=1,theta = rep(1,p),delta = rep(0,p))
start_time = Sys.time()
SpikeN <- nimbleModel(code = SpikeNorm,constants = ConstN,data = Data,inits = Inits)
mcmcN <- nimbleMCMC(model = SpikeN,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,monitors = c('theta','delta'))
TimeN <- difftime(Sys.time(),start_time,units='secs')
mcmcN$summary


### Spike continu Student

SpikeStud <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2)  
  prob ~ dbeta(a,b)                    # probabilité d'acceptation a posteriori
  for(j in 1:p)
  {
    delta[j] ~ dbern(prob)
    theta[j] ~ dt(mu = 0,sigma = sqrt((Q/nu)*((1-r)*(delta[j] == 1)+r)),df = 2*nu)  # Spike & Slab
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)}   # rajouter [1:p] à theta pour avoir la liste
})

ConstT <- list(n=n,p=p,mu=mu,a=1,b=1,Q=4,nu=5,r=10^(-4),X=X) 
Data <- list(Y=Y_data)
Inits <- list(sigma=1,theta = rep(1,p))
start_time = Sys.time()
SpikeT <- nimbleModel(code = SpikeStud,constants = ConstT,data = Data,inits = Inits)
mcmcT <- nimbleMCMC(model = SpikeT,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,monitors = c('theta','delta'))
TimeT <- difftime(Sys.time(),start_time,units='secs')
#mcmcT$summary


### Spike de Dirac (i-slab)

SpikeDirac <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2)
  prob ~ dbeta(a,b)           # probabilité d'acceptation a posteriori
  Mean[1:p] <- matrix(rep(0,p))
  VarMatr[1:p,1:p] <- sigma^2*c*diag(p)
  theta0[1:p] ~ dmnorm(mean = Mean[1:p], cov = VarMatr[1:p,1:p])
  for(j in 1:p)   
  {
    delta[j] ~ dbern(prob)
    theta[j] <- delta[j]*theta0[j]
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)} 
})

ConstD <- list(n=n,p=p,mu=mu,a=1,b=1,c=1,X=X) 
Data <- list(Y=Y_data)
Inits <- list(sigma=1,theta = rep(1,p))
start_time = Sys.time()
SpikeD <- nimbleModel(code = SpikeDirac,constants = ConstD,data = Data,inits = Inits)
mcmcD <- nimbleMCMC(model = SpikeD,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,monitors = c('theta','delta'))
TimeD <- difftime(Sys.time(),start_time,units='secs')
#mcmcD$summary 

if(situation != 'multicollinearity'){
### Spike de Dirac (g-slab)

SpikeDiracG <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2)
  prob ~ dbeta(a,b)           # probabilité d'acceptation a posteriori
  Mean[1:p] <- matrix(rep(0,p))
  VarMatr[1:p,1:p] <- sigma^2*g*inverse(t(X[,]) %*% X[,])   # g-prior
  theta0[1:p] ~ dmnorm(mean = Mean[1:p], cov = VarMatr[1:p,1:p])
  for(j in 1:p)   
  {
    delta[j] ~ dbern(prob)
    theta[j] <- delta[j]*theta0[j]
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)} 
})

ConstDg <- list(n=n,p=p,mu=mu,a=1,b=1,g=n) 
Data <- list(Y=Y_data,X=X)
Inits <- list(sigma=1,theta = rep(1,p))
start_time = Sys.time()
SpikeDg <- nimbleModel(code = SpikeDiracG,constants = ConstDg,data = Data,inits = Inits)
mcmcDg <- nimbleMCMC(model = SpikeDg,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,monitors = c('theta','delta'))
TimeDg <- difftime(Sys.time(),start_time,units='secs')
#mcmcDg$summary
}
# Donne des résultats bons, mais ne devrait pas car toute la matrice X est donnée dans dmnorm,
# alors que seules certaines colonnes devraient l'être (restriction aux indices j tel que delta_j = 1)
# Modèle non utilisable si X n'est pas de rang plein


### Horseshoe : -------------------------------------------------------

PriorHorseshoe <- nimbleCode({
  for(j in 1:p)
  {
    lambda[j] ~ T(dt(mu = 0, tau = 1, df = 1),0,10^300)   # Cauchy(0,1) = Student(df = 1), de plus tronquée sur R+ (approx)
    theta[j] ~ dnorm(0, sd = lambda[j]*tau)               # Horseshoe prior
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)}
})

ConstHS <- list(n=n,p=p,mu=mu,tau=1,sigma=sigma,X=X) 
Data <- list(Y=Y_data)
Inits <- list(theta = rep(1,p))
start_time = Sys.time()
Horseshoe <- nimbleModel(code = PriorHorseshoe,constants = ConstHS,data = Data,inits = Inits)
mcmcHS <- nimbleMCMC(model=Horseshoe,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,monitors = 'theta')
TimeHS <- difftime(Sys.time(),start_time,units='secs')
#mcmcHS$summary


### Critères de comparaison : ------------------------------------------

# Temps d'exécution :

if(situation != 'multicollinearity'){res_Time <- c(TimeN,TimeT,TimeD,TimeDg,TimeHS)
} else {res_Time <- c(TimeN,TimeT,TimeD,TimeHS)}

# Effective Sample Size (ESS, indépendance des tirages) :

essN = mean(effectiveSize(mcmc.list(mcmc(mcmcN$samples$chain1[,(p+1):(2*p)]),mcmc(mcmcN$samples$chain2[,(p+1):(2*p)]))))
essT = mean(effectiveSize(mcmc.list(mcmc(mcmcT$samples$chain1[,(p+1):(2*p)]),mcmc(mcmcT$samples$chain2[,(p+1):(2*p)]))))
essD = mean(effectiveSize(mcmc.list(mcmc(mcmcD$samples$chain1[,(p+1):(2*p)]),mcmc(mcmcD$samples$chain2[,(p+1):(2*p)]))))
essHS = mean(effectiveSize(mcmc.list(mcmc(mcmcHS$samples$chain1[,1:p]),mcmc(mcmcHS$samples$chain2[,1:p]))))

if(situation != 'multicollinearity')
{
  essDg = mean(effectiveSize(mcmc.list(mcmc(mcmcDg$samples$chain1[,(p+1):(2*p)]),mcmc(mcmcDg$samples$chain2[,(p+1):(2*p)]))))
  res_ESS = c(essN,essT,essD,essDg,essHS)  
} else{res_ESS = c(essN,essT,essD,essHS) }


# ESS/s (compromis vitesse/mixing) :

ESS_per_sec = res_ESS/as.numeric(res_Time)


# Variance Inflation Factor (VIF, Séparation des variables) :

ind_N = which(mcmcN$summary$all.chains[-(1:p),1] < seuil)
formula = parse(text = paste("V",ind_N[1],sep="","~."))[[1]]
modelN <- lm(formula = formula, data = X[,ind_N])
#vif(modelN)

ind_T = which(mcmcT$summary$all.chains[-(1:p),1] < seuil)
formula = parse(text = paste("V",ind_T[1],sep="","~."))[[1]]
modelT <- lm(formula = formula, data = X[,ind_T])
#vif(modelT)

ind_D = which(mcmcD$summary$all.chains[-(1:p),1] < seuil)
formula = parse(text = paste("V",ind_D[1],sep="","~."))[[1]]
modelD <- lm(formula = formula, data = X[,ind_D])
#vif(modelD)

ind_HS = which(mcmcHS$summary$all.chains[,1] < seuil)
formula = parse(text = paste("V",ind_HS[1],sep="","~."))[[1]]
modelHS <- lm(formula = formula, data = X[,ind_HS])
#vif(modelHS)

if(situation != 'multicollinearity')
{
  ind_Dg = which(mcmcDg$summary$all.chains[-(1:p),1] < seuil)
  formula = parse(text = paste("V",ind_Dg[1],sep="","~."))[[1]]
  modelDg <- lm(formula = formula, data = X[,ind_Dg])
  vif(modelDg)
  res_VIF <- c(mean(vif(modelN)),mean(vif(modelT)),mean(vif(modelD)),mean(vif(modelDg)),mean(vif(modelHS)))
} else {res_VIF <- c(mean(vif(modelN)),mean(vif(modelT)),mean(vif(modelD)),mean(vif(modelHS)))}


# Root Mean Squared Error (RMSE, Estimation des paramètres) :

RMSE_N = sqrt(mean((mcmcN$summary$all.chains[(p+1):(2*p),1] - theta )^2))
RMSE_T = sqrt(mean((mcmcT$summary$all.chains[(p+1):(2*p),1] - theta )^2))
RMSE_D = sqrt(mean((mcmcD$summary$all.chains[(p+1):(2*p),1] - theta )^2))
RMSE_HS = sqrt(mean((mcmcHS$summary$all.chains[1:p,1] - theta )^2))

if(situation != 'multicollinearity'){
  RMSE_Dg = sqrt(mean((mcmcDg$summary$all.chains[(p+1):(2*p),1] - theta )^2))
  res_RMSE <- c(RMSE_N,RMSE_T,RMSE_D,RMSE_Dg,RMSE_HS)
} else {res_RMSE <- c(RMSE_N,RMSE_T,RMSE_D,RMSE_HS)}


# Taux de mauvais classement (Misclassification rate) : 

Misc_N = 1-mean((theta == 0) == (mcmcN$summary$all.chains[(p+1):(2*p),1] < seuil))
Misc_T = 1-mean((theta == 0) == (mcmcT$summary$all.chains[(p+1):(2*p),1] < seuil))
Misc_D = 1-mean((theta == 0) == (mcmcD$summary$all.chains[(p+1):(2*p),1] < seuil))
Misc_HS = 1-mean((theta == 0) == (mcmcHS$summary$all.chains[1:p,1] < seuil))

if(situation != 'multicollinearity'){
  Misc_Dg = 1-mean((theta == 0) == (mcmcDg$summary$all.chains[(p+1):(2*p),1] < seuil))
  res_Misc <- c(Misc_N,Misc_T,Misc_D,Misc_Dg,Misc_HS)
} else {res_Misc <- c(Misc_N,Misc_T,Misc_D,Misc_HS)}

# Sensibilité (= VP/(VP+FN)) et Spécificité (= VN/(VN+FP)) :

tabN = table(vrai = (theta != 0), estim = (mcmcN$summary$all.chains[(p+1):(2*p),1] > seuil))
Sensi_N = tabN[2,2]/(tabN[2,2]+tabN[2,1])
Speci_N = tabN[1,1]/(tabN[1,1]+tabN[1,2])

tabT = table(vrai = (theta != 0), estim = (mcmcT$summary$all.chains[(p+1):(2*p),1] > seuil))
Sensi_T = tabT[2,2]/(tabT[2,2]+tabT[2,1])
Speci_T = tabT[1,1]/(tabT[1,1]+tabT[1,2])

tabD = table(vrai = (theta != 0), estim = (mcmcD$summary$all.chains[(p+1):(2*p),1] > seuil))
Sensi_D = tabD[2,2]/(tabD[2,2]+tabD[2,1])
Speci_D = tabD[1,1]/(tabD[1,1]+tabD[1,2])

tabHS = table(vrai = (theta != 0), estim = (mcmcHS$summary$all.chains[1:p,1] > seuil))
Sensi_HS = tabHS[2,2]/(tabHS[2,2]+tabHS[2,1])
Speci_HS = tabHS[1,1]/(tabHS[1,1]+tabHS[1,2])

if(situation != 'multicollinearity'){
  tabDg = table(vrai = (theta != 0), estim = (mcmcDg$summary$all.chains[(p+1):(2*p),1] > seuil))
  Sensi_Dg = tabDg[2,2]/(tabDg[2,2]+tabDg[2,1])
  Speci_Dg = tabDg[1,1]/(tabDg[1,1]+tabDg[1,2])
  res_Sensi <- c(Sensi_N,Sensi_T,Sensi_D,Sensi_Dg,Sensi_HS)
  res_Speci <- c(Speci_N,Speci_T,Speci_D,Speci_Dg,Speci_HS)
} else {res_Sensi <- c(Sensi_N,Sensi_T,Sensi_D,Sensi_HS)
res_Speci <- c(Speci_N,Speci_T,Speci_D,Speci_HS)}


# Continuous Ranked Probability Score (CRPS) :

list_mean <- numeric(n)
for(i in 1:n){ list_mean[i] <- mu + sum(X[i,1:p]*mcmcN$samples$chain1[i,(p+1):(2*p)]) }
crpsN <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

list_mean <- numeric(n)
for(i in 1:n){ list_mean[i] <- mu + sum(X[i,1:p]*mcmcT$samples$chain1[i,(p+1):(2*p)]) }
crpsT <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

list_mean <- numeric(n)
for(i in 1:n){ list_mean[i] <- mu + sum(X[i,1:p]*mcmcD$samples$chain1[i,(p+1):(2*p)]) }
crpsD <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

list_mean <- numeric(n)
for(i in 1:n){ list_mean[i] <- mu + sum(X[i,1:p]*mcmcHS$samples$chain1[i,1:p]) }
crpsHS <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

if(situation != 'multicollinearity'){
  list_mean <- numeric(n)
  for(i in 1:n){ list_mean[i] <- mu + sum(X[i,1:p]*mcmcDg$samples$chain1[i,(p+1):(2*p)]) }
  crpsDg <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))
  res_crps <- c(crpsN$CRPS,crpsT$CRPS,crpsD$CRPS,crpsDg$CRPS,crpsHS$CRPS)
} else {res_crps <- c(crpsN$CRPS,crpsT$CRPS,crpsD$CRPS,crpsHS$CRPS)}


# Mean Interval Score (MIS)

IS_N <- IS(alpha,mcmcN$summary$all.chains[(p+1):(2*p),4],mcmcN$summary$all.chains[(p+1):(2*p),5],theta)
IS_T <- IS(alpha,mcmcT$summary$all.chains[(p+1):(2*p),4],mcmcT$summary$all.chains[(p+1):(2*p),5],theta)
IS_D <- IS(alpha,mcmcD$summary$all.chains[(p+1):(2*p),4],mcmcD$summary$all.chains[(p+1):(2*p),5],theta)
IS_HS <- IS(alpha,mcmcHS$summary$all.chains[1:p,4],mcmcHS$summary$all.chains[1:p,5],theta)

if(situation != 'multicollinearity')
{
  IS_Dg <- IS(alpha,mcmcDg$summary$all.chains[(p+1):(2*p),4],mcmcDg$summary$all.chains[(p+1):(2*p),5],theta)
  res_MIS <- c(mean(IS_N),mean(IS_T),mean(IS_D),mean(IS_Dg),mean(IS_HS))
} else {res_MIS <- c(mean(IS_N),mean(IS_T),mean(IS_D),mean(IS_HS))}

Matr <- cbind(res_RMSE,res_crps,res_Misc,res_Sensi,res_Speci,res_MIS,res_Time,res_ESS,ESS_per_sec,res_VIF)
return(Matr)
}

### Parallel foreach loop : 

totalCores <- detectCores()
cluster <- makeCluster(totalCores-1)
registerDoParallel(cluster)

t_init = Sys.time()
Output <- foreach(k=1:q, .packages = c('nimble','coda','car','verification','MASS')) %dopar%  {
  ComputeMCMC()
}

stopCluster(cluster)
t_final = difftime(Sys.time(),t_init,units='secs')
t_final
Matr3D <- simplify2array(Output)    # Liste de matrices -> Matrice 3D
Final <- apply(Matr3D,c(1,2),mean)  # Matrice des moyennes sur les datasets (ie. sur la 3ème dimension)





### Deuxième version : pré-compilée  -------------------------------------------------------------------
# Semble être incompatible avec le calcul en parallèle à cause de la façon dont la compilation C++ de NIMBLE
# gère les modèles MCMC, dommage car cette méthode évite de recompiler à chaque fois

### Génération de données et Précompilation des modèles MCMC ---------------------------------

set.seed(1)
X = matrix(0,nrow=n,ncol=p)
if(situation == 'independant'){
  for(j in 1:p)
  {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # variables hétéroscédastiques, normales et indépendantes
}

if(situation == 'few correlations'){
  for(j in 1:(p-2))
  {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # 2 dernières variables fonctions des 2 premières
  X[,p-1] <- X[,1]^2
  X[,p] <- sin(X[,2])
}
if(situation == 'multicollinearity')
{
  for(j in 1:(p-1))
  {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}    # p_ème variable = somme de toutes les autres
  X[,p] <- sum(X[,1:(p-1)])
}
eps = rnorm(n,0,sd = sigma)
Y_data = c(mu*rep(1,n) + X %*% theta + eps)
Df_X <- as.data.frame(X)

SpikeNorm <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2) 
  prob ~ dbeta(a,b)       # probabilité d'acceptation a posteriori
  for(j in 1:p)
  {
    delta[j] ~ dbern(prob)
    theta[j] ~ dnorm(mean = 0,sd = sqrt(V*((1-r)*(delta[j] == 1)+r)))  # Spike & Slab
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)}   # rajouter [1:p] à theta pour avoir la liste
})

ConstN <- list(n=n,p=p,mu=mu,a=1,b=1,V=1,r=10^(-4)) 
Data <- list(Y=Y_data,X=X)
Inits <- list(sigma=1,theta = rep(1,p))
SpikeN <- nimbleModel(code = SpikeNorm,constants = ConstN,data = Data,inits = Inits)
CspikeN <- compileNimble(SpikeN)
confN <- configureMCMC(model = CspikeN,monitors = c('theta'))
BuildN <- buildMCMC(confN)
CmcmcN <- compileNimble(BuildN)

SpikeStud <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2)  
  prob ~ dbeta(a,b)                    # probabilité d'acceptation a posteriori
  for(j in 1:p)
  {
    delta[j] ~ dbern(prob)
    theta[j] ~ dt(mu = 0,sigma = sqrt((Q/nu)*((1-r)*(delta[j] == 1)+r)),df = 2*nu)  # Spike & Slab
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)}   # rajouter [1:p] à theta pour avoir la liste
})

ConstT <- list(n=n,p=p,mu=mu,a=1,b=1,Q=4,nu=5,r=10^(-4)) 
Data <- list(Y=Y_data,X=X)
Inits <- list(sigma=1,theta = rep(1,p))
SpikeT <- nimbleModel(code = SpikeStud,constants = ConstT,data = Data,inits = Inits)
CspikeT <- compileNimble(SpikeT)
confT <- configureMCMC(model = CspikeT,monitors = c('theta'))
BuildT <- buildMCMC(confT)
CmcmcT <- compileNimble(BuildT)

SpikeDirac <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2, scale = 1/2)
  prob ~ dbeta(a,b)           # probabilité d'acceptation a posteriori
  Mean[1:p] <- matrix(rep(0,p))
  VarMatr[1:p,1:p] <- sigma^2*c*diag(p)
  theta0[1:p] ~ dmnorm(mean = Mean[1:p], cov = VarMatr[1:p,1:p])
  for(j in 1:p)   
  {
    delta[j] ~ dbern(prob)
    theta[j] <- delta[j]*theta0[j]
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)} 
})

ConstD <- list(n=n,p=p,mu=mu,a=1,b=1,c=1) 
Data <- list(Y=Y_data,X=X)
Inits <- list(sigma=1,theta = rep(1,p))
SpikeD <- nimbleModel(code = SpikeDirac,constants = ConstD,data = Data,inits = Inits)
CspikeD <- compileNimble(SpikeD)
confD <- configureMCMC(model = CspikeD,monitors = c('theta'))
BuildD <- buildMCMC(confD)
CmcmcD <- compileNimble(BuildD)

if(situation != 'multicollinearity'){
  
  SpikeDiracG <- nimbleCode({
    sigma ~ dinvgamma(shape = 1/2, scale = 1/2)
    prob ~ dbeta(a,b)           # probabilité d'acceptation a posteriori
    Mean[1:p] <- matrix(rep(0,p))
    VarMatr[1:p,1:p] <- sigma^2*g*inverse(t(X[,]) %*% X[,])   # g-prior
    theta0[1:p] ~ dmnorm(mean = Mean[1:p], cov = VarMatr[1:p,1:p])
    for(j in 1:p)   
    {
      delta[j] ~ dbern(prob)
      theta[j] <- delta[j]*theta0[j]
    }
    for(i in 1:n)
    {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)} 
  })
  
  ConstDg <- list(n=n,p=p,mu=mu,a=1,b=1,g=n) 
  Data <- list(Y=Y_data,X=X)
  Inits <- list(sigma=1,theta = rep(1,p))
  SpikeDg <- nimbleModel(code = SpikeDiracG,constants = ConstDg,data = Data,inits = Inits)
  CspikeDg <- compileNimble(SpikeDg)
  confDg <- configureMCMC(model = CspikeDg,monitors = c('theta'))
  BuildDg <- buildMCMC(confDg)
  CmcmcDg <- compileNimble(BuildDg)
}

PriorHorseshoe <- nimbleCode({
  sigma ~ dinvgamma(shape = 1/2,scale = 1/2)
  tau ~ T(dt(mu = 0, tau = 1/sigma^4, df = 1),0,10^300)
  for(j in 1:p)
  {
    lambda[j] ~ T(dt(mu = 0, tau = 1, df = 1),0,10^300)  
    theta[j] ~ dnorm(0, sd = lambda[j]*tau)               # Horseshoe prior
  }
  for(i in 1:n)
  {Y[i] ~ dnorm(mu + sum(X[i,1:p]*theta[1:p]),sd = sigma)}
}) 

ConstHS <- list(n=n,p=p,mu=mu) 
Data <- list(Y=Y_data,X=X)
Inits <- list(sigma=1,theta = rep(1,p))
Horseshoe <- nimbleModel(code = PriorHorseshoe,constants = ConstHS,data = Data,inits = Inits)
C_HS <- compileNimble(Horseshoe)
confHS <- configureMCMC(model = C_HS,monitors = c('theta'))
BuildHS <- buildMCMC(confHS)
CmcmcHS <- compileNimble(BuildHS)


### Calcul non parallèle du MCMC précompilé -------------------

ComputeMCMCbis <- function()
{
  X = matrix(0,nrow=n,ncol=p)
  if(situation == 'independant'){
    for(j in 1:p)
    {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # variables hétéroscédastiques, normales et indépendantes
  }
  if(situation == 'few correlations'){
    for(j in 1:(p-2))
    {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # 2 dernières variables fonctions des 2 premières
    X[,p-1] <- X[,1]^2
    X[,p] <- sin(X[,2])
  }
  if(situation == 'multicollinearity')
  {
    for(j in 1:(p-1))
    {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}    # p_ème variable = somme de toutes les autres
    X[,p] <- sum(X[,1:(p-1)])
  }
  eps = rnorm(n,0,sd = sigma)
  Y_data = c(mu*rep(1,n) + X %*% theta + eps)
  
  # Boom Spike and Slab package (independant normal, SSVS prior)
  start_time = Sys.time()
  priorBoomSS <- IndependentSpikeSlabPrior(cbind(1,X),Y_data,expected.model.size = p-nbzero)
  SamplesBSS <- array(0,dim=c(niter,p,nchains))
  resBSS = matrix(0,nrow=p,ncol=nchains)
  for(z in 1:nchains){
  ChainBoomSS <- lm.spike(Y_data ~ X,niter = niter,prior = priorBoomSS,model.options = SsvsOptions(),
                           error.distribution = 'gaussian',seed = z)
  SamplesBSS[,,z] <- ChainBoomSS$beta[,-1]
  resBSS[,z] <- summary.lm.spike(ChainBoomSS,burn = nburnin,order=FALSE)$coefficients[-1,][,1]
  }
  Mean_allchains_BoomSS <- apply(resBSS,1,mean)
  SamplesBSS <- SamplesBSS[-(1:nburnin),,]
  TimeBSS <- difftime(Sys.time(),start_time,units='secs')
  
  
  Df_X <- as.data.frame(X)
  Data <- list(Y=Y_data,X=X)
  CspikeN$setData(Data)
  CspikeT$setData(Data)
  CspikeD$setData(Data)
  C_HS$setData(Data)
  if(situation != 'multicollinearity') {CspikeDg$setData(Data)}
  
  start_time = Sys.time()
  mcmcN <- runMCMC(CmcmcN,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  TimeN <- difftime(Sys.time(),start_time,units='secs')
  
  start_time = Sys.time()
  mcmcT <- runMCMC(CmcmcT,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  TimeT <- difftime(Sys.time(),start_time,units='secs')
  
  start_time = Sys.time()
  mcmcD <- runMCMC(CmcmcD,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  TimeD <- difftime(Sys.time(),start_time,units='secs')
  
  start_time = Sys.time()
  mcmcHS <- runMCMC(CmcmcHS,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  TimeHS <- difftime(Sys.time(),start_time,units='secs')
  
  if(situation != 'multicollinearity'){
  start_time = Sys.time()
  mcmcDg <- runMCMC(CmcmcDg,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  TimeDg <- difftime(Sys.time(),start_time,units='secs')
  res_Time <- c(TimeN,TimeT,TimeD,TimeDg,TimeHS,TimeBSS)
  } else {res_Time <- c(TimeN,TimeT,TimeD,TimeHS,TimeBSS)}
  
  # Effective Sample Size (ESS, indépendance des tirages) :
  
  # (à voir pour généraliser pour nchain =/= 2, mcmc$samples$chain1 = mcmc$samples[[1]])
  
  
  essN = mean(effectiveSize(mcmc.list(mcmc(mcmcN$samples$chain1[,1:p]),mcmc(mcmcN$samples$chain2[,1:p]))))
  essT = mean(effectiveSize(mcmc.list(mcmc(mcmcT$samples$chain1[,1:p]),mcmc(mcmcT$samples$chain2[,1:p]))))
  essD = mean(effectiveSize(mcmc.list(mcmc(mcmcD$samples$chain1[,1:p]),mcmc(mcmcD$samples$chain2[,1:p]))))
  essHS = mean(effectiveSize(mcmc.list(mcmc(mcmcHS$samples$chain1[,1:p]),mcmc(mcmcHS$samples$chain2[,1:p]))))
  essBSS = mean(effectiveSize(mcmc.list(mcmc(SamplesBSS[,1:p,1]),mcmc(SamplesBSS[,1:p,2]))))
  
  if(situation != 'multicollinearity')
  {
    essDg = mean(effectiveSize(mcmc.list(mcmc(mcmcDg$samples$chain1[,1:p]),mcmc(mcmcDg$samples$chain2[,1:p]))))
    res_ESS = c(essN,essT,essD,essDg,essHS,essBSS)  
  } else{res_ESS = c(essN,essT,essD,essHS,essBSS) }
  
  # ESS/s (compromis vitesse/mixing) :
  
  ESS_per_sec = res_ESS/as.numeric(res_Time)
  
  # Variance Inflation Factor (VIF, Séparation des variables) :
  
  ind_N = which(abs(mcmcN$summary$all.chains[,1]) > seuil)
  formula = parse(text = paste("V",ind_N[1],sep="","~."))[[1]]
  modelN <- lm(formula = formula, data = Df_X[,ind_N])
  
  ind_T = which(abs(mcmcT$summary$all.chains[,1]) > seuil)
  formula = parse(text = paste("V",ind_T[1],sep="","~."))[[1]]
  modelT <- lm(formula = formula, data = Df_X[,ind_T])
  
  ind_D = which(abs(mcmcD$summary$all.chains[,1]) > seuil)
  formula = parse(text = paste("V",ind_D[1],sep="","~."))[[1]]
  modelD <- lm(formula = formula, data = Df_X[,ind_D])
  
  ind_HS = which(abs(mcmcHS$summary$all.chains[,1]) > seuil)
  formula = parse(text = paste("V",ind_HS[1],sep="","~."))[[1]]
  modelHS <- lm(formula = formula, data = Df_X[,ind_HS])
  
  ind_BSS = which(abs(Mean_allchains_BoomSS) > seuil)
  formula = parse(text = paste("V",ind_BSS[1],sep="","~."))[[1]]
  modelBSS <- lm(formula = formula, data = Df_X[,ind_BSS])
  
  if(situation != 'multicollinearity')
  {
    ind_Dg = which(abs(mcmcDg$summary$all.chains[,1]) > seuil)
    formula = parse(text = paste("V",ind_Dg[1],sep="","~."))[[1]]
    modelDg <- lm(formula = formula, data = Df_X[,ind_Dg])
    res_VIF <- c(mean(vif(modelN)),mean(vif(modelT)),mean(vif(modelD)),mean(vif(modelDg)),mean(vif(modelHS)),mean(vif(modelBSS)))
  } else {res_VIF <- c(mean(vif(modelN)),mean(vif(modelT)),mean(vif(modelD)),mean(vif(modelHS)),mean(vif(modelBSS)))}
  
  
  # Root Mean Squared Error (RMSE, Estimation des paramètres) :
  
  RMSE_N = sqrt(mean((mcmcN$summary$all.chains[1:p,1] - theta )^2))
  RMSE_T = sqrt(mean((mcmcT$summary$all.chains[1:p,1] - theta )^2))
  RMSE_D = sqrt(mean((mcmcD$summary$all.chains[1:p,1] - theta )^2))
  RMSE_HS = sqrt(mean((mcmcHS$summary$all.chains[1:p,1] - theta )^2))
  RMSE_BSS = sqrt(mean((Mean_allchains_BoomSS - theta )^2))
  
  if(situation != 'multicollinearity'){
    RMSE_Dg = sqrt(mean((mcmcDg$summary$all.chains[1:p,1] - theta )^2))
    res_RMSE <- c(RMSE_N,RMSE_T,RMSE_D,RMSE_Dg,RMSE_HS,RMSE_BSS)
  } else {res_RMSE <- c(RMSE_N,RMSE_T,RMSE_D,RMSE_HS,RMSE_BSS)}
  
  
  # Taux de mauvais classement (Misclassification rate) : 
  
  Misc_N = 1-mean((theta == 0) == (abs(mcmcN$summary$all.chains[1:p,1]) < seuil))
  Misc_T = 1-mean((theta == 0) == (abs(mcmcT$summary$all.chains[1:p,1]) < seuil))
  Misc_D = 1-mean((theta == 0) == (abs(mcmcD$summary$all.chains[1:p,1]) < seuil))
  Misc_HS = 1-mean((theta == 0) == (abs(mcmcHS$summary$all.chains[1:p,1]) < seuil))
  Misc_BSS = 1-mean((theta == 0) == (abs(Mean_allchains_BoomSS) < seuil))
  
  if(situation != 'multicollinearity'){
    Misc_Dg = 1-mean((theta == 0) == (abs(mcmcDg$summary$all.chains[1:p,1]) < seuil))
    res_Misc <- c(Misc_N,Misc_T,Misc_D,Misc_Dg,Misc_HS,Misc_BSS)
  } else {res_Misc <- c(Misc_N,Misc_T,Misc_D,Misc_HS,Misc_BSS)}
  
  # Sensibilité (= VP/(VP+FN)) et Spécificité (= VN/(VN+FP)) :
  
  tabN = table(vrai = (theta != 0), estim = (abs(mcmcN$summary$all.chains[1:p,1]) > seuil))
  Sensi_N = tabN[2,2]/(tabN[2,2]+tabN[2,1])
  Speci_N = tabN[1,1]/(tabN[1,1]+tabN[1,2])
  
  tabT = table(vrai = (theta != 0), estim = (abs(mcmcT$summary$all.chains[1:p,1]) > seuil))
  Sensi_T = tabT[2,2]/(tabT[2,2]+tabT[2,1])
  Speci_T = tabT[1,1]/(tabT[1,1]+tabT[1,2])
  
  tabD = table(vrai = (theta != 0), estim = (abs(mcmcD$summary$all.chains[1:p,1]) > seuil))
  Sensi_D = tabD[2,2]/(tabD[2,2]+tabD[2,1])
  Speci_D = tabD[1,1]/(tabD[1,1]+tabD[1,2])
  
  tabHS = table(vrai = (theta != 0), estim = (abs(mcmcHS$summary$all.chains[1:p,1]) > seuil))
  Sensi_HS = tabHS[2,2]/(tabHS[2,2]+tabHS[2,1])
  Speci_HS = tabHS[1,1]/(tabHS[1,1]+tabHS[1,2])
  
  tabBSS = table(vrai = (theta != 0), estim = (abs(Mean_allchains_BoomSS) > seuil))
  Sensi_BSS = tabBSS[2,2]/(tabBSS[2,2]+tabBSS[2,1])
  Speci_BSS = tabBSS[1,1]/(tabBSS[1,1]+tabBSS[1,2])
  
  
  if(situation != 'multicollinearity'){
    tabDg = table(vrai = (theta != 0), estim = (abs(mcmcDg$summary$all.chains[1:p,1]) > seuil))
    Sensi_Dg = tabDg[2,2]/(tabDg[2,2]+tabDg[2,1])
    Speci_Dg = tabDg[1,1]/(tabDg[1,1]+tabDg[1,2])
    res_Sensi <- c(Sensi_N,Sensi_T,Sensi_D,Sensi_Dg,Sensi_HS,Sensi_BSS)
    res_Speci <- c(Speci_N,Speci_T,Speci_D,Speci_Dg,Speci_HS,Speci_BSS)
  } else {res_Sensi <- c(Sensi_N,Sensi_T,Sensi_D,Sensi_HS,Sensi_BSS)
  res_Speci <- c(Speci_N,Speci_T,Speci_D,Speci_HS,Speci_BSS)}
  
  
  # Continuous Ranked Probability Score (CRPS) :
  

  list_mean <- rep(mu,n) + X %*% mcmcN$summary$all.chains[,1]
  crpsN <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

  list_mean <- rep(mu,n) + X %*% mcmcT$summary$all.chains[,1]
  crpsT <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

  list_mean <- rep(mu,n) + X %*% mcmcD$summary$all.chains[,1]
  crpsD <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

  list_mean <- rep(mu,n) + X %*% mcmcHS$summary$all.chains[,1]
  crpsHS <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))

  list_mean <- rep(mu,n) + X %*% Mean_allchains_BoomSS
  crpsBSS <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))
  
  if(situation != 'multicollinearity'){
    list_mean <- rep(mu,n) + X %*% mcmcDg$summary$all.chains[,1]
    crpsDg <- crps(Y_data,data.frame(list_mean,rep(sigma,n)))
    res_crps <- c(crpsN$CRPS,crpsT$CRPS,crpsD$CRPS,crpsDg$CRPS,crpsHS$CRPS,crpsBSS$CRPS)
  } else {res_crps <- c(crpsN$CRPS,crpsT$CRPS,crpsD$CRPS,crpsHS$CRPS,crpsBSS$CRPS)}
  
  
  # Mean Interval Score (MIS)
  
  IS_N <- IS(alpha,mcmcN$summary$all.chains[1:p,4],mcmcN$summary$all.chains[1:p,5],theta)
  IS_T <- IS(alpha,mcmcT$summary$all.chains[1:p,4],mcmcT$summary$all.chains[1:p,5],theta)
  IS_D <- IS(alpha,mcmcD$summary$all.chains[1:p,4],mcmcD$summary$all.chains[1:p,5],theta)
  IS_HS <- IS(alpha,mcmcHS$summary$all.chains[1:p,4],mcmcHS$summary$all.chains[1:p,5],theta)
  IS_BSS <- NA
  
  if(situation != 'multicollinearity')
  {
    IS_Dg <- IS(alpha,mcmcDg$summary$all.chains[1:p,4],mcmcDg$summary$all.chains[1:p,5],theta)
    res_MIS <- c(mean(IS_N),mean(IS_T),mean(IS_D),mean(IS_Dg),mean(IS_HS),mean(IS_BSS))
  } else {res_MIS <- c(mean(IS_N),mean(IS_T),mean(IS_D),mean(IS_HS),mean(IS_BSS))}
  
  Matr <- cbind(res_RMSE,res_crps,res_Misc,res_Sensi,res_Speci,res_MIS,res_Time,res_ESS,ESS_per_sec,res_VIF)
  return(Matr)
}

set.seed(1)
Matr3D <- NULL
t_init = Sys.time()
for(k in 1:q){
  print(paste0("Itération n°",k))
  Matr <- ComputeMCMCbis()
  Matr3D <- abind(Matr3D,Matr,along = 3)
}
t_final = difftime(Sys.time(),t_init,units='secs')
t_final
Final <- apply(Matr3D,c(1,2),mean)       # Matrice des moyennes sur les datasets (ie. sur la 3ème dimension)


### Tableau de comparaison des priors : -----------------------------

df = matrix(0,nrow = nbModels,ncol=nbCrit+1)
df <- data.frame(df)
colnames(df) <- c('Model','RMSE','CRPS','Misc. rate','Sensitivity','Specificity','MIS','Execution time (s)',
                  'ESS','ESS/s','Average VIF')

if(situation != 'multicollinearity'){
  df[,1] <- c('S&S continuous normal','S&S continuous student','S&S Dirac i-slab','S&S Dirac g-slab','Horseshoe',
              'Boom S&S (indep normal)')
} else {
  df[,1] <- c('S&S continuous normal','S&S continuous student','S&S Dirac i-slab','Horseshoe','Boom S&S (indep normal)')
}

df[,2] <- round(Final[,1],4)  
df[,3] <- round(Final[,2],3)
df[,4] <- round(Final[,3],3)
df[,5] <- round(Final[,4],3)
df[,6] <- round(Final[,5],3)
df[,7] <- round(Final[,6],3)
df[,8] <- round(Final[,7],1)
df[,9] <- round(Final[,8],0)
df[,10] <- round(Final[,9],1)
df[,11] <- round(Final[,10],4)

df
myt <- ttheme_default(
  core = list(fg_params=list(cex=1.0)),
  colhead = list(fg_params=list(cex=1.0)),
  rowhead = list(fg_params=list(cex=1.0)))
grid.table(df,theme=myt,rows=NULL)

df100 <- df 
#df1000 <- df
#df5000 <- df

#grob100 <- tableGrob(df100,theme=myt,rows=NULL)
#grob1000 <- tableGrob(df1000,theme=myt,rows=NULL)
#grob5000 <- tableGrob(df5000,theme=myt,rows=NULL)
#grid.arrange(grob100,grob1000,grob5000)

### Distribution des critères selon les datasets : --------------------------------------

Titles = c('RMSE Distribution','CRPS Distribution','Misc. rate Distribution','Sensitivity Distribution',
           'Specificity Distribution','MIS Distribution','Time Distribution','ESS Distribution',
           'ESS/s Distribution','Average VIF Distribution')
xlabs = c('S&S continuous normal','S&S continuous student','S&S Dirac i-slab','S&S Dirac g-slab','Horseshoe','Boom S&S (indep normal)')


# Histogrammes :
par(mfrow=c(2,3))
for(j in 1:nbCrit){
  if(j==6) {for(i in 1:(nbModels-1)){hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i])}
    plot.new()} else {
  for(i in 1:nbModels) {
    if(j==3){
      hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i],breaks=seq(0,0.2,length=10))}
    else if(j==4 | j==5){
      hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i],breaks=seq(0.5,1,length=20))}
     else {hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i])}
  }
}
}

# ESS pour Boom S&S supérieur à nchains*(niter-nburnin) (max théorique), dû à des autocorr nég ? Comment interpréter ?

## Distributions continues :
#for(j in 1:nbCrit){
#  for(i in 1:nbModels) {
#    plot(density(Matr3D[i,j,]),main=Titles[j],xlab=xlabs[i])
#  }
#  plot.new()
#}

par(mfrow=c(1,1))
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()
#-----------------------------------------------------------------------


A = 1:4
mA = mcmc(A)
mA
mcmc.list(mA)
B = 5:8
mB = mcmc(B)
mB
mcmc.list(mB)
mcmc(A,B)
mcmc(c(A,B))
l=list(A=A,B=B)
mcmc(l)
mcmc.list(l)
L= lapply(l,FUN=mcmc)
L
mcmc.list(L)
mcmc.list(mA,mB)

mcmcN$samples[[1]][1,]
