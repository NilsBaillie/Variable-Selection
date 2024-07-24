### Comparaison de priors en linéaire (optimisé)

rm(list=objects())
graphics.off()
setwd("C:/Users/nilsb/Desktop/Fichiers R/Stage INRAE codes R") # à changer selon le répertoire utilisé

library(nimble)
library(MASS)
library(coda)
library(car)
library(verification)
library(gridExtra)
library(abind)
library(doParallel)
ncores <- detectCores() - 1        # nombre de coeurs utilisés
cluster <- makeCluster(ncores,outfile = 'output_lin.txt')   # type = 'FORK' n'est pas disponible sur Windows
registerDoParallel(cluster)

### Paramètres utilisés : ----------------------------------------------

n = 100                # nombre d'observations
p = 15                 # nombre de variables
nchains = 2            # nombre de chaînes de Markov lancées
niter = 10000          # nombre d'itérations dans le MCMC
nburnin = 5000         # Temps de chauffe du MCMC
seuil = 0.1            # Seuil en-dessous duquel le theta est estimé à zéro
nechant = 2            # nombre d'échantillons obtenus sur 1 coeur
q = nechant*ncores     # nombre de datasets

theta = rep(0,nbzero)
theta = c(theta,seq(0.1,2,length=p-nbzero))


IS <- function(alpha,lower,upper,true_val)   # fonction auxiliaire pour le calcul du MIS
{return(((upper-lower) + (2/alpha)*(lower-true_val)*(true_val < lower) + (2/alpha)*(true_val-upper)*(true_val > upper)))}
alpha <- 0.05

length <- 10^4
CRPS_custom <- function(theta,samples,length){
  P <- ecdf(samples)
  interv <- seq(theta - 10^3, beta + 10^3, length = length)
  CRPS_sing <- sum((P(interv)-(interv >= theta))^2)
  return(CRPS_sing)
}

nbModels = 5
nbCrit = 10
mu = 0               # intercept de la régression linéaire (coeff du terme constant)
sigmaX = 1           # facteur de variance sur les données générées dans X
sigma = 1            # écart-type du bruit eps (pas toujours le même)


### Génération des données
set.seed(1)
X = matrix(0,nrow=n,ncol=p)
for(j in 1:p) {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # variables hétéroscédastiques, normales et indépendantes

eps = rnorm(n,0,sd = sigma)
Y_data = c(mu*rep(1,n) + X %*% theta + eps)
Df_X <- as.data.frame(X)


Data <- list(Y=Y_data,X=X)
ConstN <- list(n=n,p=p,mu=mu,a=1,b=1,V=1,r=10^(-4))
ConstT <- list(n=n,p=p,mu=mu,a=1,b=1,Q=4,nu=5,r=10^(-4))
ConstD <- list(n=n,p=p,mu=mu,a=1,b=1,c=1)
ConstDg <- list(n=n,p=p,mu=mu,a=1,b=1,g=n)
ConstHS <- list(n=n,p=p,mu=mu,sigma=sigma)
Inits <- list(sigma=1,theta = rep(1,p))


### Calcul de tous les critères pour un modèle donné ( + génération du MCMC)
ComputeModel <- function(Cmodel,Cmcmc,Data){
  Cmodel$setData(Data)
  start_time = Sys.time()
  mcmc <- runMCMC(Cmcmc,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
  Time <- difftime(Sys.time(),start_time,units='secs')
  Summ <- mcmc$summary$all.chains
  Samples <- NULL
  for(z in 1:nchains) {Samples <- rbind(Samples,mcmc$samples[[z]])}
  ess <-  mean(effectiveSize(mcmc.list(mcmc(mcmc$samples$chain1[,1:p]),mcmc(mcmc$samples$chain2[,1:p]))))
  ESS_per_sec <- ess / as.numeric(Time)
  ind_VIF <- which(abs(Summ[,1]) > seuil)
  formula <- parse(text = paste("V",ind_VIF[1],sep="","~."))[[1]]
  modelVIF <- lm(formula = formula, data = Df_X[,ind_VIF])
  res_VIF <- mean(vif(modelVIF))
  RMSE <- sqrt(mean((Summ[1:p,1] - theta )^2))
  Misc <- 1-mean((theta == 0) == (abs(Summ[1:p,1]) < seuil))
  TP <- sum((beta_true != 0) & (abs(Summ[1:p,1]) > seuil))           # 'vrai' & 'estim'
  TN <- sum((beta_true == 0) & (abs(Summ[1:p,1]) < seuil))
  FP <- sum((beta_true == 0) & (abs(Summ[1:p,1]) > seuil))
  FN <- sum((beta_true != 0) & (abs(Summ[1:p,1]) < seuil))
  Sensi <- TP/(TP+FN)
  Speci <- TN/(TN+FP)
  list_crps <- numeric(p)
  for(l in 1:p) { list_crps[l] <- CRPS_custom(theta[l],Samples[,l],length)}
  res_crps <- mean(list_crps)
  res_MIS <- mean(IS(alpha,Summ[1:p,4],Summ[1:p,5],theta))
  result <- c(RMSE,res_crps,Misc,Sensi,Speci,res_MIS,Time,ess,ESS_per_sec,res_VIF)
  return(result)
}


### Calcul en parallèle : -------------------------------------------------

clusterExport(cluster,ls())

{
t_init <- Sys.time()
Out <- clusterEvalQ(cluster,{
  library(MASS)
  library(nimble)
  library(coda)
  library(car)
  library(verification)

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
  SpikeD <- nimbleModel(code = SpikeDirac,constants = ConstD,data = Data,inits = Inits)
  CspikeD <- compileNimble(SpikeD)
  confD <- configureMCMC(model = CspikeD,monitors = c('theta'))
  BuildD <- buildMCMC(confD)
  CmcmcD <- compileNimble(BuildD)
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
  SpikeDg <- nimbleModel(code = SpikeDiracG,constants = ConstDg,data = Data,inits = Inits)
  CspikeDg <- compileNimble(SpikeDg)
  confDg <- configureMCMC(model = CspikeDg,monitors = c('theta'))
  BuildDg <- buildMCMC(confDg)
  CmcmcDg <- compileNimble(BuildDg)
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
  Horseshoe <- nimbleModel(code = PriorHorseshoe,constants = ConstHS,data = Data,inits = Inits)
  C_HS <- compileNimble(Horseshoe)
  confHS <- configureMCMC(model = C_HS,monitors = c('theta'))
  BuildHS <- buildMCMC(confHS)
  CmcmcHS <- compileNimble(BuildHS)
  Result <- array(dim = c(nbModels,nbCrit,nechant))

  for(z in 1:nechant){
  X = matrix(0,nrow=n,ncol=p)
  for(j in 1:p) {X[,j] = rnorm(n,0,sd = sqrt(j)*sigmaX)}   # variables hétéroscédastiques, normales et indépendantes
  eps = rnorm(n,0,sd = sigma)
  Y_data = c(mu*rep(1,n) + X %*% theta + eps)
  Df_X <- as.data.frame(X)
  Data <- list(Y=Y_data,X=X)
  N_res <- ComputeModel(CspikeN,CmcmcN,Data)
  T_res <- ComputeModel(CspikeT,CmcmcT,Data)
  D_res <- ComputeModel(CspikeD,CmcmcD,Data)
  Dg_res <- ComputeModel(CspikeDg,CmcmcDg,Data)
  HS_res <- ComputeModel(C_HS,CmcmcHS,Data)
  Result[,,z] <- rbind(N_res,T_res,D_res,Dg_res,HS_res)
  }
  return(Result)
})
t_final = difftime(Sys.time(),t_init,units='secs')
stopCluster(cluster)
}
t_final

Matr3D <- NULL
for(z in 1:ncores) {Matr3D <- abind(Matr3D,Out[[z]],along=3)}
Matr3D                                   # passage d'une liste de 3d-array en un seul 3d-array
Final <- apply(Matr3D,c(1,2),mean)       # Matrice des moyennes sur les datasets (ie. sur la 3ème dimension)

### Affichage des résultats : -------------------------


df = matrix(0,nrow = nbModels,ncol=nbCrit+1)
df <- data.frame(df)
colnames(df) <- c('Model','RMSE','CRPS','Misc. rate','Sensitivity','Specificity','MIS','Execution time (s)',
                  'ESS','ESS/s','Average VIF')


df[,1] <- c('S&S continuous normal','S&S continuous student','S&S Dirac i-slab','S&S Dirac g-slab','Horseshoe')
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

### Histogrammes :

Titles = c('RMSE Distribution','CRPS Distribution','Misc. rate Distribution','Sensitivity Distribution',
           'Specificity Distribution','MIS Distribution','Time Distribution','ESS Distribution',
           'ESS/s Distribution','Average VIF Distribution')
xlabs = c('S&S continuous normal','S&S continuous student','S&S Dirac i-slab','S&S Dirac g-slab','Horseshoe')

par(mfrow=c(2,3))
for(j in 1:nbCrit){
      for(i in 1:nbModels) {
        if(j==3){
          hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i],breaks=seq(0,0.8,length=10))}
        else if(j==4 | j==5){
          hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i],breaks=seq(0,1,length=20))}
        else {hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i])}
      }
  plot.new()
}


### Temps d'exécution total selon n, p et q : -------------------

# Sur le nombre d'observations : avec p = 15, q = 14 fixés

tps_n <- c(270,333,543,860)
ns <- c(100,200,500,1000)
reglin1 <- lm(tps_n~ns)
plot(ns,tps_n,xlab = TeX('$n$'),ylab = 'Temps (en s)',
     main = TeX("Temps d'exécution total en fonction de $n$ avec $p = 15$ et $q = 14$ "),
     pch=16)
lines(x = seq(0,1000,length = 1000),reglin1$coefficients[1] + seq(0,1000,length = 1000)*reglin1$coefficients[2],col='red')
grid()

# Sur le nombre de datasets : avec n = 100, p = 15 fixés

tps_q <- c(270,315,440,612,720)
qs <- c(14,21,49,84,105)
reglin2 <- lm(tps_q ~ qs)
plot(qs,tps_q,xlab = TeX('$q$'),ylab = 'Temps (en s)',
     main = TeX("Temps d'exécution total en fonction de $q$ avec $n = 100$ et $p = 15$ "),pch=16)
lines(x = seq(0,200,length = 1000),reglin2$coefficients[1] + seq(0,200,length = 1000)*reglin2$coefficients[2],col='red')
grid()


# Sur le nombre de variables : avec n = 100, q = 14 fixés


tps_p <- c(270,423,567,973,1318)
ps <- c(15,30,50,80,100)
plot(ps,tps_p,xlab = TeX('$p$'),ylab = 'Temps (en s)',
     main = TeX("Temps d'exécution total en fonction de $p$ avec $n = 100$ et $q = 14$ "),pch=16)
reglin3 <- lm(tps_p~ps)
summary(reglin3)
lines(x = seq(0,100,length = 1000),reglin3$coefficients[1] + seq(0,100,length = 1000)*reglin3$coefficients[2],col='red')
grid()

# -> Temps
