### Simulation par un modèle non linéaire avec un prior Spike and Slab

### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4239132/

rm(list=objects())
graphics.off()
setwd("C:/Users/nilsb/Desktop/Fichiers R/Stage INRAE codes R")
library(latex2exp)
library(nimble)
library(MASS)
library(coda)
library(car)
library(verification)
library(gridExtra)
library(abind)

# Fonction du modèle :

f <- nimbleFunction(
  run = function(phi=double(0),psi1=double(0),psi2=double(0),t=double(0)) {
    result <- psi1 / (1+exp(-(t-phi)/psi2))
    return(result)
    returnType(double(0))
  })

IS <- function(alpha,lower,upper,true_val)   # fonction auxiliaire pour le calcul du MIS
{return(((upper-lower) + (2/alpha)*(lower-true_val)*(true_val < lower) + (2/alpha)*(true_val-upper)*(true_val > upper)))}
alpha <- 0.05

length <- 10^4
CRPS_custom <- function(beta,samples,length){
  P <- ecdf(samples)
  interv <- seq(beta - 10^3, beta + 10^3, length = length)
  CRPS_sing <- sum((P(interv)-(interv >= beta))^2)
  return(CRPS_sing)
}

ComputeModel <- function(Cmodel,Cmcmc,Data){
  Cmodel$setData(Data)
  start_time = Sys.time()
  MCMC <- runMCMC(Cmcmc,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
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
  Misc <- 1-mean((beta_true == 0) == (abs(Summ[1:p,1]) < seuil))
  TP <- sum((beta_true != 0) & (abs(Summ[1:p,1]) > seuil))           # 'vrai' & 'estim'
  TN <- sum((beta_true == 0) & (abs(Summ[1:p,1]) < seuil))
  FP <- sum((beta_true == 0) & (abs(Summ[1:p,1]) > seuil))
  FN <- sum((beta_true != 0) & (abs(Summ[1:p,1]) < seuil))
  Sensi <- TP/(TP+FN)
  Speci <- TN/(TN+FP)
  list_crps <- numeric(p)
  for(l in 1:p) { list_crps[l] <- CRPS_custom(beta_true[l],Samples[,l],length)}
  res_crps <- mean(list_crps)
  res_MIS <- mean(IS(alpha,Summ[1:p,4],Summ[1:p,5],beta_true))
  result <- c(RMSE,res_crps,Misc,Sensi,Speci,res_MIS,Time,ess,ESS_per_sec#,res_VIF)
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

### Paramètres fixés et données générées : ------------------------------

{
nchains = 2        # nombre de chaînes de Markov lancées
niter = 10000      # nombre d'itérations dans le MCMC
nburnin = 5000     # Temps de chauffe du MCMC
seuil = 1
nbModels = 2
nbCrit = 10

multivariateNodesAsScalars = TRUE

n = 50
p = 5
sigma = sqrt(30)
gamma = 200
psi1 = 200
psi2 = 300
J = 10
t = seq(150,3000,length=J)
mu = 1200
beta_true = c(100,50,20,rep(0,p-3))

Cas_Corr = 0     # Choix du scénario de corrélation
rho_sigma = 0.99
Matr_rho = matrix(0,nrow=p,ncol=p)
for(j in 1:p) {Matr_rho[,j] <- rho_sigma^(abs(1:p-j))}
A <- matrix(0,nrow=3,ncol=p-3)
A[3,] <- rho_sigma^(abs(3-4:p))
Corr1 <- rbind(cbind(diag(3),matrix(0,nrow=3,ncol=p-3)),cbind(matrix(0,nrow=p-3,ncol=3),Matr_rho[4:p,4:p]))
Corr2 <- rbind(cbind(diag(3),A),cbind(t(A),diag(p-3))) 
Corr3 <- rbind(cbind(Matr_rho[1:3,1:3],matrix(0,nrow=3,ncol=p-3)),cbind(matrix(0,nrow=p-3,ncol=3),diag(p-3)))
Sigma <- (Cas_Corr == 0)*diag(p)+(Cas_Corr == 1)*Corr1+(Cas_Corr == 2)*Corr2+(Cas_Corr == 3)*Corr3+(Cas_Corr == 4)*Matr_rho
#eigen(Sigma)$values

set.seed(1)
V <- mvrnorm(n,rep(0,p),Sigma)   # Vi : ième ligne de V
eps = rnorm(n*J,mean = 0,sd = sigma)
Eps <- matrix(eps,nrow=n,ncol=J,byrow=T)
ksi <- rnorm(n,0,sd = sqrt(gamma))
phi <- rep(mu,n) + apply(V,1,function(x) {t(beta_true) %*% x}) + ksi

nu0 = 1          # paramètres sur les priors
nu1 = 10^4
g = n
a = 1
b = p
}
x <- matrix(0,nrow=n,ncol=J)
for(j in 1:J){
  x[,j] <- f(phi,psi1,psi2,t[j])
}
#det(t(V)%*%V)
Df_V <- as.data.frame(V)
y_data <- x + Eps

Rnorm <- diag(p)
Rgslab <- g*inverse(t(V) %*% V)

Const <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,mu=mu,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)
Datanorm <- list(y=y_data,V=V,R=Rnorm)
Datagslab <- list(y=y_data,V=V,R=Rgslab)
Inits <- list(delta = rep(1,p))


### Modèle Spike and Slab : ----------------------------

SpikeSlab <- nimbleCode({
  alpha ~ dbeta(a,b)
  for(l in 1:p){delta[l] ~ dbern(alpha)}
  D_half[1:p,1:p] <- diag(sqrt(nu0*(1-delta[1:p]) + nu1*delta[1:p]))     # nu0 : spike, nu1 : slab
  Mean[1:p] <- rep(0,p)
  CovMatr[1:p,1:p] <-  D_half[1:p,1:p] %*% R[1:p,1:p] %*% D_half[1:p,1:p]  # R : normal ou g-slab
  beta[1:p] ~ dmnorm(mean = Mean[1:p],cov = CovMatr[1:p,1:p])
  for(i in 1:n){
    phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})

Modelnorm <- nimbleModel(code = SpikeSlab,const = Const, data = Datanorm,inits = Inits)
Cnorm <- compileNimble(Modelnorm)
confnorm <- configureMCMC(model = Cnorm, monitors = 'beta',multivariateNodesAsScalars = multivariateNodesAsScalars)
buildnorm <- buildMCMC(confnorm)
Cmcmcnorm <- compileNimble(buildnorm)

Modelgslab <- nimbleModel(code = SpikeSlab,const = Const, data = Datagslab,inits = Inits)
Cgslab <- compileNimble(Modelgslab)
confgslab <- configureMCMC(model = Cgslab, monitors = 'beta',multivariateNodesAsScalars = multivariateNodesAsScalars)
buildgslab <- buildMCMC(confgslab)
Cmcmcgslab <- compileNimble(buildgslab)


### Chaînes de Markov et distributions a posteriori : -----------------------------

set.seed(1)
MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
MCMCnorm$summary$all.chains
MCMCgslab <- runMCMC(Cmcmcgslab,nchains = nchains, niter = niter, nburnin = nburnin,summary=TRUE)
MCMCgslab$summary$all.chains

Titles <- paste("beta",1:p,sep = "")

# Distributions a posteriori
par(mfrow=c(3,3))
for(l in 1:p){
  plot(density(MCMCnorm$samples$chain1[,l]),type='l',main=Titles[l],xlab = 'Normal S&S')
  lines(density(MCMCnorm$samples$chain2[,l]),col='red')
}
par(mfrow=c(3,3))
for(l in 1:p){
  plot(density(MCMCgslab$samples$chain1[,l]),type='l',main=Titles[l],xlab = 'G-slab S&S')
  lines(density(MCMCgslab$samples$chain2[,l]),col='red')
}

# Chaînes de Markov
par(mfrow=c(3,3))
for(l in 1:p){
  plot(MCMCnorm$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'Normal S&S')
  lines(MCMCnorm$samples$chain2[,l],col='red')
}

par(mfrow=c(3,3))
for(l in 1:p){
  plot(MCMCgslab$samples$chain1[,l],type='l',main=Titles[l],ylab="",xlab = 'G-slab S&S')
  lines(MCMCgslab$samples$chain2[,l],col='red')
}

### Comparaison par critères : -----------------------

set.seed(1)
norm_res <- ComputeModel(Cnorm,Cmcmcnorm,Datanorm)
gslab_res <- ComputeModel(Cgslab,Cmcmcgslab,Datagslab)
Final <- rbind(norm_res,gslab_res)

df = matrix(0,nrow = nbModels,ncol=nbCrit+1)
df <- data.frame(df)
colnames(df) <- c('Model','RMSE','CRPS','Misc. rate','Sensitivity','Specificity','MIS','Execution time (s)',
                  'ESS','ESS/s','Average VIF')

df[,1] <- c('S&S Continuous Normal','S&S Continuous g-slab')
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
par(mfrow=c(1,1))
plot.new()
myt <- ttheme_default(
  core = list(fg_params=list(cex=1.0)),
  colhead = list(fg_params=list(cex=1.0)),
  rowhead = list(fg_params=list(cex=1.0)))
grid.table(df,theme=myt,rows=NULL)


par(mfrow=c(1,1))
# Etude des distributions a priori de beta : ----------------------------------------------
nsamples <- 10^4
library(plotly)


# 1) delta fixé : -------------------------

### Cnorm :
Cnorm$delta <- c(1,1,0,0,0)
nodesdelta <- Cnorm$getDependencies(c('delta'),self = FALSE)
Cnorm$simulate(nodesdelta)
Cnorm$delta
Cnorm$CovMatr
Cnorm$beta

set.seed(1)
Scpl <- matrix(0,nrow=nsamples,ncol=p)   ### Echantillons pour scatterplots
for(i in 1:nsamples){
  Cnorm$simulate(nodesdelta)
  Scpl[i,] <- Cnorm$beta
}
plot(density(Scpl[,1]),main = latex2exp::TeX("Distributions a priori"),xlab = latex2exp::TeX('$\\beta$'))
lines(density(Scpl[,2]),col='blue')
lines(density(Scpl[,3]),col='red')
grid()
#hist(Scpl[,1])
plot(Scpl[,1],Scpl[,2],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_2$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
grid()
den3d <- kde2d(Scpl[,1],Scpl[,2])
contour(den3d,lwd=2,add=TRUE,col = hcl.colors(10,'Spectral'))
filled.contour(den3d,xlab=latex2exp::TeX('$\\beta_1$'),ylab=latex2exp::TeX('$\\beta_2$'),
               main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
plot(Scpl[,1],Scpl[,3],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_3$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_3$'))
grid()

plot_ly(x=den3d$x,y=den3d$y,z=den3d$z) %>% add_surface()
# x : beta1, y : beta2

### Cgslab :

set.seed(1)
Cgslab$delta <- c(1,1,0,0,0)
nodesdelta <- Cgslab$getDependencies(c('delta'),self = FALSE)
Cgslab$simulate(nodesdelta)
Cgslab$delta
Cgslab$CovMatr
Cgslab$beta

set.seed(1)
nsamples <- 10000
Scpl <- matrix(0,nrow=nsamples,ncol=p)
for(i in 1:nsamples){
  Cgslab$simulate(nodesdelta)
  Scpl[i,] <- Cgslab$beta
}

pdf('Test.pdf',width = 4,height = 4)
plot(density(Scpl[,1]),main = latex2exp::TeX("Distributions a priori"),xlab = latex2exp::TeX('$\\beta$'))
lines(density(Scpl[,2]),col='blue')
lines(density(Scpl[,3]),col='red')
grid()
#legend(x=c(200,400),y=c(0.003,0.004),legend = latex2exp::TeX(c('$\\beta_1$','$\\beta_2$','$\\beta_3$')),col=c('black','blue','red'),lty = c(1,1,1),cex=1.5,lwd=1)
#hist(Scpl[,1])
plot(Scpl[,1],Scpl[,2],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_2$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
grid()
den3d <- kde2d(Scpl[,1],Scpl[,2])
contour(den3d,lwd=2,add=TRUE,col = hcl.colors(10,'Spectral'))
filled.contour(den3d,xlab=latex2exp::TeX('$\\beta_1$'),ylab=latex2exp::TeX('$\\beta_2$'))
plot(Scpl[,1],Scpl[,3],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_3$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_3$'))
grid()

plot_ly(x=den3d$x,y=den3d$y,z=den3d$z) %>% add_surface()

dev.off()
# x : beta1, y : beta2


# 2) simulations sans paramètres fixés :   --------------------------
nsamples <- 10^4

### Cnorm :

set.seed(1)
Scpl <- matrix(0,nrow=nsamples,ncol=p)
for(i in 1:nsamples){
  nodestosim <- Cnorm$getDependencies(c('alpha','delta'),self = TRUE)
  Cnorm$simulate(nodestosim)
  Scpl[i,] <- Cnorm$beta
}

plot(density(Scpl[,1]))
lines(density(Scpl[,2]),col='red')
lines(density(Scpl[,3]),col='blue')
lines(density(Scpl[,4]),col='darkgreen')
grid()
plot(Scpl[,1],Scpl[,2],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_2$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
grid()
den3d <- kde2d(Scpl[,1],Scpl[,2])
contour(den3d,lwd=2,add=TRUE,col = hcl.colors(10,'Spectral'))
filled.contour(den3d,xlab=latex2exp::TeX('$\\beta_1$'),ylab=latex2exp::TeX('$\\beta_2$'),
               main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
plot(Scpl[,1],Scpl[,3],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_3$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_3$'))
grid()

plot_ly(x=den3d$x,y=den3d$y,z=den3d$z) %>% add_surface()
# x : beta1, y : beta2

#rmnorm_chol(1,rep(0,5),chol(Cnorm$CovMatr),prec_param = FALSE)

### Cgslab :

set.seed(1)
Scpl <- matrix(0,nrow=nsamples,ncol=p)
for(i in 1:nsamples){
  nodestosim <- Cgslab$getDependencies(c('alpha','delta'),self = TRUE)
  Cgslab$simulate(nodestosim)
  Scpl[i,] <- Cgslab$beta
}

plot(density(Scpl[,1]))
lines(density(Scpl[,2]),col='red')
lines(density(Scpl[,3]),col='blue')
lines(density(Scpl[,4]),col='darkgreen')
grid()
plot(Scpl[,1],Scpl[,2],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_2$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
grid()
den3d <- kde2d(Scpl[,1],Scpl[,2])
contour(den3d,lwd=2,add=TRUE,col = hcl.colors(10,'Spectral'))
filled.contour(den3d,xlab=latex2exp::TeX('$\\beta_1$'),ylab=latex2exp::TeX('$\\beta_2$'),
               main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_2$'))
plot(Scpl[,1],Scpl[,3],xlab = latex2exp::TeX("$\\beta_1$"),ylab = latex2exp::TeX("$\\beta_3$"),
     main = latex2exp::TeX('Distribution jointe a priori de $\\beta_1$ et $\\beta_3$'))
grid()

plot_ly(x=den3d$x,y=den3d$y,z=den3d$z) %>% add_surface()
# x : beta1, y : beta2

### Samplers : -------------------------
confgslab$printSamplers('beta')  # -> RW_block sampler
#confgslab$getSamplerDefinition('alpha')


#inte = seq(-50,50,length=10000)
#plot(inte,dnorm(inte,0,sd=1),type='l',xlim=c(-50,50),col='red',main='Spike and Slab',xlab=TeX('$x$'),ylab = 'Densité')
#lines(inte,dnorm(inte,0,sd=10),col='blue')
#grid()
#legend(x=27,y=0.3,legend = c('Spike','Slab'),col=c('red','blue'),lty=c(1,1))


### Calculs en parallèle (comparaison des priors) : ----------------------------

library(doParallel)
ncores <- detectCores() - 1        # nombre de coeurs utilisés
cluster <- makeCluster(ncores,outfile = 'output_Gslab.txt')   # type = 'FORK' n'est pas disponible sur Windows
registerDoParallel(cluster)

nechant = 7
nbDataset <- ncores*nechant

clusterExport(cluster,ls())
set.seed(1)
t_init <- Sys.time()
{
  Output <- clusterEvalQ(cluster,{
    library(MASS)
    library(nimble)
    library(coda)
    library(car)
    library(verification)

    SpikeSlab <- nimbleCode({
      alpha ~ dbeta(a,b)
      for(l in 1:p){delta[l] ~ dbern(alpha)}
      D_half[1:p,1:p] <- diag(sqrt(nu0*(1-delta[1:p]) + nu1*delta[1:p]))     # nu0 : spike, nu1 : slab
      Mean[1:p] <- rep(0,p)
      CovMatr[1:p,1:p] <-  D_half[1:p,1:p] %*% R[1:p,1:p] %*% D_half[1:p,1:p]  # R : normal ou g-slab
      beta[1:p] ~ dmnorm(mean = Mean[1:p],cov = CovMatr[1:p,1:p])
      for(i in 1:n){
        phi[i] ~ dnorm(mu + sum(beta[1:p]*V[i,1:p]),var = gamma)
        for(j in 1:J){
          y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
        }}
    })
    Modelnorm <- nimbleModel(code = SpikeSlab,const = Const, data = Datanorm,inits = Inits)
    Cnorm <- compileNimble(Modelnorm)
    confnorm <- configureMCMC(model = Cnorm, monitors = 'beta',multivariateNodesAsScalars = multivariateNodesAsScalars)
    buildnorm <- buildMCMC(confnorm)
    Cmcmcnorm <- compileNimble(buildnorm)
    Modelgslab <- nimbleModel(code = SpikeSlab,const = Const, data = Datagslab,inits = Inits)
    Cgslab <- compileNimble(Modelgslab)
    confgslab <- configureMCMC(model = Cgslab, monitors = 'beta',multivariateNodesAsScalars = multivariateNodesAsScalars)
    buildgslab <- buildMCMC(confgslab)
    Cmcmcgslab <- compileNimble(buildgslab)
    Result <- array(dim = c(2,nbCrit,nechant))
    for(z in 1:nechant){
      V <- mvrnorm(n,rep(0,p),Sigma)
      y_data <- Gener_Data(V,beta_true,sigma,gamma)
      Rgslab <- g*inverse(t(V) %*% V)
      Datanorm <- list(y=y_data,V=V,R=Rnorm)
      Datagslab <- list(y=y_data,V=V,R=Rgslab)
      norm_res <- ComputeModel(Cnorm,Cmcmcnorm,Datanorm)
      gslab_res <- ComputeModel(Cgslab,Cmcmcgslab,Datagslab)
      Result[,,z] <- rbind(norm_res,gslab_res)
    }
    return(Result)
})
t_final = difftime(Sys.time(),t_init,units='secs')
stopCluster(cluster)
}
t_final
Output

Matr3D <- NULL
for(z in 1:ncores) {Matr3D <- abind(Matr3D,Output[[z]],along=3)}
Matr3D                                   # passage d'une liste de 3d-array en un seul 3d-array
Final <- apply(Matr3D,c(1,2),mean)       # Matrice des moyennes sur les datasets (ie. sur la 3ème dimension)
Final

df = matrix(0,nrow = nbModels,ncol=nbCrit+1)
df <- data.frame(df)
colnames(df) <- c('Model','RMSE','CRPS','Misc. rate','Sensitivity','Specificity','MIS','Execution time (s)',
                  'ESS','ESS/s','Average VIF')


df[,1] <- c('S&S continuous normal','S&S continuous G-slab')
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



### Histogrammes : ----------------------------

Titles = c('RMSE Distribution','CRPS Distribution','Misc. rate Distribution','Sensitivity Distribution',
           'Specificity Distribution','MIS Distribution','Time Distribution','ESS Distribution',
           'ESS/s Distribution','Average VIF Distribution')
xlabs = c('S&S continuous normal','S&S continuous G-slab')

par(mfrow=c(1,2))
for(j in 1:nbCrit){
  for(i in 1:nbModels) {
    if(j==3){
      hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i],breaks=seq(0,0.8,length=10))}
    else if(j==4 | j==5){
      hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i],breaks=seq(0,1,length=20))}
    else {hist(Matr3D[i,j,],main=Titles[j],xlab=xlabs[i])}
  }
}









