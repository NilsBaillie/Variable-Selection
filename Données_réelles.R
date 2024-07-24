### Etude sur les données réelles

rm(list=objects())
graphics.off()
setwd("G:/Mon Drive/M1 Maths Applis/Documents Stage INRAE/Stage INRAE codes R")  # à modifier 
library(latex2exp)
library(nimble)
library(MASS)
library(coda)
library(gridExtra)

######################## Reprise du code de Marion ############################################

load("Saves_Senescence/data_chr6A.Rdata")

load("Saves_Senescence/data_obs.Rdata")
sspop <- readRDS("Data_Senescence/resPCOAdf.Rds")
data_chromosome <- read.csv2("Data_Senescence/carte_Axiom-TABW420k_WGAv1.csv",header=TRUE, sep = " ")

varieties=unique(data_obs$GENOTYPE)
obs_date=unique(data_obs$Day)
n=length(unique(data_obs$GENOTYPE))
J=length(obs_date)
Id=rep(c(1:n),each=J)
id=as.matrix(Id)
Y=c()
for (i in varieties){
  for (j in obs_date){
    Y=c(Y,data_obs[(data_obs$GENOTYPE==i)&(data_obs$Day==j),]$Y)
  }
}
t=rep(unique(data_obs$Day),n)

V_tilde=data_chr6A$V_tilde
V=data_chr6A$V
nb_QTL_chr6A=data_chr6A$nb_QTL #nb de QTL de floraison présents sur chr6A après pre-processing

p=dim(V)[2] #nb de covariables

################################################################################################

# paramètres généraux

f <- nimbleFunction(
  run = function(phi=double(0),psi1=double(0),psi2=double(0),t=double(0)) {
    result <- psi1 / (1+exp(-(t-phi)/psi2))
    return(result)
    returnType(double(0))
  })

a = 1
b = p
nu0 = 10^-4
nu1 = 1
psi1 = 100

# paramètres chaînes de Markov
nchains = 2        # nombre de chaînes de Markov lancées
niter = 15000      # nombre d'itérations dans le MCMC
nburnin = 5000     # Temps de chauffe du MCMC

seuil <- 0.01

# valeurs estimées par M. Naveau : 
sigma = sqrt(15.2)
gamma = 1.1
psi2 = 2.8
V_pca <- V_tilde[,1:6]

y_data <- matrix(Y,nrow=n,ncol=J,byrow=TRUE)
Data <- list(y=y_data,V=V,V_pca=V_pca)
Inits <- vector('list',nchains)
Inits[[1]] <- list(beta = rep(0,p))
Inits[[2]] <- list(beta = rep(0.5,p))
ConstNorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)

### Spike and Slab Continu Normal :

SpikeSlab_Norm <- nimbleCode({
  alpha ~ dbeta(a,b)
  for(z in 1:6){
    lambda[z] ~ dnorm(0,var = nu1)}  # intercept et 5 premières covariables sont significatives -> slab
  for(l in 1:p) {
    delta[l] ~ dbern(alpha)
    beta[l] ~ dnorm(0,var=(1-delta[l])*nu0 + delta[l]*nu1)}    # nu0 : spike, nu1 : slab
  for(i in 1:n){
    phi[i] ~ dnorm(sum(lambda[1:6]*V_pca[i,1:6]) + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})

### Spike and Slab Dirac Normal : 

Dirac_Norm <- nimbleCode({
  alpha ~ dbeta(a,b)
  for(z in 1:6){
    lambda[z] ~ dnorm(0,var = nu1)}  # intercept et 5 premières covariables sont significatives -> slab
  for(l in 1:p) {
    delta[l] ~ dbern(alpha)
    slab[l] ~  dnorm(0,var=nu1)
    beta[l] <- delta[l]*slab[l]}   
  for(i in 1:n){
    phi[i] ~ dnorm(sum(lambda[1:6]*V_pca[i,1:6]) + sum(beta[1:p]*V[i,1:p]),var = gamma)
    for(j in 1:J){
      y[i,j] ~ dnorm(f(phi[i],psi1,psi2,t[j]),sd = sigma)
    }}
})

######### Calcul sur tout l'échantillon #########

Modelnorm <- nimbleModel(code = SpikeSlab_Norm,const = ConstNorm, data = Data,inits = Inits)
Cnorm <- compileNimble(Modelnorm)
confnorm <- configureMCMC(model = Cnorm, monitors = c('beta','lambda'))
buildnorm <- buildMCMC(confnorm)        
Cmcmcnorm <- compileNimble(buildnorm)

# set.seed(1)
# t_init <- Sys.time()  # 6.9 heures à calculer ! 
# MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
# t_final = difftime(Sys.time(),t_init,units='secs')
# t_final
# saveRDS(MCMCnorm,file = 'Senescence_Norm.rds')

MCMCnorm <- readRDS('Senescence_Norm.RDS')
Summnorm <- MCMCnorm$summary$all.chains
head(Summnorm)
tail(Summnorm)

indexes <- which(abs(Summnorm[,1]) > seuil  )
Summnorm[indexes,1]

par(mfrow=c(3,3))
for(l in indexes){
  plot(MCMCnorm$samples$chain1[,l],type = 'l')
  lines(MCMCnorm$samples$chain2[,l],col='red')
}

ModelDirac <- nimbleModel(code = Dirac_Norm,const = ConstNorm, data = Data,inits = Inits)
CDirac <- compileNimble(ModelDirac)
confDirac <- configureMCMC(model = CDirac, monitors = c('beta','lambda'))
buildDirac <- buildMCMC(confDirac)
CmcmcDirac <- compileNimble(buildDirac)

# set.seed(1)
# t_init <- Sys.time()  # ? heures à calculer ! 
# MCMCDirac <- runMCMC(CmcmcDirac,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE)
# t_final = difftime(Sys.time(),t_init,units='secs')
# t_final
# saveRDS(MCMCDirac,file = 'Senescence_Dirac.RDS')

MCMCDirac <- readRDS('Senescence_Dirac.RDS')
SummDirac <- MCMCDirac$summary$all.chains
head(SummDirac)
tail(SummDirac)

indexes <- which(abs(SummDirac[,1]) > seuil  )
SummDirac[indexes,1]

par(mfrow=c(3,3))
for(l in indexes){
  plot(MCMCDirac$samples$chain1[,l],type = 'l')
  lines(MCMCDirac$samples$chain2[,l],col='red')
}
### absence de convergence des chaînes ###


######### Calcul sur un sous-échantillon des données avec des paramètres différents #########

nu0 = 10^-4
nu1 = 1
p = 300
b = p
which(colnames(V)=='cfn2905337') # indice du marqueur conservé par Marion = 355
keep_var <- V[,355]
V = V[,1:(p-1)]
V <- cbind(V,keep_var)
Data <- list(y=y_data,V=V,V_pca=V_pca)
Inits <- vector('list',nchains)
Inits[[1]] <- list(beta = rep(0,p))
Inits[[2]] <- list(beta = rep(0.5,p))

ConstNorm <- list(n=n,p=p,J=J,psi1=psi1,psi2=psi2,a=a,b=b,nu0=nu0,nu1=nu1,sigma=sigma,gamma=gamma,t=t)
Modelnorm <- nimbleModel(code = SpikeSlab_Norm,const = ConstNorm, data = Data,inits = Inits)
Cnorm <- compileNimble(Modelnorm)
confnorm <- configureMCMC(model = Cnorm, monitors = c('beta','lambda'))
buildnorm <- buildMCMC(confnorm)        
Cmcmcnorm <- compileNimble(buildnorm)

set.seed(1)
t_init <- Sys.time()   
MCMCnorm <- runMCMC(Cmcmcnorm,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
t_final = difftime(Sys.time(),t_init,units='secs')   # ~ 35 min
t_final

Summnorm <- MCMCnorm$summary$all.chains
head(Summnorm)
tail(Summnorm)
seuil <- 0.01

indexes <- which(abs(Summnorm[,1]) > seuil  )
Summnorm[indexes,1]
par(mfrow=c(3,3))
for(l in indexes){
  plot(MCMCnorm$samples$chain1[,l],type = 'l')
  lines(MCMCnorm$samples$chain2[,l],col='red')
}

heading_QTL <- read.csv("Data_Senescence/marker_HD_INRmon13LN.csv",sep="")[,1]      #list of heading QTLs names
heading_QTL <- as.character(heading_QTL)
pos_QTL6A=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr6A"]
lenind <- length(indexes)
indexes <- indexes[-((lenind-5):lenind)]
names_marq <- colnames(V)[indexes]
names_marq[length(names_marq)] <- 'cfn2905337'
pos_selec <- data_chromosome$V2[data_chromosome$name%in%names_marq & data_chromosome$V1=="chr6A"] # marq. sélectionnés

Bool <- matrix(0,ncol=length(pos_selec),nrow = length(pos_QTL6A))
for(z in 1:length(pos_selec))
{
  Bool[,z] <- abs(pos_selec[z] - pos_QTL6A) <= 10^6
}
Bool   # seul le dernier marqueur sélectionné correspond à un marqueur de floraison


ModelDirac <- nimbleModel(code = Dirac_Norm,const = ConstNorm, data = Data,inits = Inits)
CDirac <- compileNimble(ModelDirac)
confDirac <- configureMCMC(model = CDirac, monitors = c('beta','lambda'))
buildDirac <- buildMCMC(confDirac)
CmcmcDirac <- compileNimble(buildDirac)

set.seed(1)
t_init <- Sys.time()   
MCMCDirac <- runMCMC(CmcmcDirac,nchains = nchains,niter = niter,nburnin = nburnin,summary = TRUE,inits = Inits)
t_final = difftime(Sys.time(),t_init,units='secs')  # ~ 46 min
t_final

SummDirac <- MCMCDirac$summary$all.chains
head(SummDirac)
tail(SummDirac)
seuil <- 0.01

indexes <- which(abs(SummDirac[,1]) > seuil  )
SummDirac[indexes,1]
par(mfrow=c(3,3))
for(l in indexes){
  plot(MCMCDirac$samples$chain1[,l],type = 'l')
  lines(MCMCDirac$samples$chain2[,l],col='red')
}


heading_QTL <- read.csv("Data_Senescence/marker_HD_INRmon13LN.csv",sep="")[,1]      #list of heading QTLs names
heading_QTL <- as.character(heading_QTL)
pos_QTL6A=data_chromosome$V2[data_chromosome$name%in%heading_QTL & data_chromosome$V1=="chr6A"]
lenind <- length(indexes)
indexes <- indexes[-((lenind-5):lenind)]
names_marq <- colnames(V)[indexes]
pos_selec <- data_chromosome$V2[data_chromosome$name%in%names_marq & data_chromosome$V1=="chr6A"] # marq. sélectionnés

Bool <- matrix(0,ncol=length(pos_selec),nrow = length(pos_QTL6A))
for(z in 1:length(pos_selec))
{
  Bool[,z] <- abs(pos_selec[z] - pos_QTL6A) <= 10^6
}
Bool        # aucun des marqueurs sélectionnés ne correspond à un marqueur de floraison

