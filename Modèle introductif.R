### Simulation par un modèle non linéaire (1)

rm(list=objects())
graphics.off()
library(latex2exp)

# Fonction du modèle 

Mod_NL = function(y0,yinf,tau,t)
{
  return(yinf+2*(y0-yinf)/(1+exp(t/tau)))
}

# Interprétation des paramètres -------------------

# y0 : valeur initiale de la taille de la plante
# ymax : valeur maximale théorique (pas forcément atteinte) de la taille
# yinf : valeur asymptotique de la taille 
# tau : temps caractéristique pour atteindre yinf
# sigma : écart-type du bruit gaussien centré qui est la partie aléatoire du modèle
# T : nombre de mesures
# t : ensemble des temps où les mesures ont été faites (ici, à intervalles réguliers)


# Exemple où la croissance exponentielle (concave) est visible et où la variance du bruit est faible.
# On se base sur cet exemple pour perturber les différents paramètres et faire des observations.

set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10

yinf = 8
#yinf = runif(1,0,ymax)

sigma = 0.25
#sigma = runif(1,0,0.01*yinf)
eps = rnorm(T,0,sd = sigma)

y0 = 5
#y0 = runif(1,0,yinf)

#tau = runif(1,0,0.5*t[T])
tau = 10

X_ex = Mod_NL(y0,yinf,tau,t)
Y_ex = Mod_NL(y0,yinf,tau,t) + eps  # Exemple de données
plot(t,Y_ex,type='l')
lines(t,X_ex,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2)      

#acf(Y_ex)
#pacf(Y_ex)


# yinf = 8 -> yinf = 6 ---------------
set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10
yinf = 6
sigma = 0.25
eps = rnorm(T,0,sd = sigma)
y0 = 5
tau = 10

X = Mod_NL(y0,yinf,tau,t)
Y = Mod_NL(y0,yinf,tau,t) + eps  
par(mfrow=c(1,2))
plot(t,Y,type='l',ylim=c(4.5,8.2))
lines(t,X,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 
plot(t,Y_ex,type='l',ylim=c(4.5,8.2))
lines(t,X_ex,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 

# Quand on diminue yinf avec y0 constant, la courbe semble plus aplatie dans le sens des ordonnées,
# la croissance exponentielle est plus lente et l'effet inverse apparaît si on augmente yinf sans modifier y0.
# De même, cet aplatissement apparaît si on augmente y0 sans changer yinf (effet opposé à l'augmentation de yinf)



# (y0,yinf) = (5,8) -> (y0,yinf) = (4.5,7.5) ---------------
# y0 - yinf est constant
set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10
yinf = 7.5
sigma = 0.25
eps = rnorm(T,0,sd = sigma)
y0 = 4.5
tau = 10

X = Mod_NL(y0,yinf,tau,t)
Y = Mod_NL(y0,yinf,tau,t) + eps  
par(mfrow=c(1,2))
plot(t,Y,type='l',ylim=c(4.5,8.2))
lines(t,X,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 
plot(t,Y_ex,type='l',ylim=c(4.5,8.2))
lines(t,X_ex,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 

# Ici, l'apparence du graphe ne change pas, il est seulement translaté (vers le bas), donc les variations 
# individuelles de y0 et yinf peuvent grandement impacter l'aspect des données, or quand ils varient
# simultanément de la "même manière", l'aspect général du graphe reste très similaire à celui de départ.
# Le paramètre delta_y = y0 - yinf n'est alors pas bien identifiable dans ce modèle.
# Une idée pour l'estimer serait d'utiliser la relation f'(0) = - delta_y / (2*tau),
# or le graphe des données n'est pas du tout régulier à cause du bruit gaussien, il est difficile d'apprécier sa dérivée. 


# sigma = 0.25 -> sigma = 0.75 ---------------------------
set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10
yinf = 8
sigma = 0.75
eps = rnorm(T,0,sd = sigma)
y0 = 5
tau = 10

X = Mod_NL(y0,yinf,tau,t)
Y = Mod_NL(y0,yinf,tau,t) + eps  
par(mfrow=c(1,2))
plot(t,Y,type='l',ylim=c(4.5,9))
lines(t,X,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 
plot(t,Y_ex,type='l',ylim=c(4.5,9))
lines(t,X_ex,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 

# La variance du bruit gaussien est très grande, il est difficile d'affirmer que les données suivent bien notre 
# modèle exponentiel, excepté peut-être pour les premières valeurs (t < 20). Pour que les données soient 
# similaires à la courbe théorique du modèle, il faut que sigma soit bien inférieur à y0 et yinf.
# (à l'ordre de grandeur des mesures plus généralement)


# tau = 10 -> tau = 80 ------------------------
set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10
yinf = 8
sigma = 0.25
eps = rnorm(T,0,sd = sigma)
y0 = 5
tau = 80

X = Mod_NL(y0,yinf,tau,t)
Y = Mod_NL(y0,yinf,tau,t) + eps  
par(mfrow=c(1,2))
plot(t,Y,type='l',ylim=c(4.5,9))
lines(t,X,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 
plot(t,Y_ex,type='l',ylim=c(4.5,9))
lines(t,X_ex,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 

# Quand tau est suffisamment grand par rapport à t[T] (tau > 50), la croissance exponentielle n'est pas visible 
# car l'échelle de temps n'est pas assez grande. Au lieu de cela, l'évolution apparaît comme étant affine / linéaire.


# tau = 10 -> tau = 250 ---------------
set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10
yinf = 8
sigma = 0.25
eps = rnorm(T,0,sd = sigma)
y0 = 5
tau = 250

X = Mod_NL(y0,yinf,tau,t)
Y = Mod_NL(y0,yinf,tau,t) + eps  
par(mfrow=c(1,2))
plot(t,Y,type='l',ylim=c(4.5,9))
lines(t,X,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 
plot(t,Y_ex,type='l',ylim=c(4.5,9))
lines(t,X_ex,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2) 

# Quand tau est très grand par rapport à t[T], ce qui est possible si on utilise comme prior pour tau la loi uniforme
# entre 0 et 3*t[T], l'évolution de la tendance des données semble presque stationnaire sur l'échelle de temps choisie.
# Cela est cohérent car quand tau tend vers +inf, y_t = y0 + eps_t.
# On peut aussi dire que tau est difficilement identifiable quand l'échelle de temps n'est pas adaptée, car par rapport à l'exemple
# précédent, seule l'inclinaison de la droite est légérement modifiée mais tau est passé de 80 à 250 !

# Une remarque supplémentaire est que, lorsque tau est très grand, on ne voit pas la valeur asymptotique yinf directement sur le graphe
# si l'échelle de temps est trop courte. On peut estimer y0 en regardant la première valeur des mesures, mais, comme on a vu que
# le paramètre delta_y = y0 - yinf ne donne pas d'informations particulières, c'est-à-dire qu'on ne peut pas le déduire de la 
# forme de la courbe, on ne peut pas avoir une estimation de yinf à partir de y0 et de delta_y.
# C'est une situation où le paramètre yinf est difficile à estimer, on n'a pas d'informations supplémentaires à part que 
# yinf est inclus dans [y0,ymax].


### Utilisation de NIMBLE ----------------------------------------- 
library(nimble)
library(coda)

Code <- nimbleCode({
  # Priors
  yinf ~ dunif(0,ymax)
  y0 ~ dunif(0,ymax)
  sigma ~ dunif(0,ymax)
  tau ~ dunif(0, alpha*t[T])
  # Model
  for(i in 1:T)
  {
    y[i] ~ dnorm(yinf+2*(y0-yinf)/(1+exp(t[i]/tau)), sd = sigma)
  }
})

Constants <- list(T=500, t = seq(0,100,length=T),alpha=3,ymax=10)
Data <- list(y = Y_ex) # Y_ex obtenu avec tau = 10, sigma = 0.25, y0 = 5, yinf = 8
Inits <- list(y0 = 2,yinf=6,sigma = 0.7,tau=15)
Model1 <- nimbleModel(code = Code, name = 'mod1', constants = Constants,data = Data,inits = Inits)

#Model1$getNodeNames()
#Model1$y0    
#Cplant <- compileNimble(Model1)

set.seed(1)
mcmc.out <- nimbleMCMC(code = Code, constants = Constants,data = Data, nchains = 2, inits = Inits,
                       niter = 10000,nburnin = 5000,summary=TRUE,monitors = c('y0','yinf','tau','sigma'))

#mcmc.out$samples
mcmc.out$summary

# Les moyennes a posteriori ressemblent fortement aux paramètres utilisés pour générer les données

sigma_sample1 = mcmc.out$samples$chain1[,1]
sigma_sample2 = mcmc.out$samples$chain2[,1]
tau_sample1 = mcmc.out$samples$chain1[,2]
tau_sample2 = mcmc.out$samples$chain2[,2]
y0_sample1 = mcmc.out$samples$chain1[,3]
y0_sample2 = mcmc.out$samples$chain2[,3]
yinf_sample1 = mcmc.out$samples$chain1[,4]
yinf_sample2 = mcmc.out$samples$chain2[,4]



### Distributions a posteriori --------------------------

plot(density(sigma_sample1),main=TeX("Distributions a posteriori de $\\sigma$"),xlab=TeX("$\\sigma$"),col='red')
lines(density(sigma_sample2),col='blue')
abline(v=mcmc.out$summary$all.chains[1,1])
abline(h = 1/10,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

plot(density(tau_sample1),main=TeX("Distributions a posteriori de $\\tau$"),xlab=TeX("$\\tau$"),col='red')
lines(density(tau_sample2),col='blue')
abline(v=mcmc.out$summary$all.chains[2,1])
abline(h = 1/300,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

plot(density(y0_sample1),main=TeX("Distributions a posteriori de $\\y_0$"),xlab=TeX("$\\y_0$"),col='red')
lines(density(y0_sample2),col='blue')
abline(v=mcmc.out$summary$all.chains[3,1])
abline(h = 1/10,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

plot(density(yinf_sample1),main=TeX("Distributions a posteriori de $\\y_{\\infty}$"),xlab=TeX("$\\y_{\\infty}$"),col='red')
lines(density(yinf_sample2),col='blue')
abline(v=mcmc.out$summary$all.chains[4,1])
abline(h = 1/10,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

# Bien que l'on ait utilisé des priors (très) peu informatifs, les lois a posteriori sont très concentrées
# autour des vraies valeurs.



### Traceplots ------------------------

plot(sigma_sample1,type='l',col='red',ylab = TeX("$\\sigma$"),main =TeX('Chaînes de Markov pour $\\sigma$'))
lines(sigma_sample2,col='blue')
legend('topright',c('chain1','chain2'),lty =c(1,1), col = c('red','blue'))
grid()

plot(tau_sample1,type='l',col='red',ylab = TeX("$\\tau$"),main =TeX('Chaînes de Markov pour $\\tau$'))
lines(tau_sample2,col='blue')
legend('topright',c('chain1','chain2'),lty =c(1,1), col = c('red','blue'))
grid()

plot(y0_sample1,type='l',col='red',ylab = TeX("$\\y_0$"),main =TeX('Chaînes de Markov pour $\\y_0$'))
lines(y0_sample2,col='blue')
legend('topright',c('chain1','chain2'),lty =c(1,1), col = c('red','blue'))
grid()

plot(yinf_sample1,type='l',col='red',ylab = TeX("$\\y_{\\infty}$"),main =TeX('Chaînes de Markov pour $\\y_{\\infty}$'))
lines(yinf_sample2,col='blue')
legend('topright',c('chain1','chain2'),lty =c(1,1), col = c('red','blue'))
grid()

### Effective sample size

effectiveSize(mcmc.list(mcmc(sigma_sample1),mcmc(sigma_sample2)))   # 2316.871
effectiveSize(mcmc.list(mcmc(tau_sample1),mcmc(tau_sample2)))       # 658.55
effectiveSize(mcmc.list(mcmc(y0_sample1),mcmc(y0_sample2)))         # 820.4836
effectiveSize(mcmc.list(mcmc(yinf_sample1),mcmc(yinf_sample2)))     # 1081.463


### Inférence avec très grande variance -----------------


set.seed(1)
T = 500
t = seq(0,100,length=T)
alpha = 3
ymax = 10
yinf = 8
sigma = 5  # (20 fois plus importante !)
eps = rnorm(T,0,sd = sigma)
y0 = 5
tau = 10
X = Mod_NL(y0,yinf,tau,t)
Y = Mod_NL(y0,yinf,tau,t) + eps  # Exemple de données
plot(t,Y,type='l')
lines(t,X,col='red')
grid(nx = NULL, ny = NULL,lty = 2,col = "gray", lwd = 2)   



# Même modèle, idem pour Code, juste Data qui change

set.seed(1)
Data2 = list(y = Y)
mcmc2 = nimbleMCMC(code = Code, constants = Constants,data = Data2, nchains = 2, inits = Inits,
           niter = 10000,nburnin = 5000,summary=TRUE,monitors = c('y0','yinf','tau','sigma'))

mcmc2$summary

sigma_sample1 = mcmc2$samples$chain1[,1]
sigma_sample2 = mcmc2$samples$chain2[,1]
tau_sample1 = mcmc2$samples$chain1[,2]
tau_sample2 = mcmc2$samples$chain2[,2]
y0_sample1 = mcmc2$samples$chain1[,3]
y0_sample2 = mcmc2$samples$chain2[,3]
yinf_sample1 = mcmc2$samples$chain1[,4]
yinf_sample2 = mcmc2$samples$chain2[,4]

# Graphes pour sigma_true = 5

plot(density(sigma_sample1),main=TeX("Distributions a posteriori de $\\sigma$"),xlab=TeX("$\\sigma$"),col='red',ylim=c(0,3))
lines(density(sigma_sample2),col='blue')
abline(v=sigma)
abline(h = 1/10,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

plot(density(tau_sample1),main=TeX("Distributions a posteriori de $\\tau$"),xlab=TeX("$\\tau$"),col='red')
lines(density(tau_sample2),col='blue')
abline(v=tau)
abline(h = 1/300,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

plot(density(y0_sample1),main=TeX("Distributions a posteriori de $\\y_0$"),xlab=TeX("$\\y_0$"),col='red')
lines(density(y0_sample2),col='blue') 
abline(v=y0)
abline(h = 1/10,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

plot(density(yinf_sample1),main=TeX("Distributions a posteriori de $\\y_{\\infty}$"),xlab=TeX("$\\y_{\\infty}$"),col='red')
lines(density(yinf_sample2),col='blue')
abline(v=yinf)
abline(h = 1/10,col='orange')
grid()
legend('topright',c('chain1','chain2','prior'),lty =c(1,1,1), col = c('red','blue','orange'))

# Lois a posteriori (beaucoup) moins concentrées qu'avec sigma_true = 0.25 et pas forcément autour de la vraie valeur ! 
# (voir y0)
# La différence avec le prior uniforme est moins visible


# Chaînes de Markov : 
plot(sigma_sample1,type='l',col='red',ylab = TeX("$\\sigma$"),main =TeX('Chaînes de Markov pour $\\sigma$'))
lines(sigma_sample2,col='blue')
legend('topright',c('chain1','chain2'),lty =c(1,1), col = c('red','blue'))
grid()
plot(tau_sample1,type='l',col='red',ylab = TeX("$\\tau$"),main =TeX('Chaînes de Markov pour $\\tau$'))
lines(tau_sample2,col='blue')
legend('topright',c('chain1','chain2'),lty =c(1,1), col = c('red','blue'))
grid()



effectiveSize(mcmc.list(mcmc(sigma_sample1),mcmc(sigma_sample2)))   # 2420.652 
effectiveSize(mcmc.list(mcmc(tau_sample1),mcmc(tau_sample2)))       # 191.1938
effectiveSize(mcmc.list(mcmc(y0_sample1),mcmc(y0_sample2)))         # 144.3337
effectiveSize(mcmc.list(mcmc(yinf_sample1),mcmc(yinf_sample2)))     # 557.7154

# ESS assez faibles, sigma élevé rend l'estimation des paramètres difficile



