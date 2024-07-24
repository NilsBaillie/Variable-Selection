### Nimble Tests
library(nimble)

n = 100 
p = 5  
set.seed(1)
y <- rnorm(n)

Code <- nimbleCode({ y[1:n] ~ dnorm(0,sd = 1) })
Data <- list(y=y)
Const <- list(n=n)
Model <- nimbleModel(code = Code,const = Const,data = Data)
### -> Pas possible de vectoriser une distribution à 1 dimension (nécessité de passer par une boucle for)
# Error in model$checkBasics() : 
#Dimension of 'y[1:100]' does not match required dimension for the distribution 'dnorm'. Necessary dimension is 0.


### Bon à savoir pour les priors Horseshoe : Cauchy(mu,gamma) = Student(mu, tau = 1/gamma^2, df = 1)
# la distribution de Cauchy n'étant pas de base dans Nimble

Code <- nimbleCode({
  for(j in 1:p) { delta[j] ~ dbern(0.5)}
  d <- sum(delta[1:p])
  vect[1:d] <- 5*(1:d)
  
})
Model <- nimbleModel(code = Code,const = Const,data = Data)
# -> erreur d'indice, utilisation d'indices aléatoires (d ici) interdite en indexation dynamique ?

# NB : éviter la transposée t() dans NIMBLE pour les vecteurs