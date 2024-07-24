### Graphes HS et HS+ 

rm(list=objects())
graphics.off()
library(MASS)
library(latex2exp)

I_l = seq(0.01,5,by = 0.001)
I_k = seq(0.001,0.999,by = 0.001)


######## Graphes des distributions a priori des lambda : 

HS_lambda <- function(l){
  return((2/pi)*(1/(1+l^2)))
}

HSp_lambda <- function(l){
  return((2/pi^2)*log(l^2)/(l^2-1))
} 


plot(I_l,HS_lambda(I_l),type='l',col='blue',lwd=2, main = TeX('Distributions a priori de $\\lambda$'),xlab = TeX('$\\lambda$'),ylab='Densité',
     ylim = c(0,2))
lines(I_l,HSp_lambda(I_l),col='red',lwd=2)
grid()
legend('topright',legend = c('HS','HS+'),col = c('blue','red'),lty = 1,lwd = 2)

######## Graphes des distributions a priori des kappa : 

HS_kappa <- function(k){return(dbeta(k,1/2,1/2))}

HSp_kappa <- function(k){
  beta <- dbeta(k,1/2,1/2)
  return((beta/pi)*log((1-k)/k)/(1-2*k))
}

eqscplot(I_k,HS_kappa(I_k),ratio = 1,type='l',col='blue',lwd=2,main = TeX('Distributions a priori de $\\kappa$'),xlab = TeX('$\\kappa$'),
         ylab='Densité',ylim = c(0.3,1.5))
lines(I_k,HSp_kappa(I_k),col='red',lwd=2)
#abline(v=0,lty=2,lwd=2)
#abline(v=1,lty=2,lwd=2)
grid()
legend('topright',legend = c('HS','HS+'),col = c('blue','red'),lty = 1,lwd = 2)

