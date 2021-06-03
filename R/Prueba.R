#Esto es una prueba

#Ejemplo de Valdora y Yohai (2014)
set.seed(123)

n <- 100
X <- cbind( rep(1,100), matrix(rnorm(100*5,0,1),100,5) )
beta0 <- c(2,1,0,0,0,0)

y <- rep(0,n)
for(i in 1:n){
  y[i] <- rpois(1, lambda= exp(as.numeric(beta0%*%X[i,])) )
}

#Robust estimation
library(robustbase)
sal.r <- glmrob(y~X-1, family = poisson, method="MT")
sal.r$coefficients

sal.r$fitted.values[1:5]
exp(sal.r$coefficients%*%t(X))[1:5] #Son iguales

sal.r$residuals
sal.r$linear.predictors
sal.r$optim.control
sal.r$model

#Classical estimation
sal.c <- glm(y~X-1, family=poisson)
sal.c$coefficients
sal.c$fitted.values
sal.c$residuals

#Ejemplo parcialmente lineal aditivo generalizado

#Ejemplo Poisson
set.seed(123)
n <- 100
beta <- c(3,3)
function.g1 <- function(x1) 2*sin(pi*x1)-4/pi
function.g2 <- function(x2) exp(x2)-(exp(1)-1)

x1 <- runif(n)
x2 <- runif(n)
z1 <- runif(n)
z2 <- runif(n)
Z <- cbind(z1,z2)
X <- cbind(x1,x2)
y <- rep(0,n)
for(i in 1:n){
  y[i] <- rpois(1, lambda= exp(as.numeric(beta%*%Z[i,] + function.g1(x1[i])+function.g2(x2[i]))) )
}

nknots <- 2
degree.spline <- 3
method="MT"
np.point=NULL
library(robustbase)
sal1 <- gplam.rob(y, Z, X, family=poisson, method="MT", np.point=NULL, nknots=nknots, knots=NULL, degree.spline=3, maxit=100)
sal1$coef.lin
sal1$coef.const
sal1$g.matrix

sal1$prediction[1:5]
exp(sal1$coef.const + as.vector(sal1$coef.lin%*%t(Z)+rowSums(sal1$g.matrix)))[1:5]
y[1:5]


#Ejemplo Logistico
set.seed(123)
n <- 100
beta <- c(3,3)
function.g1 <- function(x1) 2*sin(pi*x1)-4/pi
function.g2 <- function(x2) exp(x2)-(exp(1)-1)

x1 <- runif(n)
x2 <- runif(n)
z1 <- runif(n)
z2 <- runif(n)
Z <- cbind(z1,z2)
X <- cbind(x1,x2)
y <- rep(0,n)
log.dis <- function(x){
  return(exp(x)/(1+exp(x)))
}
for(i in 1:n){
  y[i] <- rbinom(1, 1, prob= log.dis(as.numeric(-2+ beta%*%Z[i,] + function.g1(x1[i])+function.g2(x2[i]))))
}

nknots <- 2
degree.spline <- 3
method="MT"
np.point=NULL
library(RobStatTM)
sal1 <- gplam.rob(y, Z, X, family=binomial, method="MT", np.point=NULL, nknots=nknots, knots=NULL, degree.spline=3, maxit=100)
sal1$coef.lin
sal1$coef.const
sal1$g.matrix

sal1$prediction[1:5]
exp(sal1$coef.const + as.vector(sal1$coef.lin%*%t(Z)+rowSums(sal1$g.matrix)))[1:5]
y[1:5]

