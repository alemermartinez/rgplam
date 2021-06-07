#' Tukey Loss Function
#' @export
tukey.loss <- function(x,k=4.685){
  n <- length(x)
  salida <- rep(0,n)
  for(i in 1:n){
    if(abs(x[i])<=k){salida[i] <- 1-(1-(x[i]/k)^2)^3 #(x[i])^6/(6*k^4)-(x[i])^4/(2*k^2)+(x[i])^2/2
    }else{
      salida[i] <- 1 #k^2/6
    }
  }
  return(salida)
}

#' Derivative of Huber's loss function.
#'
#' This function evaluates the first derivative of Huber's loss function.
#'
#' This function evaluates the first derivative of Huber's loss function.
#'
#' @param r a vector of real numbers
#' @param k a positive tuning constant.
#'
#' @return A vector of the same length as \code{x}.
#'
#' @author Matias Salibian-Barrera, \email{matias@stat.ubc.ca}, Alejandra Martinez
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.huber(r=x, k = 1.5)
#'
#' @export
psi.huber <- function(r, k=1.345)
  pmin(k, pmax(-k, r))

#Huber's weight function "Psi(r)/r"
psi.huber.w <- function(r, k=1.345)
  pmin(1, k/abs(r))


#Huber's loss function
rho.huber <- function(r, k=1.345){
  if(abs(r)<=k){
    return(r^2)
  }else{
    return(2*k*abs(r)-k^2)
  }
}


#' Derivative of CH loss function.
#'
#' This function evaluates the first derivative of CH loss function proposed by Croux and Haesbroeck (2003).
#'
#' This function evaluates the first derivative of CH loss function.
#'
#' @param r a vector of real numbers
#' @param k a positive tuning constant.
#'
#' @return A vector of the same length as \code{x}.
#'
#' @author Alejandra Martinez, \email{ale_m_martinez@hotmail.com}
#'
#' @examples
#' x <- seq(-2, 2, length=10)
#' psi.ch(r=x, k = 0.5)
#'
#' @export
psi.ch <- function(r, k=0.5)
  exp(-sqrt(pmax(r, k)))

#CH loss function
rho.ch <- function(r, k=0.5){
  #r <- abs(r)
  return( pmin(r,k)*exp(-sqrt(pmax(r,k)))*as.numeric(r<=k)+
        (-2*exp(-sqrt(pmax(r,k)))*(1+sqrt(pmax(r,k)))+exp(-sqrt(pmin(r,k)))*(2*(1+sqrt(pmin(r,k)))+pmin(r,k)))*as.numeric(r>k)  )
}

#prueba <- function(r,k=0.5){
#  r <- abs(r)
#  if(r<=k){
#    return(r*exp(-sqrt(k)))
#  }else{
#    return( -2*exp(-sqrt(r))*(1+sqrt(r))+exp(-sqrt(k))*(2*(1+sqrt(k))+k)   )
#  }
#}

#' Tukey Loss Function
#' @export
my.norm.2 <- function(x){
  return( sqrt(sum(x^2)) )
}



#' Classical knot selection
# #' @importFrom splines bs
# #' @importFrom stats lm
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
select.nknots.cl <- function(y,Z,X,degree.spline=3){
  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  lim.inf.kj <- ceiling(max(n^(1/5)/2,4))
  lim.sup.kj <- floor(8+2*n^(1/5))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  BIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]

    sal <- stats::lm(y~Z.aux+Xspline)
    betas <- as.vector(sal$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      gs.hat[,ell] <- aux - mean(aux)
    }

    regresion.hat <- stats::predict(sal) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
    BIC[nknots+1] <- log(sum((y - regresion.hat)^2))+(log(n)/(2*n))*(nbasis+q+1) #q+1 es la cantidad de lineales
  }
  posicion <- which.min(BIC)
  nknots <- posicion-1 #Decía "knots" en lugar de decir nknots... creo

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline
  salida <- list(nknots=nknots, BIC=BIC, grid.nknots=grid.nknots, nbasis = nbasis, kj=kj)
  return(salida)
}


#' Robust knot selection
#' @examples
#' x <- seq(-2, 2, length=10)
#' @importFrom splines bs
#' @importFrom robustbase glmrob
#' @importFrom RobStatTM logregWBY
#' @export
select.nknots.rob.gplam <- function(y, Z, X, family=family, method="MT", degree.spline=3, maxit=100){

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  fami <- family$family
  if (is.null(fami))
    stop(gettextf("'%s' is not a valid family (see ?family)",
                  as.character(call[["family"]])), domain = NA)
  if (!(fami %in% c("binomial", "poisson", "Gamma", "gaussian"))) {
    stop(gettextf("Robust GLM fitting not yet implemented for family %s",
                  fami), domain = NA)
  }

  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }


  lim.inf.kj <- ceiling(max(n^(1/5)/2,4))
  lim.sup.kj <- floor(8+2*n^(1/5))
  lim.sup.nknots <- lim.sup.kj - degree.spline - 1
  lim.inf.nknots <- lim.inf.kj - degree.spline - 1
  grid.nknots <- lim.inf.nknots:lim.sup.nknots

  RBIC <- rep(0,length(grid.nknots))

  for(nknots in grid.nknots){
    Mat.X <- as.list(rep(0,d))
    #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
    Xspline <- NULL
    for (ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
      Xspline <- cbind(Xspline,Mat.X[[ell]])
    }
    nMat <- dim(Mat.X[[ell]])[2]
    #dim(Xspline)[2]/4


    #Robust estimator
    if(fami=="poisson"){
      if(is.null(method)){
        cat("MT method applied")
        method <- "MT"
      }
      sal.r  <- glmrob(y ~ Z.aux+Xspline, family=family, method=method)
      perdida <- tukey.los
      residuos <- sal.r$residuals #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl
    }
    if(fami=="binomial"){
      cat("WBY method applied")

      sal.r  <- try( logregWBY(cbind(Z.aux,Xspline), y, intercept = 1))#logregBY(cbind(Z.aux,Xspline), y, intercept = 1) #glmrob(y ~ Z.aux+Xspline, family=family, method="Mqle") #logregWBY(cbind(Z.aux,Xspline), y, intercept = 1)
      if( class(sal.r) == 'try-error'){
        RBIC[ (nknots+1): length(grid.nknots)] <- +Inf
        break
      }
      perdida <- rho.ch
      #stop("No se pueden calcular los residuos para la selección automática de ventanas")
      residuos <- sal.r$residual.deviances
    }
    if(fami=="gaussian" | fami=="Gamma"){
      if(is.null(method)){
        cat("Mqle method applied")
        method <- "Mqle"
      }
      sal.r  <- glmrob(y ~ Z.aux+Xspline, family=family, method=method)
      perdida <- rho.huber
      residuos <- sal.r$residuals #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

    }


    betas <- as.vector(sal.r$coefficients)
    beta.hat <- betas[-1]
    coef.lin <- betas[2:(q+1)]
    coef.spl <- betas[(q+2):(1+q+nMat*d)]
    alpha.hat <- betas[1]

    gs.hat <- matrix(0,n,d)
    for(ell in 1:d){
      aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
      gs.hat[,ell] <- aux - mean(aux)
    }


    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline + 1)
    tuk <- perdida( abs(residuos) )
    RBIC[nknots+1] <- log( sum(tuk) )+ (log(n)/(2*n))*(nbasis+q+1) #q+1 porque q de la parte lineal y 1 de la constante. O sea, q+1 es la cantidad de lineales.
  }
  posicion <- which.min(RBIC)
  nknots <- posicion-1

  nbasis <- d*(nknots + degree.spline)
  kj <- nknots + degree.spline

  salida <- list(nknots=nknots, RBIC=RBIC, grid.nknots=grid.nknots, nbasis = nbasis, kj = kj)
  return(salida)


}

#' Classical Partial Linear Additive Model
# #' @importFrom splines bs
# #' @importFrom stats lm
#' @examples
#' x <- seq(-2, 2, length=10)
#' @export
plam.cl <- function(y, Z, X, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3){
  # y continuos response variable (n)
  # Z a discret or cathegorical vector (n) or matrix (n x q) for the linear part.
  # In case it is a cathegorical variable, class of Z should be 'factor'.
  # X a vector (n) or a matrix (n x d) for the additive part.
  # nknots number of internal knots
  # knots specific internal knots

  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  if( is.null(nknots) ){
    AUX <- select.nknots.cl(y,Z,X,degree.spline=degree.spline)
    nbasis <- AUX$nbasis
    nknots <- AUX$nknots
    kj <- AUX$kj
  }else{
    nbasis <- d*(nknots + degree.spline)
    kj <- nknots + degree.spline
  }

  Mat.X <- as.list(rep(0,d))
  #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
  Xspline <- NULL
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]

  sal <- stats::lm(y~Z.aux+Xspline)
  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]

  gs.hat <- matrix(0,n,d)
  correc <- rep(0,d)
  for(ell in 1:d){
    aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    correc[ell] <- mean(aux)
    gs.hat[,ell] <- aux - mean(aux)
  }

  regresion.hat <- as.vector(stats::predict(sal)) #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const = alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- punto <- as.matrix(np.point)
      }else{
        prediccion <- punto <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- punto <- np.point
    }
    np <- dim(punto)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X.new[[ell]] <- splines::bs( punto[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])
    }


    for(k in 1:np){
      for(ell in 1:d){
        aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        prediccion[,ell] <- aux - correc[ell]
      }
    }
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, coef.const=alpha.hat+sum(correc), g.matrix=gs.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y,X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}


#' Robust Partial Linear Additive Model
#' @examples
#' x <- seq(-2, 2, length=10)
# #' @importFrom splines bs
# #' @importFrom robustbase glmrob
#' @export
gplam.rob <- function(y, Z, X, family, method=NULL, np.point=NULL, nknots=NULL, knots=NULL, degree.spline=3, maxit=100){
  # y continuos response variable (n)
  # Z a discret or cathegorical vector (n) or matrix (n x q) for the linear part.
  # In case it is a cathegorical variable, class of Z should be 'factor'.
  # X a vector (n) or a matrix (n x d) for the additive part.
  # nknots number of internal knots
  # knots specific internal knots

  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  fami <- family$family
  if (is.null(fami))
    stop(gettextf("'%s' is not a valid family (see ?family)",
                  as.character(call[["family"]])), domain = NA)
  if (!(fami %in% c("binomial", "poisson", "Gamma", "gaussian"))) {
    stop(gettextf("Robust GLM fitting not yet implemented for family %s",
                  fami), domain = NA)
  }

  n <- length(y)
  d <- dim(X)[2]

  if(is.factor(Z)){
    q <- nlevels(as.factor(Z))-1 #Ahora son 4 las variables "discretas" porque z tiene rango 5
    lev.Z <- levels(Z)
    Z.aux <- matrix(0,n,nlevels(Z)-1)
    for(k in 1:(nlevels(Z)-1)){
      Z.aux[,k] <- as.numeric(Z == lev.Z[k+1]) #Dummies
    }
  }else{
    Z.aux <- Z
    q <- dim(Z)[2]
  }

  if( is.null(nknots) ){
    AUX <- select.nknots.rob.gplam(y, Z, X, family=family, method=method, degree.spline=degree.spline, maxit=maxit)
    nknots <- AUX$nknots
    nbasis <- AUX$nbasis
    kj <- AUX$kj
  }else{
    nbasis <- d*(nknots + degree.spline) #d*(nknots + degree.spline+1)
    kj <- (nknots + degree.spline) #(nknots + degree.spline + 1)
  }

  Mat.X <- as.list(rep(0,d))
  #nMat.X <- rep(0,d) #Esto lo tengo si los grados son distintos. Por ahora D=3
  Xspline <- NULL
  for (ell in 1:d){
    if(nknots>0){
      knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
    }else{
      knots <- NULL
    }
    Mat.X[[ell]] <- splines::bs( X[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
    #nMat.X[ell] <- dim(Mat.X[[ell]])[2]
    Xspline <- cbind(Xspline,Mat.X[[ell]])
  }
  nMat <- dim(Mat.X[[ell]])[2]

  #Robust estimator
  if(fami=="poisson"){
    if(is.null(method)){
      cat("MT method applied")
      method <- "MT"
    }
    sal  <- glmrob(y ~ Z.aux+Xspline, family=family, method=method)
  }
  if(fami=="binomial"){
    cat("WBY method applied")
    sal  <- logregWBY(cbind(Z.aux,Xspline), y, intercept = 1) #logregBY(cbind(Z.aux,Xspline), y, intercept = 1) #glmrob(y ~ Z.aux+Xspline, family=family, method="Mqle") #logregWBY(cbind(Z.aux,Xspline), y, intercept = 1)
  }
  if(fami=="gaussian" | fami=="Gamma"){
    if(is.null(method)){
      cat("Mqle method applied")
      method <- "Mqle"
    }
    sal  <- glmrob(y ~ Z.aux+Xspline, family=family, method=method)
  }

  betas <- as.vector(sal$coefficients)
  beta.hat <- betas[-1]
  coef.lin <- betas[2:(q+1)]
  coef.spl <- betas[(q+2):(1+q+nMat*d)]
  alpha.hat <- betas[1]

  gs.hat <- matrix(0,n,d)
  correc <- rep(0,d)
  for(ell in 1:d){
    aux <- as.vector( Xspline[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
    correc[ell] <- mean(aux)
    gs.hat[,ell] <- aux - mean(aux)
  }

  regresion.hat <- sal$fitted.values #alpha.hat + dummies%*%coef.lin + Xspline%*%coef.spl

  if(is.null(np.point)){
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, coef.const=alpha.hat+sum(correc), g.matrix=gs.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat, alpha.clean=alpha.hat, nbasis=nbasis, kj=kj)
    return(salida)
  }else{
    if(is.null(dim(np.point))){
      if(q==1){
        prediccion <- punto <- as.matrix(np.point)
      }else{
        prediccion <- punto <- t(as.matrix(np.point))
      }
    }else{
      prediccion <- punto <- np.point
    }
    np <- dim(punto)[1]
    Mat.X.new <- as.list(rep(0,d))
    Xspline.new <- NULL
    for(ell in 1:d){
      if(nknots>0){
        knots <- stats::quantile(X[,ell],(1:nknots)/(nknots+1))
      }else{
        knots <- NULL
      }
      Mat.X.new[[ell]] <- splines::bs( punto[,ell], knots=knots, degree=degree.spline, intercept=FALSE)
      Xspline.new <- cbind(Xspline.new,Mat.X.new[[ell]])
    }


    for(k in 1:np){
      for(ell in 1:d){
        aux <- as.vector( Xspline.new[,(nMat*(ell-1)+1):(nMat*ell)] %*% coef.spl[(nMat*(ell-1)+1):(nMat*ell)] )
        prediccion[,ell] <- aux - correc[ell]
      }
    }
    salida <- list(prediction=regresion.hat, coef.lin=coef.lin, alpha=alpha.hat+sum(correc), g.matrix=gs.hat, coef.const=alpha.hat, coef.spl=coef.spl, nknots=nknots, knots=knots, y=y, X=X, Z=Z.aux, Xspline=Xspline, nMat=nMat,alpha.clean=alpha.hat, nbasis=nbasis, kj=kj, np.prediction=prediccion)
    return(salida)
  }
}
