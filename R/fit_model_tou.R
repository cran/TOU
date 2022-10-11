

# #' @title The mean of the Transformed Ornstein-Uhlenbeck (TOU) stochastic
# #' process.
# #' @description Evaluates the mean function of the Transformed
# #' Ornstein-Uhlenbeck (TOU) stochastic process.
# #' @param x a scalar or a vector with the value for the time.
# #' @param parameters a vector with the values of qe, lambda, and a.
# #' @return the corresponding mean TOU value for x.
# #' @examples
# #' # evaluating a specific value of time
# #' mean.tou(x=1.5, parameters=c(14.75, 19.2, 66.27))
# #' # evaluating some values of time
# #' mean.tou(x=c(0.5,1.0,1.5,5.0,9.2), parameters=c(14.75, 19.2, 66.27))
mean.tou <- function(x, parameters) {
  qe <- parameters[1]
  lambda <- parameters[2]
  a <- parameters[3]

  value <- NA*x

  if(0<qe & 0<lambda & 0<=a) {

    if( a>0 ) {
      value <- qe*(1 - (1/(1+a*x)^(lambda/a)))
    } else {
      value <- qe*(1 - exp( -lambda*x ) )
    }

  } else {
    warning("At least one parameter is invalid: \n qe = ", qe, " lambda = ",
            lambda, " a = ", a, "\n")
  }
  return(value)
}



# #' @title The log-likelihood function for a Transformed Ornstein-Uhlenbeck (TOU)
# #' stochastic process.
# #' @description Evaluates the log-likelihood function for a Transformed
# #' Ornstein-Uhlenbeck (TOU) stochastic process.
# #' @param x a vector with the values of qe, lambda, a and tau.
# #' @param y a nx2 matrix with the observed times in the first column and with
# #' the observed values of the dependent variable in the second column.
# #' @return the value of the corresponding log-likelihood.
# #' @examples
# #' observed.time <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200,
# #'                    250, 300)
# #' observed.values <- c(0, 1.271, 1.693, 1.845, 1.919, 1.940, 1.999, 2.008,
# #'                      2.101, 2.064, 2.046, 2.105, 2.013, 2.040, 2.017)
# #' observed.process <- cbind(observed.time, observed.values)
# #' loglikelihood.tou(x=c(2.1, 0.4, 0.51, 0.15), y=observed.process)

loglikelihood.tou <- function(x, y) {
  n <- nrow(y)
  qe <- x[1]
  lambda <- x[2]
  a <- x[3]
  tau <- x[4]

  if(qe <= 0) {
    warning("The value of qe is invalid and will be changed \n
            from qe = ", qe, " to qe = ", mean(y[,2]), "\n")
    qe <- mean(y[,2])
  }
  if(lambda <= 0) {
    warning("The value of lambda is invalid and will be changed \n
            from lambda = ", lambda, " to lambda = ", qe/2, "\n")
    lambda <- qe/2
  }
  if(a < 0) {
    warning("The value of a is invalid and will be changed \n
            from a = ", a, " to a = ", 0, "\n")
    a <- 0
  }
  if(tau <= 0) {
    warning("The value of tau is invalid and will be changed \n
            from tau = ", tau, " to tau = ", 1e-8, "\n")
    tau <- 1e-8
  }

  t.increments <- y[2:n, 1] - y[1:(n-1), 1]
  y.increments <- y[2:n, 2] - y[1:(n-1), 2]
  if( a>0 ) {
    g.increments <- (1/a)*log(1+a*y[2:n, 1]) -
      (1/a)*log(1+a*y[1:(n-1), 1])
  } else {
    g.increments <- y[2:n, 1] - y[1:(n-1), 1]
  }
  mu.increments <- (qe - y[1:(n-1), 2])*(1-exp(-lambda*g.increments))
  sigma.increments <- tau*sqrt((1-exp(-2*lambda*g.increments)))
  value <- sum(stats::dnorm(y.increments, mean=mu.increments, sd=sigma.increments,
                            log=TRUE))
  return(value)
}



# #' @title The value of lambda based on parameters qe, kn and a.
# #' @description Calculates the value of lambda that satisfies the following
# #' equation lambda - kn*(qe^(a/lambda)) = 0.
# #' @param y a vector with the values of qe, kn and a.
# #' @return the value of lambda.
# #' @examples
# #' abc(y=c(14.75, 0.001762, 66.27))

lambda.root <- function(y) {
  a <- 0.0001
  b <- max(c(y[3],100))

  f.1 <- function(x) {
    if(a<=x & x<=b) {
      valor <- x - (y[1]^(y[3]/x))*y[2]
    } else if(x < a) {
      valor <- a - (y[1]^(y[3]/a))*y[2]
    } else {
      valor <- b - (y[1]^(y[3]/b))*y[2]
    }
    return(valor)
  }

  salida <- (a+b)/2
  if(f.1(a)*f.1(b) < 0) {
    salida <- stats::uniroot(f=f.1, interval = c(a, b))$root
  }
  return(salida)
}



# #' @title The intermediate log-likelihood function for a Transformed
# #' Ornstein-Uhlenbeck (TOU) stochastic process.
# #' @description Creates the log-likelihood function for the parameters
# #' (qe,lambda,a,tau ) of a TOU stochastic process. This function allows to fix
# #' some parameters of the pseudo-n-order (PNO) model, such as the maximum
# #' adsorption capacity (qe), the adsorption rate constant (kn) and the order of
# #' the model (n).
# #' @param y a matrix with the observed values. The first column has the
# #' observed times and the following columns have the observed values
# #' of the dependent variable, one column per experiment.
# #' @param qe an optional scalar indicating a fixed value for the parameter qe.
# #' @param kn an optional scalar indicating a fixed value for the parameter kn.
# #' @param n an optional scalar indicating a fixed value for the parameter n.
# #' @return a function which evaluates the corresponding log-likelihood.
# #' @examples
# #' observed.time <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200,
# #'                    250, 300)
# #' observed.values <- c(0, 1.271, 1.693, 1.845, 1.919, 1.940, 1.999, 2.008,
# #'                      2.101, 2.064, 2.046, 2.105, 2.013, 2.040, 2.017)
# #' observed.process <- cbind(observed.time, observed.values)
# #' # none parameters are fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process)
# #' ll.tou(x=c(2.1, 0.4, 0.51, 0.15))
# #' # the parameter qe is fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, qe=2.1)
# #' ll.tou(x=c(0.4, 0.51, 0.15))
# #' # the parameter kn is fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, kn=0.0355144)
# #' ll.tou(x=c(2.1, 0.51, 0.15))
# #' # the parameter n is fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, n=3.822317)
# #' ll.tou(x=c(2.1, 0.4, 0.15))
# #' # the parameters qe and kn are fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, qe=2.1, kn=0.0355144)
# #' ll.tou(c(0.51, 0.15))
# #' # the parameters qe and n are fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, qe=2.1, n=3.822317)
# #' ll.tou(c(0.4, 0.15))
# #' # the parameters kn and n are fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, kn=0.0355144, n=3.822317)
# #' ll.tou(c(2.1, 0.15))
# #' # the parameters qe, kn and n are fixed
# #' ll.tou <- make.loglikelihood.tou(y=observed.process, qe=2.1, kn=0.0355144,
# #'                           n=3.822317)
# #' ll.tou(c(0.15))
# #'
make.loglikelihood.tou <- function(y, qe=NULL, kn=NULL, n=NULL) {
  m <- ncol(y)
  null.size <- sum(c(is.null(qe), is.null(kn), is.null(n)))

  function(x) {
    if(length(x)==(null.size+1)) {
      parameters <- c(0,0,0,0)
      if(null.size==3) {
        # x = (qe, lambda, a, tau)
        parameters <- x
      } else if(null.size==2) {
        if(!is.null(qe)) {
          # x = (lambda, a, tau)
          parameters[c(2,3,4)] <- x
          parameters[1] <- qe
        } else if(!is.null(kn)) {
          # x = (qe, a, tau)
          parameters[c(1,3,4)] <- x
          parameters[2] <- lambda.root(y=c(parameters[1], kn, parameters[3]))
        } else {
          # x = (qe, lambda, tau)
          parameters[c(1,2,4)] <- x
          parameters[3] <- parameters[2]*(n-1)
        }
      } else if(null.size==1) {
        if(!is.null(qe) & !is.null(kn)) {
          # x = (a, tau)
          parameters[c(3,4)] <- x
          parameters[1] <- qe
          parameters[2] <- lambda.root(y=c(parameters[1], kn, parameters[3]))
        } else if(!is.null(qe) & !is.null(n)) {
          # x = (lambda, tau)
          parameters[c(2,4)] <- x
          parameters[1] <- qe
          parameters[3] <- parameters[2]*(n-1)
        } else {
          # x = (qe, tau)
          parameters[c(1,4)] <- x
          parameters[2] <- kn * (parameters[1]^(n-1))
          parameters[3] <- parameters[2]*(n-1)
        }
      } else {
        # x = (tau)
        parameters[4] <- x
        parameters[1] <- qe
        parameters[2] <- kn * (parameters[1]^(n-1))
        parameters[3] <- parameters[2]*(n-1)
      }
      value <- 0
      for(i in 2:m) {
        value <- value + loglikelihood.tou(x=parameters, y=y[,c(1,i)])
      }
    }
    return(value)
  }
}


#' @title Estimation of Transformed Ornstein-Uhlenbeck model for
#' adsorption kinetics.
#'
#' @description
#' This function finds the best values for the parameters of the pseudo-n-order
#' (PNO) model and its related Transformed Ornstein-Uhlenbeck (TOU) model. It also provides information about
#' some goodness of fit measures.
#'
#' This function allows to freely estimate the parameters of the TOU model
#' by fixing some parameters of the related pseudo-n-order (PNO) model, such as
#' the maximum adsorption capacity (qe), the adsorption rate constant (kn) and
#' the order of the model (n).
#'
#' @param w a matrix with the observed values. The first column has the
#' observed times and the following columns have the observed values of the
#' dependent variable, one column per experimental trajectory.
#' @param qe an optional scalar indicating a fixed value for the parameter qe.
#' @param kn an optional scalar indicating a fixed value for the parameter kn.
#' @param n an optional scalar indicating a fixed value for the parameter n.
#'
#' @details
#' This package is based on the methodology provided by
#' \insertCite{rodriguez.2021}{TOU} and it is designed to model the adsorption
#' kinetics of removal of contaminants such as dyes, metal ions, fluorides,
#' cyanides, arsenates, arsenites, antibiotics, hormones, etc. from an aqueous
#' phase by several materials. In this case, the only limiting factor in the
#' decrease of pollutant concentration is adsorption on the surface of the
#' adsorbent and does not exclude the existence of intraparticle and film
#' diffusion phenomena.
#'
#' This function provides the following parameters: the maximum adsorption (qe),
#' the value of the fractional-order kinetics (n), the adsorption rate constant
#' (kn), the parameter related to the reaction order n of the kinect process (a),
#' rate constant (lambda), the rate constant related to the variance of the
#' Brownian motion (sigma), the variance of the long term kinetics (tau), the
#' time needed to reach half of the maximum adsorption (t.half.qe) and the time
#' needed to reach the maximum adsorption (t.reach.qe). The goodness of fit
#' measures provided for this function are the loglikelihood of the model
#' (logLikelihood), the coefficient of determination (R2), the mean squared
#' error (MSE) and the normalized standard deviation (delta.q).
#'
#' @return a list with the estimated parameters and the goodness of fit
#' measures.
#'
#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{rodriguez.2021}{TOU}
#'
#' @examples
#'
#' # an example with one trajectory of experimental adsorption
#' observed.time <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200,
#'                    250, 300)
#' observed.values <- c(0, 1.684, 2.341, 2.581, 2.842, 2.890, 2.959, 3.042,
#'                      3.083, 3.043, 3.017, 2.954, 2.996, 2.886, 2.844)
#' observed.process <- cbind(observed.time, observed.values)
#' # fitting the model without any fixed parameters
#' result <- fit.model.tou(w=observed.process)
#' print(result)
#' # fitting the model with a fixed value for the parameter qe
#' result <- fit.model.tou(w=observed.process, qe=2.9)
#' print(result)
# #' # fitting the model with a fixed value for the parameter kn
# #' result <- fit.model.tou(w=observed.process, kn=0.07)
# #' print(result)
# #' # fitting the model with a fixed value for the parameter n
# #' result <- fit.model.tou(w=observed.process, n=1.2)
# #' print(result)
# #' # fitting the model with fixed values for the parameters qe and kn
# #' result <- fit.model.tou(w=observed.process, qe=2.9, kn=0.07)
# #' print(result)
#' # fitting the model with fixed values for the parameters qe and n
#' result <- fit.model.tou(w=observed.process, qe=2.9, n=1.2)
#' print(result)
# #' # fitting the model with fixed values for the parameters kn and n
# #' result <- fit.model.tou(w=observed.process, kn=0.07, n=1.2)
# #' print(result)
# #' # fitting the model with fixed values for the parameters qe, kn and n
# #' result <- fit.model.tou(w=observed.process, qe=2.9, kn=0.07, n=1.2)
# #' print(result)
#'
#' # an example with three trajectories of experimental adsorption
#' observed.time <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200,
#'                    250, 300)
#' observed.values.1 <- c(0, 1.684, 2.341, 2.581, 2.842, 2.890, 2.959, 3.042,
#'                        3.083, 3.043, 3.017, 2.954, 2.996, 2.886, 2.844)
#' observed.values.2 <- c(0, 1.618, 2.217, 2.571, 2.763, 2.841, 2.866, 2.898,
#'                        2.935, 2.973, 2.906, 2.919, 2.910, 3.012, 3.071)
#' observed.values.3 <- c(0, 1.596, 2.333, 2.611, 2.750, 2.731, 2.829, 2.838,
#'                        2.864, 2.884, 2.886, 2.911, 2.896, 2.877, 2.969)
#' observed.processes <- cbind(observed.time, observed.values.1,
#'                             observed.values.2, observed.values.3)
# #' # fitting the model without any fixed parameters
# #' result <- fit.model.tou(w=observed.processes)
# #' print(result)
# #' # fitting the model with a fixed value for the parameter qe
# #' result <- fit.model.tou(w=observed.processes, qe=2.95)
# #' print(result)
# #' # fitting the model with a fixed value for the parameter kn
# #' result <- fit.model.tou(w=observed.processes, kn=0.07)
# #' print(result)
#' # fitting the model with a fixed value for the parameter n
#' result <- fit.model.tou(w=observed.processes, n=1.21)
#' print(result)
# #' # fitting the model with fixed values for the parameters qe and kn
# #' result <- fit.model.tou(w=observed.processes, qe=2.95, kn=0.07)
# #' print(result)
# #' # fitting the model with fixed values for the parameters qe and n
# #' result <- fit.model.tou(w=observed.processes, qe=2.95, n=1.21)
# #' print(result)
# #' # fitting the model with fixed values for the parameters kn and n
# #' result <- fit.model.tou(w=observed.processes, kn=0.07, n=1.21)
# #' print(result)
#' # fitting the model with fixed values for the parameters qe, kn and n
#' result <- fit.model.tou(w=observed.processes, qe=2.95, kn=0.07, n=1.21)
#' print(result)
#'
#' @export
fit.model.tou <- function(w, qe=NULL, kn=NULL, n=NULL) {

  ll.process <- make.loglikelihood.tou(y=w, qe=qe, kn=kn, n=n)
  ll.process.neg <- function(x) { value <- -1*ll.process(x) }

  fixed <- c(!is.null(qe), !is.null(kn), !is.null(n), FALSE)
  null.size <- sum(c(is.null(qe), is.null(kn), is.null(n)))

  m <- ncol(w)
  qe.aux <- max(w[,2:m])

  qe.min <- 0.5*qe.aux
  lambda.min <- 1e-8
  a.min <- 0
  tau.min <- 1e-8
  aux <- c(qe.min, lambda.min, a.min, tau.min)
  lower <- aux[!fixed]

  qe.max <- 2*qe.aux
  a.max <- 1000
  lambda.max <- a.max
  tau.max <- max(c((5*stats::sd(as.matrix(w[,2:m])) / sqrt(2*lambda.max)), 100))
  aux <- c(qe.max, lambda.max, a.max, tau.max)
  upper <- aux[!fixed]

  repeticiones <- 3
  resultados <- data.frame(matrix(0, nrow=repeticiones, ncol=10))
  colnames(resultados) <- c("qe0", "lambda0", "a0", "tau0", "logvero0",
                            "qe", "lambda", "a", "tau", "logvero")
  aux <- rep(0,4)
  for(i in 1:repeticiones) {
    DE.salida <- DEoptim::DEoptim(fn=ll.process.neg, lower=lower, upper=upper,
                         control=DEoptim::DEoptim.control(trace=FALSE))
    aux[!fixed] <- as.numeric(DE.salida$optim$bestmem)
    if(null.size==3) {
      resultados$qe0[i] <- aux[1]
      resultados$lambda0[i] <- aux[2]
      resultados$a0[i] <- aux[3]
    } else if(null.size==2) {
      if(!is.null(qe)) {
        resultados$qe0[i] <- qe
        resultados$lambda0[i] <- aux[2]
        resultados$a0[i] <- aux[3]
      } else if(!is.null(kn)) {
        resultados$qe0[i] <- aux[1]
        resultados$lambda0[i] <- lambda.root(y=c(aux[1], kn, aux[3]))
        resultados$a0[i] <- aux[3]
      } else {
        resultados$qe0[i] <- aux[1]
        resultados$lambda0[i] <- aux[2]
        resultados$a0[i] <- aux[2]*(n-1)
      }
    } else if(null.size==1) {
      if(!is.null(qe) & !is.null(kn)) {
        resultados$qe0[i] <- qe
        resultados$lambda0[i] <- lambda.root(y=c(qe, kn, aux[3]))
        resultados$a0[i] <- aux[3]
      } else if(!is.null(qe) & !is.null(n)) {
        resultados$qe0[i] <- qe
        resultados$lambda0[i] <- aux[2]
        resultados$a0[i] <- aux[2]*(n-1)
      } else {
        resultados$qe0[i] <- aux[1]
        resultados$lambda0[i] <- kn * (aux[1]^(n-1))
        resultados$a0[i] <- resultados$lambda0[i]*(n-1)
      }
    } else {
      resultados$qe0[i] <- qe
      resultados$lambda0[i] <- kn * (qe^(n-1))
      resultados$a0[i] <- resultados$lambda0[i]*(n-1)
    }
    resultados$tau0[i] <- aux[4]
    resultados$logvero0[i] <- ll.process(DE.salida$optim$bestmem)
    salida <- stats::optim(par=as.numeric(DE.salida$optim$bestmem), fn=ll.process,
                    method="L-BFGS-B", lower=lower, upper=upper,
                    control=list(fnscale=-1))
    aux[!fixed] <- salida$par
    if(null.size==3) {
      resultados$qe[i] <- aux[1]
      resultados$lambda[i] <- aux[2]
      resultados$a[i] <- aux[3]
    } else if(null.size==2) {
      if(!is.null(qe)) {
        resultados$qe[i] <- qe
        resultados$lambda[i] <- aux[2]
        resultados$a[i] <- aux[3]
      } else if(!is.null(kn)) {
        resultados$qe[i] <- aux[1]
        resultados$lambda[i] <- lambda.root(y=c(aux[1], kn, aux[3]))
        resultados$a[i] <- aux[3]
      } else {
        resultados$qe[i] <- aux[1]
        resultados$lambda[i] <- aux[2]
        resultados$a[i] <- aux[2]*(n-1)
      }
    } else if(null.size==1) {
      if(!is.null(qe) & !is.null(kn)) {
        resultados$qe[i] <- qe
        resultados$lambda[i] <- lambda.root(y=c(qe, kn, aux[3]))
        resultados$a[i] <- aux[3]
      } else if(!is.null(qe) & !is.null(n)) {
        resultados$qe[i] <- qe
        resultados$lambda[i] <- aux[2]
        resultados$a[i] <- aux[2]*(n-1)
      } else {
        resultados$qe[i] <- aux[1]
        resultados$lambda[i] <- kn * (aux[1]^(n-1))
        resultados$a[i] <- resultados$lambda[i]*(n-1)
      }
    } else {
      resultados$qe[i] <- qe
      resultados$lambda[i] <- kn * (qe^(n-1))
      resultados$a[i] <- resultados$lambda[i]*(n-1)
    }
    resultados$tau[i] <- aux[4]
    resultados$logvero[i] <- salida$value
  }
  # Selection of the best estimated parameters.
  posicion <- which.max(resultados$logvero)

  qe <- resultados$qe[posicion]
  lambda <- resultados$lambda[posicion]
  a <- resultados$a[posicion]
  tau <- resultados$tau[posicion]
  logLike <- resultados$logvero[posicion]

  # when qe, kn and n are fixed, their values must be assigned
  n <- 1 + (a/lambda)
  t.reach.qe <- (1/(2*lambda)) * log(1+((qe/(tau*stats::qnorm(0.25)))^2))
  if(a==0) {
    kn <- lambda
    t.half.qe <- log(2) / lambda
  } else {
    kn <- lambda / (qe^(n-1))
    t.half.qe <- (2^(n-1)-1) / (lambda*(n-1))
    t.reach.qe <- (exp(lambda*(n-1)*t.reach.qe) - 1) / (lambda*(n-1))
  }
  s <- tau*sqrt(2*lambda)

  m <- ncol(w)
  y <- w[,2:m]
  y.media <- mean(as.matrix(y))
  y.hat <- mean.tou(x=w[,1], parameters = c(qe, lambda, a, tau))
  SST <- sum((y - y.media)^2)
  SSR <- sum((y - y.hat)^2)
  R2 <- 1 - SSR / SST
  # the first row is eliminated because to adsorption values equal to zero
  if(2 < m) {
    MSE <- SSR / (nrow(y)*(m-1))
    delta.q <- 100*sqrt( sum( ( ( y[-1,] - y.hat[-1] ) / y[-1,] )^2 ) /
                           ( ( nrow(y) - 1 )*( m - 1 ) - 1) )
  } else {
    MSE <- SSR / (length(y)*(m-1))
    delta.q <- 100*sqrt( sum( ( ( y[-1] - y.hat[-1] ) / y[-1] )^2 ) /
                           ( ( length(y) - 1 )*( m - 1 ) - 1) )
  }
  return(list(qe=qe, n=n, kn=kn, a=a, lambda=lambda, sigma=s, tau=tau,
              t.half.qe=t.half.qe, t.reach.qe=t.reach.qe, logLikelihood=logLike,
              R2=R2, MSE=MSE, delta.q=delta.q, fitted.adsorption=y.hat,
              observed.data=w))
}

