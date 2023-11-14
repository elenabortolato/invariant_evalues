# original data for the bivariate Normal example


yy=scan()
2.011697937521708
0.7921920599823331
4.059964394271654
3.849935432214579
-0.3828757347661147
-1.786989437099709
8.5493200729399
9.626909979428621
1.9912225134779111
0.9156626601532678
1.9002411714222933
2.911700390583553
3.153415988632702
3.0740093081191184
4.503721534524546
5.272374068467039
-4.302161635928462
-4.2442463064297336
1.917966000959758
-0.0838248867455969
-4.45153059816264
-4.089786111578896
-5.666003038637369
-6.438208447056342
2.8330917054120945
0.6309802519006878
0.4088046761713947
0.7905085299848531
2.3452102292002084
3.9888037705431625
1.1728355820583325
0.21972736291710904
2.5435955449112906
1.1879622915346262
-4.413827637998951
-5.161852901916431
-1.4233367585643513
-1.586462155251945
-3.765791656213326
-4.314847359807742

 
yy=matrix(yy, ncol=2, byrow = T)
cor(yy)
plot(yy)

xx=matrix(c(1.0,  3.55764,    1.80829,
            1.0 , 3.55764,    1.80829,
            1.0,  0.124464,  -0.122073,
            1.0,  0.124464,  -0.122073,
            1.0,  0.25296,   -0.104636,
            1.0,  0.25296,   -0.104636,
            1.0,  0.709986,   0.583007,
            1.0,  0.70998,   0.583007,
            1.0,  1.12412,    1.28759,
            1.0,  1.12412,    1.28759,
            1.0,   1.03211,     1.2298,
            1.0,   1.03211,     1.22982,
            1.0,   1.23229,     0.578231,
            1.0,   1.23229 ,    0.578231,
            1.0,  -0.265314,   -0.355591,
            1.0,  -0.265314,   -0.355591,
            1.0,   0.958524,    1.06946,
            1.0 ,  0.958524,    1.06946,
            1.0,   0.0377049,   0.882677,
            1.0,   0.0377049,   0.882677,
            1.0,   0.564627,   2.21928,
            1.0 ,  0.564627,   2.21928,
            1.0  , 0.474404,   1.29291,
            1.0   ,0.474404,   1.29291,
            1.0,  -0.847348,   0.968852,
            1.0,  -0.847348,   0.968852,
            1.0,   0.85870,   1.31583,
            1.0,   0.85870,   1.31583,
            1.0 ,  0.23109,   -1.16238,
            1.0  , 0.23109,   -1.16238,
            1.0,   1.73201,    0.109923,
            1.0,   1.73201,    0.109923,
            1.0,  -0.198519,  -0.0882185,
            1.0,  -0.198519,  -0.0882185,
            1.0,   2.82557,    1.70376,
            1.0,   2.82557,    1.70376,
            1.0,   1.14067,    1.14333,
            1.0,   1.14067,    1.14333,
            1.0,   2.45737,    1.14838,
            1.0,   2.45737,    1.14838), ncol=3, byrow=T)
xx
idd=rep(1:20,2)
idd


#define function for maximizing the modified likelihood:

breqn <- function(x, y, id, method = "ML", a = 0.5,
                  maxit = 1000, epsilon = 10^(-11), max_step_factor = 11, slowit = 1) {
  x <- as.matrix(x)
  p <- ncol(x)
  n <- length(table(id))
  q <- as.vector(table(id)[1])
  betas <-  solve(crossprod(x,x)) %*% (t(x) %*% y)
  
  mus <- drop(x %*% betas)
  residuals <- y - mus
  residualsmat <- matrix(residuals, nrow = n, byrow = T)
  globalmean <- mean(residuals)
  SW <- apply(residualsmat, 1, function (x) sum((x - mean(x))^2))
  SSW <- sum(SW)
  subjectmeans <- apply(residualsmat, 1, mean)
  SSB <-  sum((subjectmeans - globalmean)^2)
  sigma2 <- (SSW + q * (SSB))/(n * q)
  rho <- (1 - q/(q - 1)*SSW/(SSW + q*SSB) )
  
  if (rho <= -1/(q-1)) rho <- 0; flag = 1
  pars <- c(sigma2, rho) 
  
  ## score function of sigma2 and rho
  score <- function(pars, y, mus) {
    sigma2 <- pars[1]
    rho <- pars[2]
    omega_s2 <- matrix( -rho/((rho - 1) * (q * rho - rho + 1) * sigma2^2), 
                        ncol = q, nrow = q)
    diag(omega_s2) <-  (q * rho - 2 * rho + 1)/((rho - 1) * (q * rho - rho + 1) * sigma2^2)
    omega_rho <- matrix(-(q * rho^2 - rho^2 + 1)/((rho - 1)^2 * (q * rho - rho + 1)^2 * sigma2) ,
                        ncol = q, nrow = q)
    diag(omega_rho) <- ((q - 1) * rho * (q * rho - 2 * rho + 2))/((rho - 1)^2 * (q * rho - rho + 1)^2 * sigma2)
    
    c( ( -n*q/(2 * sigma2) -1/2 * sum(t(y - mus) %*% ( diag(n)%x%omega_s2) %*% (y - mus)) ),
       (-(n * (q - 1) * q * rho)/(2 * (rho - 1) * (q * rho - rho + 1))-1/2 * sum(t(y - mus) %*% ( diag(n)%x%omega_rho) %*% (y - mus)) )
    )
  }
  
  ## Fisher information
  # information <- function(pars, level = 0, inverse = FALSE) {
  #   sigma2 <- pars[1]
  #   rho <- pars[2]
  #   
  #   wx <- sqrt(q/((q * rho - rho + 1) * sigma2)) * x
  #   qr_decomposition <- qr(wx)
  #   if (level == 0) {
  #     R_matrix <- qr.R(qr_decomposition)
  #     if (inverse) {
  #       return(chol2inv(R_matrix))
  #     }
  #     else {
  #       return(crossprod(R_matrix))
  #     }
  #   }
  #   if (level == 1) {
  #     ans <- matrix(NA, 2, 2)
  #     if (inverse) {
  #       ans[1,1] <- ( (2 * (q * rho^2 - rho^2 + 1) * sigma2^2)/(n * q) )
  #       ans[1,2] <- ans[2,1] <- ( -(2 * (rho - 1) * rho * (q * rho - rho + 1) * sigma2)/(n * q) ) 
  #       ans[2,2] <- ( (2 * (rho - 1)^2 * (q * rho - rho + 1)^2)/(n * (q - 1) * q) )
  #       return(ans)
  #     }
  #     else {
  #       ans[1,1] <- n*q/(2 * sigma2^2)
  #       ans[1,2] <- ans[2,1] <- ( (n * q * (q - 1) * rho)/(2 * (rho - 1) * (q * rho - rho + 1) * sigma2) )
  #       ans[2,2] <- ( n * (q - 1) * q * (q * rho^2 - rho^2 + 1)/(2 * (rho - 1)^2 * (q * rho - rho + 1)^2) )
  #       return(ans) 
  #     }
  #   }
  # }
  
  information <- function(pars, level = 0, inverse = FALSE) {
    sigma2 <- pars[1]
    rho <- pars[2]
    
    if (level == 0) {
      Iq <- diag(q)
      ones <- matrix(1, q, q)
      omega <- (1/(sigma2 * (1 - rho)))*(Iq - ones * rho/(1 + rho * (q - 1)))
      if (inverse) {
        return(solve(t(x) %*%  (diag(n)%x%omega) %*% x))
      }
      else {
        return((t(x) %*%  (diag(n)%x%omega) %*% x))
      }
    }
    if (level == 1) {
      ans <- matrix(NA, 2, 2)
      if (inverse) {
        ans[1,1] <- ( (2 * (q * rho^2 - rho^2 + 1) * sigma2^2)/(n * q) )
        ans[1,2] <- ans[2,1] <- ( -(2 * (rho - 1) * rho * (q * rho - rho + 1) * sigma2)/(n * q) ) 
        ans[2,2] <- ( (2 * (rho - 1)^2 * (q * rho - rho + 1)^2)/(n * (q - 1) * q) )
        return(ans)
      }
      else {
        ans[1,1] <- n*q/(2 * sigma2^2)
        ans[1,2] <- ans[2,1] <- ( (n * q * (q - 1) * rho)/(2 * (rho - 1) * (q * rho - rho + 1) * sigma2) )
        ans[2,2] <- ( n * (q - 1) * q * (q * rho^2 - rho^2 + 1)/(2 * (rho - 1)^2 * (q * rho - rho + 1)^2) )
        return(ans) 
      }
    }
  }
  
  # invInfo2 <- function(pars)
  # {
  #   sigma2 <- pars[1]
  #   rho <- pars[2]
  #   Iq <- diag(q)
  #   ones <- matrix(1, q, q)
  #   omega <- (1/(sigma2 * (1 - rho)))*(Iq - ones * rho/(1 + rho * (q - 1)))
  #   solve(t(x) %*%  (diag(n)%x%omega) %*% x)
  # }
  
  ## Adjustment functions
  ## adjustment with Jeffreys
  as_Jeffreys_adjustment <- function(pars) {
    sigma2 <- pars[1]
    rho <- pars[2]
    a*c(
      c(-(2 + p)/sigma2,
        -((4 + p) *q * rho - (4 + p) * rho -(2 + p) * q + (4 + p))/((rho - 1) * (q * rho - rho + 1)))
    )
  }
  ## mean adjustment
  as_mean_adjustment <- function(pars) {
    
    sigma2 <- pars[1]
    rho <- pars[2]
    
    w0od=rho/((rho-1)*(q*rho-rho+1)*sigma2)
    w0d=-(q*rho-2*rho+1)/((rho-1)*(q*rho-rho+1)*sigma2)
    Omega=matrix(w0od,ncol=q, nrow=q)-diag(w0od-w0d, nrow = q, ncol=q)
    w_rho_d=((q-1)*rho*(q*rho-2*rho+2))/((rho-1)^2*(q*rho-rho+1)^2*sigma2)
    w_rho_od=-(q*rho^2-rho^2+1)/((rho-1)^2*(q*rho-rho+1)^2*sigma2)
    
    Omega_rho=matrix( w_rho_od,ncol=q,nrow=q)-diag(w_rho_od-w_rho_d, ncol=q, nrow=q)
    x.omega.rho.x=  (t(x) %*%  (diag(n)%x%Omega_rho) %*% x)
    x.omega.x=  solve(t(x) %*%  (diag(n)%x%Omega) %*% x)
    xprodx=x.omega.rho.x%*%x.omega.x
    c( 
      -(2 * q * rho^2 - 2 * rho^2 - p)/(2 * sigma2),
      #-(((q - 1) * (2 * q * rho^3 - 2 * rho^3 - p * rho + 2 * rho + p))/(2 * (rho - 1) * (q * rho - rho + 1)))
      -(q-1)*rho*(q*rho^2-rho^2+1)/((rho-1)*(q*rho-rho+1))-0.5*sum(diag(xprodx))
      
    )
  }
  
  ## median adjustment
  as_median_adjustment <- function(pars) {
    sigma2 <- pars[1]
    rho <- pars[2]
    # c(  
    #  ((3 * p * q * rho^2 + 6 * q * rho^2 - 3 * p * rho^2 - 6 * rho^2 - 2 * q * rho + 4 * rho + 3 * p + 2)/(6 * (q * rho^2 - rho^2 + 1) * sigma2)),
    #   ((3 * p * q^2 * rho^3 + 6 * q^2 * rho^3 -6 * p * q * rho^3 - 12 * q * rho^3 + 3 * p * rho^3 + 6 * rho^3 - 3 * p * q^2 * rho^2 - 4 * q^2 * rho^2 +  
    #      6 * p * q * rho^2 + 12 * q * rho^2 - 3 * p * rho^2 - 8 * rho^2 + 3 * p * q * rho + 2 * q * rho - 3 * p * rho - 2 * rho - 3 * p * q - 2 * q + 3 * p + 4)/
    #     (6 * (rho - 1) * (q * rho - rho + 1) * (q * rho^2 - rho^2 + 1)))
    #)
    as_mean_adjustment(pars)-c( 
      -(3*q^2*rho^4-6*q*rho^4+3*rho^4+6*q*rho^2-6*rho^2-q*rho+2*rho+1)/(3*(q*rho^2-rho^2+1)*sigma2),
      -(3*q^3*rho^5-9*q^2*rho^5+9*q*rho^5-3*rho^5+9*q^2*rho^3-18*q*rho^3+9*rho^3-2*q^2*rho^2+6*q*rho^2
        -4*rho^2+4*q*rho-4*rho-q+2)/(3*(rho-1)*(q*rho-rho+1)*(q*rho^2-rho^2+1))
    )
  }
  
  ## Bias ###
  
  bias <- function(pars) {
    sigma2 <- pars[1]
    rho <- pars[2]
    c(
      -(p * (q * rho - rho + 1) * sigma2)/(n * q),
      ((rho - 1) * (2 * rho + p) * (q * rho - rho + 1))/(n * q)
    )
  }
  
  ## adjustment term function ##
  
  adjustment_function <- switch(method, 
                                AS_mean = as_mean_adjustment, 
                                AS_median = as_median_adjustment,
                                MPL_Jeffreys = as_Jeffreys_adjustment,
                                correction = function(pars, ...) 0,
                                ML = function(pars, ...) 0)
  
  inverse_info <- information(pars, level = 1, inverse = TRUE)
  adjusted_grad <- score(pars, y, mus) + adjustment_function(pars)
  step <- drop(inverse_info %*% adjusted_grad)
  if (maxit == 0) {
    iter <- 0
    failed <- FALSE
  }
  else 
  {
    for (iter in seq.int(maxit)) {
      step_factor <- 0
      testhalf <- TRUE
      while (testhalf & step_factor < max_step_factor) {
        step_previous <- step
        pars <- pars + slowit * 2^(-step_factor) * step
        
        inverse_info <- information(pars, level = 1, inverse = TRUE)
        adjusted_grad <- score(pars, y, mus) + adjustment_function(pars)
        step <- drop(inverse_info %*% adjusted_grad)
        
        if (step_factor == 0 & iter == 1) {
          testhalf <- TRUE
        }
        else {
          testhalf <- sum(abs(step), na.rm = TRUE) > sum(abs(step_previous),  na.rm = TRUE)
        }
        step_factor <- step_factor + 1
      }
      if ( all(abs( adjusted_grad ) < epsilon)) {
        break
      }
    }
  }
  if(iter >= maxit) {
    convergence <- FALSE
    warning("optimization failed to converge")
  }
  else
  {
    convergence <- TRUE
  }
  if (method == "correction") 
  {
    pars <- pars - bias(pars)
  }
  
  inverse_info_betas <- information(pars, level = 0, inverse = TRUE)
  inverse_info_sigma2rho <- information(pars, level = 1, inverse = TRUE)
  se <- c(sqrt(diag(inverse_info_betas)),
          sqrt(diag(inverse_info_sigma2rho)))
  
  list(estimate = c(drop(betas),pars[1],pars[2]), sigma2 = pars[1], rho = pars[2], se = se, converged = convergence,
       iter =iter, grad = c(rep(0,p),adjusted_grad ), method = method)
  
}


# estimates
breqn(x = xx, as.vector(t(yy)),idd, method="AS_median",  maxit = 1000, epsilon = 1e-19
)
 

breqn(x = xx, as.vector(t(yy)),idd, method="AS_mean",  maxit = 1000, epsilon = 1e-19
)
