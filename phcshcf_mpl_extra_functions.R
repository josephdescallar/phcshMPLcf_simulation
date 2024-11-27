#Function which generates cure fraction data
phcshcf.gen.data <- function(no.of.sims=1000, cureMin = 0.5, cureMax = 4.5,
                             pi.rå.z.coef = c(3, 1, -0.5), n.obs = 1000, gammaL = 0.5, gammaR = 2.5, 
                             beta = list(c(0.5,-2), c(1,-1)), prob_event){
  simdata <- list()
  for(s in 1:no.of.sims){
    lambda = c(1, 0.5)
    rho = 3
    seed <- s
    set.seed(seed)
    #Generate probability of cured
    z1 <- rnorm(n.obs)
    z2 <- rnorm(n.obs)
    pi.U <- runif(n.obs)
    pi.z = exp(pi.r.z.coef[1] + pi.r.z.coef[2]*z1 + pi.r.z.coef[3]*z2) /
      (1 + exp(pi.r.z.coef[1] + pi.r.z.coef[2]*z1 + pi.r.z.coef[3]*z2))
    noncure <- pi.U < pi.z
    #Make Z and Z covariates the same
    #x1 <- z1
    #x2 <- z2
    #make Z and X difference
    x1 <- rnorm(n.obs)
    x2 <- rnorm(n.obs)
    z <- data.matrix(data.frame(z1,z2,noncure))
    x <- data.matrix(data.frame(x1,x2,noncure))
    #split z and x covariates into cured and non-cured
    z.noncured <- z[noncure==1,]
    z.cured <- z[noncure==0,]
    n.cured <- nrow(z.cured)
    n.noncured <- nrow(z.noncured)
    x.noncured <- x[noncure==1,]
    x.cured <- x[noncure==0,]
    #for cured, generate censored observation times
    Ul.cured <- runif(n.cured)
    Ur.cured <- runif(n.cured,Ul.cured,1)
    # Uc <- gammaR*Ur.cured
    #Uc <- rexp(n.cured, 1/2)
    Uc <- runif(n.cured, cureMin, cureMax)
    Uc.time2 <- NA
    risk <- NA
    event <- 0
    cured <- cbind(z.cured, "time"=Uc, "time2"=Uc.time2, risk, event, 
                   x.cured[, -ncol(x.cured)])
    #for noncured, generate observation times
    u.noncured <- runif(n.noncured)
    t.noncured <- (-log(u.noncured) / 
                     (lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]])
                      + lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])))^(1/rho)
    Ul.noncured <- runif(n.noncured)
    Ur.noncured <- runif(n.noncured,Ul.noncured,1)
    Ue.noncured <- runif(n.noncured)
    event <- ifelse(Ue.noncured < prob_event,1,
                    ifelse(gammaL*Ul.noncured <= t.noncured &
                             t.noncured <= gammaR*Ur.noncured,3,ifelse(t.noncured < gammaL*Ul.noncured,
                                                                       2,0)))
    p1 <- (lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]])) / 
      ((lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]]) +
          lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])))
    p2 <- (lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])) / 
      ((lambda[1]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[1]]) +
          lambda[2]*exp(x.noncured[, -ncol(x.noncured)] %*% beta[[2]])))
    v <- runif(n.noncured)
    time1 <- ifelse(event==0,gammaR*Ur.noncured,ifelse(event==1,t.noncured,
                                                       ifelse(event==2,gammaL*Ul.noncured,gammaL*Ul.noncured)))
    time2 <- ifelse(event==0,NA,ifelse(event==1,t.noncured,ifelse(event==2,NA,
                                                                  gammaR*Ur.noncured)))
    risk.noncured <- ifelse(event==0,NA,ifelse(v <= p1, 1, 2))
    noncured <- data.frame(cbind(z.noncured, "time"=time1, time2, 
                                 "risk"=risk.noncured, event, x.noncured[, -ncol(x.noncured)]))
    #combine cured and noncured data
    simdata[[s]] <- data.frame(rbind(cured, noncured))
    simdata[[s]]$tmid <- ifelse(simdata[[s]]$event==0, simdata[[s]]$time,
                                ifelse(simdata[[s]]$event==1, simdata[[s]]$time,
                                       ifelse(simdata[[s]]$event==2, simdata[[s]]$time/2,
                                              (simdata[[s]]$time +simdata[[s]]$time2) / 2)))
    simdata[[s]]$risk_1 <- ifelse(is.na(simdata[[s]]$risk), 0,
                                  ifelse(simdata[[s]]$risk==1, 1, 0))
    simdata[[s]]$risk_2 <- ifelse(is.na(simdata[[s]]$risk), 0,
                                  ifelse(simdata[[s]]$risk==2, 1, 0))
  }
  simdata
}

#Function which summarizes MPL cure results from simulation
cf.mpl.results <- function(object, true.gamma, true.beta){
  #calculate bias
  all.gamma.list <- all.beta1.list <- all.beta2.list <- list()
  for(i in 1:n.sim){
    all.gamma.list[[i]] <- object[i,]$gamma
    all.beta1.list[[i]] <- object[i,]$beta[[1]]
    all.beta2.list[[i]] <- object[i,]$beta[[2]]
  }
  all.gamma <- Reduce("rbind", all.gamma.list)
  all.beta1 <- Reduce("rbind", all.beta1.list)
  all.beta2 <- Reduce("rbind", all.beta2.list)
  gamma.bias <- colMeans(all.gamma, na.rm=TRUE) - true.gamma
  beta1.bias <- colMeans(all.beta1, na.rm=TRUE) - true.beta[[1]]
  beta2.bias <- colMeans(all.beta2, na.rm=TRUE) - true.beta[[2]]
  all.gamma.std.list <- all.beta1.std.list <- all.beta2.std.list <- list()
  for(i in 1:n.sim){
    all.gamma.std.list[[i]] <- object[i,]$seG
    all.beta1.std.list[[i]] <- object[i,]$seB[[1]]
    all.beta2.std.list[[i]] <- object[i,]$seB[[2]]
  }
  all.gamma.std <- Reduce("rbind", all.gamma.std.list)
  all.beta1.std <- Reduce("rbind", all.beta1.std.list)
  all.beta2.std <- Reduce("rbind", all.beta2.std.list)
  gamma.std.mc <- apply(all.gamma, 2, sd)
  gamma.std.asymp <- colMeans(all.gamma.std, na.rm=TRUE)
  beta1.std.mc <- apply(all.beta1, 2, sd)
  beta1.std.asymp <- colMeans(all.beta1.std, na.rm=TRUE)
  beta2.std.mc <- apply(all.beta2, 2, sd)
  beta2.std.asymp <- colMeans(all.beta2.std, na.rm=TRUE)
  #coverage probability
  all.gamma.lower.list <- all.gamma.upper.list <- all.gamma.cp.list <- list()
  all.beta1.lower.list <- all.beta1.upper.list <- all.beta1.cp.list <- list()
  all.beta2.lower.list <- all.beta2.upper.list <- all.beta2.cp.list <- list()
  for(i in 1:n.sim){
    all.gamma.lower.list[[i]] <- object[i,]$gamma - 1.96*object[i,]$seG
    all.gamma.upper.list[[i]] <- object[i,]$gamma + 1.96*object[i,]$seG
    all.gamma.cp.list[[i]] <- ifelse(all.gamma.lower.list[[i]] < true.gamma & 
                                       true.gamma < all.gamma.upper.list[[i]], 1, 0)
    
    all.beta1.lower.list[[i]] <- object[i,]$beta[[1]] - 1.96*object[i,]$seB[[1]]
    all.beta1.upper.list[[i]] <- object[i,]$beta[[1]] + 1.96*object[i,]$seB[[1]]
    all.beta1.cp.list[[i]] <- ifelse(all.beta1.lower.list[[i]] < true.beta[[1]] & 
                                       true.beta[[1]] < all.beta1.upper.list[[i]], 1, 0)
    
    all.beta2.lower.list[[i]] <- object[i,]$beta[[2]] - 1.96*object[i,]$seB[[2]]
    all.beta2.upper.list[[i]] <- object[i,]$beta[[2]] + 1.96*object[i,]$seB[[2]]
    all.beta2.cp.list[[i]] <- ifelse(all.beta2.lower.list[[i]] < true.beta[[2]] & 
                                       true.beta[[2]] < all.beta2.upper.list[[i]], 1, 0)
  }
  all.gamma.lower <- Reduce("rbind", all.gamma.lower.list)
  all.gamma.upper <- Reduce("rbind", all.gamma.upper.list)
  all.gamma.cp <- Reduce("rbind", all.gamma.cp.list) 
  all.beta1.lower <- Reduce("rbind", all.beta1.lower.list)
  all.beta1.upper <- Reduce("rbind", all.beta1.upper.list)
  all.beta1.cp <- Reduce("rbind", all.beta1.cp.list)
  all.beta2.lower <- Reduce("rbind", all.beta2.lower.list)
  all.beta2.upper <- Reduce("rbind", all.beta2.upper.list)
  all.beta2.cp <- Reduce("rbind", all.beta2.cp.list)
  gamma.cp <- colMeans(all.gamma.cp, na.rm=TRUE)
  beta1.cp <- colMeans(all.beta1.cp, na.rm=TRUE)
  beta2.cp <- colMeans(all.beta2.cp, na.rm=TRUE)
  gamma.list <- list("bias"=gamma.bias, "std.asymp"=gamma.std.asymp,
                     "std.mc"=gamma.std.mc, "gamma.cp"=gamma.cp)
  beta1.list <- list("bias"=beta1.bias, "std.asymp"=beta1.std.asymp, 
                     "std.mc"=beta1.std.mc, "cp"=beta1.cp)
  beta2.list <- list("bias"=beta2.bias, "std.asymp"=beta2.std.asymp, 
                     "std.mc"=beta2.std.mc, "cp"=beta2.cp) 
  rlist <- list("gamma"=gamma.list, "beta1"=beta1.list, "beta2"=beta2.list)
  return(rlist)
}

#Function to display bias, std, cp for cox regression results
sim.results.cox <- function(object, beta, simdata){
  betalist = seBlist = biaslist = mselist = lowerlist = upperlist = list()
  for(sim in 1:simdata){
    betalist[[sim]] = object[[sim]]$coefficients
    biaslist[[sim]] = beta - betalist[[sim]]
    mselist[[sim]] = biaslist[[sim]]^2
    seBlist[[sim]] = sqrt(diag(object[[sim]]$var))
    lowerlist[[sim]] = betalist[[sim]] - 1.96*seBlist[[sim]]
    upperlist[[sim]] = betalist[[sim]] + 1.96*seBlist[[sim]]
  }
  lower = Reduce("rbind", lowerlist)
  upper = Reduce("rbind", upperlist)
  betamat = matrix(beta, nrow = simdata, ncol = length(beta), byrow = TRUE)
  fit = list()
  fit$bias = colMeans(Reduce("rbind",biaslist))
  fit$seB = colMeans(Reduce("rbind",seBlist))
  fit$seBmc = apply(Reduce("rbind",betalist), 2, sd)
  fit$mse = colMeans(Reduce("rbind",mselist))
  fit$cov.prob = colMeans(ifelse(betamat >= lower & betamat <= upper,1,0))
  fit
}

# Function which provides predictions for baseline hazards
pred.bh <- function(object, t.points, r, sand = FALSE, rho=0, lambda=0, 
         var.max=99999, coef){
  psif <- function(x, bknots, iknots){
    mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
            intercept = object$basis.intercept)}
  pred.psi = psif(t.points, object$b.knots, object$i.knots[[r]])
  pred.h0r = as.vector(pred.psi %*% object$"theta"[[r]]) 
  if(sand == FALSE){
    VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                       object$theta.index[[r]]]
  }
  else if(sand == TRUE){ #up to here
    VarCovMat.theta = object$sand[object$theta.index[[r]], object$theta.index[[r]]]
  }
  pred.h0r.var = diag(pred.psi %*% VarCovMat.theta %*% t(pred.psi))
  if(object$pos.def==1){
    pred.var = pred.h0r.var / pred.h0r^2
    pred.var[pred.var > var.max] <- var.max
    pred.h0r.lower = as.vector(pred.h0r) * exp(-1.96*sqrt(pred.var))
    pred.h0r.upper = as.vector(pred.h0r) * exp(1.96*sqrt(pred.var))
  }
  true.h0r <- lambda*rho*t.points^(rho-1)
  if(object$valid==1){
    rlist <- list("pred.h0r"=pred.h0r, "pred.h0r.lower"=pred.h0r.lower,
                  "pred.h0r.upper"=pred.h0r.upper, "true.h0r"=true.h0r)
  }
  else{
    rlist <- list("pred.h0r"=pred.h0r, "true.h0r"=true.h0r)
  }
  rlist
}

# Functions used for integrand of cumulative incidence functions plots
integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho)) 
}