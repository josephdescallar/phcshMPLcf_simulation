library(parallel)
library(MASS)
library(foreach)
library(doParallel)
library(doSNOW)
library(survival)
library(splines2)
library(gbm)
library(statmod)
n.sim = 1000

# Run phcshcf_mpl.R file

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

#sim11, n = 300, cm = 12
true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim11.data.n300 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, 
cureMax = 4.378516,pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
                                        beta = true.beta)

allt <- unlist(lapply(cf.sim11.data.n300, 
    function(a) na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), 
      ]$time2))))
max(allt)

#Determine censoring proportions
total.noncure <- lapply(cf.sim11.data.n300, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim11.data.n300, function(a) a$event))) /
  length(unlist(lapply(cf.sim11.data.n300, function(a) a$event)))

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim11n300.results.time <- system.time({cf.sim11n300.results <- foreach(
  i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim11.data.n300[[i]]$risk, data = cf.sim11.data.n300[[i]],
  z = as.matrix(cbind(cf.sim11.data.n300[[i]]$z1, cf.sim11.data.n300[[i]]$z2)), 
  max.outer = 10, max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, 
  iknots.pos = NULL,aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, 
  iter.disp = FALSE, gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim11n300.results, true.gamma, true.beta)

#figure 1
#caluclate values for plots
cf.sim11.data.n300.min <- max(sapply(cf.sim11.data.n300, function(a)
  min(na.omit(c(a$time, a$time2)))))
cf.sim11.data.n300.max <- min(sapply(cf.sim11.data.n300, function(a)
  max(na.omit(c(a$time, a$time2)))))
cf1.scen1.plot.t <- seq(cf.sim11.data.n300.min, cf.sim11.data.n300.max, 
                        length.out=1000) 
cf1.plot.bh.r1 <- cf1.plot.bh.r2 <- list()
for(sim in 1:1000){
  cf1.plot.bh.r1[[sim]] <- pred.bh(cf.sim11n300.results[sim,],
  cf1.scen1.plot.t, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.bh.r2[[sim]] <- pred.bh(cf.sim11n300.results[sim,],
  cf1.scen1.plot.t, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
}

cf1.plot.bh.r1.all <- lapply(cf1.plot.bh.r1, function(a) a$pred.h0r)
cf1.plot.bh.r1.h0r.mean <- colMeans(Reduce("rbind", cf1.plot.bh.r1.all))
cf1.plot.bh.r1.all.lower <- lapply(cf1.plot.bh.r1, function(a) a$pred.h0r.lower)
cf1.plot.bh.r1.h0r.mean.lower<- colMeans(Reduce("rbind", 
                                                cf1.plot.bh.r1.all.lower))
cf1.plot.bh.r1.all.upper <- lapply(cf1.plot.bh.r1, function(a) a$pred.h0r.upper)
cf1.plot.bh.r1.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                cf1.plot.bh.r1.all.upper)) 
cf1.bh.cp.r1 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r1, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,0
         )))) / sum(temp.valid)

cf1.plot.bh.r2.all <- lapply(cf1.plot.bh.r2, function(a) a$pred.h0r)
cf1.plot.bh.r2.h0r.mean <- colMeans(Reduce("rbind", cf1.plot.bh.r2.all))
cf1.plot.bh.r2.all.lower <- lapply(cf1.plot.bh.r2, function(a) a$pred.h0r.lower)
cf1.plot.bh.r2.h0r.mean.lower<- colMeans(Reduce("rbind", 
                                                cf1.plot.bh.r2.all.lower))
cf1.plot.bh.r2.all.upper <- lapply(cf1.plot.bh.r2, function(a) a$pred.h0r.upper)
cf1.plot.bh.r2.h0r.mean.upper<- colMeans(Reduce("rbind", 
                                                cf1.plot.bh.r2.all.upper))

cf1.bh.cp.r2 <- colSums(Reduce("rbind", lapply(cf1.plot.bh.r2, function(a) 
  ifelse((a$pred.h0r.lower < a$true.h0r) & (a$true.h0r < a$pred.h0r.upper),1,
         0)))) / sum(temp.valid)


par(mfrow=c(2,2))
plot(cf1.scen1.plot.t, cf1.plot.bh.r1.h0r.mean, type='l', col='black', 
     xlab='t', ylab='Baseline hazard', main="(A) Risk 1", ylim=c(0,100))
lines(cf1.scen1.plot.t, cf1.plot.bh.r1.h0r.mean.lower, col='black',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r1.h0r.mean.upper, col='black',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r1[[1]]$true.h0r, col = "grey", cex=10)
legend("topleft", legend=c("MPL", "True"),
       col=c("black","grey"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t, cf1.plot.bh.r2.h0r.mean, type='l', col='black',
     xlab='t', ylab='Baseline hazard', main="(B) Risk 2", ylim=c(0,100))
lines(cf1.scen1.plot.t, cf1.plot.bh.r2.h0r.mean.lower, col='black',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r2.h0r.mean.upper, col='black',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.bh.r2[[1]]$true.h0r, col = "grey", cex=10) 
legend("topleft", legend=c("MPL", "True"),
       col=c("black","grey"), lty=1, bty = 'n')


plot(cf1.scen1.plot.t, cf1.bh.cp.r1, ylab="Coverage probability", xlab="t", 
     col='black', type='l', ylim=c(0.5,1), main="(C) Risk 1")
lines(cf1.scen1.plot.t, rep(0.95, 1000), lty=2, col='grey')
legend("bottomleft", legend=c("MPL", "0.95"), col=c("black","grey"), lty=c(1,2),
       bty = 'n')

plot(cf1.scen1.plot.t, cf1.bh.cp.r2, ylab="Coverage probability", xlab="t", 
     col='black', type='l', ylim=c(0.5,1), main="(D) Risk 2")
lines(cf1.scen1.plot.t, rep(0.95, 1000), lty=2, col='grey')
legend("bottomleft", legend=c("MPL", "0.95"), col=c("black","grey"), lty=c(1,2),
       bty = 'n')

#CIF graph
cf1.plot.cif.r1 <- cf1.plot.cif.r2 <- list()
for(sim in 1:1000){
  if(cf.sim11n300.results[sim,]$valid==1){ 
  cf1.plot.cif.r1[[sim]] <- pred.CIF(cf.sim11n300.results[sim,],
  cf1.scen1.plot.t, r=1, sand = FALSE, rho=3, lambda=1, var.max=5)
  cf1.plot.cif.r2[[sim]] <- pred.CIF(cf.sim11n300.results[sim,],
  cf1.scen1.plot.t, r=2, sand = FALSE, rho=3, lambda=0.5, var.max=5)
  }
}
cf1.plot.cif.r1.all <- lapply(cf1.plot.cif.r1, function(a) a$pred.F0r)
cf1.plot.cif.r1.F0r.mean <- colMeans(Reduce("rbind", cf1.plot.cif.r1.all))
cf1.plot.cif.r1.all.lower <- lapply(cf1.plot.cif.r1, function(a) 
  a$pred.F0r.lower)
cf1.plot.cif.r1.F0r.mean.lower<- colMeans(Reduce("rbind", 
                                                 cf1.plot.cif.r1.all.lower))
cf1.plot.cif.r1.all.upper <- lapply(cf1.plot.cif.r1, 
                                    function(a) a$pred.F0r.upper)
cf1.plot.cif.r1.F0r.mean.upper<- colMeans(Reduce("rbind", 
                                                 cf1.plot.cif.r1.all.upper)) 

cf1.plot.cif.r2.all <- lapply(cf1.plot.cif.r2, function(a) a$pred.F0r)
cf1.plot.cif.r2.F0r.mean <- colMeans(Reduce("rbind", cf1.plot.cif.r2.all))
cf1.plot.cif.r2.all.lower <- lapply(cf1.plot.cif.r2, function(a) 
  a$pred.F0r.lower)
cf1.plot.cif.r2.F0r.mean.lower<- colMeans(Reduce("rbind", 
                                                 cf1.plot.cif.r2.all.lower))
cf1.plot.cif.r2.all.upper <- lapply(cf1.plot.cif.r2, 
                                    function(a) a$pred.F0r.upper)
cf1.plot.cif.r2.F0r.mean.upper<- colMeans(Reduce("rbind", 
                                                 cf1.plot.cif.r2.all.upper))

integrand.F01 <- function(x, lambda1, lambda2, rho){
  (lambda1*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho))
}

integrand.F02 <- function(x, lambda1, lambda2, rho){
  (lambda2*rho*x^(rho-1))*exp(-lambda1*(x^rho))*exp(-lambda2*(x^rho)) 
}
true.F0r.r1 <- sapply(cf1.scen1.plot.t, function(a) integrate(integrand.F01, 0, 
                                  a, rho=3,lambda1=1, lambda2 = 0.5)$value)
true.F0r.r2 <- sapply(cf1.scen1.plot.t, function(a) integrate(integrand.F02, 0, 
                                  a, rho=3, lambda1=1, lambda2=0.5)$value)
F0r.plot.mat <- sapply(cf1.plot.cif.r1, function(a) a$pred.F0r)
#calculate length of each list
lengthr1 <- lengthr2 <- list()
for(sim in 1:1000){
  lengthr1[[sim]] <- length(cf1.plot.cif.r1[[sim]]$pred.F0r)
  lengthr2[[sim]] <- length(cf1.plot.cif.r2[[sim]]$pred.F0r)
}

F0r.plot.mat.list <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r1[[sim]]$pred.F0r)==1000){ 
    F0r.plot.mat.list[[count]] <- cf1.plot.cif.r1[[sim]]$pred.F0r
  }
  count=count+1
}
F0r.plot.mat <- Reduce("cbind", F0r.plot.mat.list)
F0r.plot.mean <- apply(F0r.plot.mat, 1, mean)
F0r.plot.l95 <- apply(F0r.plot.mat, 1, quantile, 0.025)
F0r.plot.u95 <- apply(F0r.plot.mat, 1, quantile, 0.975)
F0r.plot.r2.mat.list <- list()
count=1
for(sim in 1:1000){
  if(length(cf1.plot.cif.r2[[sim]]$pred.F0r)==1000){
    F0r.plot.r2.mat.list[[count]] <- cf1.plot.cif.r2[[sim]]$pred.F0r
  }
  count=count+1 
}
F0r.plot.r2.mat <- Reduce("cbind", F0r.plot.r2.mat.list)
F0r.plot.r2.mean <- apply(F0r.plot.r2.mat, 1, mean)
F0r.plot.r2.l95 <- apply(F0r.plot.r2.mat, 1, quantile, 0.025)
F0r.plot.r2.u95 <- apply(F0r.plot.r2.mat, 1, quantile, 0.975)

par(mfrow=c(2,1))
plot(cf1.scen1.plot.t, cf1.plot.cif.r1.F0r.mean, type='l', col='grey',
     xlab='t', ylab='CIF', ylim=c(0,1), main="(A) Risk 1", xlim=c(0, 2))
lines(cf1.scen1.plot.t, cf1.plot.cif.r1.F0r.mean.lower, col='grey',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.cif.r1.F0r.mean.upper, col='grey',
      lty='dashed')
lines(cf1.scen1.plot.t, true.F0r.r1, col = "black", cex=10)
lines(cf1.scen1.plot.t, F0r.plot.l95, col='black',
      lty=2)
lines(cf1.scen1.plot.t, F0r.plot.u95, col='black', 
      lty=2)
legend("topleft", legend=c("MPL", "True"),
       col=c("grey","black"), lty=1, bty = 'n')

plot(cf1.scen1.plot.t, cf1.plot.cif.r2.F0r.mean, type='l', col='grey',
     xlab='t', ylab='CIF', ylim=c(0,1), main="(B) Risk 2", xlim=c(0, 2))
lines(cf1.scen1.plot.t, cf1.plot.cif.r2.F0r.mean.lower, col='grey',
      lty='dashed')
lines(cf1.scen1.plot.t, cf1.plot.cif.r2.F0r.mean.upper, col='grey',
      lty='dashed')
lines(cf1.scen1.plot.t, true.F0r.r2, col = "black", cex=10)
lines(cf1.scen1.plot.t, F0r.plot.r2.l95, col='black',
      lty=2)
lines(cf1.scen1.plot.t, F0r.plot.r2.u95, col='black',
      lty=2)
legend("topleft", legend=c("MPL", "True"), 
       col=c("grey","black"), lty=1, bty = 'n')



#create variables for Cox Regression application
for(i in 1:n.sim){
  cf.sim11.data.n300[[i]]$tmid <- ifelse(cf.sim11.data.n300[[i]]$event==0, 
  cf.sim11.data.n300[[i]]$time, ifelse(cf.sim11.data.n300[[i]]$event==1, 
  cf.sim11.data.n300[[i]]$time,ifelse(cf.sim11.data.n300[[i]]$event==2, 
  cf.sim11.data.n300[[i]]$time/2,
  (cf.sim11.data.n300[[i]]$time +cf.sim11.data.n300[[i]]$time2) / 2)))
  cf.sim11.data.n300[[i]]$risk_1 <- ifelse(is.na(cf.sim11.data.n300[[i]]$risk), 
  0,ifelse(cf.sim11.data.n300[[i]]$risk==1, 1, 0))
  cf.sim11.data.n300[[i]]$risk_2 <- ifelse(is.na(cf.sim11.data.n300[[i]]$risk), 
  0, ifelse(cf.sim11.data.n300[[i]]$risk==2, 1, 0)) 
}

cf.sim1.r1.n300.cox = cf.sim1.r2.n300.cox = list()
for(sim in 1:n.sim){
  cf.sim1.r1.n300.cox[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                     data = cf.sim11.data.n300[[sim]])
  cf.sim1.r2.n300.cox[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                     data = cf.sim11.data.n300[[sim]])
}
sim.results.cox(cf.sim1.r1.n300.cox, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.r2.n300.cox, beta = true.beta[[2]], simdata = n.sim)

#sim11, n = 1000, cm = 12
true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim11.data.n1000 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, 
cureMax = 4.939527,pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
beta = true.beta) 
total.noncure <- lapply(cf.sim11.data.n1000, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim11.data.n1000, function(a) a$event))) /
  length(unlist(lapply(cf.sim11.data.n1000, function(a) a$event)))

# allt <- unlist(lapply(cf.sim11.data.n1000, 
# function(a) na.omit(c(a[which(a$event != 0), ]$time,a[which(a$event != 0),
#                                                       ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim11n1000.results.time <- system.time({cf.sim11n1000.results <- foreach( 
  i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim11.data.n1000[[i]]$risk, data = cf.sim11.data.n1000[[i]],
  z = as.matrix(cbind(cf.sim11.data.n1000[[i]]$z1, 
  cf.sim11.data.n1000[[i]]$z2)), max.outer = 10, max.iter = 500, 
  lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
}) 
close(pb)
stopCluster(cl)

true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim11n1000.results, true.gamma, true.beta) 

cf.sim1.r1.n1000.cox = cf.sim1.r2.n1000.cox = list()
for(sim in 1:n.sim){
  cf.sim1.r1.n1000.cox[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                      data = cf.sim11.data.n1000[[sim]])
  cf.sim1.r2.n1000.cox[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                      data = cf.sim11.data.n1000[[sim]])
}
sim.results.cox(cf.sim1.r1.n1000.cox, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.r2.n1000.cox, beta = true.beta[[2]], simdata = n.sim)


#sim11, n = 300, cm = 12, 40% cured
true.gamma <- c(0.5, 1, -0.5) 
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n300.c40 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, 
   cureMax = 4.504286,pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
                                       beta = true.beta)
total.noncure <- lapply(cf.sim.data.n300.c40, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n300.c40, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n300.c40, function(a) a$event)))

allt <- unlist(lapply(cf.sim.data.n300.c40, 
function(a) na.omit(c(a[which(a$event != 0), ]$time,a[which(a$event != 0),
                                                      ]$time2))))

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3) #up to here
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n300.c40.time <- system.time({
cf.sim.results.n300.c40 <- foreach(i = 1:n.sim, .combine = rbind, 
.options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
   risk = cf.sim.data.n300.c40[[i]]$risk, data = cf.sim.data.n300.c40[[i]],
   z = as.matrix(cbind(cf.sim.data.n300.c40[[i]]$z1, 
   cf.sim.data.n300.c40[[i]]$z2)), max.outer = 10,  max.iter = 500, 
   lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE, 
   knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
   gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n300.c40, true.gamma, true.beta)

cf.sim1.c40.n300.cox.r1 = cf.sim1.c40.n300.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n300.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                  data = cf.sim.data.n300.c40[[sim]])
  cf.sim1.c40.n300.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                  data = cf.sim.data.n300.c40[[sim]])
}
sim.results.cox(cf.sim1.c40.n300.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c40.n300.cox.r2, beta = true.beta[[2]], simdata = n.sim)

#sim11, n = 1000, cm = 12, 40% cured
true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n1000.c40 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, 
cureMax = 5.507274, pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
 beta = true.beta)
total.noncure <- lapply(cf.sim.data.n1000.c40, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n1000.c40, function(a) a$event))) /
length(unlist(lapply(cf.sim.data.n1000.c40, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n1000.c40, 
# function(a) na.omit(c(a[which(a$event != 0), ]$time,a[which(a$event != 0), 
#                                                         ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n1000.c40.time <- system.time({cf.sim.results.n1000.c40 <- 
foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n1000.c40[[i]]$risk, data = cf.sim.data.n1000.c40[[i]],
  z = as.matrix(cbind(cf.sim.data.n1000.c40[[i]]$z1, 
 cf.sim.data.n1000.c40[[i]]$z2)), max.outer = 10, max.iter = 500, 
 lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL, aps = TRUE, 
 knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
               gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n1000.c40, true.gamma, true.beta)

cf.sim1.c40.n1000.cox.r1 = cf.sim1.c40.n1000.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n1000.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                          data = cf.sim.data.n1000.c40[[sim]])
  cf.sim1.c40.n1000.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                          data = cf.sim.data.n1000.c40[[sim]])
}
sim.results.cox(cf.sim1.c40.n1000.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c40.n1000.cox.r2, beta = true.beta[[2]], simdata = n.sim)

#sim11, n = 300, cm = 12, 60% cured
true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n300.c60 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, 
    cureMax = 4.672316, pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
    beta = true.beta)
total.noncure <- lapply(cf.sim.data.n300.c60, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n300.c60, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n300.c60, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n300.c60, 
#         function(a) na.omit(c(a[which(a$event != 0), ]$time,
#         a[which(a$event != 0), ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n300.c60.time <- system.time({cf.sim.results.n300.c60 <- 
foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n300.c60[[i]]$risk, data = cf.sim.data.n300.c60[[i]],
  z = as.matrix(cbind(cf.sim.data.n300.c60[[i]]$z1, 
  cf.sim.data.n300.c60[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n300.c60, true.gamma, true.beta)

cf.sim1.c60.n300.cox.r1 = cf.sim1.c60.n300.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n300.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                         data = cf.sim.data.n300.c60[[sim]])
  cf.sim1.c60.n300.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                         data = cf.sim.data.n300.c60[[sim]])
}
sim.results.cox(cf.sim1.c60.n300.cox.r1, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.c60.n300.cox.r2, beta = true.beta[[2]], simdata = n.sim)

#sim11, n = 1000, cm = 12, 60% cured
true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n1000.c60 <- phcshcf.gen.data(no.of.sims=n.sim, cureMin = 0.5, 
    cureMax = 4.553649,pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
    beta = true.beta)
total.noncure <- lapply(cf.sim.data.n1000.c60, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n1000.c60, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n1000.c60, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n1000.c60, 
#   function(a) na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), 
#                                                          ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores() - 1
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n1000.c60.time <- system.time({cf.sim.results.n1000.c60 <- 
  foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n1000.c60[[i]]$risk, data = cf.sim.data.n1000.c60[[i]],
  z = as.matrix(cbind(cf.sim.data.n1000.c60[[i]]$z1, 
  cf.sim.data.n1000.c60[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n1000.c60, true.gamma, true.beta)

cf.sim1.c60.n1000.cox.r1 = cf.sim1.c60.n1000.cox.r2 = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n1000.cox.r1[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                          data = cf.sim.data.n1000.c60[[sim]])
  cf.sim1.c60.n1000.cox.r2[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                          data = cf.sim.data.n1000.c60[[sim]])
}
sim.results.cox(cf.sim1.c60.n1000.cox.r1, beta = true.beta[[1]],simdata = n.sim)
sim.results.cox(cf.sim1.c60.n1000.cox.r2, beta = true.beta[[2]],simdata = n.sim)

#Binary variable simulations
#sim11, n = 300, cm = 12
true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim11.data.n300.bin <- phcshcf.gen.data.bin(no.of.sims=n.sim, cureMin = 0.5, 
cureMax = 3.765278, pi.r.z.coef = true.gamma, n.obs = 300, prob_event = 0.2,
beta = true.beta, gammaL = 0.5, gammaR = 2.8)
total.noncure <- lapply(cf.sim11.data.n300.bin, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim11.data.n300.bin, function(a) a$event))) /
  length(unlist(lapply(cf.sim11.data.n300.bin, function(a) a$event)))

# allt <- unlist(lapply(cf.sim11.data.n300.bin, 
#   function(a) na.omit(c(a[which(a$event != 0), ]$time,a[which(a$event != 0), 
#                                                         ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores()-1
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim11n300.results.time.bin <- system.time({cf.sim11n300.results.bin <- 
  foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim11.data.n300.bin[[i]]$risk, data = cf.sim11.data.n300.bin[[i]],
  z = as.matrix(cbind(cf.sim11.data.n300.bin[[i]]$z1, 
  cf.sim11.data.n300.bin[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

cf.mpl.results(cf.sim11n300.results.bin, true.gamma, true.beta)

#create variables for Cox Regression application
for(i in 1:n.sim){
cf.sim11.data.n300.bin[[i]]$tmid <- ifelse(
cf.sim11.data.n300.bin[[i]]$event==0, cf.sim11.data.n300.bin[[i]]$time,
ifelse(cf.sim11.data.n300.bin[[i]]$event==1, cf.sim11.data.n300.bin[[i]]$time,
ifelse(cf.sim11.data.n300.bin[[i]]$event==2, cf.sim11.data.n300.bin[[i]]$time/2,
(cf.sim11.data.n300.bin[[i]]$time +cf.sim11.data.n300.bin[[i]]$time2) / 2)))
cf.sim11.data.n300.bin[[i]]$risk_1 <- 
ifelse(is.na(cf.sim11.data.n300.bin[[i]]$risk), 0,
ifelse(cf.sim11.data.n300.bin[[i]]$risk==1, 1, 0))
  cf.sim11.data.n300.bin[[i]]$risk_2 <- 
  ifelse(is.na(cf.sim11.data.n300.bin[[i]]$risk), 0,
 ifelse(cf.sim11.data.n300.bin[[i]]$risk==2, 1, 0))
}

cf.sim1.r1.n300.cox.bin = cf.sim1.r2.n300.cox.bin = list()
for(sim in 1:n.sim){
  cf.sim1.r1.n300.cox.bin[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                         data = cf.sim11.data.n300.bin[[sim]])
  cf.sim1.r2.n300.cox.bin[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                         data = cf.sim11.data.n300.bin[[sim]])
}
sim.results.cox(cf.sim1.r1.n300.cox.bin, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.r2.n300.cox.bin, beta = true.beta[[2]], simdata = n.sim)


#sim11, n = 1000, cm = 12
true.gamma <- c(1.75, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim11.data.n1000.bin <- phcshcf.gen.data.bin(no.of.sims=n.sim, cureMin = 0.5, 
    cureMax = 4.323413, pi.r.z.coef = true.gamma, n.obs = 1000, prob_event = 0.2,
    beta = true.beta, gammaL = 0.5, gammaR = 2.8)
total.noncure <- lapply(cf.sim11.data.n1000.bin, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim11.data.n1000.bin, function(a) a$event))) /
  length(unlist(lapply(cf.sim11.data.n1000.bin, function(a) a$event)))

# allt <- unlist(lapply(cf.sim11.data.n1000.bin, function(a) 
#   na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim11n1000.results.time.bin <- system.time({cf.sim11n1000.results.bin <- 
foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim11.data.n1000.bin[[i]]$risk, data = cf.sim11.data.n1000.bin[[i]],
  z = as.matrix(cbind(cf.sim11.data.n1000.bin[[i]]$z1, 
  cf.sim11.data.n1000.bin[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

cf.mpl.results(cf.sim11n1000.results.bin, true.gamma, true.beta)

#create variables for Cox Regression application
for(i in 1:n.sim){
  cf.sim11.data.n1000.bin[[i]]$tmid <- 
ifelse(cf.sim11.data.n1000.bin[[i]]$event==0, cf.sim11.data.n1000.bin[[i]]$time,
ifelse(cf.sim11.data.n1000.bin[[i]]$event==1, cf.sim11.data.n1000.bin[[i]]$time,
ifelse(cf.sim11.data.n1000.bin[[i]]$event==2, 
cf.sim11.data.n1000.bin[[i]]$time/2, 
(cf.sim11.data.n1000.bin[[i]]$time + cf.sim11.data.n1000.bin[[i]]$time2) / 2)))
  cf.sim11.data.n1000.bin[[i]]$risk_1 <- 
ifelse(is.na(cf.sim11.data.n1000.bin[[i]]$risk), 0,
 ifelse(cf.sim11.data.n1000.bin[[i]]$risk==1, 1, 0))
  cf.sim11.data.n1000.bin[[i]]$risk_2 <- 
  ifelse(is.na(cf.sim11.data.n1000.bin[[i]]$risk), 0,
  ifelse(cf.sim11.data.n1000.bin[[i]]$risk==2, 1, 0))
}

cf.sim1.r1.n1000.cox.bin = cf.sim1.r2.n1000.cox.bin = list()
for(sim in 1:n.sim){
  cf.sim1.r1.n1000.cox.bin[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                          data = cf.sim11.data.n1000.bin[[sim]])
  cf.sim1.r2.n1000.cox.bin[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                          data = cf.sim11.data.n1000.bin[[sim]])
}
sim.results.cox(cf.sim1.r1.n1000.cox.bin, beta = true.beta[[1]], simdata = n.sim)
sim.results.cox(cf.sim1.r2.n1000.cox.bin, beta = true.beta[[2]], simdata = n.sim)

# sim11, n = 300, cm = 12, 40% cured
true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n300.c40.bin <- phcshcf.gen.data.bin(no.of.sims=n.sim, 
  cureMin = 0.5, cureMax = 4.047462, pi.r.z.coef = true.gamma, n.obs = 300, 
  prob_event = 0.2, beta = true.beta, gammaL = 0.5, gammaR = 2.8)
total.noncure <- lapply(cf.sim.data.n300.c40.bin, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n300.c40.bin, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n300.c40.bin, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n300.c40.bin, 
# function(a) na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), 
#                                                        ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores() - 1
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n300.c40.time.bin <- system.time({cf.sim.results.n300.c40.bin <- 
 foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n300.c40.bin[[i]]$risk, 
  data = cf.sim.data.n300.c40.bin[[i]],
  z = as.matrix(cbind(cf.sim.data.n300.c40.bin[[i]]$z1, 
  cf.sim.data.n300.c40.bin[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
 aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
               gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n300.c40.bin, true.gamma, true.beta)

cf.sim1.c40.n300.cox.r1.bin = cf.sim1.c40.n300.cox.r2.bin = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n300.cox.r1.bin[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                        data = cf.sim.data.n300.c40.bin[[sim]])
  cf.sim1.c40.n300.cox.r2.bin[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                        data = cf.sim.data.n300.c40.bin[[sim]])
}
sim.results.cox(cf.sim1.c40.n300.cox.r1.bin, beta = true.beta[[1]], 
                simdata = n.sim)
sim.results.cox(cf.sim1.c40.n300.cox.r2.bin, beta = true.beta[[2]], 
                simdata = n.sim)


# sim11, n = 1000, cm = 12, 40% cured
true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n1000.c40.bin <- phcshcf.gen.data.bin(no.of.sims=n.sim, 
  cureMin = 0.5, cureMax = 3.663038, pi.r.z.coef = true.gamma, n.obs = 1000, 
 prob_event = 0.2, beta = true.beta, gammaL = 0.5, gammaR = 2.8)
total.noncure <- lapply(cf.sim.data.n1000.c40.bin, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n1000.c40.bin, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n1000.c40.bin, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n1000.c40.bin, 
#   function(a) na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), 
#                                                          ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores() - 1
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n1000.c40.time.bin <- system.time({
  cf.sim.results.n1000.c40.bin <- foreach(i = 1:n.sim, .combine = rbind, 
  .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n1000.c40.bin[[i]]$risk, 
  data = cf.sim.data.n1000.c40.bin[[i]],
  z = as.matrix(cbind(cf.sim.data.n1000.c40.bin[[i]]$z1, 
  cf.sim.data.n1000.c40.bin[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n1000.c40.bin, true.gamma, true.beta)

cf.sim1.c40.n1000.cox.r1.bin = cf.sim1.c40.n1000.cox.r2.bin = list()
for(sim in 1:n.sim){
  cf.sim1.c40.n1000.cox.r1.bin[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                        data = cf.sim.data.n1000.c40.bin[[sim]])
  cf.sim1.c40.n1000.cox.r2.bin[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                        data = cf.sim.data.n1000.c40.bin[[sim]])
}
sim.results.cox(cf.sim1.c40.n1000.cox.r1.bin, beta = true.beta[[1]], 
                simdata = n.sim)
sim.results.cox(cf.sim1.c40.n1000.cox.r2.bin, beta = true.beta[[2]], 
                simdata = n.sim)



#sim11, n = 300, cm = 12, 60% cured
true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n300.c60.bin <- phcshcf.gen.data.bin(no.of.sims=n.sim, 
  cureMin = 0.5, cureMax = 3.710337, pi.r.z.coef = true.gamma, n.obs = 300, 
  prob_event = 0.2, beta = true.beta, gammaL = 0.5, gammaR = 2.8)
total.noncure <- lapply(cf.sim.data.n300.c60.bin, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n300.c60.bin, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n300.c60.bin, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n300.c60.bin, function(a) 
#   na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores() - 1
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n300.c60.time.bin <- system.time({cf.sim.results.n300.c60.bin <- 
  foreach(i = 1:n.sim, .combine = rbind, .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n300.c60.bin[[i]]$risk, 
  data = cf.sim.data.n300.c60.bin[[i]],
  z = as.matrix(cbind(cf.sim.data.n300.c60.bin[[i]]$z1, 
  cf.sim.data.n300.c60.bin[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n300.c60.bin, true.gamma, true.beta)

cf.sim1.c60.n300.cox.r1.bin = cf.sim1.c60.n300.cox.r2.bin = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n300.cox.r1.bin[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                      data = cf.sim.data.n300.c60.bin[[sim]])
  cf.sim1.c60.n300.cox.r2.bin[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                  data = cf.sim.data.n300.c60.bin[[sim]])
}
sim.results.cox(cf.sim1.c60.n300.cox.r1.bin, beta = true.beta[[1]], 
                simdata = n.sim)
sim.results.cox(cf.sim1.c60.n300.cox.r2.bin, beta = true.beta[[2]], 
                simdata = n.sim)

#sim11, n = 1000, cm = 12, 60% cured
true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.sim.data.n1000.c60.bin <- phcshcf.gen.data.bin(no.of.sims=n.sim, 
  cureMin = 0.5, cureMax = 3.829909, pi.r.z.coef = true.gamma, n.obs = 1000, 
  prob_event = 0.2, beta = true.beta, gammaL = 0.5, gammaR = 2.8)
total.noncure <- lapply(cf.sim.data.n1000.c60.bin, function(a) a$noncure)
table(unlist(total.noncure)) / length(unlist(total.noncure))

table(unlist(lapply(cf.sim.data.n1000.c60.bin, function(a) a$event))) /
  length(unlist(lapply(cf.sim.data.n1000.c60.bin, function(a) a$event)))

# allt <- unlist(lapply(cf.sim.data.n1000.c60.bin, 
#     function(a) na.omit(c(a[which(a$event != 0), ]$time, a[which(a$event != 0), 
#                                                            ]$time2))))
# max(allt)

#Use parallel processing
cores = detectCores()
cl = makeCluster(cores[1])
registerDoSNOW(cl)
pb <- txtProgressBar(max = n.sim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cf.sim.results.n1000.c60.time.bin <- system.time({
cf.sim.results.n1000.c60.bin <- foreach(i = 1:n.sim, .combine = rbind, 
                                        .options.snow = opts) %dopar% {
  phcshcf_mpl(Surv(time, time2, event=event, "interval") ~ x1 + x2,
  risk = cf.sim.data.n1000.c60.bin[[i]]$risk, 
  data = cf.sim.data.n1000.c60.bin[[i]],
  z = as.matrix(cbind(cf.sim.data.n1000.c60.bin[[i]]$z1, 
  cf.sim.data.n1000.c60.bin[[i]]$z2)), max.outer = 10,
  max.iter = 500, lambda = c(1e-5,1e-5), n.basis = NULL, iknots.pos = NULL,
  aps = TRUE, knots.perc.limit = c(0.05,0.95), tmid = TRUE, iter.disp = FALSE,
  gq.points = 50, z_names=c("intercept", "z1", "z2"))
}
})
close(pb)
stopCluster(cl)

true.gamma <- c(-0.5, 1, -0.5)
true.beta <- list(c(0.5, -2), c(0.5, -0.5))
cf.mpl.results(cf.sim.results.n1000.c60.bin, true.gamma, true.beta)

cf.sim1.c60.n1000.cox.r1.bin = cf.sim1.c60.n1000.cox.r2.bin = list()
for(sim in 1:n.sim){
  cf.sim1.c60.n1000.cox.r1.bin[[sim]] = coxph(Surv(tmid, risk_1==1) ~ x1 + x2, 
                                      data = cf.sim.data.n1000.c60.bin[[sim]])
  cf.sim1.c60.n1000.cox.r2.bin[[sim]] = coxph(Surv(tmid, risk_2==1) ~ x1 + x2, 
                                    data = cf.sim.data.n1000.c60.bin[[sim]])
}
sim.results.cox(cf.sim1.c60.n1000.cox.r1.bin, beta = true.beta[[1]], 
                simdata = n.sim)
sim.results.cox(cf.sim1.c60.n1000.cox.r2.bin, beta = true.beta[[2]], 
                simdata = n.sim)
