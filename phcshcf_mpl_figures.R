#Figure 1 - Baseline hazard and their coverage probabilities

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

# Figure 2 - Cumulative Incidence Function plots
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