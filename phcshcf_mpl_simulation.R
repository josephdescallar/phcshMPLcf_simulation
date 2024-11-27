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

# Table 1, n = 300, cure fraction = 20%, int cens = 58%, right cens = 26
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

# Table 1, n = 1000, cure fraction = 20%, int cens = 58%, right cens = 26
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

# Table 1, n = 300, cure fraction = 40%, int cens = 43%, right cens = 45
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

# Table 1, n = 1000, cure fraction = 40%, int cens = 43%, right cens = 45
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

# Table 1, n = 300, cure fraction = 60%, int cens = 29%, right cens = 63
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

# Table 1, n = 1000, cure fraction = 60%, int cens = 29%, right cens = 63
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

#Binary variable simulations - Supplementary Materials
#Supplementary Table 1, n = 300, cure fraction = 20%, int cens = 58%, right cens = 26%
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


#Supplementary Table 1, n = 1000, cure fraction = 20%, int cens = 58%, right cens = 26%
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

#Supplementary Table 1, n = 300, cure fraction = 40%, int cens = 43%, right cens = 45%
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

#Supplementary Table 1, n = 1000, cure fraction = 40%, int cens = 43%, right cens = 45%
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

#Supplementary Table 1, n = 300, cure fraction = 60%, int cens = 29%, right cens = 63%
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

#Supplementary Table 1, n = 1000, cure fraction = 60%, int cens = 29%, right cens = 63%
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
