#This file fits the same 6 models to Acadian redfish data in Miller and Hyun (in press at CJFAS)

load("red_SSM_0.RData")
library(TMB)
#dyn.unload(dynlib("SSM_0"))
compile("SSM_0.cpp","-O0 -g")
dyn.load(dynlib("SSM_0"))
#gdbsource("SSM_0.r", TRUE)

#statistical catch-at-age model, fixed M constant wrt age
temp = red
temp$map = list(M_pars1 = factor(rep(NA,length(temp$par$M_pars1))), 
  log_NAA = factor(rep(NA,length(temp$par$log_NAA))), log_NAA_sigma = factor(rep(NA, length(temp$par$log_NAA_sigma))), log_b = factor(NA))
red.M1 <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = "log_R", map = temp$map)
red.M1$opt <- nlminb(red.M1$par,red.M1$fn,red.M1$gr, control = list(iter.max = 1000, eval.max = 1000))
red.M1$rep = red.M1$report()
x = red.M1
red.M1$sdrep = sdreport(x)

#statistical catch-at-age model, estimate M constant wrt age
temp = red
temp$map = list(#M_pars1 = factor(rep(NA,length(temp$par$M_pars1))), 
  log_NAA = factor(rep(NA,length(temp$par$log_NAA))), log_NAA_sigma = factor(rep(NA, length(temp$par$log_NAA_sigma))), log_b = factor(NA))
red.M2 <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = "log_R", map = temp$map)
red.M2$opt<-nlminb(red.M2$par,red.M2$fn,red.M2$gr, control = list(iter.max = 1000, eval.max = 1000))
x = red.M2
red.M2$sdrep = sdreport(x)
red.M2$rep = red.M2$report()

#statistical catch-at-age, Lorenzen M with fixed b
temp = red
temp$dat$M_model = 1
temp$par$M_pars1 = log(0.05)
temp$map = list(#M_pars1 = factor(rep(NA,length(temp$par$M_pars1))), 
  log_NAA = factor(rep(NA,length(temp$par$log_NAA))), log_NAA_sigma = factor(rep(NA, length(temp$par$log_NAA_sigma))), log_b = factor(NA))
red.M3 <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = "log_R", map = temp$map)
red.M3$opt <-nlminb(red.M3$par,red.M3$fn,red.M3$gr, control = list(iter.max = 1000, eval.max = 1000))
red.M3$rep = red.M3$report()
x = red.M3
red.M3$sdrep = sdreport(x)

#state-space model, fixed M constant wrt age (abundance is the state vector)
temp = red
temp$dat$use_NAA_re = 1
temp$dat$random_recruitment = 0
temp$map = list(M_pars1 = factor(rep(NA,length(temp$par$M_pars1))), 
  log_R = factor(rep(NA,length(temp$par$log_R))), log_R_sigma = factor(rep(NA, length(temp$par$log_R_sigma))), log_b = factor(NA))
red.M4 <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = "log_NAA", map = temp$map)
red.M4$opt <- nlminb(red.M4$par,red.M4$fn,red.M4$gr, control = list(iter.max = 1000, eval.max = 1000))
red.M4$rep = red.M4$report()
x = red.M4
red.M4$sdrep = sdreport(x)

#state-space model, estimate M constant wrt age (abundance is the state vector)
temp = red
temp$dat$use_NAA_re = 1
temp$dat$random_recruitment = 0
temp$map = list(#M_pars1 = factor(rep(NA,temp$dat$n_M_re)), 
  log_R = factor(rep(NA,length(temp$par$log_R))), log_R_sigma = factor(rep(NA, length(temp$par$log_R_sigma))), log_b = factor(NA))
red.M5 <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = "log_NAA", map = temp$map)
red.M5$opt <- nlminb(red.M5$par,red.M5$fn,red.M5$gr, control = list(iter.max = 1000, eval.max = 1000))
red.M5$rep = red.M5$report()
x = red.M5
red.M5$sdrep = sdreport(x)

#state-space model  (abundance is the state vector), Lorenzen M with fixed b
temp = red
temp$dat$use_NAA_re = 1
temp$dat$random_recruitment = 0
temp$dat$M_model = 1
temp$par$M_pars1 = log(0.05)
temp$map = list(#M_pars1 = factor(rep(NA,temp$dat$n_M_re)), 
  log_R = factor(rep(NA,length(temp$par$log_R))), log_R_sigma = factor(rep(NA, length(temp$par$log_R_sigma))), log_b = factor(NA))
red.M6 <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = "log_NAA", map = temp$map)
red.M6$opt <- nlminb(red.M6$par, red.M6$fn,red.M6$gr, control = list(iter.max = 1000, eval.max = 1000))
red.M6$rep = red.M6$report()
x = red.M6
red.M6$sdrep = sdreport(x)

peel.fit.fn = function(peel, model = red.M1)
{
  print(peel)
  print(unique(names(model$env$par[model$env$random])))
  temp = list(dat = model$env$data, par = model$env$parList(par = model$env$last.par.best), map = model$env$map, 
    random = unique(names(model$env$par[model$env$random])))
  temp$dat$n_years = temp$dat$n_years - peel
  ind = numeric()
  for(i in 1:temp$dat$n_indices) ind = c(ind, 1:temp$dat$n_years + (i-1)*model$env$data$n_years)
  temp$dat$index_paa = temp$dat$index_paa[ind,]
  ind = numeric()
  for(i in 1:temp$dat$n_fleets) ind = c(ind, 1:temp$dat$n_years + (i-1)*model$env$data$n_years)
  temp$dat$catch_paa = temp$dat$catch_paa[ind,]
  log_NAA_na_ind = rbind(matrix(1:(temp$dat$n_ages*(temp$dat$n_years-1)), temp$dat$n_years-1), matrix(rep(NA, peel*temp$dat$n_ages), peel))
  F_devs_na_ind = rbind(matrix(1:(temp$dat$n_fleets * (temp$dat$n_years-1)), temp$dat$n_years-1), matrix(rep(NA, peel * temp$dat$n_fleets)))
  log_R_na_ind = c(1:(temp$dat$n_years-1), rep(NA, peel))
  if("log_R" %in% temp$random) temp$map$log_R = factor(log_R_na_ind)
  if("log_NAA" %in% temp$random) temp$map$log_NAA = factor(log_NAA_na_ind)
  temp$map$F_devs = factor(F_devs_na_ind)
  temp.mod <- MakeADFun(temp$dat,temp$par,DLL="SSM_0", random = temp$random, map = temp$map)
  temp.opt = nlminb(temp.mod$par,temp.mod$fn,temp.mod$gr, control = list(iter.max = 1000, eval.max = 1000))
  return(list(opt = temp.opt, rep = temp.mod$report()))
}
red.M1$peels = list(peel.fit.fn(1, model = red.M1))
red.M1$peels[[2]] = peel.fit.fn(2, model = red.M1)
red.M1$peels[[3]] = peel.fit.fn(3, model = red.M1)
red.M1$peels[[4]] = peel.fit.fn(4, model = red.M1)
red.M1$peels[[5]] = peel.fit.fn(5, model = red.M1)

red.M2$peels = list(peel.fit.fn(1, model = red.M2))
red.M2$peels[[2]] = peel.fit.fn(2, model = red.M2)
red.M2$peels[[3]] = peel.fit.fn(3, model = red.M2)
red.M2$peels[[4]] = peel.fit.fn(4, model = red.M2)
red.M2$peels[[5]] = peel.fit.fn(5, model = red.M2)

red.M3$peels = list(peel.fit.fn(1, model = red.M3))
red.M3$peels[[2]] = peel.fit.fn(2, model = red.M3)
red.M3$peels[[3]] = peel.fit.fn(3, model = red.M3)
red.M3$peels[[4]] = peel.fit.fn(4, model = red.M3)
red.M3$peels[[5]] = peel.fit.fn(5, model = red.M3)

red.M4$peels = list(peel.fit.fn(1, model = red.M4))
red.M4$peels[[2]] = peel.fit.fn(2, model = red.M4)
red.M4$peels[[3]] = peel.fit.fn(3, model = red.M4)
red.M4$peels[[4]] = peel.fit.fn(4, model = red.M4)
red.M4$peels[[5]] = peel.fit.fn(5, model = red.M4)

red.M5$peels = list(peel.fit.fn(1, model = red.M5))
red.M5$peels[[2]] = peel.fit.fn(2, model = red.M5)
red.M5$peels[[3]] = peel.fit.fn(3, model = red.M5)
red.M5$peels[[4]] = peel.fit.fn(4, model = red.M5)
red.M5$peels[[5]] = peel.fit.fn(5, model = red.M5)

red.M6$peels = list(peel.fit.fn(1, model = red.M6))
red.M6$peels[[2]] = peel.fit.fn(2, model = red.M6)
red.M6$peels[[3]] = peel.fit.fn(3, model = red.M6)
red.M6$peels[[4]] = peel.fit.fn(4, model = red.M6)
red.M6$peels[[5]] = peel.fit.fn(5, model = red.M6)

#Mohn's rho
temp.fn = function(mod = red.M1) mean(sapply(1:5, function(x) mod$peels[[x]]$rep$SSB[mod$env$data$n_years-x]/mod$rep$SSB[mod$env$data$n_years-x] - 1))
x = cbind(SSB = sapply(paste0("red.M", 1:6), function(x) temp.fn(eval(parse(text = x))))
temp.fn = function(mod = red.M1) mean(sapply(1:5, function(x) mod$peels[[x]]$rep$F[mod$env$data$n_years-x]/mod$rep$F[mod$env$data$n_years-x] - 1))
x = cbind(x, "F" = sapply(paste0("red.M", 1:6), function(x) temp.fn(eval(parse(text = x))))
x
