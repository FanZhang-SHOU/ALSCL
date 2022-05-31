# stock assessment using ACL for M10
library(TMB)
setwd("E:\\simulation study\\assessment models\\ACL")
# compile("ACL.cpp")
compile(file = "ACL.cpp", "&> /tmp/logfile.log")

for(iter in 4:100){
  setwd("E:\\simulation study\\operating models\\M10")
  load(paste0("sim_rep",iter))
  
  tmb.data=list(
    logN_at_len = t(log(sim.data$SN_at_len)),
    log_q = log(sim.data$q_surv),
    len_border = (sim.data$len_mid + 2.5)[1:(sim.data$nlen-1)],
    age = sim.data$ages,
    Y = sim.data$nstep,
    A = sim.data$nage,
    L = sim.data$nlen,
    weight = sim.data$weight,
    mat=sim.data$mat,
    M=0.2
  )
  
  parameters = list(
    log_init_Z = 0.5,
    log_std_log_N0 = log(0.5),
    
    mean_log_R = 5,
    log_std_log_R = log(0.2),
    logit_log_R = log(0.75/0.25),
    
    mean_log_F = log(0.3),
    log_std_log_F = log(0.3),
    logit_log_F_y = log(0.75/0.25),
    logit_log_F_a = log(0.75/0.25),
    
    log_vbk = log(0.38),
    log_Linf = log(152),
    log_t0 = log(1/152),
    log_cv_len = log(0.3),
    
    log_std_index = log(0.1),
    
    # random effects
    dev_log_R = rep(0,sim.data$nstep),
    dev_log_F = array(0,c(sim.data$nage,sim.data$nstep)),
    dev_log_N0 = rep(0,(sim.data$nage-1))
  )
  
  parameters.L = list(
    # fixed effects
    log_init_Z = log(0.01),
    log_std_log_N0 = -Inf,
    
    mean_log_R = log(10),
    log_std_log_R = log(0.01),
    logit_log_R = -30,
    
    mean_log_F = log(0.01),
    #log_std_log_F = log(0.01),
    #logit_log_F_y = -10,
    #logit_log_F_a = -10,
    
    log_vbk = log(0.1),
    log_Linf = log(100),
    #log_t0 = -20,
    log_cv_len = log(0.01),
    
    log_std_index = log(0.001)
  )
  
  parameters.U = list(
    # fixed effects
    log_init_Z = log(1),
    log_std_log_N0 = log(5),
    
    mean_log_R = 10,
    log_std_log_R = log(1),
    logit_log_R = 20,
    
    mean_log_F = log(1),
    #log_std_log_F = log(2),
    #logit_log_F_y = 20,
    #logit_log_F_a = 10,
    
    log_vbk = log(1),
    log_Linf = log(200),
    #log_t0 = 0,
    log_cv_len = log(1),
    
    log_std_index = log(1)
  )
  
  lower=unlist(parameters.L)
  upper=unlist(parameters.U)
  
  map = list(
    log_std_log_F=factor(NA),
    logit_log_F_y = factor(NA),
    logit_log_F_a = factor(NA),
    #log_vbk=factor(NA),
    #log_Linf=factor(NA),
    log_t0=factor(NA)
  )
  
  rnames=c("dev_log_R","dev_log_F","dev_log_N0")
  
  setwd("E:\\simulation study\\assessment models\\ACL")
  dyn.load("ACL")
  obj<-MakeADFun(tmb.data,parameters,random=rnames,map=map,DLL="ACL",inner.control=list(trace=F, maxit=500))
  opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))
  # opt1<-nlminb(opt$par,obj$fn,obj$gr,lower=lower,upper=upper,control=list(trace=0,iter.max=2000,eval.max=10000))
  # obj$gr(opt$par)
  # cbind(opt$par,lower,upper)
  report<-obj$report()
  bound_check<-c((as.vector(opt$par)-as.vector(lower)),(as.vector(upper)-as.vector(opt$par)))
  bound_hit<-min(bound_check)==0
  result<-list(obj=obj,opt=opt,report=report,bound_hit=bound_hit,bound_check=bound_check,converge=opt$message)
  dyn.unload("ACL")
  
  setwd("E:\\simulation study\\assessment models\\ACL\\A10")
  save(result,file=paste0("result_rep_",iter))
  rm(obj,opt,report,result)
}
