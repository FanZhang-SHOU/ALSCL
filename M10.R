#
# operating model M10: Tuna, Flat, Age-dependent
#-----------------------------------------------------
#---------------------------------------------------

# define functions
VB_func<-function(Linf,k,t0,age)
{
  Lt = Linf*(1-exp(-k*(age-t0))) # po is the length of age 0 fish.
  return(Lt)
}

mat_func<-function(L50,L95,length)
{
  b1=log(0.95/0.05)/(L95-L50)
  bo = -L50*b1
  logit_pt =  bo+b1*length
  matp= exp(logit_pt)/(1+exp(logit_pt))
  return(matp)
}

#---------------------------------------------------------------------------------
# define biological and fishing variables

nyear<-100 # number of years
rec.age<-0.25 # the age of recruitment to fishery (one quater)
first.year<-2020 # the first year of projection

# biological variables
max_age<-5 # number of ages

ages<-seq(rec.age,max_age,0.25) # list of ages in the simulation
years<-c(first.year:(first.year+nyear-1)) # list of years in the simulation

nage=length(ages)

len_mid<-seq(12.5,122.5,5)
len_border<-seq(10,125,5)
len_border[1]=-Inf
len_border[length(len_border)]=Inf
nlen<-length(len_mid)
len_lower=len_border[1:nlen]
len_upper=len_border[2:(nlen+1)]

M=0.2 # quaterly natural mortality
init_Z=0.1 # initial total mortality lead to equilibrium age structure

vbk=0.38
Linf=152
t0=1/Linf

LatA=VB_func(Linf,vbk,t0,ages)
plot(ages,LatA,type="l",ylim=c(0,200),xlim=c(0,10))

# weigth at length
a=2e-5; b=3
W_at_len = a*(len_mid^b)
plot(len_mid,W_at_len,type="l")

# maturation at length
mat_L50=100; mat_L95=120
mat=mat_func(mat_L50,mat_L95,len_mid)

plot(len_mid,mat,type="l")

cv_L=0.2 # the coefficient of variance in length at age at the beginning year
cv_inc=0.3 # the coefficient of variance in growth increment

std_logR=0.5 # standard deviation of recruitment
std_logN0=0.3 # standar deviation of initial number at age

# SR relationship (BH model) aS/(b+S)
alpha=15000
beta=150

ssb_test=seq(0,4000,1)
rec_test=alpha*ssb_test/(beta+ssb_test)

plot(ssb_test,rec_test,ylim=c(0,20000),type="l")

# age-length tansition matrix
pla=matrix(NA,nrow=nage,ncol=nlen) # transfer age to length
ml=VB_func(Linf,vbk,t0,ages)
sl=cv_L*ml
for(i in 1:nlen){
  pla[,i]=pnorm(len_border[i+1],ml,sl)-pnorm(len_border[i],ml,sl)
}

# survey variables

std_SN=0.2 # survey measurement error

q_surv_L50=30 # survey catchability at length
q_surv_L95=50
q_surv=mat_func(q_surv_L50,q_surv_L95,len_mid)

# fishing variable

# sel_L50=30 # fishing selectivity at length
# sel_L95=35
# sel=mat_func(sel_L50,sel_L95,len_mid)

sel=rep(1,nage)

#-------------------------------------------------------------------------------------------------
# start the simulation
setwd("E:\\simulation study\\operating models\\M10")

for(iter in 4:100){
  
  # to mimic continuous recruitment, tuna population dynamics is modelled as quater time step
  nstep=nyear*4
  
  # recruitment
  Rec<-rep(NA,nstep)
  R_init=10000 # initial number of recruitment(unit in '000)
  set.seed(iter)
  dev_logR = arima.sim(list(order=c(1,0,0), ar=0.1), n=nstep)*std_logR # std of rec is 0.2
  for(i in 1:nstep){
    Rec[i]<-exp(log(R_init)+dev_logR[i])
  }
  
  # mortality
  Z_at_age<-matrix(NA,nrow=nstep,ncol=nage) # total mortality at age
  F_at_age<-matrix(NA,nrow=nstep,ncol=nage) # fishing mortality at age
  M_at_age<-matrix(M,nrow=nstep,ncol=nage) # natural mortality at age
  set.seed(iter)
  F_yr<-0.2*exp(arima.sim(list(order=c(1,0,0), ar=0.75), n=nstep)*0.2) # annual fishing mortality
  for(i in 1:nstep){
    F_at_age[i,]<-F_yr[i]*sel
    Z_at_age[i,]<-F_at_age[i,]+M_at_age[i,]
  }
  
  # first year age structure
  set.seed(iter)
  N0_at_age = rep(NA,nage)
  dev_logN0 = rnorm((nage-1),0,std_logN0) # std of N0 is 0.3
  N0_at_age[1]<-Rec[1]
  for(i in 2:nage){
    N0_at_age[i]=N0_at_age[i-1]*exp(-init_Z+dev_logN0[i-1])
  }
  
  # age-based cohort dynamics
  
  # define variables
  N_at_age<-matrix(NA,nrow=nstep,ncol=nage) # abundance at age
  N_at_len<-matrix(NA,nrow=nstep,ncol=nlen) # abundance at length
  B_at_len<-matrix(NA,nrow=nstep,ncol=nlen) # biomass at length
  SB_at_len<-matrix(NA,nrow=nstep,ncol=nlen) # spawning biomass at length
  CN_at_age<-matrix(NA,nrow=nstep,ncol=nage) # catch number at age
  CN_at_len<-matrix(NA,nrow=nstep,ncol=nlen) # catch number at length
  CB_at_len<-matrix(NA,nrow=nstep,ncol=nlen) # catch biomass at length
  SSB=rep(NA,nstep)
  TN=rep(NA,nstep)
  TB=rep(NA,nstep)
  CN=rep(NA,nstep)
  CB=rep(NA,nstep)
  
  # Initialization
  N_at_age[1,]=N0_at_age
  N_at_len[1,]=N_at_age[1,]%*%pla
  B_at_len[1,]=N_at_len[1,]*W_at_len
  SB_at_len[1,]=B_at_len[1,]*mat
  SSB[1]=sum(SB_at_len[1,])
  
  # forward dynamics
  for(i in 2:nstep)
  {
    # recruitment process
    logR=log(alpha*SSB[i-1]/(beta + SSB[i-1])) # BH model
    Rec[i]<-exp(logR+dev_logR[i]) # add autocorrelated error
    N_at_age[i,1]<-Rec[i]
    
    for(j in 2:(nage-1))
    {
      # survival process
      N_at_age[i,j]<-N_at_age[i-1,j-1]*exp(-Z_at_age[i-1,j-1])
    }
    # plus group
    N_at_age[i,nage]<-N_at_age[i-1,nage-1]*exp(-Z_at_age[i-1,nage-1]) + N_at_age[i-1,nage]*exp(-Z_at_age[i-1,nage])
    
    N_at_len[i,]=N_at_age[i,]%*%pla
    B_at_len[i,]=N_at_len[i,]*W_at_len
    SB_at_len[i,]=B_at_len[i,]*mat
    SSB[i]=sum(SB_at_len[i,])
    
  }
  
  # derive other quantities
  for(i in 1:nstep)
  {
    for(j in 1:nage)
    {
      CN_at_age[i,j]=N_at_age[i,j]*(1-exp(-Z_at_age[i,j]))*(F_at_age[i,j]/Z_at_age[i,j])
    }
    CN_at_len[i,]=CN_at_age[i,]%*%pla
    CB_at_len[i,]=CN_at_len[i,]*W_at_len
    TN[i]=sum(N_at_len[i,])
    TB[i]=sum(B_at_len[i,])
    CN[i]=sum(CN_at_len[i,])
    CB[i]=sum(CB_at_len[i,])
  }
  
  # calculate research vessel survey index
  
  RVN_at_len=matrix(NA,nrow=nstep,ncol=nlen)
  RVB_at_len=matrix(NA,nrow=nstep,ncol=nlen)
  
  set.seed(iter)
  surv_error<-rnorm(nstep,0,std_SN)
  for(i in 1:nstep)
  {
    RVN_at_len[i,]=N_at_len[i,]*q_surv*exp(surv_error[i])
    RVB_at_len[i,]=RVN_at_len[i,]*W_at_len
  }
  
  for(i in 1:nstep)
  {
    for(j in 1:nlen)
    {
      if(RVN_at_len[i,j]<1e-6){RVN_at_len[i,j]=1e-6}
    }
  }
  
  RVN=rowSums(RVN_at_len)
  RVB=rowSums(RVB_at_len)
  
  sim.data=list(
    SN_at_len = RVN_at_len[381:400,],
    q_surv = q_surv,
    len_mid = len_mid,
    nstep=(nyear-95)*4,
    nage=nage,
    nlen=nlen,
    ages=seq(0.25,5,0.25),
    weight=matrix(rep(W_at_len,(nyear-95)*4),nrow=nlen,ncol=(nyear-95)*4),
    mat=matrix(rep(mat,(nyear-95)*4),nrow=nlen,ncol=(nyear-95)*4),
    
    sel = sel,
    F_at_age = F_at_age[381:400,],
    M_at_age = M_at_age[381:400,],
    N_at_len = N_at_len[381:400,],
    N_at_age = N_at_age[381:400,],
    B_at_len = B_at_len[381:400,],
    SB_at_len = SB_at_len[381:400,],
    CN_at_len = CN_at_len[381:400,],
    CN_at_age = CN_at_age[381:400,],
    CB_at_len = CB_at_len[381:400,],
    RVN_at_len = RVN_at_len[381:400,],
    RVB_at_len = RVB_at_len[381:400,],
    Rec = N_at_age[381:400,1],
    TN = TN[381:400],
    TB = TB[381:400],
    CN = CN[381:400],
    CB = CB[381:400],
    SSB = SSB[381:400]
  )
  
  save(sim.data,file=paste0("sim_rep",iter))
}
