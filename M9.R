#
# operating model M9: Flatfish, Flat, Dome shape
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
rec.age<-1 # the age of recruitment to fishery
first.year<-2020 # the first year of projection

# biological variables
nage<-15 # number of ages

ages<-c(rec.age:(rec.age+nage-1)) # list of ages in the simulation
years<-c(first.year:(first.year+nyear-1)) # list of years in the simulation

len_mid<-seq(6,50,2)
len_border<-seq(5,51,2)
len_border[1]=-Inf
len_border[length(len_border)]=Inf
nlen<-length(len_mid)
len_lower=len_border[1:nlen]
len_upper=len_border[2:(nlen+1)]

M=0.2 # natural mortality
init_Z=0.5 # initial total mortality lead to equilibrium age structure

vbk=0.2
Linf=60
t0=1/60

LatA=VB_func(Linf,vbk,t0,ages)
plot(ages,LatA,type="l",ylim=c(0,80),xlim=c(0,25))

# weigth at length
a=exp(-12); b=3
W_at_len = a*(len_mid^b)

#plot(len_mid,W_at_len,type="l")

# maturation at length
mat_L50=35; mat_L95=40
mat=mat_func(mat_L50,mat_L95,len_mid)

#plot(len_mid,mat,type="l")

cv_L=0.2 # the coefficient of variance in length at age at the beginning year
cv_inc=0.2 # the coefficient of variance in growth increment

std_logR=0.3 # standard deviation of recruitment
std_logN0=0.2 # standar deviation of initial number at age

# SR relationship (BH model) aS/(b+S)
alpha=400
beta=10

ssb_test=seq(0,300,1)
rec_test=alpha*ssb_test/(beta+ssb_test)

#plot(ssb_test,rec_test,xlim=c(0,300),type="l")

# age-length tansition matrix
pla=matrix(NA,nrow=nage,ncol=nlen) # transfer age to length
ml=VB_func(Linf,vbk,t0,ages)
sl=cv_L*ml
for(i in 1:nlen){
  pla[,i]=pnorm(len_border[i+1],ml,sl)-pnorm(len_border[i],ml,sl)
}

# survey variables

std_SN=0.1 # survey measurement error

q_surv_L50=15 # survey catchability at length
q_surv_L95=20
q_surv=mat_func(q_surv_L50,q_surv_L95,len_mid)

# fishing variable

sel=c(mat_func(3,5,ages[1:10]),mat_func(15,10,ages[11:15])) # dome shape double logistic equation
plot(sel,type="l")

#-------------------------------------------------------------------------------------------------
# start the simulation
setwd("E:\\simulation study\\operating models\\M9")

for(iter in 4:100){
  
  # recruitment
  Rec<-rep(NA,nyear)
  R_init=500 # initial number of recruitment
  set.seed(iter)
  dev_logR = arima.sim(list(order=c(1,0,0), ar=0.1), n=nyear)*std_logR # std of rec is 0.2
  for(i in 1:nyear){
    Rec[i]<-exp(log(R_init)+dev_logR[i])
  }
  
  # mortality
  Z_at_age<-matrix(NA,nrow=nyear,ncol=nage) # total mortality at age
  F_at_age<-matrix(NA,nrow=nyear,ncol=nage) # fishing mortality at age
  M_at_age<-matrix(M,nrow=nyear,ncol=nage) # natural mortality at age
  set.seed(iter)
  F_yr<-0.3*exp(arima.sim(list(order=c(1,0,0), ar=0.75), n=nyear)*0.2) # annual fishing mortality
  for(i in 1:nyear){
    F_at_age[i,]<-F_yr[i]*sel
    Z_at_age[i,]<-F_at_age[i,]+M_at_age[i,]
  }
  
  N_at_age<-matrix(NA,nrow=nyear,ncol=nage) # abundance at age
  N_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # abundance at length
 
  # first year age structure
  set.seed(iter)
  N0_at_age = rep(NA,nage)
  dev_logN0 = rnorm((nage-1),0,std_logN0) 
  N0_at_age[1]<-Rec[1]
  for(i in 2:nage){
    N0_at_age[i]=N0_at_age[i-1]*exp(-init_Z+dev_logN0[i-1])
  }
  
  # age-based cohort dynamics
  
  # define variables
  N_at_age<-matrix(NA,nrow=nyear,ncol=nage) # number at age
  N_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # number at length
  #B_at_age<-matrix(NA,nrow=nyear,ncol=nage) # biomass at age
  B_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # biomass at length
  #SB_at_age<-matrix(NA,nrow=nyear,ncol=nage) # spawning biomass at age
  SB_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # spawning biomass at length
  CN_at_age<-matrix(NA,nrow=nyear,ncol=nage) # catch number at age
  CN_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # catch number at length
  #CB_at_age<-matrix(NA,nrow=nyear,ncol=nage) # catch biomass at age
  CB_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # catch biomass at length
  SSB=rep(NA,nyear)
  TN=rep(NA,nyear)
  TB=rep(NA,nyear)
  CN=rep(NA,nyear)
  CB=rep(NA,nyear)
  
  # Initialization
  N_at_age[1,]=N0_at_age
  N_at_len[1,]=N_at_age[1,]%*%pla
  B_at_len[1,]=N_at_len[1,]*W_at_len
  SB_at_len[1,]=B_at_len[1,]*mat
  SSB[1]=sum(SB_at_len[1,])
  
  # forward dynamics
  for(i in 2:nyear)
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
  for(i in 1:nyear)
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
  
  RVN_at_len=matrix(NA,nrow=nyear,ncol=nlen)
  RVB_at_len=matrix(NA,nrow=nyear,ncol=nlen)
  
  set.seed(iter)
  surv_error<-rnorm(nyear,0,std_SN)
  for(i in 1:nyear)
  {
    RVN_at_len[i,]=N_at_len[i,]*q_surv*exp(surv_error[i])
    RVB_at_len[i,]=RVN_at_len[i,]*W_at_len
  }
  
  for(i in 1:nyear)
  {
    for(j in 1:nlen)
    {
      if(RVN_at_len[i,j]<1e-6){RVN_at_len[i,j]=1e-6}
    }
  }
    
  RVN=rowSums(RVN_at_len)
  RVB=rowSums(RVB_at_len)
  
  sim.data=list(
    SN_at_len = RVN_at_len[81:100,],
    q_surv = q_surv,
    len_mid = len_mid,
    nyear=nyear-80,
    nage=nage,
    nlen=nlen,
    ages=c(1:nage),
    weight=matrix(rep(W_at_len,nyear-80),nrow=nlen,ncol=nyear-80),
    mat=matrix(rep(mat,nyear-80),nrow=nlen,ncol=nyear-80),
    
    sel = sel,
    F_at_age = F_at_age[81:100,],
    M_at_age = M_at_age[81:100,],
    N_at_len = N_at_len[81:100,],
    N_at_age = N_at_age[81:100,],
    B_at_len = B_at_len[81:100,],
    SB_at_len = SB_at_len[81:100,],
    CN_at_len = CN_at_len[81:100,],
    CN_at_age = CN_at_age[81:100,],
    CB_at_len = CB_at_len[81:100,],
    RVN_at_len = RVN_at_len[81:100,],
    RVB_at_len = RVB_at_len[81:100,],
    Rec = N_at_age[81:100,1],
    TN = TN[81:100],
    TB = TB[81:100],
    CN = CN[81:100],
    CB = CB[81:100],
    SSB = SSB[81:100]
  )
  
  save(sim.data,file=paste0("sim_rep",iter))
}
