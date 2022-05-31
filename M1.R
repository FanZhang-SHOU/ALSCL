#
# operating model M1: Flatfish, Flat, Length-dependent
#-----------------------------------------------------
#---------------------------------------------------
# define functions
VB_func<-function(Linf,k,t0,age)
{
  Lt = Linf*(1-exp(-k*(age-t0))) 
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
#plot(ages,LatA,type="l",ylim=c(0,80),xlim=c(0,25))

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

# growth transition matrix
Gij=matrix(NA,nrow=nlen,ncol=nlen) # probability from length bin j to lengh bin i.
for(j in 1:nlen){
  
  # ml_inc = (1-exp(-vbk))*(Linf-len_mid[j]) # expected growth increment based on VB model
  
  # decreasing logistic function
  delta_max=(1-exp(-vbk))*Linf # maximum expected size increment
  l50 = 0.5 * Linf # initial size having 50% maximum size increment
  l95 = 0.05 * Linf # initial size having 95% maximum size increment
  ml_inc = delta_max / (1 + exp(-log(19)*(len_mid[j]-l50)/(l95-l50)))
  
  sl_inc = ml_inc*cv_inc # variability of growth increment
  
  for(i in 1:nlen){
    if(i<j){
      Gij[i,j]=0
    }
    if(i==j){
      Gij[i,j]=pnorm((len_upper[i]-len_mid[j]),ml_inc,sl_inc)
    }
    if(i>j){
      Gij[i,j]= pnorm((len_upper[i]-len_mid[j]),ml_inc,sl_inc)-pnorm((len_lower[i]-len_mid[j]),ml_inc,sl_inc)
    }
  }
}
colSums(Gij)

# survey variables

std_SN=0.1 # survey measurement error

q_surv_L50=15 # survey catchability at length
q_surv_L95=20
q_surv=mat_func(q_surv_L50,q_surv_L95,len_mid)

# fishing variable

# sel_L50=30 # fishing selectivity at length
# sel_L95=35
# sel=mat_func(sel_L50,sel_L95,len_mid)

sel=rep(1,length(len_mid))

#-------------------------------------------------------------------------------------------------
# start the simulation
setwd("E:\\simulation study\\operating models\\M1")

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
  Z_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # total mortality at length
  F_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # fishing mortality at length
  M_at_len<-matrix(M,nrow=nyear,ncol=nlen) # natural mortality at length
  set.seed(iter)
  F_yr<-0.2*exp(arima.sim(list(order=c(1,0,0), ar=0.75), n=nyear)*0.2) # annual fishing mortality
  for(i in 1:nyear){
    F_at_len[i,]<-F_yr[i]*sel
    Z_at_len[i,]<-F_at_len[i,]+M_at_len[i,]
  }
  
  N_at_age<-matrix(NA,nrow=nyear,ncol=nage) # abundance at age
  N_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # abundance at length
  
  # first year age structure
  set.seed(iter)
  N0_at_age = rep(NA,nage)
  dev_logN0 = rnorm((nage-1),0,std_logN0) # std of N0 is 0.3
  N0_at_age[1]<-Rec[1]
  for(i in 2:nage){
    N0_at_age[i]=N0_at_age[i-1]*exp(-init_Z+dev_logN0[i-1])
  }
  
  # three-dimensional cohort dynamics
  nla_list=list()
  bla_list=list()
  sbla_list=list()
  cnla_list=list()
  cbla_list=list()
  
  nla=matrix(NA,nrow=nlen,ncol=nage) # number of length at age
  bla=matrix(NA,nrow=nlen,ncol=nage) # biomass of length at age
  sbla=matrix(NA,nrow=nlen,ncol=nage) # spawning biomass of length at age
  cnla=matrix(NA,nrow=nlen,ncol=nage) # catch number of length at age
  cbla=matrix(NA,nrow=nlen,ncol=nage) # catch biomass of length at age
  
  SSB=rep(NA,nyear)
  TN=rep(NA,nyear)
  TB=rep(NA,nyear)
  CN=rep(NA,nyear)
  CB=rep(NA,nyear)
  
  pla=matrix(NA,nrow=nlen,ncol=nage) # probability of length at age
  ml=VB_func(Linf,vbk,t0,ages)
  sl=cv_L*ml
  for(i in 1:nlen){
    pla[i,]=pnorm(len_border[i+1],ml,sl)-pnorm(len_border[i],ml,sl)
  }
  
  for(i in 1:nage){
    nla[,i]=N0_at_age[i]*pla[,i]
    bla[,i]=nla[,i]*W_at_len
    sbla[,i]=bla[,i]*mat
    cnla[,i]=nla[,i]*(1-exp(-Z_at_len[1,]))*(F_at_len[1,]/Z_at_len[1,])
    cbla[,i]=cnla[,i]*W_at_len
  }
  
  nla_list[[1]]=nla
  bla_list[[1]]=bla
  sbla_list[[1]]=sbla
  cnla_list[[1]]=cnla
  cbla_list[[1]]=cbla
  
  TN[1]=sum(nla)
  TB[1]=sum(bla)
  SSB[1]=sum(sbla)
  CN[1]=sum(cnla)
  CB[1]=sum(cbla)
  
  for(i in 2:nyear)
  {
    # recruitment process
    logR=log(alpha*SSB[i-1]/(beta + SSB[i-1])) # BH model
    Rec[i]<-exp(logR+dev_logR[i]) # add autocorrelated error
    nla[,1]=Rec[i]*pla[,1]
    bla[,1]=nla[,1]*W_at_len
    sbla[,1]=bla[,1]*mat
    cnla[,1]=nla[,1]*(1-exp(-Z_at_len[i,]))*(F_at_len[i,]/Z_at_len[i,])
    cbla[,1]=cnla[,1]*W_at_len
    
    for(j in 2:(nage-1))
    {
      # survival process
      nla_survival=nla_list[[i-1]][,j-1]*exp(-Z_at_len[i-1,])
      
      # growth process
      nla[,j]=Gij %*% nla_survival
      bla[,j]=nla[,j]*W_at_len
      sbla[,j]=bla[,j]*mat
      cnla[,j]=nla[,j]*(1-exp(-Z_at_len[i,]))*(F_at_len[i,]/Z_at_len[i,])
      cbla[,j]=cnla[,j]*W_at_len
    }
    # plus group
    # survival process
    nla_survival=nla_list[[i-1]][,nage-1]*exp(-Z_at_len[i-1,]) + nla_list[[i-1]][,nage]*exp(-Z_at_len[i-1,])
    
    # growth process
    nla[,nage]=Gij %*% nla_survival
    bla[,nage]=nla[,nage]*W_at_len
    sbla[,nage]=bla[,nage]*mat
    cnla[,nage]=nla[,nage]*(1-exp(-Z_at_len[i,]))*(F_at_len[i,]/Z_at_len[i,])
    cbla[,nage]=cnla[,nage]*W_at_len
    
    nla_list[[i]]=nla
    bla_list[[i]]=bla
    sbla_list[[i]]=sbla
    cnla_list[[i]]=cnla
    cbla_list[[i]]=cbla
    
    TN[i]=sum(nla)
    TB[i]=sum(bla)
    SSB[i]=sum(sbla)
    CN[i]=sum(cnla)
    CB[i]=sum(cbla)
  }
  
  # derive two-dimensional number, biomass and catch at length and age
  N_at_age<-matrix(NA,nrow=nyear,ncol=nage) # number at age
  N_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # number at length
  B_at_age<-matrix(NA,nrow=nyear,ncol=nage) # biomass at age
  B_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # biomass at length
  SB_at_age<-matrix(NA,nrow=nyear,ncol=nage) # spawning biomass at age
  SB_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # spawning biomass at length
  CN_at_age<-matrix(NA,nrow=nyear,ncol=nage) # catch number at age
  CN_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # catch number at length
  CB_at_age<-matrix(NA,nrow=nyear,ncol=nage) # catch biomass at age
  CB_at_len<-matrix(NA,nrow=nyear,ncol=nlen) # catch biomass at length
  
  for(i in 1:nyear){
    N_at_age[i,]<-colSums(nla_list[[i]])
    N_at_len[i,]<-rowSums(nla_list[[i]])
    B_at_age[i,]<-colSums(bla_list[[i]])
    B_at_len[i,]<-rowSums(bla_list[[i]])
    SB_at_age[i,]<-colSums(sbla_list[[i]])
    SB_at_len[i,]<-rowSums(sbla_list[[i]])
    CN_at_age[i,]<-colSums(cnla_list[[i]])
    CN_at_len[i,]<-rowSums(cnla_list[[i]])
    CB_at_age[i,]<-colSums(cbla_list[[i]])
    CB_at_len[i,]<-rowSums(cbla_list[[i]])
  }
  
  # derive one dimensional number, biomass and catch
  TN<-rep(NA,nyear)
  TB<-rep(NA,nyear)
  SSB<-rep(NA,nyear)
  CN<-rep(NA,nyear)
  CB<-rep(NA,nyear)
  
  for(i in 1:nyear){
    TN[i]<-sum(N_at_len[i,])
    TB[i]<-sum(B_at_len[i,])
    SSB[i]<-sum(SB_at_len[i,])
    CN[i]<-sum(CN_at_len[i,])
    CB[i]<-sum(CB_at_len[i,])
  }
  
  # derive age-based mortality
  Z_at_age<-matrix(NA,nrow=nyear-1,ncol=nage-1)
  F_at_age<-matrix(NA,nrow=nyear-1,ncol=nage-1)
  M_at_age<-matrix(NA,nrow=nyear-1,ncol=nage-1)
  
  for(i in 1:(nyear-1)){
    for(j in 1:(nage-1)){
      Z_at_age[i,j]=log(N_at_age[i,j])-log(N_at_age[i+1,j+1])
      #F_at_age[i,j]=log(N_at_age[i,j]) - log(N_at_age[i,j]-CN_at_age[i,j])
      F_at_age[i,j]=(CN_at_age[i,j]*Z_at_age[i,j])/(N_at_age[i,j]*(1-exp(-Z_at_age[i,j])))
      M_at_age[i,j]=Z_at_age[i,j]-F_at_age[i,j]
    }
  }
  
  # calculate research vessel survey index
  rvnla_list=list()
  rvbla_list=list()
  rvnla=matrix(NA,nrow=nlen,ncol=nage) # number of length at age
  rvbla=matrix(NA,nrow=nlen,ncol=nage) # biomass of length at age
  RVN_at_age=matrix(NA,nrow=nyear,ncol=nage)
  RVN_at_len=matrix(NA,nrow=nyear,ncol=nlen)
  RVB_at_age=matrix(NA,nrow=nyear,ncol=nage)
  RVB_at_len=matrix(NA,nrow=nyear,ncol=nlen)
  
  set.seed(iter)
  surv_error<-rnorm(nyear,0,std_SN)
  for(i in 1:nyear){
    temp=nla_list[[i]]
    for(j in 1:nage){
      rvnla[,j]=temp[,j]*q_surv*exp(surv_error[i])
      rvbla[,j]=rvnla[,j]*W_at_len
    }
    rvnla_list[[i]]=rvnla
    rvbla_list[[i]]=rvbla
    
    RVN_at_age[i,]<-colSums(rvnla_list[[i]])
    RVN_at_len[i,]<-rowSums(rvnla_list[[i]])
    RVB_at_age[i,]<-colSums(rvbla_list[[i]])
    RVB_at_len[i,]<-rowSums(rvbla_list[[i]])
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
    vbk=vbk,
    Linf=Linf,
    
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
    F_at_len = F_at_len[81:100,],
    F_at_age = F_at_age[81:99,],
    M_at_len = M_at_len[81:100,],
    M_at_age = M_at_age[81:99,],
    N_at_len = N_at_len[81:100,],
    N_at_age = N_at_age[81:100,],
    B_at_len = B_at_len[81:100,],
    B_at_age = B_at_age[81:100,],
    SB_at_len = SB_at_len[81:100,],
    SB_at_age = SB_at_age[81:100,],
    CN_at_len = CN_at_len[81:100,],
    CN_at_age = CN_at_age[81:100,],
    CB_at_len = CB_at_len[81:100,],
    CB_at_age = CB_at_age[81:100,],
    RVN_at_len = RVN_at_len[81:100,],
    RVN_at_age = RVN_at_age[81:100,],
    RVB_at_len = RVB_at_len[81:100,],
    RVB_at_age = RVB_at_age[81:100,],
    Rec = N_at_age[81:100,1],
    TN = TN[81:100],
    TB = TB[81:100],
    CN = CN[81:100],
    CB = CB[81:100],
    SSB = SSB[81:100]
  )
  
  save(sim.data,file=paste0("sim_rep",iter))
}
