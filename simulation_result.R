
library(survival)
library(glrt)
library(Icens)
library(interval)
#n=100 b=0.8 q=0.5 m=1000
m=1000
p_value1=rep(0,m)
p_value2=rep(0,m)
grange=(1:m)
for (g in grange){
  set.seed(g)
  n=100# 50
  b=0.8#,0.4 0 
  q=0.5#,0.7 0.2
  t_i1=rexp(n/2,0.2)
  t_i2=rexp(n/2,0.2*exp(b))
  ti=c(t_i1,t_i2)
  Li=rep(0,n)
  Ri=rep(0,n)
  for (i in 1:n){
    k=rbinom(1,9,1-q)
    ki=sample(1:9,k,replace =F)
    v=sort(ki)
    a=findInterval(ti[i],v)
    if (a==0){
      Li[i]=0
      Ri[i]=v[a+1]
    }
    else if(a>0&&a<k){
      Li[i]=v[a]
      Ri[i]=v[a+1]
    }
    else {
      Li[i]=v[a]
      Ri[i]= 10000
    }
  }

  Ri[is.na(Ri)] = 5
  Ii=c(rep(0,n/2),rep(1,n/2))
  
  mat <- matrix(c(Li,Ri,Ii),nrow=length(Li))
  
  p_tem1=gLRT1(mat)
  p_value1[g]=as.vector(p_tem1$p)
  p_tem2=gLRT4(mat)
  p_value2[g]=as.vector(p_tem2$p)
}

#Fay's method

n=100 
b=0
q=0.5
set.seed(1)
t_i1=rexp(n/2,0.2)
t_i2=rexp(n/2,0.2*exp(b))
ti=c(t_i1,t_i2)
Li=rep(0,n)
Ri=rep(0,n)
for (i in 1:n){
  k=rbinom(1,9,1-q)
  ki=sample(1:9,k,replace =F)
  v=sort(ki)
  a=findInterval(ti[i],v)
  if (a==0){
    Li[i]=0
    Ri[i]=v[a+1]
  }
  else if(a>0&&a<k){
    Li[i]=v[a]
    Ri[i]=v[a+1]
  }
  else {
    Li[i]=v[a]
    Ri[i]= Inf
  }
}
Ri[is.na(Ri)] = 5
Ii=c(rep(0,n/2),rep(1,n/2))

df = data.frame("left" = Li, "right" = Ri, "Treat" = Ii)

icout=ictest(Surv(left,right,type = 'interval2')~Treat,data=df)

fit1=ictest(Surv(left,right,type = 'interval2')~Treat,data=df)#Sun(1996)
fit1
##############################

fit2 =ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,scores='logrank2')#Finkelstein's
fit2
############################
L=with(df,left)
R=with(df,right)
trt=with(df,Treat)
fit3=ictest(L,R,trt,initfit=icout$fit,scores='wmw')#Peto,Peto Wilcoxon-type
fit3
##########################
fit4=ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,method='wsr.HLY',mcontrol=mControl(nwsr = 99),scores='logrank1')#Huang's method
fit4
###########################
fit5=ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,scores='general',dqfunc=function(x){dlogis(qlogis(x))})
fit5

#######################

m=1
p_value1=rep(0,m)
p_value2=rep(0,m)
p_value3=rep(0,m)
p_value4=rep(0,m)
p_value5=rep(0,m)
p_value6=rep(0,m)
#p_value7=rep(0,m)
grange=(1:m)
for (g in 1:m){
 set.seed(1691)
  n=50
  b=0.4
  q=0.2
  t_i1=rexp(n/2,0.2)
  t_i2=rexp(n/2,0.2*exp(b))
  ti=c(t_i1,t_i2)
  Li=rep(0,n)
  Ri=rep(0,n)
  for (i in 1:n){
    k=rbinom(1,9,1-q)
    ki=sample(1:9,k,replace =F)
    v=sort(ki)
    a=findInterval(ti[i],v)
    if (a==0){
      Li[i]=0
      Ri[i]=v[a+1]
    }
    else if(a>0&&a<k){
      Li[i]=v[a]
      Ri[i]=v[a+1]
    }
    else {
      Li[i]=v[a]
      Ri[i]= 10000
    }
  }
  
  Ri[is.na(Ri)] = 5
  Ii=c(rep(0,n/2),rep(1,n/2))
  
  mat <- matrix(c(Li,Ri,Ii),nrow=length(Li))
  
  p_tem1=gLRT1(mat)#Zhao and Sun(2003)
  p_value1[g]=as.vector(p_tem1$p)
  p_tem2=gLRT4(mat)#Zhao,Duan and Sun
  p_value2[g]=as.vector(p_tem2$p)
  
  Ri[which(Ri==10000)]=Inf
  df = data.frame("left" = Li, "right" = Ri, "Treat" = Ii)
  icout=ictest(Surv(left,right,type = 'interval2')~Treat,data=df)
  fit1=ictest(Surv(left,right,type = 'interval2')~Treat,data=df)#Sun(1996)
  p_value3[g]=fit1$p.value
  ##############################
  
  fit2 =ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,scores='logrank2')#Finkelstein's
  p_value4[g]=fit2$p.value
  ############################
  L=with(df,left)
  R=with(df,right)
  trt=with(df,Treat)
  fit3=ictest(L,R,trt,initfit=icout$fit,scores='wmw')#Peto,Peto Wilcoxon-type
  p_value5[g]=fit3$p.value
  ##########################
  fit4=ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,method='wsr.HLY',mcontrol=mControl(nwsr = 99),scores='logrank1')#Huang's method
  p_value6[g]=fit4$p.value
  ###########################
 # fit5=ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,scores='general',dqfunc=function(x){dlogis(qlogis(x))})
 # p_value7[g]=fit5$p.value
}
 
print(length(which(p_value1<=0.05))/length(p_value1))
print(length(which(p_value2<=0.05))/length(p_value2))
print(length(which(p_value3<=0.05))/length(p_value3))
print(length(which(p_value4<=0.05))/length(p_value4))
print(length(which(p_value5<=0.05))/length(p_value5))
print(length(which(p_value6<=0.05))/length(p_value6))
