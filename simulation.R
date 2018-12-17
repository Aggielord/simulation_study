
library(survival)
library(glrt)
library(Icens)
library(interval)
library(FHtest)
library(doParallel)
library(KMsurv)
cl=makeCluster(16)
registerDoParallel(cl)

rep_m=1000
grange=(1:rep_m)#[-c(84,185,238,239,421,422,447,462,501)]

getX<- function(x1,x2,y1,y2,yn){
  slope <- (y2-y1)/(x2-x1)
  intercept <- y2-slope*x2
  (yn-intercept)/slope
}

sfunction<- function(left,right1 ){
  #r6_left= round(left,6)
  #r6_right= round(right1,6)
  r6_left= left
  r6_right= right1
  
  turn_tau= sort(unique(c(r6_left,r6_right)))
  len_tt= length(turn_tau)
  survobj= survfit(Surv(turn_tau[1:(len_tt-1)],rep(1,len_tt-1))~1)
  init_s= c(1, survobj$surv)
  turn_p= -diff(init_s)
  
  turn_tau12= cbind(turn_tau[-len_tt],turn_tau[-1])
  alpha_fn= function(xx,infi,supr) ifelse(xx[1]>=infi & xx[2]<=supr, 1, 0)
  alpha_mat= apply(turn_tau12,1,alpha_fn,infi=left,supr=right1)
  
  
  id_zero= which(apply(alpha_mat==0,1,all))
  if(length(id_zero)>0) alpha_mat= alpha_mat[-id_zero,]
  
  new_s= init_s
  
  turn_err= 1
  while(turn_err>0.01)
  {
    old_s= new_s
    turn_p= -diff(old_s)
    dj_temp= t(t(alpha_mat)*turn_p)/as.vector(alpha_mat%*%turn_p)
    turn_dj= apply(dj_temp,2,sum)
    turn_yj= rev(cumsum(rev(turn_dj)))
    new_s= c(1,cumprod(1-turn_dj/turn_yj))
    new_s[which(is.na(new_s)==TRUE)]= 0
    turn_err= sum(abs((new_s[-len_tt]-old_s[-len_tt])))
  }  
  new_s[len_tt]=0.5*new_s[len_tt-1]
  return(new_s)
}

#for (g in grange){
oper= foreach(g=grange, .combine='rbind', .packages=c("glrt","survival","Icens","interval","FHtest")) %dopar%{
  set.seed(g)
#  n=30
#  b=0.0
#  q=0.7
#  t_i1=rexp(n/2,0.2)
#  t_i2=rexp(n/2,0.2*exp(b))
#  ti=c(t_i1,t_i2)
#  Li=rep(0,n)
#  Ri=rep(0,n)
#  for (i in 1:n){
#    k=rbinom(1,9,1-q)
#    ki=sample(1:9,k,replace =F)
#    v=sort(ki)
#    a=findInterval(ti[i],v)
#    if (a==0){
#      Li[i]=0
#      Ri[i]=v[a+1]
#    }
#    else if(a>0&&a<k){
#      Li[i]=v[a]
#      Ri[i]=v[a+1]
#    }
#    else {
#      Li[i]=v[a]
#      Ri[i]= 10000
#    }
#  }
#  Ri[which(Ri==10000)]=Inf
#  Ri[is.na(Ri)] = 5
  n=100 #100 50
  M=10
  rho0=0
  rho5=0.5
  rho1=1
  #Wa=0.5  #1 0.5
  #Wb=1.73 #2 1.73
  Ca=8
  Wk=5
  len=1#0.25
  
  #t_i1=(-log(runif(n/2)))^(1/Wa)
  #t_i2=(-log(runif(n/2))/Wb)^(1/Wa)
  
  #b=0.8
  #t_i1=rexp(n/2,2/3)
  #t_i2=rweibull(n/2,2.5,0.7)
  
  t_i1=rexp(n/2,0.2)
  t_i2=rexp(n/2,0.2+1/sqrt(n))
  
  t_i1=round(t_i1,5)
  t_i2=round(t_i2,5)
  
  
  ti=c(t_i1,t_i2)
  Li=rep(0,n)
  Ri=rep(0,n)
  Ci=runif(n,0,Ca)
  
  for (i in 1:n){
    Wtao=c(Ci[i]+((1:Wk)-1)*len,10000)
    prob_miss=c(0,rbinom(Wk-1,1,prob = c(0.3,0.3,0.4,0.4)),0)
    avaWtao=Wtao*(1-prob_miss)
    avaWtao1=c(0,avaWtao[avaWtao!=0])
    
    if(avaWtao1[1]>=ti[i]){
      leftpoint=1
    }else{
    leftpoint=max(which(avaWtao1-ti[i]<0))
    }
    Li[i]=avaWtao1[leftpoint]
    Ri[i]=avaWtao1[leftpoint+1]
  }
  

  
  Li=round(Li,5)
  Ri=round(Ri,5)

############################################################### 
  li1=Li[1:(n/2)]
  ri1=Ri[1:(n/2)]
  li2=Li[(n/2+1):n]
  ri2=Ri[(n/2+1):n]
  
  G0=rep(0,M)
  V0=rep(0,M)
  G1=rep(0,M)
  V1=rep(0,M)
  G5=rep(0,M)
  V5=rep(0,M)
  
  for (m in 1:M){
    
    df1 = data.frame(li1,ri1)
    df1=df1[sample((1:nrow(df1)),(n/2),replace=TRUE), ]
    df2 = data.frame(li2,ri2)
    df2=df2[sample((1:nrow(df2)),(n/2),replace=TRUE), ]#bootstrap from each (U,V)
    newli1=df1[,"li1"]
    newri1=df1[,"ri1"]
    newli2=df2[,"li2"]
    newri2=df2[,"ri2"]
    
    s1=sfunction(newli1,newri1)
    s2=sfunction(newli2,newri2)
    
    
    L1=sort(unique(c(newli1,newri1)))
    L2=sort(unique(c(newli2,newri2)))
    
    x1=rep(0,(n/2))
    x2=rep(0,(n/2))
    statu1=rep(1,(n/2))
    statu2=rep(1,(n/2))
      
    for (i in 1:(n/2)) {   
    ax=c(0,L1); ay=c(1,s1)
    lpoint=findInterval(li1[i],ax); rpoint=findInterval(ri1[i],ax)
    if(ri1[i]>=tail(ax, n=1)) rpoint=rpoint-1
    xl=ax[lpoint];xr=ax[rpoint+1]
    ps1=rpoint-lpoint+1
    deno=ay[lpoint]-ay[rpoint+1]
    p=runif(1)
    ry=ay[lpoint]-deno*p
    if(ri1[i]==10000|li1[i]>=tail(ax, n=1)){
    x1[i]=li1[i]
    statu1[i]=0
    }else{
    for (k in 1:ps1){
    x1.tem=ax[lpoint+k-1];x2.tem=ax[lpoint+k]
    y1.tem=ay[lpoint+k-1];y2.tem=ay[lpoint+k]
    if(y1.tem==y2.tem){
    x1[i]=runif(1,min=x1.tem,max=x2.tem)
    }
    else{
    a=getX(x1.tem,x2.tem,y1.tem,y2.tem,ry)
    if(x1.tem<=a&x2.tem>=a) x1[i]=a
    }
  }
  }
  }

for (i in 1:(n/2)) {   
  ax=c(0,L2); ay=c(1,s2)
  lpoint=findInterval(li2[i],ax); rpoint=findInterval(ri2[i],ax)
  if(ri2[i]>=tail(ax, n=1)) rpoint=rpoint-1
  xl=ax[lpoint];xr=ax[rpoint+1]
  ps1=rpoint-lpoint+1
  deno=ay[lpoint]-ay[rpoint+1]
  p=runif(1)
  ry=ay[lpoint]-deno*p
  if(ri2[i]==10000|li2[i]>=tail(ax, n=1)) {
    x2[i]=li2[i]
    statu2[i]=0
    }else{
  for (k in 1:ps1){
    x1.tem=ax[lpoint+k-1];x2.tem=ax[lpoint+k]
    y1.tem=ay[lpoint+k-1];y2.tem=ay[lpoint+k]
    if(y1.tem==y2.tem){
    x2[i]=runif(1,min=x1.tem,max=x2.tem)
    }
    else{
    a=getX(x1.tem,x2.tem,y1.tem,y2.tem,ry)
    if(x1.tem<=a&x2.tem>=a) x2[i]=a
    }
  }
  }
}
   
    x1=round(x1,5)
    x2=round(x2,5)
    #x2[is.na(x2)] <- 5
    u=sort(unique(c(x1[which(statu1==1)] ,x2[which(statu2==1)])))
    
    n1=rep(0,length(u))
    y1=rep(0,length(u))
    n2=rep(0,length(u))
    y2=rep(0,length(u))
    dn1=rep(0,length(u))
    dn2=rep(0,length(u))
    for (i in 1:length(u)){
      for (k in 1:length(x1)){
        if(x1[k]<=u[i]&&statu1[k]==1){
          n1[i]=n1[i]+1
        }
        else if(x1[k]>=u[i]){
          y1[i]=y1[i]+1
        }
      }
      
    }
    
    for (i in 1:length(u)){
      for (k in 1:length(x2)){
        if(x2[k]<=u[i]&&statu2[k]==1){
          n2[i]=n2[i]+1
        }
        else if(x2[k]>=u[i]){
          y2[i]=y2[i]+1
        }
      }
    }
    dn1=c(n1[1],diff(n1))
    dn2=c(n2[1],diff(n2))
    
    my.surv <- Surv(c(x1,x2), c(statu1,statu2),type="right")
    my.fit <- survfit(my.surv~1,type="kaplan-meier")
    su=summary(my.fit)$surv
    
    G_tem0=su^rho0*(y1*y2/(y1+y2))*(dn1/y1-dn2/y2)
    G_tem0[is.nan(G_tem0)] <- 0
    G0[m]=sum(G_tem0)
    
    G_tem1=su^rho1*(y1*y2/(y1+y2))*(dn1/y1-dn2/y2)
    G_tem1[is.nan(G_tem1)] <- 0
    G1[m]=sum(G_tem1)
    
    G_tem5=su^rho5*(y1*y2/(y1+y2))*(dn1/y1-dn2/y2)
    G_tem5[is.nan(G_tem5)] <- 0
    G5[m]=sum(G_tem5)
    
    V_tem01=su^(2*rho0)/y1*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
    V_tem02=su^(2*rho0)/y2*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
    V_tem01[is.nan(V_tem01)] <- 0
    V_tem02[is.nan(V_tem02)] <- 0
    V0[m]=sum(V_tem01)+sum(V_tem02)
    
    V_tem11=su^(2*rho1)/y1*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
    V_tem12=su^(2*rho1)/y2*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
    V_tem11[is.nan(V_tem11)] <- 0
    V_tem12[is.nan(V_tem12)] <- 0
    V1[m]=sum(V_tem11)+sum(V_tem12)
    
    V_tem51=su^(2*rho5)/y1*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
    V_tem52=su^(2*rho5)/y2*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
    V_tem51[is.nan(V_tem51)] <- 0
    V_tem52[is.nan(V_tem52)] <- 0
    V5[m]=sum(V_tem51)+sum(V_tem52)
  }
  
  
  finalG0=sum(G0)/M
  finalv0=sum(V0)/M+(1+1/M)*var(G0)
  finalG1=sum(G1)/M
  finalv1=sum(V1)/M+(1+1/M)*var(G1)
  finalG5=sum(G5)/M
  finalv5=sum(V5)/M+(1+1/M)*var(G5)
  
  p_value13=2*pnorm(abs(finalG0),mean=0,sd = sqrt(finalv0),lower.tail = FALSE)
  p_value14=2*pnorm(abs(finalG5),mean=0,sd = sqrt(finalv5),lower.tail = FALSE)
  p_value15=2*pnorm(abs(finalG1),mean=0,sd = sqrt(finalv1),lower.tail = FALSE)
   
  
#############################################################  
  Ii=c(rep(0,n/2),rep(1,n/2))
  #Li[which(Li==0)]=0.5
  #Ri[which(Ri<Li)]=0.6
  mat <- matrix(c(Li,Ri,Ii),nrow=length(Li))
  
  p_tem1=gLRT1(mat)#Zhao and Sun(2003)
  #  p_value1[g]=as.vector(p_tem1$p)
  p_value1=as.vector(p_tem1$p)
  
  #p_tem2=gLRT4(mat)#Zhao,Duan and Sun
  #  p_value2[g]=as.vector(p_tem2$p)
  #p_value2=as.vector(p_tem2$p) 
  
  df = data.frame("left" = Li, "right" = Ri, "Treat" = Ii)
  icout=ictest(Surv(left,right,type = 'interval2')~Treat,data=df)
  fit3=ictest(Surv(left,right,type = 'interval2')~Treat,data=df,scores='logrank1')#Sun(1996)
  #  p_value3[g]=fit1$p.value
  p_value3=fit3$p.value
  ##############################
  
  fit4 =ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,scores='logrank2')#Finkelstein's
  #  p_value4[g]=fit2$p.value
  p_value4=fit4$p.value
  ############################
  L=with(df,left)
  R=with(df,right)
  trt=with(df,Treat)
  fit5=ictest(L,R,trt,initfit=icout$fit,scores='wmw')#Peto,Peto Wilcoxon-type
  #  p_value5[g]=fit3$p.value
  p_value5=fit5$p.value
  ##########################
  fit6=ictest(Surv(left,right,type = 'interval2')~Treat,data=df,initfit=icout$fit,method='wsr.HLY',mcontrol=mControl(nwsr = 99),scores='logrank1')#Huang's method
  #  p_value6[g]=fit4$p.value
  p_value6=fit6$p.value
  ###########################
  
  fit7=FHtesticp(Surv(left,right,type = 'interval2')~Treat,data=df,rho = 0, lambda = 0)
  p_value7=fit7$pvalue
 ###########################
  fit8=FHtesticp(Surv(left,right,type = 'interval2')~Treat,data=df,rho = 1, lambda = 0)#early diff
  p_value8=fit8$pvalue
  
  fit9=FHtesticp(Surv(left,right,type = 'interval2')~Treat,data=df,rho = 0, lambda = 1)#midd or late diff
  p_value9=fit9$pvalue
  
  fit10=FHtesticp(Surv(left,right,type = 'interval2')~Treat,data=df,rho = 1, lambda = 1)#midd or late diff
  p_value10=fit10$pvalue
  
  #fit11=FHtestics(Surv(left,right,type = 'interval2')~Treat,data=df,rho = 0, lambda = 0)#midd or late diff
  #p_value11=fit11$pvalue
  
  #fit12=FHtestics(Surv(left,right,type = 'interval2')~Treat,data=df,rho = 1, lambda = 0)#early diff
  #p_value12=fit12$pvalue
  #c(p_value1,p_value2,p_value3,p_value4,p_value5,p_value6,p_value7)
  c(p_value1,p_value3,p_value4,p_value5,p_value6,p_value7,p_value8,p_value9,p_value10,p_value13,p_value14,p_value15)
  #c(p_value1,p_value3,p_value4,p_value5,p_value6)
}

write.table(oper,file = "local_100_10.txt")

