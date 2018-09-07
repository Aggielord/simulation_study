#Pan's method
library(interval)
library(FHtest)
library(Icens)
library(survival)
library(KMsurv)

sfunction<- function(left,right1 ){
  r6_left= round(left,3)
  r6_right= round(right1,3)
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

M=30
n=100 
b=0.4
q=0.2
rho=0
p_value=rep(0,10000)
for (p in 1:10000){
#p=1
set.seed(p)
t_i1=rexp(n/2,0.2)
t_i2=rexp(n/2,0.2*exp(b))
#t_i1=rweibull(n/2,2,scale=1)
#t_i2=rweibull(n/2,2,scale=2.25)
ti=c(t_i1,t_i2)
Li=rep(0,n)
Ri=rep(0,n)
for (i in 1:n){
  k=rbinom(1,9,1-q)
  ki=sample(1:9,k,replace =FALSE)
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
li1=Li[1:(n/2)]
ri1=Ri[1:(n/2)]
li2=Li[(n/2+1):n]
ri2=Ri[(n/2+1):n]

G=rep(0,M)
V=rep(0,M)

for (m in 1:M){
  
  df1 = data.frame(li1,ri1)
  df1=df1[sample((1:nrow(df1)),(n/2),replace=TRUE), ]
  df2 = data.frame(li2,ri2)
  df2=df2[sample((1:nrow(df2)),(n/2),replace=TRUE), ]
  newli1=df1[,"li1"]
  newri1=df1[,"ri1"]
  newli2=df2[,"li2"]
  newri2=df2[,"ri2"]
  
 
  
  s1=sfunction(li1,ri1)
  s2=sfunction(li2,ri2)
  L1=sort(unique(c(li1,ri1)))
  L2=sort(unique(c(li2,ri2)))
  
  x1=rep(0,(n/2))
  x2=rep(0,(n/2))
  statu1=rep(1,(n/2))
  statu2=rep(1,(n/2))
  for (i in 1:(n/2)){
    
    point=match(ri1[i],L1)-match(li1[i],L1)
    prob1=rep(0,point)
    if (ri1[i]==10000) {
      x1[i]=li1[i]
      statu1[i]=0
    }
    else if (point==1) x1[i]=ri1[i]
    else{
      fenmu1=s1[match(li1[i],L1)]-s1[match(ri1[i],L1)]
      fenzi1=rep(0,point)
      value1=rep(0,point)
      ma1=match(li1[i],L1)
      for (k in 1:point){
        #ma=match(li1[i],L1)
        #fenzi[k]=s1[match(li1[i]+k-1,L1)]-s1[match(li1[i]+k-1,L1)+1]#change
        fenzi1[k]=s1[ma1+k-1]-s1[ma1+k]
        prob1[k]=fenzi1[k]/fenmu1
        value1[k]=L1[ma1+k]
      }
      x1[i]=sample(value1,1,replace = FALSE,prob = prob1)
    }
  }
  for (i in 1:(n/2)){
    point=match(ri2[i],L2)-match(li2[i],L2)
    prob2=rep(0,point)
    if (ri2[i]==10000) {
      x2[i]=li2[i]
      statu2[i]=0
    }
    else if (point==1) x2[i]=ri2[i]
    else{
      fenmu2=s2[match(li2[i],L2)]-s2[match(ri2[i],L2)]
      fenzi2=rep(0,point)
      value2=rep(0,point)
      ma2=match(li2[i],L2)
      for (k in 1:point){
        fenzi2[k]=s2[ma2+k-1]-s2[ma2+k]
        prob2[k]=fenzi2[k]/fenmu2
        value2[k]=L2[ma2+k]
      }
      x2[i]=sample(value2,1,replace = FALSE,prob = prob2)
    }
  }

#  x1=c(10,12,8,7,6,5)
#  x2=c(13,9,5,2,3,8)
#  statu1=c(0,1,1,1,0,1) 
#  statu2=c(1,1,1,1,0,0)#1 is observed, 0 is right cencored.
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

G_tem=su^rho*(y1*y2/(y1+y2))*(dn1/y1-dn2/y2)
G_tem[is.nan(G_tem)] <- 0
G[m]=sum(G_tem)

V_tem1=su^(2*rho)/y1*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
V_tem2=su^(2*rho)/y2*(y1*y2/(y1+y2))^2*(1-(dn1+dn2-1)/(y1+y2-1))*(dn1+dn2)/(y1+y2)
V_tem1[is.nan(V_tem1)] <- 0
V_tem2[is.nan(V_tem2)] <- 0
V[m]=sum(V_tem1)+sum(V_tem2)
}


finalG=sum(G)/M
finalv=sum(V)/M+(1+1/M)*var(G)

p_value[p]=2*pnorm(abs(finalG),mean=0,sd = sqrt(finalv),lower.tail = FALSE)

}

print(b)
print(q)
print(length(p_value[which(p_value<=0.05)])/length(p_value))
