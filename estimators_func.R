#functions for IPW, G-computation, TMLE
#dat must be a data frame created by datagen (or with the same structure)

ipwfunc<-function(dat,k){

A1<-dat$A1
A2<-dat$A2
A3<-dat$A3
A4<-dat$A4
D1<-dat$D1
D2<-dat$D2
D3<-dat$D3

W1<-dat$W1
W2<-dat$W2
W3<-dat$W3
W4<-dat$W4
Y<-dat$Y

Amat<-cbind(A1,A2,A3,A4)
Dmat<-cbind(D1,D2,D3,1)
Wmat<-cbind(W1,W2,W3,W4)

ssize=length(W1)
abar<-rep(0,4); if(k<5){abar[k:4]<-1}


TS<-rep(1,ssize)
TS<-(Amat[,1]==abar[1])
for (t in 2:4){ TS<-(Amat[,t]==abar[t]|Dmat[,t-1]==1)&TS==1 } #checked

Pmat<-matrix(0,nrow=ssize,ncol=4)
Pmat[,1]<-predict( glm(A1~W1,family=binomial) , newdata=as.data.frame(cbind(W1=W1)) , type="response") #probability of exposure first time point
weight<-1/(Pmat[,1]*abar[1]+(1-Pmat[,1])*(1-abar[1]))

for (t in 2:min(k,4)){ #only need to multiply probabilities until a(k)=1 (since after it's systematically a=1)
Wt<-Wmat[,t]
Pmat[Dmat[,t-1]==0,t]<-predict( glm(Amat[,t]~Wt,family=binomial,subset=(Amat[,t-1]==0&Dmat[,t-1]==0)) , newdata=as.data.frame(cbind(Wt=(Wt[Dmat[,t-1]==0]))) , type="response") #probability of exposure subsequent for those eligible
weight<-weight*1/( ( Pmat[,t]*abar[t]+(1-Pmat[,t])*(1-abar[t]) )*(1-Dmat[,t-1]) + Dmat[,t-1] ) #i.e. no multiplier if delivered
}

result<-(lm(Y~1,weights=(TS*weight))$coef)

return(result)

}

gcompfunc<-function(dat,k){


A1<-dat$A1
A2<-dat$A2
A3<-dat$A3
A4<-dat$A4
D1<-dat$D1
D2<-dat$D2
D3<-dat$D3

W1<-dat$W1
W2<-dat$W2
W3<-dat$W3
W4<-dat$W4
Y<-dat$Y

Amat<-cbind(A1,A2,A3,A4)
Dmat<-cbind(D1,D2,D3,1)
Wmat<-cbind(W1,W2,W3,W4)

ssize=length(W1)
abar<-rep(0,4); if(k<5){abar[k:4]<-1}


Qmod<-list()
Qbart<-Y #initialize at Y

Qmod[[4]]<-glm(Qbart~A1+A2+A3+A4+W1+W2+W3+W4,subset=(D3==0),family=binomial())
Q11n<-predict(Qmod[[4]],newdata=as.data.frame(cbind(A1=abar[1],A2=abar[2],A3=(abar[3]),A4=(abar[4]),W1)),type="response")
 Qbart[D3==0]<-Q11n[D3==0]

Qmod[[3]]<-glm(Qbart~A1+A2+A3+W1+W2+W3,family=binomial(),subset=(D2==0))
Q21n<-predict(Qmod[[3]],type="response",newdata=as.data.frame(cbind(A1=abar[1],A2=abar[2],A3=(abar[3]),W1)))
 Qbart[D2==0]<-Q21n[D2==0]

Qmod[[2]]<-glm(Qbart~A1+A2+W1+W2,family=binomial())
Qbart<-predict(Qmod[[2]],type="response",newdata=as.data.frame(cbind(A1=abar[1],A2=abar[2],W1)))
Qmod[[1]]<-glm(Qbart~A1+W1,family=binomial())
Qbart<-predict(Qmod[[1]],type="response",newdata=as.data.frame(cbind(A1=abar[1],W1)))

result<-mean(Qbart)

return(result)

}


tmlefunc<-function(dat,k){


A1<-dat$A1
A2<-dat$A2
A3<-dat$A3
A4<-dat$A4
D1<-dat$D1
D2<-dat$D2
D3<-dat$D3

W1<-dat$W1
W2<-dat$W2
W3<-dat$W3
W4<-dat$W4
Y<-dat$Y

Amat<-cbind(A1,A2,A3,A4)
Dmat<-cbind(D1,D2,D3,1)
Wmat<-cbind(W1,W2,W3,W4)

ssize=length(W1)
abar<-rep(0,4); if(k<5){abar[k:4]<-1}


weightmat<-matrix(nrow=ssize,ncol=4)

TS<-matrix(nrow=ssize,ncol=4);
TS[,1]<-(Amat[,1]==abar[1])
for (t in 2:4){ TS[,t]<-(Amat[,t]==abar[t]|Dmat[,t-1]==1)&TS[,t-1]==1 } #checked

Pmat<-matrix(0,nrow=ssize,ncol=4)
Pmat[,1]<-predict( glm(A1~W1,family=binomial) , newdata=as.data.frame(cbind(W1=W1)) , type="response") #probability of exposure first time point
weightmat[,1]<-weight<-TS[,1]/(Pmat[,1]*abar[1]+(1-Pmat[,1])*(1-abar[1]))

for (t in 2:min(k,4)){ #only need to multiply probabilities until a(k)=1 (since after it's systematically a=1)
Wt<-Wmat[,t]
Pmat[Dmat[,t-1]==0,t]<-predict( glm(Amat[,t]~Wt,family=binomial,subset=(Amat[,t-1]==0&Dmat[,t-1]==0)) , newdata=as.data.frame(cbind(Wt=(Wt[Dmat[,t-1]==0]))) , type="response") #probability of exposure subsequent for those eligible
weightmat[,t]<-weight<-weight*TS[,t]/( ( Pmat[,t]*abar[t]+(1-Pmat[,t])*(1-abar[t]) )*(1-Dmat[,t-1]) + Dmat[,t-1] ) #i.e. no multiplier if delivered
}
if(k<4){for (t in k:3){ #if no additional weight multiplier, fill in last weights
weightmat[,t+1]<-weightmat[,t]
}}

Qmod<-list()
Qbart<-Y #initialize at Y

Qmod[[4]]<-glm(Qbart~A1+A2+A3+A4+W1+W2+W3+W4,subset=(D3==0),family=binomial)
Q41n<-predict(Qmod[[4]],newdata=as.data.frame(cbind(A1=abar[1],A2=abar[2],A3=(abar[3]),A4=(abar[4]),W1)),type="response")
ep4<-coef(glm(Qbart~offset(qlogis(Q41n)),weights=weightmat[,4],subset=(D3==0),family=binomial))
Q41n_star<-plogis( qlogis(Q41n) + ep4)
 Qbart[D3==0]<-Q41n_star[D3==0]
#check sum(weightmat[,4]*(Y-Qbart))

Qmod[[3]]<-glm(Qbart~A1+A2+A3+W1+W2+W3,family=binomial(),subset=(D2==0))
Q31n<-predict(Qmod[[3]],type="response",newdata=as.data.frame(cbind(A1=abar[1],A2=abar[2],A3=(abar[3]),W1)))
ep3<-coef(glm(Qbart~offset(qlogis(Q31n)),weights=weightmat[,3],subset=(D2==0),family=binomial))
Q31n_star<-plogis( qlogis(Q31n) + ep3)
 Qbart[D2==0]<-Q31n_star[D2==0]

Qmod[[2]]<-glm(Qbart~A1+A2+W1+W2,family=binomial())
Q21n<-predict(Qmod[[2]],type="response",newdata=as.data.frame(cbind(A1=abar[1],A2=abar[2],W1)))
ep2<-coef(glm(Qbart~offset(qlogis(Q21n)),weights=weightmat[,2],family=binomial))
Qbart<-plogis( qlogis(Q21n) + ep2)

Qmod[[1]]<-glm(Qbart~A1+W1,family=binomial())
Q11n<-predict(Qmod[[1]],type="response",newdata=as.data.frame(cbind(A1=abar[1],W1)))
ep1<-coef(glm(Qbart~offset(qlogis(Q11n)),weights=weightmat[,1],family=binomial))
Qbart<-plogis( qlogis(Q11n) + ep1)

result<-mean(Qbart)
return(result)

}


#IPW where you ignore delivery time but adjust for other time-dep covariates
ipw_standfunc<-function(dat,k){

A1<-dat$A1
A2<-dat$A2
A3<-dat$A3
A4<-dat$A4
D1<-dat$D1
D2<-dat$D2
D3<-dat$D3

W1<-dat$W1
W2<-dat$W2
W3<-dat$W3
W4<-dat$W4
Y<-dat$Y

Amat<-cbind(A1,A2,A3,A4)
Dmat<-cbind(D1,D2,D3,1)
Wmat<-cbind(W1,W2,W3,W4)

ssize=length(W1)
abar<-rep(0,4); if(k<5){abar[k:4]<-1}


TS<-rep(1,ssize)
#Is the first exposure according to the regime? Ignoring delivery
for (t in 1:min(4,k)){ 
TS<- ( Amat[,t]==abar[t]&TS==1 )
} #


Pmat<-matrix(0,nrow=ssize,ncol=4)
Pmat[,1]<-predict( glm(A1~W1,family=binomial) , newdata=as.data.frame(cbind(W1=W1)) , type="response") #probability of exposure first time point
weight<-1/(Pmat[,1]*abar[1]+(1-Pmat[,1])*(1-abar[1]))
for (t in 2:min(4,k)){
Wt<-Wmat[,t]
Pmat[Dmat[,t-1]==0,t]<-predict( glm(Amat[,t]~Wt,family=binomial,subset=(Amat[,t-1]==0&Dmat[,t-1]==0)) , newdata=as.data.frame(cbind(Wt=(Wt[Dmat[,t-1]==0]))) , type="response") #probability of exposure subsequent for those eligible
weight<-weight*1/( ( Pmat[,t]*abar[t]+(1-Pmat[,t])*(1-abar[t]) ) ) #

}


result<-(lm(Y~1,weights=(TS*weight*(1-D3)))$coef)

return(result)

}


