#data generation function for the simulation study

#parameters are 
#beta, the strength of the effect of D(2) and D(3) on Y (note that beta is on the odds ratio scale in the paper but here it is on the log-odds scale)
#theta1, the strength of the effect of A(2) on Y, log-odds scale (dashed arrow in Figure 3)
#theta2, the strength of the effect of A(4) on Y, log-odds scale (dotted arrow in Figure 3)

#scenario 1: beta=theta1=theta2=0, 
#scenario 2:  beta=log(2) theta1=theta2=0
#scenario 3:  beta=log(2); theta1=0.4; theta2=0 #only scenario where there is an effect of A(2) on Y
#scenario 4:  beta=log(2); theta1=0; theta2=0.4 #only scenario where there is an effect of A(4) on Y

datagen<-function(ssize=20000,seed=sample(size=1,1:1000000),beta,theta1,theta2,counterfactual=F,k=NA){

set.seed(seed)

D1=D2=rep(0,ssize)
U<-rbinom(p=0.5,n=ssize,size=1)

alpha1=alpha2=2

if(!counterfactual){

#pregnancy onset

W1<-rbinom(ssize, size=1, p=0.3)
A1<-rbinom(size=1,p=(-0.1*W1+0.15),n=ssize)

#during pregnancy, can't deliver

W2<-rbinom(ssize, size=1, p=(0.2+0.1*A1+0.2*W1))
A2<-rbinom(size=1,p=(0.2-0.1*W2),n=ssize)

A2[A1==1]<-1

#during pregnancy, can deliver #but low risk

D2<-rbinom(size=1,p=(0.1)*(alpha1^(A2+U)),n=ssize) #Delivery risk mult for U=1 and A2=1
W3<-rbinom(ssize, size=1, p=(0.2+0.1*A2+0.2*W2))
A3<-rbinom(size=1,p=(0.2-0.1*W3),n=ssize)

A3[A2==1]<-1
A3[D2==1]<-0
W3[D2==1]<-0

D3<-rbinom(size=1,p=0.20*(alpha2^(U+A3)),n=ssize) #delivery risk mult by alpha2 for U=1 and also for A3=1
W4<-rbinom(ssize, size=1, p=(0.2+0.1*A3+0.2*W3))
A4<-rbinom(size=1,p=(0.2-0.1*W4),n=ssize)

D3[D2==1]<-1
A4[A3==1]<-1
A4[D3==1]<-0
W4[D3==1]<-0

Y<-rbinom(size=1,p=plogis(-2+0.1*(W1+W2)+beta*D2+beta*D3+theta1*A2+theta2*A4+ 0.2*U),n=ssize) 


return(as.data.frame(cbind(A1,A2,A3,A4,D1,D2,D3,W1,W2,W3,W4,Y)))

}else{

#what is the counterfactual treatment pattern?
abar<-rep(0,4); if(k<5){abar[k:4]<-1}
cat("Generating a counterfactual data set with treatment pattern", abar, "\n")

#pregnancy onset

W1<-rbinom(ssize, size=1, p=0.3)
A1<-rep(abar[1],ssize)

#during pregnancy, can't deliver

W2<-rbinom(ssize, size=1, p=(0.2+0.1*A1+0.2*W1))
A2<-rep(abar[2],ssize)

A2[A1==1]<-1

#during pregnancy, can deliver #but low risk

D2<-rbinom(size=1,p=(0.1)*(alpha1^(A2+U)),n=ssize) #Delivery risk mult for U=1 and A2=1
W3<-rbinom(ssize, size=1, p=(0.2+0.1*A2+0.2*W2))
A3<-rep(abar[3],ssize)

A3[A2==1]<-1
A3[D2==1]<-0
W3[D2==1]<-0

D3<-rbinom(size=1,p=0.20*(alpha2^(U+A3)),n=ssize) #delivery risk mult by alpha2 for U=1 and also for A3=1
W4<-rbinom(ssize, size=1, p=(0.2+0.1*A3+0.2*W3))
A4<-rep(abar[4],ssize)

D3[D2==1]<-1
A4[A3==1]<-1
A4[D3==1]<-0
W4[D3==1]<-0

Y<-rbinom(size=1,p=plogis(-2+0.1*(W1+W2)+beta*D2+beta*D3+theta1*A2+theta2*A4+ 0.2*U),n=ssize) 


return(as.data.frame(cbind(A1,A2,A3,A4,D1,D2,D3,W1,W2,W3,W4,Y)))
}
}


