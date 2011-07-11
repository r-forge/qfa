funcMCurveVar_LG<-function(sim,QFA){
x<-QFA$x
samp<-QFA$samp
PO<-QFA$PO
M<-QFA$M
N<-QFA$N
for (l in 1:nrow(samp)){
PO=samp[l,(M+2*N+3)]
KK=exp(rnorm(sim,samp[l,1],1/samp[l,(N+2)]^0.5))
rr=exp(rnorm(sim,samp[l,(M+2*N+4)],1/samp[l,(M+3*N+5)]^0.5))
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)),xlimmin, xlimmax,add=TRUE,col=5) 
}
}



funcICurveVar_LG<-function(sim,QFA){
x<-QFA$x
samp<-QFA$samp
PO<-QFA$PO
M<-QFA$M
N<-QFA$N
for (i in 1:N){

 plot(-1,-1,main=paste(gene[i],"Curve",i),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])

 for (l in 1:nrow(samp)){
  PO=samp[l,(M+2*N+3)]
  A=1/samp[l,M+N+3+i-1]^0.5
  B=1/samp[l,2*M+3*N+6+i-1]^0.5
  KK=exp(rnorm(sim,samp[l,(2+i-1)],A))
  rr=exp(rnorm(sim,samp[l,(M+2*N+5+i-1)],B))
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)),xlimmin, xlimmax,add=TRUE,col=5) 
 }

KK=exp(K_i[i])
rr=exp(r_i[i])

curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=5)

plot(-1,-1,main=paste(gene[i],"Repeat Curves"),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])
for (j in 1:NoORF[i])
{
KK=K_ij[(j+NoSum[i])]
rr=r_ij[(j+NoSum[i])]
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1+j) 
} 
}

}





funcModelVarPost_LG<-function(QFA,i){
K<-QFA$K
alpha<-QFA$alpha
K_i<-QFA$K_i
k_tau<-QFA$k_tau
K_ij<-QFA$K_ij
r<-QFA$r
gamma<-QFA$gamma
r_i<-QFA$r_i
gamma_i<-QFA$gamma_i
r_ij<-QFA$r_ij
NoSum<-QFA$NoSum

len<-length(K_ij[((1+NoSum[i]):NoSum[i+1])])
list(
A=rnorm(3*len,K,alpha),
B=rnorm(2*len,K_i[i],1/k_tau[i]^0.5),
C=log(K_ij[((1+NoSum[i]):NoSum[i+1])]),
D=rnorm(3*len,r,gamma),
E=rnorm(2*len,r_i[i],gamma_i),
F=log(r_ij[((1+NoSum[i]):NoSum[i+1])])
)
}

funcDen_LG<-function(sampsize,QFA){
K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta
alpha_i_tau<-QFA$alpha_i_tau
gamma_i_tau<-QFA$gamma_i_tau
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij
alpha_ij_tau<-QFA$alpha_ij_tau
gamma_ij_tau<-QFA$gamma_ij_tau

den<-matrix(0,sampsize,13)
den[,1]<-(rnorm(sampsize,K_s,alpha))
den[,2]<-(rnorm(sampsize,K_s,alpha))
den[,3]<-(rnorm(sampsize,alpha_i,1/alpha_i_tau^0.5))
den[,4]<-exp(rnorm(sampsize,K_s,alpha))
den[,5]<-exp(rnorm(sampsize,alpha_ij,1/alpha_ij_tau^0.5))
den[,6]<-exp(rnorm(sampsize,PO_s,beta))
den[,7]<-(rnorm(sampsize,r_s,gamma))
den[,8]<-(rnorm(sampsize,r_s,gamma))
den[,9]<-exp(rnorm(sampsize,gamma_i,1/gamma_i_tau^0.5))
den[,10]<-exp(rnorm(sampsize,r_s,gamma))
den[,11]<-exp(rnorm(sampsize,gamma_ij,1/gamma_ij_tau^0.5))
den[,12]<-rnorm(sampsize,tau_s,1/delta^0.5)
den[,13]<-exp(rnorm(sampsize,tau_s,1/delta^0.5))
den
}

funcPostPred_LG<-function(iter,QFA){

N<-QFA$N
M<-QFA$M
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime

K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta
alpha_i_tau<-QFA$alpha_i_tau
gamma_i_tau<-QFA$gamma_i_tau
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij
alpha_ij_tau<-QFA$alpha_ij_tau
gamma_ij_tau<-QFA$gamma_ij_tau

namesamp<-QFA$namesamp
K<-QFA$K
K_i<-QFA$K_i
K_ij<-QFA$K_ij
PO<-QFA$PO
k_tau<-QFA$k_tau
r<-QFA$r
r_i<-QFA$r_i
r_ij<-QFA$r_ij
r_tau<-QFA$r_tau
taui<-QFA$taui
tau<-QFA$tau





postpred<-matrix(0,iter,length(namesamp))

postpred[,1]<-rnorm(iter,K_s,alpha)

for (i in 1:length(K_i)){
j=(2)+i-1
postpred[,j]<-rnorm(iter,K,alpha_i)
}

t=1
for (i in (N+2):(M+N+1) ){
if(((i-(N+2))==c(NoSum[-1])[t])) t=t+1
postpred[,i]<-exp(rnorm(iter,K_i[t],1/k_tau[t]^0.5))
}
postpred[,(M+N+2)]<-runif(iter,PO_s,beta)

for (i in 1:length(k_tau)){
j=(M+N+3)+i-1
postpred[,j]<-rgamma(iter,(alpha_ij^2)/alpha_ij_sd^2,(alpha_ij^2)/alpha_ij_sd^2)
}

postpred[,(M+2*N+3)]<-rnorm(iter,r_s,gamma)

for (i in 1:length(r_i)){
j=(M+2*N+4)+i-1
postpred[,j]<-rnorm(iter,r,gamma_i)
}

t=1
for (i in (M+2*N+4):(M+3*N+3) ){
if(((i-(M+2*N+4))==c(NoSum[-1])[t])) t=t+1	
postpred[,i]<-exp(rnorm(iter,r_i[t],gamma_ij))
}

for (i in 1:length(r_tau)){
j=(2*M+3*N+4)+i-1
postpred[,j]<-rgamma(iter,(gamma_ij^2)/gamma_ij_sd^2,(gamma_ij_sd^2)/alpha_ij_sd^2)
}

for (i in 1:length(taui)){
j=(2*M+4*N+4)+i-1
postpred[,j]<-rgamma(iter,(tau^2)/delta^2,tau/delta^2)
}

postpred[,(2*M+5*N+4)]<-rgamma(iter,(tau_s^2)/delta^2,tau_s/delta^2)
postpred
}
