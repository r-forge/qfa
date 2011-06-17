#### Interaction Model ####

fun5<-function(x){
c(defa[(aNoSum[x]+1):(aNoSum[x+1])],rep(NA,dimr-length(defa[(aNoSum[x]+1):(aNoSum[x+1])])),defb[(bNoSum[x]+1):(bNoSum[x+1])],rep(NA,dimr-length(defb[(bNoSum[x]+1):(bNoSum[x+1])])))
}

fun6<-function(x){
x$r/log(2*(x$K-x$g)/(x$K-2*x$g)) #MDR
}

fun7<-function(x){
log(x$K/x$g)/log(2) #MDP
}

fun2<-function(x){
x[is.na(x)]=-Inf
x[x<0]=NA
x
}

fun8<-function(x){
vecMDRa<-x$r_ij/log(2*(x$K_ij-x$PO)/(x$K_ij-2*x$PO)) #MDR
vecMDPa<-log(x$K_ij/x$PO)/log(2) #MDP
as.numeric(lapply(vecMDPa,fun2))*as.numeric(lapply(vecMDRa,fun2))
}

fun9<-function(x){CUSTOMDEF(x)}

funcPosterior_I<-function(samp,N,M,iter,thin,upd){
if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-samp}
if(nrow(samp)>1) {deltagamma<-colMeans(samp[,(N+3):(2*N+2)]*samp[,(3*N+3):(4*N+2)])} else {deltagamma<-(samp[,(N+3):(2*N+2)]*samp[,(3*N+3):(4*N+2)])}
delta=vecsamp[(N+3):(2*N+2)]
sig=sum(rep(1,N)[delta>0.5])
list(
vecsamp=vecsamp,
namesamp=names(vecsamp),
A1=vecsamp[1],
A2=vecsamp[2],
delta=delta,
gamma=vecsamp[(3*N+3):(4*N+2)],
mu_i=vecsamp[(4*N+4):(5*N+3)],
mu=vecsamp[4*N+3],
nu=vecsamp[(5*N+4)],
nuj=vecsamp[(5*N+5):(5*N+6)],
tau=vecsamp[(5*N+7)],
taui=vecsamp[(5*N+8):(6*N+7)],
deltagamma=deltagamma,
samp=samp,
iter=iter,
thin=thin,
burnandupd=(1000+upd),
sig=sig,
order=order(1-delta),
vecorder=order(1-delta)[1:sig]
)
}

qfa.Interaction<-function
(a,b,iter,upd,thin,PlotOutput=TRUE,work,CustomModel=FALSE,Priors=FALSE){
QFAuni<-a$QFAuni
aNoSum<-a$NoSum
bNoSum<-b$NoSum
N<-a$N
NoORF<-cbind(a$NoORF,b$NoORF)
dimr<-max(NoORF)
defa<-fun8(a)
defb<-fun8(b)
vec<-unlist(lapply(1:N,fun5))
y=array(vec,dim=c(dimr,2,N))

funcMODELInteraction()

if (Priors==FALSE){
p=0.10
mu_a=mean(na.omit(c(y)))
mu_b=max(na.omit(c(y)))/2
mu_b=1/(mu_b)^2
alpha_a=1
alpha_b=max(na.omit(c(y)))
gam_b=max(na.omit(c(y)))/2
gam_b=1/(gam_b)^2
tau_a=1
tau_b=max(na.omit(c(y)))
}
QFA.P<-list(p=p,mu_a=mu_a,mu_b=mu_b,alpha_a=alpha_a,alpha_b=alpha_b,gam_b=gam_b,tau_a=tau_a,tau_b=tau_b)

library("rjags")
jags <- jags.model('model1.bug',data = list('y'=y,'NoORF'=NoORF,'p'=p,'N' = N,'mu_a'=mu_a,'mu_b'=mu_b,'alpha_a'=alpha_a,'alpha_b'=alpha_b,'gam_b'=gam_b,'tau_a' = tau_a,'tau_b' = tau_b),n.chains = 1,n.adapt = 100)

TimeC<-(iter+upd)*system.time(update(jags,990))[1]

print(paste("Time till completion",TimeC/(60*60*990),"(hours)",TimeC/(60*990),"(minutes)"))
update(jags, upd)
samp<-coda.samples(jags,
          c('mui','gam','delt','tau','nu','alpha','mu','taui','nuj'),
            iter,thin=thin)
samp<-samp[[1]]

QFA.S<-funcPosterior_I(samp,N,M,iter,thin,upd)
QFA.l<-list(N=N,gene=a$gene,treat="27",y=y,NoORF=NoORF)
QFA<-c(QFA.S,QFA.P,QFA.l)
return(QFA)
}


### Interaction Model Plots to Pdf###
###########################################
print("plot fitted with Conditioning on delta=1")
###########################################
QFA.I.Plots<-function(work,QFA){

vecsamp=QFA$vecsamp
namesamp=QFA$namesamp
A1=QFA$A1
A2=QFA$A2
delta=QFA$delta
gamma=QFA$gamma
mu_i=QFA$mu_i
mu=QFA$mu
nu=QFA$nu
nuj=QFA$nuj
tau=QFA$tau
taui=QFA$taui
deltagamma=QFA$deltagamma
samp=QFA$samp
iter=QFA$iter
thin=QFA$thin
burnandupd=QFA$burnandupd
sig=QFA$sig
order=QFA$order
vecorder=QFA$vecorder
N=QFA$N
gene=QFA$gene
treat=QFA$treat
y=QFA$y
NoORF=QFA$NoORF

funplot1<-function(){
limmin<-min(A1*mu_i, A2*(mu_i+gamma))
limmax<-max(A1*mu_i, A2*(mu_i+gamma))
i=1:N
plot(1,type="n",main=paste("Treatment",treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control (=Alpha1*mu_i)",ylab="Query (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
lines(A1*c(0,1000),A2*c(0,1000),lwd=2)
points(A1*mu_i[i], A2*(mu_i[i]+gamma[i]),main=paste("Treatment",treat,"Degrees","(Conditioning on deltas=1)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[gamma[vecorder]>0]
points(A1*mu_i[i],A2*(mu_i[i]+gamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[gamma[vecorder]<=0]  
points(A1*mu_i[i],A2*(mu_i[i]+gamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(A1*mu_i[i],A2*(mu_i[i]+gamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}
###########################################
print("plot fitted NO conditioning")
###########################################
funplot2<-function(){
limmin<-min(A1*mu_i, A2*(mu_i+deltagamma))
limmax<-max(A1*mu_i, A2*(mu_i+deltagamma))
i=1:N
plot(1,type="n",main=paste("Treatment",treat,"Degrees","(delta=Posterior Expectations)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control (=Alpha1*mu_i)",ylab="Query (=Alpha2*(mu_i+delta_i*gamma_i))",col=8,pch=19,cex=0.5)
lines(A1*c(0,1000),A2*c(0,1000),lwd=2)
points(A1*mu_i[i], A2*(mu_i[i]+deltagamma[i]),main=paste("Treatment",treat,"Degrees","(delta=Posterior Expectations)"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single (=Alpha1*mu_i)",ylab="Double (=Alpha2*(mu_i+delta_i*gamma_i))",col=8,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]>0]
points(A1*mu_i[i],A2*(mu_i[i]+deltagamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]<=0]  
points(A1*mu_i[i],A2*(mu_i[i]+deltagamma[i]),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(A1*mu_i[i],A2*(mu_i[i]+deltagamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}
##########################################
print("Data plot with highlighted Interactions")
##########################################
funplot3<-function(){
#################################### VECORDER LENGTH>2
if(length(vecorder)>1){
limmin<-min(y,na.rm=TRUE)
limmax<-max(y,na.rm=TRUE)*1.1
plot(1,type="n",main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control",ylab="Query",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),A2*c(-1000,10000),col="cadetblue",lwd=2,)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=2)
abline(lm(colMeans(y[,2,],na.rm=TRUE)~0+colMeans(y[,1,],na.rm=TRUE)),col="grey",lty=3)
lines(c(mean(defa[gene=="HIS3"]),mean(defa[gene=="HIS3"])),c(-1000,1000),lwd=2)
lines(c(-1000,1000),c(mean(defb[gene=="HIS3"]),mean(defb[gene=="HIS3"])),lwd=2)
points(colMeans(y[,1,],na.rm=TRUE),colMeans(y[,2,],na.rm=TRUE),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,col=8,cex=0.5)
i=vecorder[deltagamma[vecorder]>0]
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]<=0]  
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder
text(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),gene[i],pos=4,offset=0.1,cex=0.4)
legend(1,limmax, c("Model Fit (y=Alpha2*x)","1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("cadetblue","grey","grey","black"), lty=c(1,2,3,1))
plot(1,type="n",main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Control",ylab="Query",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),A2*c(-1000,10000),col="cadetblue",lwd=2,)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=2)
abline(lm(colMeans(y[,2,],na.rm=TRUE)~0+colMeans(y[,1,],na.rm=TRUE)),col="grey",lty=3)
lines(c(mean(defa[gene=="HIS3"]),mean(defa[gene=="HIS3"])),c(-1000,1000),lwd=2)
lines(c(-1000,1000),c(mean(defb[gene=="HIS3"]),mean(defb[gene=="HIS3"])),lwd=2)
points(colMeans(y[,1,],na.rm=TRUE),colMeans(y[,2,],na.rm=TRUE),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,col=8,cex=0.5)
i=vecorder[deltagamma[vecorder]>0]
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),ylim=c(0,5),xlim=c(0,5),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder[deltagamma[vecorder]<=0]  
points(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),xlim=c(0,5),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=1:N
text(colMeans(y[,1,i],na.rm=TRUE),colMeans(y[,2,i],na.rm=TRUE),gene[i],pos=4,offset=0.1,cex=0.4)
legend(1,limmax, c("Model Fit (y=Alpha2*x)","1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("cadetblue","grey","grey","black"), lty=c(1,2,3,1))
}
}
###########################################
print("Diagnostics trace acf density")
###########################################
funplot4<-function(){
par(mfrow=c(4,1))
den<-matrix(N,2000,8)
den[,1]<-rgamma(2000,(alpha_a^2)/(alpha_b^2),alpha_a/(alpha_b^2))
den[,2]<-rnorm(2000,0,gam_b^(-0.5))
den[,3]<-rnorm(2000,mu_a,mu_b^(-0.5))
den[,4]<-rnorm(2000,mu_a,mu_b^(-0.5))
den[,5]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
den[,6]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
den[,7]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
den[,8]<-rgamma(2000,(tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
namesampden<-unique(substring(namesamp,1,4))
names<-namesampden[-2]
for (i in 1:ncol(den))
{
plot(density(den[,i]),paste(names[i],"Prior Density"))
}
par(mfrow=c(4,3))
for (i in c(2,(3*N+3):c(6*N+7)))
{
plot(as.numeric(samp[,i]),main=paste("Trace Top",namesamp[i],"(Treatment",treat,"Degrees)"),type="l")
plot(density(as.numeric(samp[,i])),main=paste("Density",namesamp[i]))
acf(as.numeric(samp[,i]),main=paste("ACF",namesamp[i]))
}
par(mfrow=c(4,2))
for (i in (N+3):(2*N+2)){
plot(as.numeric(samp[,i]),main=paste("Trace Top",namesamp[i],"(Treatment",treat,"Degrees)"),type="l")
plot(density(as.numeric(samp[,i])),main=paste("Density",namesamp[i]))
}
par(mfrow=c(1,1))
}
#####################################################################
print("Number of Repeats for (Single+Double)Versus Order of Sig Inter")
#####################################################################
funplot5<-function(){
plot(rowSums(NoORF)[order],main="No. of Rep for (Single+Double) Versus Order of Sig Inter",ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",pch=19,cex=0.5)
points(rowSums(NoORF)[order][1:sig],ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",col=2,pch=19,cex=0.5)
plot(rowSums(NoORF)[order],main="No. of Rep for (Single+Double) Versus Order of Sig Inter",ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",pch=19,cex=0.5)
points(rowSums(NoORF)[order][1:sig],ylab="Total Number of Repeats (Single+Double)",xlab="Order of Sig Interacting",col=2,pch=19,cex=0.5)
i=order
text(rowSums(NoORF)[order],gene[i],pos=4,offset=0.1,cex=0.4)
}
##########################################
print("Individual plots")
##########################################
funplot6<-function(){
limmin<-min(y,na.rm=TRUE)
limmax<-max(y,na.rm=TRUE)
vec<-order(1-delta)[1:sig];col=1
for (i in vec)
{
plot(y[,1,i],y[,2,i],main=paste(gene[i],"Interacting","Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,cex=0.5,col=i)
lines(c(-1000,10000),A2*c(-1000,10000),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",cex=0.5,col=i,lwd=2)
}
vec<-order(1-delta)[N:(N-sig)];col=8
for (i in vec)
{
plot(y[,1,i],y[,2,i],main=paste(gene[i],"Non Interacting","Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",pch=19,cex=0.5,col=i)
lines(c(-1000,10000),A2*c(-1000,10000),main=paste("Treatment",treat,"Degrees"),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",lwd=2,cex=0.5,col=i)
}}

pdf(paste("Plots_",work,".pdf",sep=""))
print("1")
funplot1()
print("2")
funplot2()
print("3")
funplot3()
dev.off()
pdf(paste("Plots_Diag_",work,".pdf",sep=""))
print("4")
funplot4()
print("5")
funplot5()
dev.off()
pdf(paste("Plots_Indiv_",work,".pdf",sep=""))
print("6")
funplot6()
dev.off()
print(paste("Plots_",work,".pdf",sep=""))
print(paste("Plots_Diag_",work,".pdf",sep=""))
print(paste("Plots_Indiv_",work,".pdf",sep=""))
}
