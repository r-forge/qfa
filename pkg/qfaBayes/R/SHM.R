#### Hierachical Logistic Curve Model ####
qfa.Hierachical<-function(experiment,Scaling,iter,upd,thin,PlotOutput=TRUE,work,CustomModel=FALSE){
a<-experiment
a<-funcIDORDER(a)
IDuni<-unique(a$ID)
ORFuni<-unique(a$ORF)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))
N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))
dimr<-max(NoORF_a);dimc<-max(NoTime_a)
y<-funcXY(a$Growth,M,N,NoTime_a,NoSum_a,dimr,dimc)
x<-funcXY(a$Expt.Time,M,N,NoTime_a,NoSum_a,dimr,dimc)
QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),"N"=N,"M"=M,"gene"=gene)
if (Scaling==TRUE){y<-funcSCALING(a,y)}
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)
if (!(CustomModel==FALSE)){source(CustomModel)} else {funcMODELHierarchical()}
QFA.P<-funcPRIORS(CustomModel)
samp<-funcFITandUPDATE(QFA.I,QFA.D,QFA.P)
QFA.O<-funcPosterior(samp,N,M,iter,thin,upd)
QFA<-c(QFA.O,QFA.I,QFA.D,QFA.P)
if(PlotOutput==TRUE){qfaplots.H(QFA,work)}
return(QFA)
}

### Hierachical Logistic Curve Model Plots to Pdf###
qfaplots.H<-function(QFA,work,LinearGaussian=FALSE){

samp<-QFA$samp
iter<-QFA$iter
thin<-QFA$thin

y<-QFA$y
x<-QFA$x

N<-QFA$N
M<-QFA$M
NoSum<-QFA$NoSum
NoORF<-QFA$NoORF
NoTime<-QFA$NoTime
gene<-QFA$gene

K_s<-QFA$K_s
r_s<-QFA$r_s
PO_s<-QFA$PO_s
beta<-QFA$beta
tau_s<-QFA$tau_s
delta<-QFA$delta
alpha_ij_sd<-QFA$alpha_ij_sd
gamma_ij_sd<-QFA$gamma_ij_sd
alpha<-QFA$alpha
gamma<-QFA$gamma
alpha_i<-QFA$alpha_i
gamma_i<-QFA$gamma_i
alpha_ij<-QFA$alpha_ij
gamma_ij<-QFA$gamma_ij

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


if(LinearGaussian==TRUE){
K<-exp(QFA$K)
K_i<-exp(QFA$K_i)
K_ij<-QFA$K_ij
r<-exp(QFA$r)
r_i<-exp(QFA$r_i)
}

################################################
print("Plots")
################################################

ylimmin<-0
ylimmax<-max(na.omit(as.numeric(y)))
xlimmin<-0
xlimmax<-max(na.omit(as.numeric(x)))

pdf(paste("Plots_M",work,".pdf",sep=""))
################################################
print("Master Curve")
################################################
plot(x,y,main="Master Curve",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
curve((K*PO*exp(r*x))/(K+PO*(exp(r*x)-1)), 0, 8,add=TRUE,col=1)
################################################
print("ORF Curves")
################################################
plot(-1,-1,main="ORF Curves",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
for (i in 1:N){
points(x[,,i],y[,,i],main="ORF Curves",xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax),col=i)
}
for (i in 1:N)
{
curve((K_i[i]*PO*exp(r_i[i]*x))/(K_i[i]+PO*(exp(r_i[i]*x)-1)), 0, 8,add=TRUE,col=i) 
}
################################################
print("Repeat Curves")
################################################
plot(x,y,main="Repeat Curves", xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
for (i in 1:M)
{
curve((K_ij[i]*PO*exp(r_ij[i]*x))/(K_ij[i]+PO*(exp(r_ij[i]*x)-1)), 0, 8,add=TRUE,col=i) 
}

################################################
print("Model Variation tau")#fix
################################################
plot(x,y,main="Curve variation tau_m", xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
KK=K
rr=r
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)),xlimmin, xlimmax,add=TRUE,col=1) 
if (LinearGaussian==FALSE){KK=funcCurveVarK(K,tau,1)} else {KK=funcCurveVarK_LG(K,tau,1)}
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=3) 
if (LinearGaussian==FALSE){KK=funcCurveVarK(K,tau,-1)} else {KK=funcCurveVarK_LG(K,tau,-1)}
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=3)

################################################
print("Model Variation posterior")
################################################

par(mfrow=c(4,2))
for (i in 1:N)
{
if (LinearGaussian==FALSE){MVP<-funcModelVarPost(QFA)} else {MVP<-funcModelVarPost_LG(QFA)}
plot(density(MVP[,1]),main="Density Curve variation", xlab="Time (days)", ylab="Culture Density (AU)")
lines(density(MVP[,2]),main="Master Curve variation", xlab="Time (days)", ylab="Culture Density (AU)",col=2)
lines(density(MVP[,3]),main="Master Curve variation", xlab="Time (days)", ylab="Culture Density (AU)",col=3)
plot(density(MVP[,4]),main="Density Curve variation", xlab="Time (days)", ylab="Culture Density (AU)")
lines(density(MVP[,5]),main="Master Curve variation", xlab="Time (days)", ylab="Culture Density (AU)",col=2)
lines(density(MVP[,6]),main="Master Curve variation", xlab="Time (days)", ylab="Culture Density (AU)",col=3)
}
dev.off()

pdf(paste("Plots_M_indiv",work,".pdf",sep=""))
###########################################
print("plots for individual Logistic curve fits")#fix
###########################################
for (i in 1:N)
{
plot(-1,-1,main=paste(gene[i],"Curve"),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])
KK=K_i[i]
rr=r_i[i]
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1) 
if (LinearGaussian==FALSE){KK=funcCurveVarK(K_i[i],k_tau[i],1)} else {KK=funcCurveVarK_LG(K_i[i],tau,1)}
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmax-0.8*xlimmax, xlimmax,add=TRUE,col=3) 
if (LinearGaussian==FALSE){KK=funcCurveVarK(K_i[i],k_tau[i],-1)} else {KK=funcCurveVarK_LG(K_i[i],tau,-1)}
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmax-0.8*xlimmax, xlimmax,add=TRUE,col=3) 
plot(-1,-1,main=paste(gene[i],"Repeat Curves"),xlab="Time (days)", ylab="Culture Density (AU)",xlim=c(xlimmin,xlimmax),ylim=c(ylimmin,ylimmax))
points(x[,,i],y[,,i])
for (j in 1:NoORF[i])
{
KK=K_ij[(j+NoSum[i])]
rr=r_ij[(j+NoSum[i])]
curve((KK*PO*exp(rr*x))/(KK+PO*(exp(rr*x)-1)), xlimmin, xlimmax,add=TRUE,col=1+j) 
}
}
dev.off()

pdf(paste("Plots_M_diag",work,".pdf",sep=""))
###########################################
print("Prior density")
###########################################
par(mfrow=c(4,2))
sampsize<-round(iter/thin)
if (LinearGaussian==FALSE){den<-funcDen(sampsize,QFA)}else den<-funcDen_LG(sampsize,QFA)
namesampden<-unique(substring(namesamp,1,4))
for (i in 1:ncol(den))
{
plot(density(den[,i]),paste(namesampden[i],"Prior Density"))
}

###########################################
print("Diagnostics trace acf density")
###########################################
if (LinearGaussian==FALSE){postpred<-funcPostPred(sampsize,QFA)} else postpred<-funcPostPred_LG(sampsize,QFA)


par(mfrow=c(4,4))
for (i in 1:length(namesamp))
{
post<-density(as.numeric(samp[,i]))
pred<-density(postpred[,i])
plot(as.numeric(samp[,i]),main=paste(namesamp[i],"Trace Top"),type="l")
t<-(1:ncol(den))[namesampden==substring(namesamp[i],1,4)]
pri<-density(den[,t])
plot(post,main=paste(namesamp[i],"Density"),xlim=c(min(pri$x),max(pri$x)),
ylim=c(min(post$y,pri$y),max(post$y,pri$y)))
lines(pri,col=2)
acf(as.numeric(samp[,i]),main=paste(namesamp[i],"ACF"))
plot(post,main=paste(namesamp[i],"Density"),xlim=c(min(post$x,pred$x),max(post$x,pred$x)),
ylim=c(min(post$y,postpred),max(post$y,postpred)))
lines(pred,main=paste(namesamp[i],"Density PostPred"),col=2)
}
dev.off()
}
