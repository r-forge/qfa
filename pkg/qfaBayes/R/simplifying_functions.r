
SHM_postpro<-function(a,Treat,Screen,MPlate)
{
a<-funcREMOVE(a,Screen,Treat,MPlate)
a<-a[!a$Row==1,]
a<-a[!a$Row==16,]
a<-a[!a$Col==1,]
a<-a[!a$Col==24,]

Row<-paste(a$Row)
Col<-paste(a$Col)
for (i in 1:nrow(a)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

a$ID<-paste(a$Barcode,a$MasterPlate.Number,Row,Col,sep="")

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]

IDuni<-unique(a$ID)
ORFuni=unique(a$ORF)

gene<-unlist(lapply(ORFuni,funcGENE,data=a))

N<-length(ORFuni);M<-length(IDuni)
NoORF_a<-unlist(lapply(ORFuni,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))
dimr<-max(NoORF_a);dimc<-max(NoTime_a)
y<-funcXY(a$Growth,M,N,NoTime_a,NoSum_a,dimr,dimc)
x<-funcXY(a$Expt.Time,M,N,NoTime_a,NoSum_a,dimr,dimc)

QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),"N"=N,"M"=M,"gene"=gene)
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)
x[is.na(x)]=-999
y[is.na(y)]=-999
xx<-aperm(x,c(2,1,3))
yy<-aperm(y,c(2,1,3))

QFA.x=as.double(xx)
QFA.y=as.double(yy)
QFA.NoORF=as.integer(NoORF_a)
QFA.NoTIME=as.integer(c(NoTime_a)[-1])
QFA.I=as.integer(c(N,max(NoORF_a),max(NoTime_a),length(y),length(NoTime_a[-1])))
list(QFA.I=QFA.I,QFA.y=QFA.y,QFA.x=QFA.x,QFA.NoORF=QFA.NoORF,QFA.NoTIME=QFA.NoTIME,gene=gene)
}

SHM_main <- function(burn,iters,thin,CAPL,QFA.I,QFA.y,QFA.x,QFA.NoORF,QFA.NoTIME,PRIORS) {
L=min(CAPL,length(QFA.NoORF))
LM<-sum(QFA.NoORF[1:L])
NCOL=
LM+
L+
L+
1+
1+
1+
LM+
L+
L+
1+
1+
L+
1+
1+
1+
1+
1+
1
tmp <- .C("main", as.integer(burn),as.integer(iters),as.integer(thin),as.integer(L),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAI=as.integer(QFA.I),QFAy=as.double(QFA.y),QFAx=as.double(QFA.x),QFANoORF=as.integer(QFA.NoORF),QFANoTIME=as.integer(QFA.NoTIME),
PRIORS=as.double(PRIORS),PACKAGE="qfaBayes"
)
mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

IHM_main <- function(burn,iters,thin,QFA.IA,QFA.yA,QFA.NoORFA,QFA.IB,QFA.yB,QFA.NoORFB,PRIORS) {
aa<-QFA.NoORFA
bb<-QFA.NoORFB
if(!(length(aa)==length(bb))){stop()}
L=length(aa)
NCOL=
L+
1+
1+
L+
1+
1+
1+
L+
L+
1
tmp <- .C("main_IHM", as.integer(burn),as.integer(iters),as.integer(thin),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAIA=QFA.IA,QFAyA=as.double(QFA.yA),QFANoORFA=as.integer(QFA.NoORFA),QFAIB=as.integer(QFA.IB),QFAyB=as.double(QFA.yB),QFANoORFB=as.integer(QFA.NoORFB),
PRIORS=PRIORS,PACKAGE="qfaBayes"
)
mat=matrix((tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

plot_IHM_simple<-function(IHM_output,SHM){
N=SHM$QFA.I[1]
gene=SHM$gene
samp=IHM_output
if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-as.numeric(samp)}
namesamp<-names(vecsamp)
Z_l<-exp(vecsamp[1:(N)])
sigma_Z<-exp(vecsamp[N+1])
Z<-exp(vecsamp[N+2])
nu_l<-exp(vecsamp[(N+3):(2*N+2)])
sigma_nu<-exp(vecsamp[2*N+3])
nu<-exp(vecsamp[(2*N+4)])
A1<-exp(0)
A2<-exp(vecsamp[2*N+5])
delta<-vecsamp[(2*N+6):(3*N+5)]
gamma<-vecsamp[(3*N+6):(4*N+5)]
sigma_gamma<-exp(vecsamp[(4*N+6)])
if(nrow(samp)>1) {delta_gamma<-colMeans(samp[,(2*N+6):(3*N+5)]*samp[,(3*N+6):(4*N+5)])} else {delta_gamma<-colMeans(samp[,(2*N+6):(3*N+5)]*(samp[,(3*N+6):(4*N+5)]))}
delta_gamma=exp(delta_gamma)

sig<-sum(rep(1,N)[delta>0.5])
order<-order(1-delta)
vecorder<-order(1-delta)[1:sig]

limmin<-0
limmax<-max(A2*Z_l*delta_gamma)
limmaxx<-max(A1*Z_l)
i=1:N
plot(1,type="n",main=expression(paste("Treatment",Treat,degree,"C"," (delta=Posterior Expectations)")),ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),xlab="Control Fitness (=exp(Z_l))",ylab="Query Fitness (=exp(alpha+Z_l+delta_l*gamma_l))",col=8,pch=19,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
lines(A1*c(-1000,10000),A2*c(-1000,10000),col="grey",lwd=2)
points(A1*Z_l[i], A2*(Z_l[i]*delta_gamma[i]),col=8,pch=19,cex=0.5)
i=vecorder[log(delta_gamma)[vecorder]>0]
points(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),col=2,pch=19,cex=0.5)
i=vecorder[log(delta_gamma)[vecorder]<=0]  
points(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),col=3,pch=19,cex=0.5)
i=vecorder
text(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}



JHM_postpro<-function(a,TreatA,Screen_a,MPlate_a,b,TreatB,Screen_b,MPlate_b)
{
a<-funcREMOVE(a,Screen_a,TreatA,MPlate_a)

a<-a[!a$Row==1,]
a<-a[!a$Row==16,]
a<-a[!a$Col==1,]
a<-a[!a$Col==24,]

Row<-a$Row
Col<-a$Col
for (i in 1:nrow(a)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

a$ID<-paste(a$Barcode,a$MasterPlate.Number,Row,Col,sep="")

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]
ORFuni=unique(a$ORF)########
ORFuni_a<-unique(a$ORF)

b<-funcREMOVE(b,Screen_b,TreatB,MPlate_b)
b<-b[!b$Row==1,]
b<-b[!b$Row==16,]
b<-b[!b$Col==1,]
b<-b[!b$Col==24,]

Row<-b$Row
Col<-b$Col
for (i in 1:nrow(b)){
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

b$ID<-paste(b$Barcode,b$MasterPlate.Number,Row,Col,sep="")
b<-b[order(b$ORF,b$ID,b$Expt.Time), ]
ORFuni_b<-unique(b$ORF)

ORFuni<-unique(b$ORF)

IDuni<-unique(a$ID)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))

N<-length(ORFuni);M=Ma=length(IDuni)
NoORF_a<-unlist(lapply(ORFuni_a,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))

IDuni<-unique(b$ID)

N<-length(ORFuni);M=Mb=length(IDuni)
NoORF_b<-unlist(lapply(ORFuni,funcNoORF,data=b))#no of repeats each orf
NoTime_b<-c(0,unlist(lapply(IDuni,funcNoTime,data=b)))# 0+ no of time each repeat
NoSum_b<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_b)))

dimr<-max(NoORF_a,NoORF_b);dimc<-max(NoTime_a,NoTime_b)

y<-funcXY_J(a$Growth,b$Growth,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)
x<-funcXY_J(a$Expt.Time,b$Expt.Time,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)

QFA.I<-list("NoORF"=cbind(NoORF_a,NoORF_b),"NoTime_a"=NoTime_a[-1],"NoTime_b"=NoTime_b[-1],"NoSum"=cbind(NoSum_a,NoSum_b),"N"=N,"Ma"=Ma,"Mb"=Mb,"gene"=gene,SHIFT=c(0,max(NoSum_a,NoSum_b))
)
QFA.D<-list(x=x,y=y)

x[is.na(x)]=-999
y[is.na(y)]=-999
xx_a<-aperm(x[,,,1],c(2,1,3))
yy_a<-aperm(y[,,,1],c(2,1,3))

xx_b<-aperm(x[,,,2],c(2,1,3))
yy_b<-aperm(y[,,,2],c(2,1,3))

list(QFA.IA=c(N,max(NoORF_a),max(NoTime_a),length(y)/2,length(NoTime_a[-1])),
QFA.yA=c(yy_a),QFA.xA=c(xx_a), QFA.NoORFA=c(NoORF_a),QFA.NoTIMEA=c(NoTime_a)[-1],
QFA.IB=c(N,max(NoORF_b),max(NoTime_b),length(y)/2,length(NoTime_b[-1])),
QFA.yB=c(yy_b),QFA.xB=c(xx_b), QFA.NoORFB=c(NoORF_b),QFA.NoTIMEB=c(NoTime_b)[-1],gene=gene
)
}

JHM_main<- function(burn,iters,thin,QFA.IA,QFA.yA,QFA.xA,QFA.NoORFA,QFA.NoTIMEA,QFA.IB,QFA.yB,QFA.xB,QFA.NoORFB,QFA.NoTIMEB,PRIORS) {
aa<-QFA.NoORFA
bb<-QFA.NoORFB
if(!(length(aa)==length(bb))){stop()}
L=length(aa)
LMa<-sum(aa)
LMb<-sum(bb)
NCOL=
LMa+LMb+
2*L+
L+
1+
1+
1+
LMa+LMb+
2*L+
L+
1+
1+
L+
1+
1+
1+
1+
L+
L+
1+
L+
1+
2*2+
2*2
tmp <- .C("main_JHM", as.integer(burn),as.integer(iters),as.integer(thin),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAIA=as.integer(QFA.IA),QFAy=as.double(QFA.yA),QFAxA=as.double(QFA.xA),QFANoORFA=as.integer(QFA.NoORFA),QFANoTIMEA=as.integer(QFA.NoTIMEA),
QFAIB=as.integer(QFA.IB),QFAy=as.double(QFA.yB),QFAxB=as.double(QFA.xB),QFANoORFB=as.integer(QFA.NoORFB),QFANoTIMEB=as.integer(QFA.NoTIMEB),
PRIORS=PRIORS,PACKAGE="qfaBayes"
)
mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

plot_JHM_simple<-function(JHM_output,JHM){
samp<-JHM_output
gene=JHM$gene
L<-JHM$QFA.IA[1]
M=sum(c(JHM$QFA.NoORFA,c(JHM$QFA.NoORFB)))

K_clm=tau_K_cl=K_o_l=sigma_K_o=K_p=P=r_clm=tau_r_cl=r_o_l=sigma_r_o=r_p=nu_l=sigma_nu=nu_p=alpha_c=beta_c=delta_l=gamma_cl=sigma_gamma=omega_cl=sigma_omega=upsilon_c=sigma_upsilon=0
####
t=1
#K_clm
for (i in 1:c(M))
{
j=i
K_clm[t]=mean(samp[,j]);t=t+1
}

t=1
#tau_K_cl
j=M+1
for (i in (2*M+9*L+15):(2*M+11*L+14))
{
tau_K_cl[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#K_o_l
j=M+2*L+1
for (i in (M+1):(M+L))
{
K_o_l[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#sigma_K_o
i=2*M+9*L+9
j=M+3*L+1
sigma_K_o=mean(samp[,j])

t=1
#K_p
i=M+L+1
j=M+3*L+2
K_p=mean(samp[,j])

t=1
#P
i=M+L+2
j=M+3*L+3
P=mean(samp[,j])

t=1
#r_clm
j=M+3*L+4
for (i in (M+8*L+8):(2*M+8*L+7))
{
r_clm[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#tau_r_cl
j=2*M+3*L+4
for (i in (2*M+11*L+15):(2*M+13*L+14))
{
tau_r_cl[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#r_o_l
j=2*M+5*L+4
for (i in (2*M+8*L+8):(2*M+9*L+7))
{
r_o_l[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#sigma_r_o
i=2*M+9*L+13
j=2*M+6*L+4
sigma_r_o=mean(samp[,j])

t=1
#r_p
i=2*M+9*L+8
j=2*M+6*L+5
r_p=mean(samp[,j])


t=1
#nu_l
j=2*M+6*L+6
for (i in (M+5*L+7):(M+6*L+6))
{
nu_l[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#sigma_nu
i=2*M+9*L+11
j=2*M+7*L+6
sigma_nu=mean(samp[,j])

t=1
#nu_p
i=M+6*L+7
j=2*M+7*L+7
nu_p=mean(samp[,j])

t=1
#alpha_c
i=M+L+4
j=2*M+7*L+8
alpha_c=mean(samp[,j])

t=1
#beta_c
i=M+L+6
j=2*M+7*L+9
beta_c=mean(samp[,j])

t=1
#delta_l
j=2*M+7*L+10
for (i in (M+2*L+7):(M+3*L+6))
{
delta_l[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#gamma_cl
j=2*M+8*L+10
for (i in (M+4*L+7):(M+5*L+6))
{
gamma_cl[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#sigma_gamma
i=2*M+9*L+10
j=2*M+9*L+10
sigma_gamma=mean(samp[,j])

t=1
#omega_cl
j=2*M+9*L+11
for (i in (M+7*L+8):(M+8*L+7))
{
omega_cl[t]=mean(samp[,j]);t=t+1
j=j+1
}

t=1
#sigma_omega
i=2*M+9*L+12
j=2*M+10*L+11
sigma_omega=mean(samp[,j])

#upsilon_c
i=2*M+13*L+16
j=2*M+10*L+12
upsilon_c=mean(samp[,j])


t=1
#sigma_upsilon
i=2*M+9*L+14
j=2*M+10*L+13
sigma_upsilon=mean(samp[,j])



K<-exp(K_p)
K_i<-exp(K_o_l)
K_ij<-exp(K_clm)
PO<-exp(P)
r<-exp(r_p)
r_i<-exp(r_o_l)
r_ij<-exp(r_clm)

taui<-exp(nu_l)
tau<-exp(nu_p)
gam<-gamma_cl
omega<-omega_cl
nuc<-exp(upsilon_c)

gamdelt=0
j=2*M+7*L+10
jj=2*M+8*L+10
ii=M+4*L+7
t=1
for (i in (M+2*L+7):(M+3*L+6))
{
gamdelt[t]=mean(samp[,j]*samp[,jj]);t=t+1
j=j+1
ii=i+1
jj=jj+1
}

omegadelt=0
j=2*M+7*L+10
jj=2*M+9*L+11
ii=M+7*L+8
t=1
for (i in (M+2*L+7):(M+3*L+6))
{
omegadelt[t]=mean(samp[,j]*samp[,jj]);t=t+1
j=j+1
ii=i+1
jj=jj+1
}

delta<-delta_l

A1<-1
A2<-exp(alpha_c)
B1<-1
B2<-exp(beta_c)
sig<-sum(rep(1,L)[delta>0.5])
order<-order(1-delta)
vecorder<-order(1-delta)[1:sig]

K_ij=vecK=c(exp(K_o_l),A2*exp(K_o_l+gamdelt))
r_ij=c(exp(r_o_l),B2*exp(r_o_l+omegadelt))

K_ij[K_ij<2*PO]=PO+0.001

vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
vecMDPa<-log(K_ij/PO)/log(2) #MDP
vecMDRPMDR<-vecMDRa*vecMDPa
vecMDRPMDR[vecK<2*PO]=0
mu_a=(vecMDRPMDR)[1:L]
mu_b=(vecMDRPMDR)[(1+L):(2*L)]

limmin<-0
limmax<-max(na.omit(c(mu_b)))
limmaxx<-max(na.omit(c(mu_a)))

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)

vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
vecMDPa<-log(K_ij/PO)/log(2) #MDP
vecMDPa*vecMDRa

i=1:L
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]>0]#######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)

mu_a=exp(K_o_l)#####
mu_b=A2*exp(K_o_l+gamdelt)#####
limmin<-0
limmax<-max(na.omit(c(mu_b)))
limmaxx<-max(na.omit(c(mu_a)))

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
lines(c(-1000,1000),A2*c(-1000,1000),col="grey",lty=2)########
i=1:L
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]>0]####
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)

mu_a=exp(r_o_l)#####
mu_b=B2*exp(r_o_l+omegadelt)#####
limmin<-0
limmax<-max(na.omit(c(mu_b)))
limmaxx<-max(na.omit(c(mu_a)))

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
lines(c(-1000,1000),B2*c(-1000,1000),col="grey",lty=2)#######
i=1:L
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]>0]#######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}
