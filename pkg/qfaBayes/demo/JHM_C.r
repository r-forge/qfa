#a=read.delim("~/QFADatasets/SHM/URA3_Raw/data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("URA3_Raw_trim")#Control a 

#b=read.delim("~/QFADatasets/SHM/YKU70_Raw/data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("CDC13-1_Raw_trim")#Query b 

TreatA=27
TreatB=27

filename=paste("M_JHM_demo","_",TreatA,"_",TreatB,sep="")

qfa.variables(a)

Screen<-as.character(unique(a$Screen.Name))
MPlate<-as.character(unique(a$MasterPlate.Number))

a<-funcREMOVE(a,Screen,TreatA,MPlate)

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

qfa.variables(b)
Screen<-as.character(unique(b$Screen.Name))
MPlate<-as.character(unique(b$MasterPlate.Number))
b<-funcREMOVE(b,Screen,TreatB,MPlate)

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

####


if(!(sum(rep(1,length(ORFuni))[ORFuni==ORFuni_b])/length(ORFuni)==1)){
print("ORF names differ!")
print(ORFuni[!(ORFuni_b==ORFuni)])
print("ORF names differ!")
stop()
}

if(max(a$Growth,b$Growth)>1){
print("Data not scaled appropriately")
stop()
}
####




ORFuni<-unique(b$ORF)

IDuni<-unique(a$ID)
gene<-unlist(lapply(ORFuni,funcGENE,data=a))
###
if(sum(gene=="0")>0){#Data Correction
gene<-as.character(gene)
gene[gene=="0"]=ORFuni[gene=="0"]
}
###

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
xx<-aperm(x[,,,1],c(2,1,3))
yy<-aperm(y[,,,1],c(2,1,3))
write.table(file="xdata_A2.txt",c(xx))
write.table(file="ydata_A2.txt",c(yy))

write.table(file="NoORFdata_A2.txt",c(NoORF_a))
write.table(file="NoTIMEdata_A2.txt",c(NoTime_a)[-1])
write.table(file="LMNmaxdata_A2.txt",c(N,max(NoORF_a),max(NoTime_a),length(y)/2,length(NoTime_a[-1])))

xx<-aperm(x[,,,2],c(2,1,3))
yy<-aperm(y[,,,2],c(2,1,3))
write.table(file="xdata_B2.txt",c(xx))
write.table(file="ydata_B2.txt",c(yy))

write.table(file="NoORFdata_B2.txt",c(NoORF_b))
write.table(file="NoTIMEdata_B2.txt",c(NoTime_b)[-1])
write.table(file="LMNmaxdata_B2.txt",c(N,max(NoORF_b),max(NoTime_b),length(y)/2,length(NoTime_b[-1])))

save.image(paste(filename,".RData",sep=""))

#################################################
#You may use standalone C code for SHM from here
#################################################
main_JHM <- function(burn,iters,thin) {
aa<-read.table("LMNmaxdata_A2.txt",header=T)
QFA.IA=as.integer((aa)[[1]])
aa<-read.table("ydata_A2.txt",header=T)
QFA.yA=as.double((aa)[[1]])
aa<-read.table("xdata_A2.txt",header=T)
QFA.xA=as.double((aa)[[1]])
aa<-read.table("NoORFdata_A2.txt",header=T)
QFA.NoORFA=as.integer((aa)[[1]])
aa<-read.table("NoTIMEdata_A2.txt",header=T)
QFA.NoTIMEA=as.integer((aa)[[1]])
aa<-read.table("LMNmaxdata_B2.txt",header=T)
QFA.IB=as.integer((aa)[[1]])
aa<-read.table("ydata_B2.txt",header=T)
QFA.yB=as.double((aa)[[1]])
aa<-read.table("xdata_B2.txt",header=T)
QFA.xB=as.double((aa)[[1]])
aa<-read.table("NoORFdata_B2.txt",header=T)
QFA.NoORFB=as.integer((aa)[[1]])
aa<-read.table("NoTIMEdata_B2.txt",header=T)
QFA.NoTIMEB=as.integer((aa)[[1]])
#priors_JHM=read.table("priors.txt",header=T)
data("priors_JHM")
PRIORS=as.double((priors_JHM)[[1]])
PRIORS[19]=0##
aa<-read.table("NoORFdata_A2.txt",header=T)
bb<-read.table("NoORFdata_B2.txt",header=T)
if(!(nrow(aa)==nrow(bb))){stop()}
L=nrow(aa)
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
QFAIA=QFA.IA,QFAy=QFA.yA,QFAxA=QFA.xA,QFANoORFA=QFA.NoORFA,QFANoTIMEA=QFA.NoTIMEA,
QFAIB=QFA.IB,QFAy=QFA.yB,QFAxB=QFA.xB,QFANoORFB=QFA.NoORFB,QFANoTIMEB=QFA.NoTIMEB,
PRIORS=PRIORS
)
mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

#Change the following variables
burn=1#Burn in period
iters=1# sample iterations
thin=1# thining for sample
D<-main_JHM(burn,iters,thin)

plotYN=0
while(plotYN < 1 ){
  n<-readline("do you wish to plot? Y or N: ")
if(n=="Y"){plotYN=1}
if(n=="N"){stop()}
}


################################
load("M_JHM_demo_27_27.RData")
QFA.P<-read.table("priors.txt",header=T)
QFA<-c(QFA.I,QFA.P,QFA.D)
samp<-D
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
SHIFT<-QFA$SHIFT
###
L=N #4294
M=sum(NoORF)
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


########strip 
#strip=TRUE
#if(strip==TRUE){
#strip_ORF<-read.delim("~/strip_list.txt",header=T,sep="\t")$orf
#K_o_l<-K_o_l[!ORFuni%in%strip_ORF]
#r_o_l<-r_o_l[!ORFuni%in%strip_ORF]
#delta<-delta[!ORFuni%in%strip_ORF]
#omegadelt<-omegadelt[!ORFuni%in%strip_ORF]
#gamdelt<-gamdelt[!ORFuni%in%strip_ORF]
#gene<-gene[!ORFuni%in%strip_ORF]
#ORFuni<-ORFuni[!ORFuni%in%strip_ORF]
#N<-length(ORFuni)}
############

A1<-1
A2<-exp(alpha_c)
B1<-1
B2<-exp(beta_c)
sig<-sum(rep(1,N)[delta>0.5])
order<-order(1-delta)
vecorder<-order(1-delta)[1:sig]

K_ij=vecK=c(exp(K_o_l),A2*exp(K_o_l+gamdelt))####
r_ij=c(exp(r_o_l),B2*exp(r_o_l+omegadelt))####

K_ij[K_ij<2*PO]=PO+0.001

vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
vecMDPa<-log(K_ij/PO)/log(2) #MDP
vecMDRPMDR<-vecMDRa*vecMDPa
vecMDRPMDR[vecK<2*PO]=0
mu_a=(vecMDRPMDR)[1:N]
mu_b=(vecMDRPMDR)[(1+N):(2*N)]

limmin<-0
limmax<-max(na.omit(c(mu_b)))
limmaxx<-max(na.omit(c(mu_a)))

Treat=TreatA=TreatB=27

gene[gene==0]=ORFuni[gene==0]#correction

file="DEMO"
pdf(paste("JHM_plot_",file,".pdf",sep=""),useDingbats=F)

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)

#lines(lm(mu_b~0+mu_a),col="grey",lty=2)#############



vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
vecMDPa<-log(K_ij/PO)/log(2) #MDP
vecMDPa*vecMDRa
#lines(c(-1000,1000),A2*c(-1000,1000),col="grey",lty=2)########



lines(c(mu_a[gene=="HIS3"],mu_a[gene=="HIS3"]),c(-1000,1000),lwd=2,col="lightblue")
lines(c(-1000,1000),c(mu_b[gene=="HIS3"],mu_b[gene=="HIS3"]),lwd=2,col="lightblue")
i=1:N
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]>0]#######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
#legend(1,limmax, c("1-1","simple Lin Reg"), cex=0.5,col=c("grey","grey","black"), lty=c(1,2,3,1))
#if (sum((1:N)[gene=="HIS3"])==1){
#legend(1,limmax, c("1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("grey","grey","black"), lty=c(2,3,1))
#}


mu_a=exp(K_o_l)#####
mu_b=A2*exp(K_o_l+gamdelt)#####
limmin<-0
#limmax<-max(na.omit(Mu))
limmax<-max(na.omit(c(mu_b)))
limmaxx<-max(na.omit(c(mu_a)))

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
lines(c(-1000,1000),A2*c(-1000,1000),col="grey",lty=2)########
lines(c(mu_a[gene=="HIS3"],mu_a[gene=="HIS3"]),c(-1000,1000),lwd=2,col="lightblue")
lines(c(-1000,1000),c(mu_b[gene=="HIS3"],mu_b[gene=="HIS3"]),lwd=2,col="lightblue")
i=1:N
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]>0]####
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[gamdelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
#legend(1,limmax, c("1-1","simple Lin Reg"), cex=0.5,col=c("grey","grey","black"), lty=c(1,2,3,1))
#if (sum((1:N)[gene=="HIS3"])==1){
#legend(1,limmax, c("1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("grey","grey","black"), lty=c(2,3,1))
#}

mu_a=exp(r_o_l)#####
mu_b=B2*exp(r_o_l+omegadelt)#####
limmin<-0
#limmax<-max(na.omit(Mu))
limmax<-max(na.omit(c(mu_b)))
limmaxx<-max(na.omit(c(mu_a)))

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
lines(c(-1000,1000),B2*c(-1000,1000),col="grey",lty=2)#######
lines(c(mu_a[gene=="HIS3"],mu_a[gene=="HIS3"]),c(-1000,1000),lwd=2,col="lightblue")
lines(c(-1000,1000),c(mu_b[gene=="HIS3"],mu_b[gene=="HIS3"]),lwd=2,col="lightblue")
i=1:N
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]>0]#######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
#legend(1,limmax, c("1-1","simple Lin Reg"), cex=0.5,col=c("grey","grey","black"), lty=c(1,2,3,1))
#if (sum((1:N)[gene=="HIS3"])==1){
#legend(1,limmax, c("1-1","simple Lin Reg","HIS3 Fit"), cex=0.5,col=c("grey","grey","black"), lty=c(2,3,1))
#}
dev.off()

stop()

gene[gene==0]=ORFuni[gene==0]
setwd("~/")
list<-read.table("Addinall_all.txt",header=T)
list[,1]<-as.character(list[,1])
list[,1][list[,1]=="YMR169c"]="YMR169C"
list[,1][list[,1]=="YMR175w"]="YMR175W"
list[,1][list[,1]=="YML009c"]="YML009C"
list$qvalue[is.na(list$qvalue)]=1
strip_ORF<-read.delim("strip_list.txt",header=T,sep="\t")$orf
list<-list[!(list[,1]%in%strip_ORF),]
list<-list[order(list[,1]),]
list$gene<-as.character(list$gene)
list[,5][is.na(list[,5])]=1
list[,6][is.na(list[,6])]=1
list$gene[is.na(list$gene)]<-list[,1][is.na(list$gene)]
#list<-list[abs(list[,2])>0.5,]##########
list2<-list
list<-unique(as.character(list[list[,6]<0.05,1]))

delta_gamma<-exp(gamdelt)
delta_omega<-exp(omegadelt)

lORF<-ORFuni[vecorder]
llORF<-lORF[lORF%in%list]
llORF_not<-lORF[!(lORF%in%list)]
lgene<-cbind(gene[ORFuni%in%llORF],as.numeric(delta[ORFuni%in%llORF]),as.numeric(delta_gamma[ORFuni%in%llORF]))
lgene_not<-cbind(gene[ORFuni%in%llORF_not],as.numeric(delta[ORFuni%in%llORF_not]),as.numeric(delta_gamma[ORFuni%in%llORF_not]))


###
K=A2*exp(K_o_l)*delta_gamma
r=B2*exp(r_o_l)*delta_omega
P<-PO
r[K<2*P]=0
K[K<2*P]=2*PO+0.001
FIT<-(r/log(2*((K-P)/(K-2*P))))*log(K/P)/log(2)
K=A2*exp(K_o_l)
r=B2*exp(r_o_l)
r[K<2*P]=0
K[K<2*P]=2*PO+0.001
FIT<-FIT-((r/log(2*((K-P)/(K-2*P))))*log(K/P)/log(2))
###

#
list2<-list2[order(abs(list2[,4]),decreasing=T),]
list2<-cbind(list2,1:nrow(list2))
#list2<-list2[order(list2[,1]),]
ADD_position<-list2[,8]




ORDER<-cbind(ORFuni[vecorder],gene[vecorder],as.numeric(delta[vecorder]),as.numeric(delta_gamma[vecorder]),as.numeric(delta_omega[vecorder]),FIT[vecorder],ADD_position[vecorder])
write.table(file="JHM_interactions.txt",ORDER)




ORDER[,4][order(abs(as.numeric(ORDER[,3])),decreasing=T)]
#
l<-list[!(list%in%lORF)]
l<-gene[ORFuni%in%l]
write.table(l,"JHM_not_interactions.txt")
write.table(cbind(ORFuni[order],gene[order],ORFuni[order],delta[order],delta_gamma[order],delta_omega[order],FIT[order],ADD_position[order]),"JHM_all.txt")

#Percent interactors
#
length(llORF)/(length(unique(c(as.character(list),lORF))))

length(listI[listI%in%lORF])/(length(unique(c(as.character(listI),lORF))))



