setwd(tempdir())
#a=read.delim("URA3_Raw.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("URA3_Raw_trim")

qfa.variables(a)
Treat=27
Screen<-unique(a$Screen.Name)
MPlate<-(unique(a$MasterPlate.Number))
filename=paste("SHM_demo","_",Treat,sep="")

a<-funcREMOVE(a,Screen,Treat,MPlate)
a<-a[!a$Row==1,]
a<-a[!a$Row==16,]
a<-a[!a$Col==1,]
a<-a[!a$Col==24,]

Row<-paste(a$Row)
Col<-paste(a$Col)
for (i in 1:nrow(a)){
  if (nchar(Row[i])<2){
    Row[i]=paste(0,Row[i],sep="")
  }
  if (nchar(Col[i])<2){
    Col[i]=paste(0,Col[i],sep="")
  }
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

QFA.I<-list("NoORF"=c(NoORF_a),"NoTime"=c(NoTime_a)[-1],"NoSum"=c(NoSum_a),
  "N"=N,"M"=M,"gene"=gene)
QFA.D<-list(y=y,x=x,ORFuni=ORFuni)

x[is.na(x)]=-999
y[is.na(y)]=-999
xx<-aperm(x,c(2,1,3))
yy<-aperm(y,c(2,1,3))
write.table(file="xdata.txt",c(xx))
write.table(file="ydata.txt",c(yy))

write.table(file="NoORFdata.txt",c(NoORF_a))
write.table(file="NoTIMEdata.txt",c(NoTime_a)[-1])
write.table(file="LMNmaxdata.txt",c(N,max(NoORF_a),max(NoTime_a),length(y),
  length(NoTime_a[-1])))

save.image(paste(filename,".RData",sep=""))
#################################################
#You may use standalone C code for SHM from here
#################################################
aa<-read.table("LMNmaxdata.txt",header=T)
QFA.I=as.integer((aa)[[1]])
aa<-read.table("ydata.txt",header=T)
QFA.y=as.double((aa)[[1]])
aa<-read.table("xdata.txt",header=T)
QFA.x=as.double((aa)[[1]])
aa<-read.table("NoORFdata.txt",header=T)
QFA.NoORF=as.integer((aa)[[1]])
aa<-read.table("NoTIMEdata.txt",header=T)
QFA.NoTIME=as.integer((aa)[[1]])
data("priors_SHM")
#priors_SHM=read.table("priors.txt",header=T)
PRIORS=as.double((priors_SHM)[[1]])[1:18]
data("tuning_SHM")
TUNING=as.double((tuning_SHM))[1:8]

main <- function(burn,iters,thin,CAPL) {
aa<-read.table("NoORFdata.txt",header=T)
L=min(CAPL,nrow(aa))
LM<-sum(aa[1:L,])
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
tmp <- .C("main", as.integer(burn),as.integer(iters),as.integer(thin),
  as.integer(L),OUT=as.double(1:(NCOL*iters)),
  HEADER=as.character(rep("NULLNULL",NCOL)),QFAI=QFA.I,QFAy=QFA.y,QFAx=QFA.x,
  QFANoORF=QFA.NoORF,QFANoTIME=QFA.NoTIME,PRIORS=PRIORS,TUNING=TUNING)
mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

#Change the following variables
burn=1#Burn in period
iters=1# sample iterations
thin=1# thining for sample
CAPL=1#maximum no. of ORF's
A<-main(burn,iters,thin,CAPL)


plotYN=0
while(plotYN < 1 ){
  n<-readline("do you wish to plot? Y or N: ")
  if(n=="Y"){
    plotYN=1
  }
  if(n=="N"){
    stop()
  }
}

####################

filename="DEMO"
#pdf(paste("SHM_plot_",filename,".pdf",sep=""),useDingbats=F)

load("M_SHM_FULL_27.RData")

samp<-A

#########
y<-QFA.D$y
x<-QFA.D$x
L=N=CAPL
NoSum<-NoSum_a
NoORF<-NoORF_a
NoTime<-NoTime_a[-1]

M=sum(NoORF[1:CAPL])

K_lm=tau_K_l=K_o_l=sigma_K_o=K_p=P_l=r_lm=tau_r_l=r_o_l=sigma_r_o=r_p=
  nu_l=nu_p=sigma_nu=0
aa<-samp
#K_lm[%i]
t=1
for (i in 1:M){
  j=i
  K_lm[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#tau_K_l[%i]
j=M+1
for (i in (2*M+3*N+8):(2*M+4*N+7)){
  tau_K_l[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#"K_o_l[%i] 
j=M+N+1
for (i in (M+1):(M+N)){
  K_o_l[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#sigma_K_o ");
i=2*M+3*N+5
j=M+2*N+1
sigma_K_o=mean(samp[,j])

t=1
#K_p ");
i=M+1+N
j=M+2*N+2
K_p=mean(samp[,j])

t=1
#"P_l ");
i=(M+N+2)
j=M+2*N+3
P_l=mean(samp[,j])

t=1
#r_lm[%i] 
j=M+2*N+4
for (i in (M+2*N+4):(2*M+2*N+3)){
  r_lm[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#tau_r_l[%i] ",l);
j=2*M+2*N+4
for (i in (2*M+4*N+8):(2*M+5*N+7)){
  tau_r_l[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#r_o_l[%i] ",l);
j=2*M+3*N+4
for (i in (2*M+2*N+4):(2*M+3*N+3)){
  r_o_l[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#sigma_r_o ");
i=2*M+3*N+7
j=2*M+4*N+4
sigma_r_o=mean(samp[,j]);

t=1
#r_p ");
i=2*M+3*N+4
j=2*M+4*N+5
r_p=mean(samp[,j]);

t=1
#"nu_l[%i] ",l);
j=2*M+4*N+6
for (i in (M+N+3):(M+2*N+2)){
  nu_l[t]=mean(samp[,j]);t=t+1
  j=j+1
}

t=1
#sigma_nu ");
i=2*M+3*N+6
j=2*M+5*N+6
sigma_nu=mean(samp[,j]);

t=1
#nu_p ");
i=M+2*N+3
j=2*M+5*N+7
nu_p=mean(samp[,j]);

###
K<-exp(K_p)
K_i<-exp(K_o_l)
K_ij<-exp(K_lm)
P<-exp(P_l)
r<-exp(r_p)
r_i<-exp(r_o_l)
r_ij<-exp(r_lm)
taui<-exp(nu_l)
tau<-exp(nu_p)
K_i_tau<-exp(sigma_K_o)
r_i_tau<-exp(sigma_r_o)
K_ij_tau<-exp(tau_K_l)
r_ij_tau<-exp(tau_r_l)

for (i in 1:N){
  ylimmax=max(y[,,i][!is.na(y[,,i])])
  xlimmax=max(x[,,i][!is.na(y[,,i])])
  plot(-1,-1,main=paste(gene[i],"Repeat Curves"),xlab="Time (days)",
    ylab="Culture Density (AU)",xlim=c(0,xlimmax),ylim=c(0,ylimmax))
  for (j in 1:NoORF[i]){
    points(x[j,,i],y[j,,i])
    KK=K_ij[(j+NoSum[i])]
    rr=r_ij[(j+NoSum[i])]
    curve((KK*P*exp(rr*x))/(KK+P*(exp(rr*x)-1)), 0, xlimmax,add=TRUE) 
  }
  K=exp(K_o_l[i])
  r=exp(r_o_l[i])
  curve((K*P*exp(r*x))/(K+P*(exp(r*x)-1)),lwd=3,col="red",add=T)
}
#dev.off()


