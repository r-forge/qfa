setwd(tempdir())
#a=read.delim("data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("URA3_Raw_trim")
qfa.variables(a)
Treat=27
Screen<-unique(a$Screen.Name)
MPlate<-(unique(a$MasterPlate.Number))
filename=paste("M_IHM_demo_a","_",Treat,sep="")

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

ORFuni=unique(a$ORF)

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]

IDuni<-unique(a$ID)
ORFuni=ORFuni_a=unique(a$ORF)

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
write.table(file="xdata_A.txt",c(xx))
write.table(file="ydata_A.txt",c(yy))

write.table(file="NoORFdata_A.txt",c(NoORF_a))
write.table(file="NoTIMEdata_A.txt",c(NoTime_a)[-1])
write.table(file="LMNmaxdata_A.txt",c(N,max(NoORF_a),max(NoTime_a),length(y),length(NoTime_a[-1])))

save.image(paste(filename,".RData",sep=""))
#################################################
#You may use standalone C code for SHM from here
#################################################
#a=read.delim("data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("CDC13-1_Raw_trim")
a<-b
qfa.variables(a)
Treat=27
Screen<-unique(a$Screen.Name)
MPlate<-(unique(a$MasterPlate.Number))
filename=paste("M_IHM_demo_b","_",Treat,sep="")

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

ORFuni=unique(a$ORF)

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]

IDuni<-unique(a$ID)
ORFuni=ORFuni_b=unique(a$ORF)

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
write.table(file="xdata_B.txt",c(xx))
write.table(file="ydata_B.txt",c(yy))

write.table(file="NoORFdata_B.txt",c(NoORF_a))
write.table(file="NoTIMEdata_B.txt",c(NoTime_a)[-1])
write.table(file="LMNmaxdata_B.txt",c(N,max(NoORF_a),max(NoTime_a),length(y),length(NoTime_a[-1])))

save.image(paste(filename,".RData",sep=""))
#################################################
#You may use standalone C code for SHM from here
#################################################

aa<-read.table("LMNmaxdata_A.txt",header=T)
QFA.I=as.integer((aa)[[1]])
aa<-read.table("ydata_A.txt",header=T)
QFA.y=as.double((aa)[[1]])
aa<-read.table("xdata_A.txt",header=T)
QFA.x=as.double((aa)[[1]])
aa<-read.table("NoORFdata_A.txt",header=T)
QFA.NoORF=as.integer((aa)[[1]])
aa<-read.table("NoTIMEdata_A.txt",header=T)
QFA.NoTIME=as.integer((aa)[[1]])
data("priors_SHM")
#priors_SHM=read.table("priors.txt",header=T)
PRIORS=as.double((priors_SHM)[[1]])[1:18]

main <- function(burn,iters,thin,CAPL) {
aa<-read.table("NoORFdata_A.txt",header=T)
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
tmp <- .C("main_SHM", as.integer(burn),as.integer(iters),as.integer(thin),as.integer(L),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAI=QFA.I,QFAy=QFA.y,QFAx=QFA.x,QFANoORF=QFA.NoORF,QFANoTIME=QFA.NoTIME,
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
CAPL=4294#maximum no. of ORF's
A<-main(burn,iters,thin,CAPL)
write.table(colMeans(A),"data_A.txt",col.names=FALSE,row.names=FALSE)
#################################################
#You may use standalone C code for SHM from here
#################################################
main <- function(burn,iters,thin,CAPL) {
aa<-read.table("LMNmaxdata_B.txt",header=T)
QFA.I=as.integer((aa)[[1]])
aa<-read.table("ydata_B.txt",header=T)
QFA.y=as.double((aa)[[1]])
aa<-read.table("xdata_B.txt",header=T)
QFA.x=as.double((aa)[[1]])
aa<-read.table("NoORFdata_B.txt",header=T)
QFA.NoORF=as.integer((aa)[[1]])
aa<-read.table("NoTIMEdata_B.txt",header=T)
QFA.NoTIME=as.integer((aa)[[1]])
data("priors_SHM")
#priors_SHM=read.table("priors.txt",header=T)
PRIORS=as.double((priors_SHM)[[1]])[1:18]
aa<-read.table("NoORFdata_B.txt",header=T)
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
tmp <- .C("main_SHM", as.integer(burn),as.integer(iters),as.integer(thin),as.integer(L),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAI=QFA.I,QFAy=QFA.y,QFAx=QFA.x,QFANoORF=QFA.NoORF,QFANoTIME=QFA.NoTIME,
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
CAPL=4294#maximum no. of ORF's
B<-main(burn,iters,thin,CAPL)
write.table(colMeans(B),"data_B.txt",col.names=FALSE,row.names=FALSE)
#################################################
#You may use standalone C code for SHM from here
#################################################
main_IHM <- function(burn,iters,thin) {
aa<-read.table("LMNmaxdata_A.txt",header=T)
QFA.IA=as.integer((aa)[[1]])
aa<-read.table("data_A.txt",header=F)
QFA.yA=as.double((aa)[[1]])
aa<-read.table("NoORFdata_A.txt",header=T)
QFA.NoORFA=as.integer((aa)[[1]])
aa<-read.table("LMNmaxdata_B.txt",header=T)
QFA.IB=as.integer((aa)[[1]])
aa<-read.table("data_B.txt",header=F)
QFA.yB=as.double((aa)[[1]])
aa<-read.table("NoORFdata_B.txt",header=T)
QFA.NoORFB=as.integer((aa)[[1]])
#priors_IHM=read.table("priors_IHM.txt",header=T)
data("priors_IHM")
PRIORS=as.double((priors_IHM)[[1]])
aa<-read.table("NoORFdata_A.txt",header=T)
bb<-read.table("NoORFdata_B.txt",header=T)
if(!(nrow(aa)==nrow(bb))){stop()}
L=nrow(aa)
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
QFAIA=QFA.IA,QFAyA=QFA.yA,QFANoORFA=QFA.NoORFA,
QFAIB=QFA.IB,QFAyB=QFA.yB,QFANoORFB=QFA.NoORFB,
PRIORS=PRIORS
)
mat=matrix((tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

#Change the following variables
burn=1#Burn in period
iters=1# sample iterations
thin=1# thining for sample
C<-main_IHM(burn,iters,thin)

plotYN=0
while(plotYN < 1 ){
  n<-readline("do you wish to plot? Y or N: ")
if(n=="Y"){plotYN=1}
if(n=="N"){stop()}
}

###
load("M_IHM_demo_a_27.RData")
samp=C
if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-as.numeric(samp)}
namesamp<-names(vecsamp)
#write.table(samp,"backup.txt")
#write.table(vecsamp,"backup2.txt")
Z_l<-exp(vecsamp[1:(N)])
sigma_Z<-exp(vecsamp[N+1])
Z<-exp(vecsamp[N+2])
nu_cl<-exp(vecsamp[(N+3):(3*N+2)])
sigma_nu<-exp(vecsamp[3*N+3])
nu<-exp(vecsamp[(3*N+4)])
A1<-exp(0)
A2<-exp(vecsamp[3*N+5])
delta<-vecsamp[(3*N+6):(4*N+5)]
gamma<-vecsamp[(4*N+6):(5*N+5)]
sigma_gamma<-exp(vecsamp[(5*N+6)])
if(nrow(samp)>1) {delta_gamma<-colMeans(samp[,(3*N+6):(4*N+5)]*samp[,(4*N+6):(5*N+5)])} else {delta_gamma<-colMeans(samp[,(3*N+6):(4*N+5)]*(samp[,(4*N+6):(5*N+5)]))}
delta_gamma=exp(delta_gamma)

sig<-sum(rep(1,N)[delta>0.5])
order<-order(1-delta)
vecorder<-order(1-delta)[1:sig]

file="DEMO"#

#pdf(paste("IHM_plot_",file,".pdf",sep=""),useDingbats=F)

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
 plot(density((Z_l)))
 plot(density(A2*(Z_l*delta_gamma)))
#dev.off()


stop()

#
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


lORF<-ORFuni[vecorder]
llORF<-lORF[lORF%in%list]
llORF_not<-lORF[!(lORF%in%list)]
lgene<-cbind(gene[ORFuni%in%llORF],as.numeric(delta[ORFuni%in%llORF]),as.numeric(delta_gamma[ORFuni%in%llORF]))
lgene_not<-cbind(gene[ORFuni%in%llORF_not],as.numeric(delta[ORFuni%in%llORF_not]),as.numeric(delta_gamma[ORFuni%in%llORF_not]))


list2<-list2[order(abs(list2[,4]),decreasing=T),]
list2<-cbind(list2,1:nrow(list2))
list2<-list2[order(list2[,1]),]
ADD_position<-list2[,8]



ORDER<-cbind(gene[vecorder],as.numeric(delta[vecorder]),as.numeric(delta_gamma[vecorder]),ADD_position[vecorder])
write.table(file="IHM_interactions.txt",ORDER)

ORDER[,4][order(abs(as.numeric(ORDER[,3])),decreasing=T)]

l<-list[!(list%in%lORF)]
l<-gene[ORFuni%in%l]
write.table(l,"IHM_not_interactions.txt")

write.table(cbind(gene[order],ORFuni[order],delta[order],delta_gamma[order],ADD_position[order]),"IHM_all.txt")

save.image(file="IHM_4oct_8k_8k_8k.RData")
#Percent interactors
#
length(llORF)/(length(unique(c(as.character(list),lORF))))





