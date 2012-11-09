#a=read.delim("data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
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
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

a$ID<-paste(a$Barcode,a$MasterPlate.Number,Row,Col,sep="")

ORFuni=unique(a$ORF)

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]

IDuni<-unique(a$ID)
ORFuni=unique(a$ORF)

gene<-unlist(lapply(ORFuni,funcGENE,data=a))
if(sum(gene=="0")>0){#Data Correction
gene[gene=="0"]=ORFuni[gene=="0"]
}
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
if (nchar(Row[i])<2){Row[i]=paste(0,Row[i],sep="")}
if (nchar(Col[i])<2){Col[i]=paste(0,Col[i],sep="")}
}

a$ID<-paste(a$Barcode,a$MasterPlate.Number,Row,Col,sep="")

ORFuni=unique(a$ORF)

a<-a[order(a$ORF,a$ID,a$Expt.Time), ]

IDuni<-unique(a$ID)
ORFuni=unique(a$ORF)

gene<-unlist(lapply(ORFuni,funcGENE,data=a))
if(sum(gene=="0")>0){#Data Correction
gene[gene=="0"]=ORFuni[gene=="0"]
}
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
tmp <- .C("main", as.integer(burn),as.integer(iters),as.integer(thin),as.integer(CAPL),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
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
CAPL=1#maximum no. of ORF's
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
tmp <- .C("main", as.integer(burn),as.integer(iters),as.integer(thin),as.integer(CAPL),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
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
CAPL=1#maximum no. of ORF's
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
QFAIA=QFA.IA,QFAy=QFA.yA,QFANoORFA=QFA.NoORFA,
QFAIB=QFA.IB,QFAy=QFA.yB,QFANoORFB=QFA.NoORFB,
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







