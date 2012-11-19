
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
list(QFA.I=QFA.I,QFA.y=QFA.y,QFA.x=QFA.x,QFA.NoORF=QFA.NoORF,QFA.NoTIME=QFA.NoTIME)
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


