### main SHM function (saved data) ###
main_SHM<-function(burn,iters,thin,CAPL,
LMNmaxdata="LMNmaxdata.txt",
ydata="ydata.txt",
xdata="xdata.txt",
NoORFdata="NoORFdata.txt",
NoTIMEdata="NoTIMEdata.txt",
priors="priors.txt"
) {
aa<-read.table(LMNmaxdata,header=T)
QFA.I=as.integer((aa)[[1]])
aa<-read.table(ydata,header=T)
QFA.y=as.double((aa)[[1]])
aa<-read.table(xdata,header=T)
QFA.x=as.double((aa)[[1]])
aa<-read.table(NoTIMEdata,header=T)
QFA.NoTIME=as.integer((aa)[[1]])
aa=read.table(priors,header=T)
PRIORS=as.double((aa)[[1]])[1:18]
aa<-read.table(NoORFdata,header=T)
QFA.NoORF=as.integer((aa)[[1]])
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
tmp <- .C("main",PACKAGE="qfaBayes", as.integer(burn),as.integer(iters),as.integer(thin),as.integer(CAPL),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAI=QFA.I,QFAy=QFA.y,QFAx=QFA.x,QFANoORF=QFA.NoORF,QFANoTIME=QFA.NoTIME,
PRIORS=PRIORS
)
mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

### main IHM function ###

main_IHM<- function(burn,iters,thin,
LMNmaxdataA1="LMNmaxdataA1.txt",
dataA2="dataA2.txt",
NoORFdataA1="NoORFdataA1.txt",
LMNmaxdataB1="LMNmaxdataB1.txt",
dataB2="dataB2.txt",
NoORFdataB1="NoORFdataB1.txt",
priors_IHM="priors_IHM.txt"
) {
aa<-read.table(LMNmaxdataA1,header=T)
QFA.IA=as.integer((aa)[[1]])
aa<-read.table(dataA2,header=T)
QFA.yA=as.double((aa)[[1]])

aa<-read.table(LMNmaxdataB1,header=T)
QFA.IB=as.integer((aa)[[1]])
aa<-read.table(dataB2,header=T)
QFA.yB=as.double((aa)[[1]])

aa=read.table(priors_IHM,header=T)
PRIORS=as.double((aa)[[1]])

aa<-read.table(NoORFdataA1,header=T)
QFA.NoORFA=as.integer((aa)[[1]])
bb<-read.table(NoORFdataB1,header=T)
QFA.NoORFB=as.integer((bb)[[1]])
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
1+
1+
1
tmp <- .C("main_IHM",PACKAGE="qfaBayes", as.integer(burn),as.integer(iters),as.integer(thin),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAIA=QFA.IA,QFAy=QFA.yA,QFANoORFA=QFA.NoORFA,
QFAIB=QFA.IB,QFAy=QFA.yB,QFANoORFB=QFA.NoORFB,
PRIORS=PRIORS
)
mat=matrix((tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

### main JHM function ###

main_JHM<- function(burn,iters,thin,
LMNmaxdataA1="LMNmaxdataA1.txt",
ydataA1="ydataA1.txt",
xdataA1="xdataA1.txt",
NoORFdataA1="NoORFdataA1.txt",
NoTIMEdataA1="NoTIMEdataA1.txt",
LMNmaxdataB1="LMNmaxdataB1.txt",
ydataB1="ydataB1.txt",
xdataB1="xdataB1.txt",
NoORFdataB1="NoORFdataB1.txt",
NoTIMEdataB1="NoTIMEdataB1.txt",
priors="priors.txt"
) {
aa<-read.table("LMNmaxdataA1.txt",header=T)
QFA.IA=as.integer((aa)[[1]])
aa<-read.table(ydataA1,header=T)
QFA.yA=as.double((aa)[[1]])
aa<-read.table(xdataA1,header=T)
QFA.xA=as.double((aa)[[1]])
aa<-read.table(NoTIMEdataA1,header=T)
QFA.NoTIMEA=as.integer((aa)[[1]])
aa<-read.table(LMNmaxdataB1,header=T)
QFA.IB=as.integer((aa)[[1]])
aa<-read.table(ydataB1,header=T)
QFA.yB=as.double((aa)[[1]])
aa<-read.table(xdataB1,header=T)
QFA.xB=as.double((aa)[[1]])
aa<-read.table(NoTIMEdataB1,header=T)
QFA.NoTIMEB=as.integer((aa)[[1]])
aa=read.table(priors,header=T)
PRIORS=as.double((aa)[[1]])
PRIORS[19]=0

aa<-read.table(NoORFdataA1,header=T)
QFA.NoORFA=as.integer((aa)[[1]])
bb<-read.table(NoORFdataB1,header=T)
QFA.NoORFB=as.integer((bb)[[1]])
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
1+
1+
2*2+
2*2
tmp <- .C("main_JHM",PACKAGE="qfaBayes", as.integer(burn),as.integer(iters),as.integer(thin),OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
QFAIA=QFA.IA,QFAy=QFA.yA,QFAxA=QFA.xA,QFANoORFA=QFA.NoORFA,QFANoTIMEA=QFA.NoTIMEA,
QFAIB=QFA.IB,QFAy=QFA.yB,QFAxB=QFA.xB,QFANoORFB=QFA.NoORFB,QFANoTIMEB=QFA.NoTIMEB,
PRIORS=PRIORS
)
mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
mat=data.frame(mat)
names(mat)=tmp$HEADER
mat
}

