
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

