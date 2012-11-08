

#b=read.delim("~/QFADatasets/SHM/YKU70_Raw/data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("CDC13-1_Raw_trim")#Query b 
TreatB=27
b<-a#query
#a=read.delim("~/QFADatasets/SHM/URA3_Raw/data.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data("URA3_Raw_trim")#Control a 
a$ORF<-as.character(a$ORF)#corrections
a$ORF[a$ORF=="YML009c"]=rep("YML009C",length(a$ORF[a$ORF=="YML009c"]))#corrections
a$ORF[a$ORF=="YMR169c"]=rep("YMR169C",length(a$ORF[a$ORF=="YMR169c"]))#corrections
a$ORF[a$ORF=="YMR175w"]=rep("YMR175W",length(a$ORF[a$ORF=="YMR175w"]))#corrections

TreatA=27

filename=paste("M_JHM_FULL","_",TreatA,"_",TreatB,sep="")

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
####
ORFuni<-unique(b$ORF)
#a<-funcIDORDER(a)
IDuni<-unique(a$ID)
gene<-unlist(lapply(ORFuni,funcGENE,data=b))#?

gene[gene=="0"]=ORFuni[gene=="0"]  #correction


N<-length(ORFuni);M=Ma=length(IDuni)
NoORF_a<-unlist(lapply(ORFuni_a,funcNoORF,data=a))#no of repeats each orf
NoTime_a<-c(0,unlist(lapply(IDuni,funcNoTime,data=a)))# 0+ no of time each repeat
NoSum_a<-c(0,unlist(lapply(1:N,funcNoSum,NoORF_vec=NoORF_a)))

#b<-funcIDORDER(b)
IDuni<-unique(b$ID)
###

#gene<-unlist(lapply(ORFuni,funcGENE,data=a))#?
N<-length(ORFuni);M=Mb=length(IDuni)#?
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
D<-main_JHM(3,1,2)


