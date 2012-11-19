data("URA3_Raw_extratrim_15")
a<-a_15
data("CDC13-1_Raw_extratrim_15")
b<-b_15
qfa.variables(a)
qfa.variables(b)

data("priors_SHM")
PRIORS=as.double((priors_SHM)[[1]])[1:18]

Screen_a<-unique(a$Screen.Name)
SHM_a<-SHM_postpro(a=a,Treat=27,Screen=Screen_a,MPlate=15)
SHM_output_a<-SHM_main(burn=1,iters=1,thin=1,CAPL=50,QFA.I=SHM_a$QFA.I,QFA.y=SHM_a$QFA.y,QFA.x=SHM_a$QFA.x,QFA.NoORF=SHM_a$QFA.NoORF,QFA.NoTIME=SHM_a$QFA.NoTIME,PRIORS=PRIORS)

Screen_b<-unique(b$Screen.Name)
SHM_b<-SHM_postpro(a=b,Treat=27,Screen=Screen_b,MPlate=15)
SHM_output_b<-SHM_main(burn=1,iters=1,thin=1,CAPL=50,QFA.I=SHM_b$QFA.I,QFA.y=SHM_b$QFA.y,QFA.x=SHM_b$QFA.x,QFA.NoORF=SHM_b$QFA.NoORF,QFA.NoTIME=SHM_b$QFA.NoTIME,PRIORS=PRIORS)

SHM_a$QFA.yA=colMeans(SHM_output_a)
SHM_b$QFA.yB=colMeans(SHM_output_b)

IHM_output=IHM_main(burn=1,iters=1,thin=1,QFA.IA=SHM_a$QFA.I,QFA.yA=SHM_a$QFA.yA,QFA.NoORFA=SHM_a$QFA.NoORF,QFA.IB=SHM_b$QFA.I,QFA.yB=SHM_b$QFA.yB,QFA.NoORFB=SHM_b$QFA.NoORF,PRIORS)

##########################################################
plotYN=0
while(plotYN < 1 ){
  n<-readline("do you wish to plot? Y or N: ")
if(n=="Y"){plotYN=1}
if(n=="N"){stop()}
}
N=SHM_a$QFA.I[1]
gene=SHM_a$gene
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

file="DEMO"#

pdf(paste("IHM_plot_",file,".pdf",sep=""),useDingbats=F)

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
dev.off()


stop() 
