data("URA3_Raw_extratrim_15")#Control a 
a<-a_15
data("CDC13-1_Raw_extratrim_15")#Query b 
b<-b_15
qfa.variables(a)
qfa.variables(b)

data("priors_JHM")
PRIORS=as.double((priors_JHM)[[1]])
PRIORS[19]=0##

TreatA=27
Screen_a<-as.character(unique(a$Screen.Name))
MPlate_a<-as.character(unique(a$MasterPlate.Number))

TreatB=27
Screen_b<-as.character(unique(b$Screen.Name))
MPlate_b<-as.character(unique(b$MasterPlate.Number))

JHM<-JHM_postpro(a,TreatA=TreatA,Screen_a=Screen_a,MPlate_a,b,TreatB=TreatB,Screen_b,MPlate_b)
JHM_output<-JHM_main(burn=1,iters=1,thin=1,QFA.IA=JHM$QFA.IA,QFA.yA=JHM$QFA.yA,QFA.xA=JHM$QFA.xA,QFA.NoORFA=JHM$QFA.NoORFA,QFA.NoTIMEA=JHM$QFA.NoTIMEA,QFA.IB=JHM$QFA.IB,QFA.yB=JHM$QFA.yB,QFA.xB=JHM$QFA.xB,QFA.NoORFB=JHM$QFA.NoORFB,QFA.NoTIMEB=JHM$QFA.NoTIMEB,PRIORS)

#####################################################
plotYN=0
while(plotYN < 1 ){
  n<-readline("do you wish to plot? Y or N: ")
if(n=="Y"){plotYN=1}
if(n=="N"){stop()}
}

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

file="DEMO"
pdf(paste("JHM_plot_",file,".pdf",sep=""),useDingbats=F)

plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",ylab="",pch=19,col=8,cex=0.5)
lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)

vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
vecMDPa<-log(K_ij/PO)/log(2) #MDP
vecMDPa*vecMDRa

i=1:N
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
i=1:N
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
i=1:N
points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]>0]#######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
i=vecorder[omegadelt[order][1:sig]<=0]  #######
points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
i=vecorder
text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
dev.off()

stop()
