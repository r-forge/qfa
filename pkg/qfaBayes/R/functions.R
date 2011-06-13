
funcREMOVE<-function(data,uni,treat,MPlate){
data=data[data$Screen.Name%in%uni,]
data=data[data$Treatment%in%treat,]
data=data[data$MasterPlate.Number%in%MPlate,]
data
}
funcIDORDER<-function(data){
data$ID<-paste(data$Barcode,data$MasterPlate.Number,formatC(data$Row,digits=2),formatC(data$Col,digits=2),sep="")
data<-data[order(paste(data$ORF,data$ID),data$Timeseries.order), ]
}
fun1<-function(x,data){
data$Gene[data$ORF%in%x][1]
}
fun2<-function(x,data){
length(unique((data$ID[data$ORF==x])))
}
fun3<-function(x,data){
length((data$ID[data$ID==x]))
}
fun4<-function(x,NoORF_vec){
sum(NoORF_vec[1:x])
}
fun5<-function(x,NoTime_vec,data_vec,dimr,dimc){
c(data_vec[sum(1,NoTime_vec[1:x]):sum(NoTime_vec[1:(x+1)])],rep(NA,dimc-length(data_vec[sum(1,NoTime_vec[1:x]):sum(NoTime_vec[1:(x+1)])])))
}
fun6<-function(x,NoSum_vec,data_vec,dimr,dimc){
c(data_vec[(dimc*NoSum_vec[x]+1):(dimc*NoSum_vec[x+1])],rep(NA,dimr*dimc-length(data_vec[(dimc*NoSum_vec[x]+1):(dimc*NoSum_vec[x+1])])))
}
funcARRAYTRANS<-function(data_vec,dim){
vec<-array(c(data_vec),dim=dim)
vec<-aperm(vec, c(2,1,3))
vec
}

funcXY<-function(data,M,N,NoTime_vec,NoSum_vec,dimr,dimc){
XY<-unlist(lapply(1:M,fun5,NoTime_vec=NoTime_vec,data_vec=data,dimr,dimc))
XY<-unlist(lapply(1:N,fun6,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
dim<-c(dimc,dimr,N)
XY<-funcARRAYTRANS(XY,dim)
XY
}

funcARRAYTRANS_J<-function(data_vec,dim){
vec<-array(c(data_vec),dim=dim)
vec<-aperm(vec, c(2,1,3,4))
vec
}

funcXY_J<-function(data,data_b,M,N,NoTime_vec,NoSum_vec,NoTime_vec_b,NoSum_vec_b,dimr,dimc){
XY<-unlist(lapply(1:M,fun5,NoTime_vec=NoTime_vec,data_vec=data,dimr,dimc))
XY<-unlist(lapply(1:N,fun6,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
XY_b<-unlist(lapply(1:M,fun5,NoTime_vec=NoTime_vec_b,data_vec=data_b,dimr,dimc))
XY_b<-unlist(lapply(1:N,fun6,NoSum_vec=NoSum_vec_b,data_vec=XY_b,dimr,dimc))
dim<-c(dimc,dimr,N,2)
XY<-funcARRAYTRANS_J(c(XY,XY_b),dim)
XY
}

funcSCALING<-function(data,vec){
lim<-max(data$Tile.Dimensions.Y)*max(data$Tile.Dimensions.X)*255
vec<-vec/lim
vec
}


funcMODELHierarchical<-function(){
write("
model {
	for (i in 1:N){
		for (j in 1:NoORF[i]){
			for (l in 1:NoTime[(NoSum[i]+j)]){
				y[j,l,i] ~ dnorm(y.hat[j,l,i], tau[i])
                		y.hat[j,l,i] <- (K_ij[(NoSum[i]+j)]*PO*exp(r_ij[(NoSum[i]+j)]*x[j,l,i]))/(K_ij[(NoSum[i]+j)]+PO*(exp(r_ij[(NoSum[i]+j)]*x[j,l,i])-1))
        		}
        		K_ij[(NoSum[i]+j)] ~ dgamma((K_i[i]^2)*k_tau[i],K_i[i]*k_tau[i])
        		r_ij[(NoSum[i]+j)] ~ dgamma((r_i[i]^2)*r_tau[i],r_i[i]*r_tau[i])
		}
		K_i[i] ~ dgamma((K^2)/alpha_i^2,K/alpha_i^2)
		r_i[i] ~ dgamma((r^2)/gamma_i^2,r/gamma_i^2)
		tau[i] ~ dgamma((tau_m^2)/delta^2,tau_m/delta^2)
		k_tau[i] ~ dgamma((alpha_ij^2)/alpha_ij_sd^2,alpha_ij/alpha_ij_sd^2)
                r_tau[i] ~ dgamma((gamma_ij^2)/gamma_ij_sd^2,gamma_ij/gamma_ij_sd^2)
	}
	PO ~ dunif(PO_s,beta)   
	K ~ dgamma((K_s^2)/alpha^2,K_s/alpha^2)
	r ~ dgamma((r_s^2)/gamma^2,r_s/gamma^2)
	tau_m ~ dgamma((tau_s^2)/delta^2,tau_s/delta^2)}
","model1.bug")
}

funcMODELJoint<-function(){
write("
model {
for (i in 1:N){
for (c in 1:2){
	for (j in 1:NoORF[i,c]){
		for (l in 1:NoTime[(NoSum[i,c]+j),c]){
			y[j,l,i,c] ~ dnorm(y.hat[j,l,i,c], tau[i]*nuc[c])
			y.hat[j,l,i,c] <- (K_ij[(SHIFT[c]+NoSum[i,c]+j)]*PO*exp(r_ij[(SHIFT[c]+NoSum[i,c]+j)]*x[j,l,i,c]))/(K_ij[(SHIFT[c]+NoSum[i,c]+j)]+PO*(exp(r_ij[(SHIFT[c]+NoSum[i,c]+j)]*x[j,l,i,c])-1))
			}
			K_ij[(SHIFT[c]+NoSum[i,c]+j)] ~ dgamma(((alph[c]*(K_i[i]+delt[i,c]*gam[i,c]))^2)*k_tau[i],(alph[c]*(K_i[i]+delt[i,c]*gam[i,c]))*k_tau[i])
			r_ij[(SHIFT[c]+NoSum[i,c]+j)] ~ dgamma(((bet[c]*(r_i[i]+delt[i,c]*omega[i,c]))^2)*r_tau[i],(bet[c]*(r_i[i]+delt[i,c]*omega[i,c]))*r_tau[i])
			}
		}
		gam[i,1]<-0
		gam[i,2]~dnorm(0,gam_b)
		omega[i,1]<-0
		omega[i,2]~dnorm(0,omega_b)
		delt[i,1]<-0
		delt[i,2]~dbern(p)
		K_i[i] ~ dgamma((K^2)/alpha_i^2,K/alpha_i^2)
		r_i[i] ~ dgamma((r^2)/gamma_i^2,r/gamma_i^2)
		tau[i] ~ dgamma((tau_m^2)/delta^2,tau_m/delta^2)
		k_tau[i]~dgamma((alpha_ij^2)/alpha_ij_sd^2,alpha_ij/alpha_ij_sd^2)
 		r_tau[i]~dgamma((gamma_ij^2)/gamma_ij_sd^2,gamma_ij/gamma_ij_sd^2)
	}
	alph[1]<-1
	alph[2]~dgamma((alpha_a^2)/(alpha_b^2),(alpha_a)/(alpha_b^2))
	bet[1]<-1
	bet[2]~dgamma((alpha_a^2)/(alpha_b^2),(alpha_a)/(alpha_b^2))
	nuc[1]~dgamma((nu^2)/(delta^2),(nu)/(delta^2))
	nuc[2]~dgamma((nu^2)/(delta^2),(nu)/(delta^2))
	nu~dgamma((tau_s^2)/(delta^2),(tau_s)/(delta^2))
	PO ~ dunif(PO_s,beta)	
      K ~ dgamma((K_s^2)/alpha^2,K_s/alpha^2)
	r ~ dgamma((r_s^2)/gamma^2,r_s/gamma^2)
	tau_m~dgamma((tau_s^2)/delta^2,tau_s/delta^2)
}","model1.bug")
}

funcPRIORS<-function(CustomModel){
if (!(CustomModel==FALSE)){Priors<-read.delim(paste(CustomModel,"Priors",sep="."),header=F)} else {data(PriorsH)}
list("K_s"=Priors[1,],
"r_s"=Priors[2,],
"PO_s"=Priors[3,],
"beta"=Priors[4,],
"tau_s"=Priors[5,],
"delta"=Priors[6,],
"alpha"=Priors[7,],
"gamma"=Priors[8,],
"alpha_i"=Priors[9,],
"gamma_i"=Priors[10,],
"alpha_ij"=Priors[11,],
"gamma_ij"=Priors[12,],
"alpha_ij_sd"=Priors[13,],
"gamma_ij_sd"=Priors[14,]
)
}

funcPRIORS_J<-function(CustomModel){

if (!(CustomModel==FALSE)){Priors<-read.delim(paste(CustomModel,"Priors",sep="."),header=F)} else {Priors<-read.delim("Priors",header=F)}
list(
K_s=Priors[1,],
r_s=Priors[2,],
PO_s=Priors[3,],
beta=Priors[4,],
tau_s=Priors[5,],
delta=Priors[6,],
alpha=Priors[7,],
gamma=Priors[8,],
alpha_i=Priors[9,],
gamma_i=Priors[10,],
alpha_ij=Priors[11,],
gamma_ij=Priors[12,],
alpha_ij_sd=Priors[13,],
gamma_ij_sd=Priors[14,],
p=Priors[15,],
alpha_a=Priors[16,],
alpha_b=Priors[17,],
gam_b=Priors[18,],
omega_b=Priors[19,]
)
}

funcFITandUPDATE<-function(QFA.I,QFA.D,QFA.P){
jags <- jags.model('model1.bug',
                   data = list('x' = QFA.D$x,
                               'y' = QFA.D$y,
                               'N' = QFA.I$N,'alpha_ij_sd'=QFA.P$alpha_ij_sd,'gamma_ij_sd'=QFA.P$gamma_ij_sd,'NoTime' = QFA.I$NoTime,'NoORF' = QFA.I$NoORF,'NoSum' = QFA.I$NoSum, 'K_s' = QFA.P$K_s,'PO_s' = QFA.P$PO_s,'r_s' = QFA.P$r_s,'tau_s' = QFA.P$tau_s,'delta'=QFA.P$delta,'alpha' = QFA.P$alpha,'beta' = QFA.P$beta,'gamma' = QFA.P$gamma,'alpha_i' = QFA.P$alpha_i,'gamma_i' = QFA.P$gamma_i,'alpha_ij' = QFA.P$alpha_ij,'gamma_ij' = QFA.P$gamma_ij),
                   n.chains = 1,
                   n.adapt = 100)
TimeC<-(iter+upd)*system.time(update(jags,900))[1]
print(paste("Time till completion",TimeC/(60*60*900),"(hours)",TimeC/(60*900),"(minutes)"))
update(jags, upd)
samp<-coda.samples(jags,
          c('K_ij','r_ij','K_i','r_i','K','PO','r','tau','tau_m','k_tau','r_tau'),
             iter,thin=thin)
samp<-samp[[1]]
samp
}

funcFITandUPDATE_J<-function(QFA.I,QFA.D,QFA.P){
SHIFT<-c(0,max(QFA.I$NoSum[,1]))#####!!!!!!!!!!!!!!!!!!!!!!
jags <- jags.model('model1.bug',
                   data = list('x' = QFA.D$x,
                               'y' = QFA.D$y,'SHIFT'=SHIFT,'p'=QFA.P$p,'alpha_a'=QFA.P$alpha_a,'alpha_b'=QFA.P$alpha_b,'gam_b'=QFA.P$gam_b,'omega_b'=QFA.P$omega_b,
                               'N' = QFA.I$N,'alpha_ij_sd'=QFA.P$alpha_ij_sd,'gamma_ij_sd'=QFA.P$gamma_ij_sd,'NoTime' = QFA.I$NoTime,'NoORF' = QFA.I$NoORF,'NoSum' = QFA.I$NoSum, 'K_s' = QFA.P$K_s,'PO_s' = QFA.P$PO_s,'r_s' = QFA.P$r_s,'tau_s' = QFA.P$tau_s,'delta'=QFA.P$delta,'alpha' = QFA.P$alpha,'beta' = QFA.P$beta,'gamma' = QFA.P$gamma,'alpha_i' = QFA.P$alpha_i,'gamma_i' = QFA.P$gamma_i,'alpha_ij' = QFA.P$alpha_ij,'gamma_ij' = QFA.P$gamma_ij),
                   n.chains = 1,
                   n.adapt = 100)
TimeC<-(iter+upd)*system.time(update(jags,900))[1]
print(paste("Time till completion",TimeC/(60*60*900),"(hours)",TimeC/(60*900),"(minutes)"))
update(jags, upd)
samp<-coda.samples(jags,c('K','K_i','K_ij','PO','alph','bet','delt','gam','k_tau','r_tau','nu','nuc','omega','r','r_i','r_ij','tau','tau_m'),iter,thin=thin)
samp<-samp[[1]]
samp
}

funcSAVE<-function(data,samp,N,M,iter,thin,upd){
vecsamp=colMeans(samp)
list(
vecsamp=vecsamp,
namesamp=names(vecsamp),
K=vecsamp[1],
K_i=vecsamp[2:(N+1)],
K_ij=vecsamp[(N+2):(M+N+1)],
PO=vecsamp[(M+N+2)],
k_tau=vecsamp[(M+N+3):(M+2*N+2)],
r=vecsamp[(M+2*N+3)],
r_i=vecsamp[(M+2*N+4):(M+3*N+3)],
r_ij=vecsamp[(M+3*N+4):(2*M+3*N+3)],
r_tau=vecsamp[(2*M+3*N+4):(2*M+4*N+3)],
taui=vecsamp[(2*M+4*N+4):(2*M+5*N+3)],
tau=vecsamp[(2*M+5*N+4)],
samp=samp,
iter=iter,
thin=thin,
burnandupd=(1000+upd)
)
}

funcSAVE_J<-function(data,data_b,samp,N,M,iter,thin,upd){
vecsamp=colMeans(samp)
list(
vecsamp=vecsamp,
namesamp=names(vecsamp),
K=vecsamp[1],
K_i=vecsamp[2:(N+1)],
K_ij=vecsamp[(N+2):(2*M+N+1)],
PO=vecsamp[(2*M+N+2)],
k_tau=vecsamp[(2*M+5*N+7):(2*M+6*N+6)],
r=vecsamp[(2*M+8*N+10)],
r_i=vecsamp[(2*M+8*N+11):(2*M+9*N+10)],
r_ij=vecsamp[(2*M+9*N+11):(4*M+9*N+10)],
r_tau=vecsamp[(4*M+9*N+11):(4*M+10*N+10)],
taui=vecsamp[(4*M+10*N+11):(4*M+11*N+10)],
tau=vecsamp[(4*M+12*N+11)],
gam=vecsamp[(2*M+4*N+7):(2*M+5*N+6)],
omega=vecsamp[(2*M+7*N+10):(2*M+8*N+9)],
nu=vecsamp[(2*M+6*N+7)],
nuc=vecsamp[(2*M+6*N+8):(2*M+6*N+9)],
gamdelt=colMeans(samp[,(2*M+4*N+7):(2*M+5*N+6)]*samp[,(2*M+2*N+7):(2*M+3*N+6)]),
omegadelt=colMeans(samp[,(2*M+7*N+10):(2*M+8*N+9)]*samp[,(2*M+2*N+7):(2*M+3*N+6)]),
delta=vecsamp[,(2*M+2*N+7):(2*M+3*N+6)],
samp=samp,
iter=iter,
thin=thin,
burnandupd=(1000+upd))
}

#####################################################################
print("Preprocessing")
#####################################################################

qfa.variables<-function(data){
uni<-as.character(unique(data$Screen.Name))
treat<-as.character(unique(data$Treatment))
MPlate<-unique(data$MasterPlate.Number)
list(uni,treat,MPlate)
}
