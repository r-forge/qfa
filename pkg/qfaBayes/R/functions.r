### Filter by Screen Name, Temprature and Master Plate Number ###
funcREMOVE<-function(data,Screen,Treat,MPlate){
data=data[data$Screen.Name%in%Screen,]
data=data[data$Treatment%in%Treat,]
data=data[data$MasterPlate.Number%in%MPlate,]
data
}

### Orders dataset ###
funcIDORDER<-function(data){
data$ID<-paste(data$Barcode,data$MasterPlate.Number,formatC(data$Row,digits=2),formatC(data$Col,digits=2),sep="")
data<-data[order(paste(data$ORF,data$ID),data$Timeseries.order), ]
}

### Gives gene names ###
funcGENE<-function(x,data){
data$Gene[data$ORF%in%x][1]
}

### Gives number of repeats for each ORF ###
funcNoORF<-function(x,data){
length(unique((data$ID[data$ORF==x])))
}

### Gives number of time points for each repeat ###
funcNoTime<-function(x,data){
length((data$ID[data$ID==x]))
}

### Gives running total of number of number of repeats for each ORF ###
funcNoSum<-function(x,NoORF_vec){
sum(NoORF_vec[1:x])
}

### Adds NA values at the end of a repeat to give consistent row length for an array###
funcRowRep<-function(x,NoTime_vec,data_vec,dimr,dimc){
c(data_vec[sum(1,NoTime_vec[1:x]):sum(NoTime_vec[1:(x+1)])],rep(NA,dimc-length(data_vec[sum(1,NoTime_vec[1:x]):sum(NoTime_vec[1:(x+1)])])))
}

### Adds NA values at the end of an ORF to give consistent Col length for an array ###
funcColORF<-function(x,NoSum_vec,data_vec,dimr,dimc){
c(data_vec[(dimc*NoSum_vec[x]+1):(dimc*NoSum_vec[x+1])],rep(NA,dimr*dimc-length(data_vec[(dimc*NoSum_vec[x]+1):(dimc*NoSum_vec[x+1])])))
}

### Creates and transposes an Array ###
funcARRAYTRANS<-function(data_vec,dim){
vec<-array(c(data_vec),dim=dim)
vec<-aperm(vec, c(2,1,3))
vec
}

### Sorts data into array with correct dimensions ###
funcXY<-function(data,M,N,NoTime_vec,NoSum_vec,dimr,dimc){
XY<-unlist(lapply(1:M,funcRowRep,NoTime_vec=NoTime_vec,data_vec=data,dimr,dimc))
XY<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
dim<-c(dimc,dimr,N)
XY<-funcARRAYTRANS(XY,dim)
XY
}

### Creates and transposes an Array ###
funcARRAYTRANS_J<-function(data_vec,dim){
vec<-array(c(data_vec),dim=dim)
vec<-aperm(vec, c(2,1,3,4))
vec
}

### Sorts data into array with correct dimensions (Joint Model Specific) ###
funcXY_J<-function(data,data_b,Ma,Mb,N,NoTime_vec,NoSum_vec,NoTime_vec_b,NoSum_vec_b,dimr,dimc){
XY<-unlist(lapply(1:Ma,funcRowRep,NoTime_vec=NoTime_vec,data_vec=data,dimr,dimc))
XY<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec,data_vec=XY,dimr,dimc))
XY_b<-unlist(lapply(1:Mb,funcRowRep,NoTime_vec=NoTime_vec_b,data_vec=data_b,dimr,dimc))
XY_b<-unlist(lapply(1:N,funcColORF,NoSum_vec=NoSum_vec_b,data_vec=XY_b,dimr,dimc))
dim<-c(dimc,dimr,N,2)
XY<-funcARRAYTRANS_J(c(XY,XY_b),dim)
XY
}


### Creates and transposes an Array (Joint Model Specific) ###
funcSCALING<-function(data,vec){
lim<-max(data$Tile.Dimensions.Y)*max(data$Tile.Dimensions.X)*255
vec<-vec/lim
vec
}

### Saves hierarchical model to disk space for rjags to load ###
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


### Saves joint model to disk space for rJags to load ###
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

### Saves interaction model to disk space for rJags to load ###
funcMODELInteraction<-function(){
write("
model {
	for (i in 1:N){
		for (j in 1:2){
			for (k in 1:NoORF[i,j]){
				y[k,j,i]~ dnorm(alpha[j]*(mui[i]+delt[i,j]*gam[i,j]),nuj[j]*taui[i])
			}
		}
		mui[i]~dnorm(mu,mu_b)
		gam[i,1]<-0
		gam[i,2]~dnorm(0,gam_b)
		delt[i,1]<-0
		delt[i,2]~dbern(p)
		taui[i]~dgamma((tau^2)/(tau_b^2),(tau)/(tau_b^2))
	}
	mu~dnorm(mu_a,mu_b)
	alpha[1]<-1
	alpha[2]~dgamma((alpha_a^2)/(alpha_b^2),(alpha_a)/(alpha_b^2))
	tau~dgamma((tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
	nuj[1]~dgamma((nu^2)/(tau_b^2),(nu)/(tau_b^2))
	nuj[2]~dgamma((nu^2)/(tau_b^2),(nu)/(tau_b^2))
	nu~dgamma((tau_a^2)/(tau_b^2),(tau_a)/(tau_b^2))
}
","model1.bug")
}

### Load priors ###
funcPRIORS<-function(CustomModel){
if (!(CustomModel==FALSE)){source(paste(CustomModel,"Priors",sep="."))} else {data(PriorsH)}
list(
K_s=Priors$K_s,
r_s=Priors$r_s,
PO_s=Priors$PO_s,
beta=Priors$beta,
tau_s=Priors$tau_s,
delta=Priors$delta,
alpha=Priors$alpha,
gamma=Priors$gamma,
alpha_i=Priors$alpha_i,
gamma_i=Priors$gamma_i,
alpha_i_tau=Priors$alpha_i_tau,
gamma_i_tau=Priors$gamma_i_tau,
alpha_ij=Priors$alpha_ij,
gamma_ij=Priors$gamma_ij,
alpha_ij_tau=Priors$alpha_ij_tau,
gamma_ij_tau=Priors$gamma_ij_tau,
delta_mu=Priors$delta_mu,
delta_sd=Priors$delta_sd
)
}


### Load priors (Joint Model Specific) ###
funcPRIORS_J<-function(CustomModel){

if (!(CustomModel==FALSE)){source(paste(CustomModel,"Priors",sep="."))} else {data(PriorsJ)}
list(
K_s=Priors$K_s,
r_s=Priors$r_s,
PO_s=Priors$PO_s,
beta=Priors$beta,
tau_s=Priors$tau_s,
delta=Priors$delta,
alpha=Priors$alpha,
gamma=Priors$gamma,
alpha_i=Priors$alpha_i,
gamma_i=Priors$gamma_i,
alpha_ij=Priors$alpha_ij,
gamma_ij=Priors$gamma_ij,
alpha_ij_sd=Priors$alpha_ij_sd,
gamma_ij_sd=Priors$gamma_ij_sd,
p=Priors$p,
alpha_a=Priors$alpha_a,
alpha_b=Priors$alpha_b,
gam_b=Priors$gam_b,
omega_b=Priors$omega_b
)
}

funcJagsTime<-function(iter,upd,jags){
#TimeC<-(iter+upd)*system.time(update(jags,900))[2]
#print(paste("Time till completion",TimeC/(60*60*900),"(hours)",TimeC/(60*900),"(minutes)"))
}

### Fit, update and sample from the rjags model ###
funcFITandUPDATE<-function(QFA.I,QFA.D,QFA.P,inits){
jags <- jags.model('model1.bug',
                   data = list('x' = QFA.D$x,
                               'y' = QFA.D$y,
                               'N' = QFA.I$N,
'NoTime' = QFA.I$NoTime,
'NoORF' = QFA.I$NoORF,
'NoSum' = QFA.I$NoSum, 
'K_s' = QFA.P$K_s,
'r_s' = QFA.P$r_s,
'PO_s' = QFA.P$PO_s,
'beta' = QFA.P$beta,
'tau_s' = QFA.P$tau_s,
'delta'=QFA.P$delta,
'alpha' = QFA.P$alpha,
'alpha_i' = QFA.P$alpha_i,
'alpha_i_tau'=QFA.P$alpha_i_tau,
'alpha_ij' = QFA.P$alpha_ij,
'alpha_ij_tau'=QFA.P$alpha_ij_tau,
'gamma' = QFA.P$gamma,
'gamma_i' = QFA.P$gamma_i,
'gamma_i_tau'=QFA.P$gamma_i_tau,
'gamma_ij' = QFA.P$gamma_ij,
'gamma_ij_tau'=QFA.P$gamma_ij_tau,
'delta_mu'=QFA.P$delta_mu,
'delta_sd'=QFA.P$delta_sd
),
                   n.chains = 1,
                   n.adapt = 100,
			 inits=inits
)
funcJagsTime(iter,upd,jags)
update(jags, upd)
samp<-coda.samples(jags,
c('K_ij',
'r_ij',
'K_i',
'r_i',
'K',
'PO',
'r',
'tau',
'tau_i',
'K_ij_tau',
'r_ij_tau',
'K_i_tau',
'r_i_tau',
'delta_tau'),
             iter,thin=thin)
samp<-samp[[1]]
samp
}

### Fit, update and sample from the rjags model (Joint Model Specific) ###
funcFITandUPDATE_J<-function(QFA.I,QFA.D,QFA.P){
jags <- jags.model('model1.bug',
                   data = list('x' = QFA.D$x,
                               'y' = QFA.D$y,'SHIFT'=QFA.I$SHIFT,'p'=QFA.P$p,'alpha_a'=QFA.P$alpha_a,'alpha_b'=QFA.P$alpha_b,'gam_b'=QFA.P$gam_b,'omega_b'=QFA.P$omega_b,
                               'N' = QFA.I$N,'alpha_ij_sd'=QFA.P$alpha_ij_sd,'gamma_ij_sd'=QFA.P$gamma_ij_sd,'NoTime' = QFA.I$NoTime,'NoORF' = QFA.I$NoORF,'NoSum' = QFA.I$NoSum, 'K_s' = QFA.P$K_s,'PO_s' = QFA.P$PO_s,'r_s' = QFA.P$r_s,'tau_s' = QFA.P$tau_s,'delta'=QFA.P$delta,'alpha' = QFA.P$alpha,'beta' = QFA.P$beta,'gamma' = QFA.P$gamma,'alpha_i' = QFA.P$alpha_i,'gamma_i' = QFA.P$gamma_i,'alpha_ij' = QFA.P$alpha_ij,'gamma_ij' = QFA.P$gamma_ij),
                   n.chains = 1,
                   n.adapt = 100)
funcJagsTime(iter,upd,jags)
update(jags, upd)
samp<-coda.samples(jags,c('K','K_i','K_ij','PO','alph','bet','delt','gam','k_tau','r_tau','nu','nuc','omega','r','r_i','r_ij','tau','tau_m'),iter,thin=thin)
samp<-samp[[1]]
samp
}

### Outputs posterior sample in a named list ###
funcPosterior<-function(samp,N,M,iter,thin,upd){
if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-samp}
list(
vecsamp=vecsamp,
namesamp=names(vecsamp),
K=vecsamp[1],
K_i=vecsamp[2:(N+1)],
K_i_tau=vecsamp[(N+2)],
K_ij=vecsamp[(N+3):(M+N+2)],
K_ij_tau=vecsamp[(M+N+3):(M+2*N+2)],
PO=vecsamp[(M+2*N+3)],
delta_tau=vecsamp[(M+2*N+4)],
r=vecsamp[(M+2*N+5)],
r_i=vecsamp[(M+2*N+6):(M+3*N+5)],
r_i_tau=vecsamp[(M+3*N+6)],
r_ij=vecsamp[(M+3*N+7):(2*M+3*N+6)],
r_ij_tau=vecsamp[(2*M+3*N+7):(2*M+4*N+6)],
tau=vecsamp[(2*M+4*N+7)],
taui=vecsamp[(2*M+4*N+8):(2*M+5*N+7)],

samp=samp,
iter=iter,
thin=thin,
burnandupd=(1000+upd)
)
}

### Outputs posterior sample in a named list (Joint Model Specific) ###
funcPosterior_J<-function(samp,N,M,iter,thin,upd){
if(nrow(samp)>1) {vecsamp<-colMeans(samp)} else {vecsamp<-samp}
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
delta=vecsamp[(2*M+2*N+7):(2*M+3*N+6)],
samp=samp,
iter=iter,
thin=thin,
burnandupd=(1000+upd))
}

###  Gives experiment variables from ROD output###
qfa.variables<-function(data){
Screen<-as.character(unique(data$Screen.Name))
Treat<-as.character(unique(data$Treatment))
MPlate<-unique(data$MasterPlate.Number)
list(Screen=Screen,Treat=Treat,MPlate=MPlate)
}