############################### Functions common to Bayesian and Likelihoood #######################################

################ Read and bind rod output into suitable form for likelihood & bayesian #################
rod.read<-function(path=".",files=c(),inoctimes="BarcodeTimes.txt",background="",treatments=c(),barcodes=c(),
master.plates=c(),screen.names=c()){
# Create environment with inoculation times; checks if path specified
if (length(strsplit(path,".")[[1]])>1){pathT<-1
	inocenv<-bctimes(paste(path,inoctimes,sep="/"))} else {
	inocenv<-bctimes(inoctimes); pathT<-0}
# If no files specified, use all .txt in working directory
if (length(files)==0){fs<-list.files(path=path,pattern=".txt")} else {fs<-files}
fs<-fs[(fs!='BarcodeTimes.txt')&(fs!='ORF2GENE.txt')&(fs!="model.txt")&(fs!="test.txt")]
print("List of data files to read:")
print(fs)
if (pathT==1){fs<-paste(path,fs,sep="/")}
bigd<-data.frame()
for (f in fs){
	print(paste("Reading",as.character(f)))
	# Find out which ROD version and print
	rodversion<-as.numeric(scan(f,nlines=1,what="character")[4])
	print(paste("ROD Export Version",rodversion))
	if (rodversion==5){
		# Read in tab deliminated rod output
		roddata<-read.delim(f,skip=1)
		# Only take specified treatments, or all if none specified
		# and do the same for other arguments
		if (length(treatments)>0){
			roddata<-roddata[(roddata$Treatments%in%treatments),]}
		# Only take specified barcodes, or all if none specified
		if (length(barcodes)>0){
			roddata<-roddata[(roddata$Barcode%in%barcodes),]}
		if (length(master.plates)>0){
			roddata<-roddata[(roddata$MasterPlate.Number%in%master.plates),]}
		if (length(screen.names)>0){
			roddata<-roddata[(roddata$Screen.Name%in%screen.names),]}
		# Add column for inoculation times
		roddata$Inoc.Time<-""; bcs<-unique(roddata$Barcode)
		for (b in bcs){
		roddata[roddata$Barcode==as.character(b),'Inoc.Time']<-get(as.character(b),envir=inocenv)}
		# Correct column names
		names(roddata)[3]<-"Row"; names(roddata)[5]<-"Col"; names(roddata)[14]<-"ORF"
		names(roddata)[4]<-"Growth"; names(roddata)[15]<-"Date.Time"
		# Bind 
		bigd<-rbind(bigd,roddata)} else {stop("Unkown ROD version")} #ifROD5
	} #fs
# Add column for genetic background
bigd$Background<-background
# Print checks so people are sure they're using the correct data #
print(paste("Number of Barcodes =",length(unique(bigd$Barcode))))
print("Treatments :")
print(unique(bigd$Treatments))
print("Media:")
print(unique(bigd$Medium))
print(paste("Number of ORFs =",length(unique(bigd$ORF))))
print("Screens:")
print(unique(bigd$Screen.Name))
platesize<-max(as.numeric(bigd$Row))*max(as.numeric(bigd$Col))
if (length(bigd$Date.Time)%%platesize!=0){
	warning("Number of photos not multiple of plate size")}
print("Inoculation DateTimes:")
print(unique(bigd$Inoc.Time))
return(bigd)
}

#### Define Phenotype for fitness ####
mdrmdp<-function(K,r,g){MDR<-sapply(r/((2*(K-g))/(K-2*g)),na2zero)
MDP<-sapply(log(K/g)/log(2),na2zero)
return(MDR*MDP)
}

####### Convert data datetime to time from start in days ########
tconv<-function(tstring,startt,fmt){
t<-as.POSIXlt(as.character(tstring),format=fmt)
return(as.numeric(difftime(t,startt,units="days")))
}

# Timing function (internal) #
trep<-function(thing,tobj){
print(paste("Finished",thing,"in",round(tobj[3],2),"seconds"))
print(paste(round(tobj[3]/60,2),"minutes"))
}

# Converts NAs to zeros
na2zero<-function(mdr){if (is.na(mdr)){mdr<-0}
return(mdr)
}

###### Create inoculation time environment #####
bctimes<-function(bctimef){
startframe=read.delim(bctimef,colClasses=c("character"))
starts<-new.env(hash=TRUE)
for (k in 1:length(startframe$Barcode)){
st=startframe$Start.Time[k]
assign(startframe$Barcode[k],st,envir=starts)}
return(starts)
} # bctimes

# Convert barcode to start time #
bc2st<-function(bc,inocenv){get(as.character(bc),envir=inocenv)}

### Create ORF2Gene dictionary ###
orf2gdict<-function(dictionary){
cache=new.env(hash = TRUE)
# Read in dictionary txt file 
orf2g=read.delim(dictionary,header=FALSE,
colClasses=c("character"),col.names=c("ORF","Gene"))
z=apply(orf2g,1,orfun,cache)
return(cache)
}

# Function called in orf2gdict #
orfun<-function(row,environ){
orf<-as.character(row[1])
gene<-as.character(row[2])
return(assign(orf,gene,envir=environ))
} #orfun

### Function that converts orf 2 gene ###
orf2g<-function(orf,dictenv){get(orf,envir=dictenv)}


################################################## Epistasis Function ###########################################################
qfa.epi<-function(double,control,qthresh,rjags=FALSE,orfdict="ORF2GENE.txt",
GISthresh=0.5,plot=TRUE,modcheck=TRUE,fitfunct=mdrmdp){
###### Get ORF mean fitnesses for control & double #######
print("Calculating mean fitness for each ORF")
## Bayesian ##
if (is.null(dim(double))){bayes=1
# Get initial colony sizes
dg<-double$g; cg<-control$g
### Get orfs in double posterior ###
orfs<-dimnames(double[["k1"]])[2][[1]]
## Make posterior organized by ORFs ##
double<-lapply(orfs,posmake,orfs,double)
names(double)<-orfs
### Get orfs in control posterior ###
corfs<-dimnames(control[["k1"]])[2][[1]]
## Make posterior organized by ORFs ##
control<-lapply(corfs,posmake,corfs,control)
names(control)<-corfs
# Calculate fitness posteriors with common orfs
orfs<-orfs[orfs%in%corfs]
dfits<-lapply(orfs,fitmake,double,mdrmdp,dg)
cfits<-lapply(orfs,fitmake,control,mdrmdp,cg)
dFms<-sapply(dfits,mean)
names(dFms)<-orfs
cFms<-sapply(cfits,mean)
names(cFms)<-orfs} else {bayes=0
## LIK ##
# Get orfs in question
orfs<-as.character(double$ORF)
orfs<-orfs[orfs%in%as.character(control$ORF)]
# Get lists with fitnesses for each repeat
cFstats<-lapply(orfs,orfstat,control,fitfunct)
names(cFstats)<-orfs
dFstats<-lapply(orfs,orfstat,double,fitfunct)
names(dFstats)<-orfs
# Get means for each ORF
cFms<-sapply(cFstats,mean)
names(cFms)<-orfs
dFms<-sapply(dFstats,mean)
names(dFms)<-orfs}
### Fit genetic independence model ###
if (rjags==FALSE){m<-lm.epi(dFms,cFms,modcheck)} else {
print("Initializing genetic independence model")
mod<-epimod(dFms,cFms)
print("Simulating posterior")
update(mod,2*10^3)
m<-coda.samples(mod,c("m"),n.iter=2*10^3,thin=2)
m<-m[[1]][,'m']
acf(m,lag.max=10^3)
mdens<-density(m,n=10^4)
mdens<-data.frame(x=mdens$x,y=mdens$y)
mdens<-mdens[order(-mdens$y),]
m<-mdens[1,'x']}
print(paste("Ratio of background mutant fitness to wildtype fitness =",round(m,4)))
###### Estimate probability of interaction #######
print("Calculating interaction probabilities")
if (bayes==0){
p<-sapply(orfs,pmake,m,cFstats,dFstats,cFms,dFms)
} else {
p<-sapply(orfs,bayesp,m,dfits,cfits)}
# Adjust for multiple comparisons
q<-p.adjust(p,"fdr")
# Get gene names if needed
if (bayes==1){# Create orf-gene dictionary
orfdict<-orf2gdict(orfdict)
genes<-sapply(orfs,orf2g,orfdict)} else {# lik
genes<-as.character(double$Gene); genes<-genes[genes%in%as.character(control$Gene)]
# If gene names missing, find them with orf2g
if (is.null(genes)){# Create orf-gene dictionary
orfdict<-orf2gdict(orfdict)
genes<-sapply(orfs,orf2g,orfdict)}}
# Get genetic interaction scores
gis<-dFms/mean(dFms)-cFms/mean(cFms)
# Put into data.frame
results<-data.frame(ORF=orfs,Gene=genes,P=p,Q=q,GIS=gis,Mean.Double=dFms,Mean.Control=cFms)
results$Type<-apply(results,1,typemake,m)
results<-results[order(-abs(results$GIS),results$Q,results$Type),]
# Plot results
if (plot==TRUE){qfa.epiplot(results,qthresh,m)}
final<-list(Results=results,
Enhancers=gethits(results,qthresh,type="E",GISthresh=GISthresh),
Suppressors=gethits(results,qthresh,type="S",GISthresh=GISthresh))
return(final)
}

############### Epistasis Functions ##################
# Makes epistasis plot for a given fdr level #
qfa.epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W"){
enhancers<-gethits(results,qthresh,type="E")
suppressors<-gethits(results,qthresh,type="S")
others<-results[(!results$ORF%in%enhancers$ORF)&(!results$ORF%in%suppressors$ORF),]
# Get plot parameters
ymax=1.1*max(suppressors$Mean.Double,enhancers$Mean.Double)
#ymin=0.9*min(suppressors$Mean.Double,enhancers$Mean.Double)
ymin=0
xmax=1.1*max(suppressors$Mean.Control,enhancers$Mean.Control)
#xmin=0.9*min(suppressors$Mean.Control,enhancers$Mean.Control)
xmin=0
plot(NULL,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
xlab="Control Fitness",ylab="Query Fitness",main="Epistasis Plot",
col=8,pch=19,cex=0.5)
# Add line for genetic independence
if (fitratio!=FALSE){abline(0,fitratio,lwd=2,col=8)} else {
abline(0,lm.epi(results$Mean.Double,results$Mean.Control,modcheck=FALSE),lwd=2,col=8)}
# Add 1:1 fitness line
abline(0,1,lwd=2,lty=4,col=8)
# Add reference ORF fitnesses lines
if (ref.orf!=FALSE){reforf<-results[results$ORF==ref.orf,]
abline(v=reforf$Mean.Control,col="lightblue",lwd=2)
abline(h=reforf$Mean.Double,col="lightblue",lwd=2)}
# Add points for non-suppressors & non-enhancers
points(others$Mean.Control,others$Mean.Double,col="grey",cex=0.5,pch=19)
# Add suppressors & enhancers
points(enhancers$Mean.Control,enhancers$Mean.Double,col='green',pch=19,cex=0.5)
text(enhancers$Mean.Control,(enhancers$Mean.Double+ymax/150),enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
points(suppressors$Mean.Control,suppressors$Mean.Double,col='red',pch=19,cex=0.5)
text(suppressors$Mean.Control,(suppressors$Mean.Double+ymax/150),suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
}

## Extract hits from epistasis results object ##
gethits<-function(results,qthresh,type="S",all.types=FALSE,GISthresh=0){
results<-results[results$Q<qthresh,]
if (all.types==FALSE){results<-results[results$Type==type,]}
results<-results[abs(results$GIS)>=GISthresh,]
return(results[order(-abs(results$GIS),results$Q,results$Type),])
}

# Function to initialize genetic ind. model
epimod<-function(doubles,controls){require(rjags)
write("
model {
for (i in 1:N){
	Fd[i]~dnorm(m*Fc[i],sigma^(-2))
	}
m~dweib(1.9,0.7)
sigma~dlnorm(log(4),1)
}","model.txt")
jagd<-list(Fd=doubles,Fc=controls,N=length(doubles))
return(jags.model(data=jagd,file="model.txt"))
}

# Function to calculate fitnesses for repeats
orfstat<-function(orf,fitframe,fitfunct){
orfd<-fitframe[fitframe$ORF==orf,]
return(fitfunct(orfd$K,orfd$r,orfd$g))
}

# Linear regression genetic ind. model fit
lm.epi<-function(doubles,controls,modcheck){
indmod<-lm(doubles~0+controls)
m<-as.numeric(indmod$coefficients)
# Check if linear regression OK
if (modcheck==TRUE){
pdf("ModelCheck.pdf")
resids<-indmod$residuals
hist(resids,xlab="Fitness Residuals",main="Histogram of Residuals")
qqnorm(resids/sd(resids),main="Normal QQ Plot")
abline(0,1,lwd=2,col="red")
dev.off()}
return(m)
}

# Estimates probability of interaction for max lik method #
pmake<-function(orf,m,cFs,dFs,cFms,dFms){
if (dFms[orf]>m*cFms[orf]){
	p<-wilcox.test(dFs[[orf]],m*cFs[[orf]],alternative="greater")$p.value} else {
	p<-wilcox.test(dFs[[orf]],m*cFs[[orf]],alternative="less")$p.value}
return(p)
}

# Calculate probability of interaction for Bayesian #
bayesp<-function(orf,m,dfits,cfits){
d<-dfits[[orf]]; c<-cfits[[orf]]; edf<-m*c
z<-(mean(d)-mean(edf))/(sd(d)+sd(edf))
if (z<0){p<-pnorm(z)} else {p<-1-pnorm(z)}
return(p)
}

# Function to make fitnesses
fitmake<-function(orf,poslist,fitfunct,g){
orfpos<-poslist[[orf]]
return(fitfunct(orfpos$k1,orfpos$r1,g))
}

# Get type of interaction
typemake<-function(row,m){md<-as.numeric(row['Mean.Double']); c<-as.numeric(row['Mean.Control'])
if (m<1){if (md>m*c){type<-"S"} else {type<-"E"}} else {
if (md>m*c){type<-"E"} else {type<-"S"}}
return(type)
}

####### Grabs posteriors for each variable for an ORF #######
posmake<-function(orf,orfs,varpos){orfn<-match(orf,orfs)
norfs<-length(orfs)
return(lapply(varpos,varposget,orfn,norfs))
}

#### Returns posterior if variable repeated for all ORFs ####
varposget<-function(var,orfn,norfs){
if (!is.null(dim(var))){z<-var[,orfn]} else {z<-NULL}
return(z)
}


####################################################### Bayesian Functions ###############################################################

################################### Function to put data into JAGS format ############################################
qfa.data<-function(d,split=c(),fmt="%Y-%m-%d_%H-%M-%S",orflist=c()){
# Restrict to variable of interest in data.frame
d<-data.frame(Barcode=d$Barcode,ORF=d$ORF,Row=d$Row,Col=d$Col,Growth=d$Growth,
Date.Time=d$Date.Time,Inoc.Time=d$Inoc.Time)
stime<-proc.time()
print("Selecting ORFs to put in model")
# Use specified ORFs or all if none specified
if (length(orflist)==0){orfs<-unique(d$ORF)} else {orfs<-orflist}
# Split ORFs and only use one part if specified
if (length(split)>0){
nparts<-split[1]; part<-split[2]
splitsize<-floor(length(orfs)/nparts)
barriers<-c(0,1:(nparts-1)*splitsize,length(orfs))
orfparts<-list()
for (j in 1:nparts){
orfparts[[j]]<-orfs[(1+barriers[j]):barriers[j+1]]}
orfs<-orfparts[[part]]} # ifsplit
norf<-length(orfs)
print(paste(norf,"ORFs in model"))
# Find out all the barcode-positions for each ORF #
print("Finding information about each repeat for each ORF") 
poslist<-lapply(orfs,positget,d) 
names(poslist)<-orfs
# Find max number of repeats & Photos #
replist<-sapply(poslist,nrepget)
maxnphlist<-sapply(poslist,nphget)
maxnph<-max(maxnphlist); maxrep<-max(replist)
# Set up arrays for growth and time & nph #
growth<-array(dim=c(norf,maxrep,maxnph))
time<-array(dim=c(norf,maxrep,maxnph))
nph<-matrix(nrow=norf,ncol=maxrep)
# Fill in arrays
print("Collecting and structuring dynamic data for JAGS")
for (k in 1:norf){orf<-orfs[k]; rep<-0
	# Print a report at quartiles of completeness
	if (k%in%(floor(1:3*norf*(1/4)))){print(paste("Done for ",k,"/",norf," ORFs",sep=""))}
	# Get list of positions and data for ORF
	orfpos<-poslist[[as.character(orf)]]; orfd<-d[d$ORF==orf,]
	# Get list of barcodes for this ORF
	bcodes<-unique(orfd$Barcode) 
	# Fill in arrays 
	for (b in bcodes){
		# Positions and data for this ORF on this Barcode
		bcpos<-orfpos[[as.character(b)]]; dbc<-orfd[orfd$Barcode==b,]
		# Get inoculation time for this barcode
		startt<-as.POSIXlt(as.character(dbc[1,'Inoc.Time']),format=fmt)
		# Get growth and time data for each spot on this barcode
		for (pos in bcpos){row<-pos[1]; col<-pos[2]; rep<-rep+1
			drep<-dbc[(dbc$Row==row) & (dbc$Col==col),]
			# Add number of photos for this repeat to matrix #
			nphotos<-length(drep[,1]); nph[k,rep]<-nphotos
			# Get times in days since inoculation #
			drep$t<-sapply(drep$Date.Time,tconv,startt,fmt)
			# Order data by increasing time from inoculation #
			drep<-drep[order(drep$t,drep$Growth),]
			# Get growth data in correct order with weirdness ironed out #
			growth[k,rep,1:nphotos]<-monofix(sapply(drep$Growth,nozero))
			# Put corresponding times in corresponding time array (no negatives) #
			time[k,rep,1:nphotos]<-sapply(drep$t,noneg)
					} #pos
							}#bcode
									}#k-orf			
dimnames(growth)[[1]]<-orfs; dimnames(time)[[1]]<-orfs; dimnames(nph)[[1]]<-orfs
ftime<-proc.time()-stime; trep("Building data structure for JAGS",ftime)
return(list(og=growth,time=time,E=norf,M=replist,N=nph,ORF=orfs))
} #jagdata	

##### Gets a list of all positions #####
positget<-function(orf,d){
dor<-d[d$ORF==orf,]; nreps<-0
positions<-list(); bcodes<-unique(dor$Barcode); nphmax<-0
for (j in 1:length(bcodes)){dorbc<-dor[dor$Barcode==bcodes[j],]
	pos2<-list(); for (k in 1:length(dorbc[,1])){
		pos2[[k]]<-c(dorbc[k,'Row'],dorbc[k,'Col'])}#k
	pos2<-unique(pos2); positions[[j]]<-pos2; nreps<-nreps+length(pos2)
	for (pos in pos2){drep<-dorbc[(dorbc$Row==pos[1])&(dorbc$Col==pos[2]),]
	nph<-length(drep[,1]); if (nph>nphmax){nphmax<-nph} } #pos
} #j
names(positions)<-bcodes; positions$nrep<-nreps; positions$nphmax<-nphmax
return(positions)
}

# To get max photos etc #
nphget<-function(orfpos){orfpos$nphmax}
nrepget<-function(orfpos){orfpos$nrep}

# Stops too implausible data that breaks updating #
monofix<-function(tc){
# Stop anomalous first point as it easily breaks the updating
if (length(tc)>2){if (tc[1]>(tc[2]+tc[3])){tc[1]<-max(1,0.9*min(tc[2],tc[3]))}}
# Stop flat at 1
if (length(tc[tc==1])==length(tc)){tc[length(tc)]<-2}
return(tc)
}

# Function to get rid of zero growth data #
nozero<-function(z){if (z<1){z=1} else {z=z}}

## Remove negative times ##
noneg<-function(z){if (z<0){z=0} else {z=z}}


######################### Build Bayesian Model with user input for priors #################################
qfa.model<-function(jagd,Kloc,gloc,Ksd="(2/3)*Kloc",g.range=10,...){require(rjags)
# Write model with user inputted guess for g and typical final colony size
write(paste("
################# Bayesian QFA Model ####################
model {
### User specified prior parameters ###
# Typical Inoculum level 
gloc<-",gloc,"
grange<-",g.range," 
# Typical final colony size
Kloc<-",Kloc,"
Ksd<-",Ksd," 

# For each ORF (mutant)
for (orf in 1:E){
      # For each repeat for that ORF
	for (i in 1:M[orf]){
            # For each photo for that repeat
		for (j in 1:N[orf,i]) {
			# Assume lognormal error
			og[orf,i,j] ~ dnorm(logm[orf,i,j],(sigma[orf]^(-2)))
			# Model growth as logistic
			logm[orf,i,j]<-(K[orf,i]*g*exp(r[orf,i]*time[orf,i,j]))/(K[orf,i]+g*(exp(r[orf,i]*time[orf,i,j])-1))
	} #j	
	
	## Priors for parameters for each repeat
	# Specific to an ORF (except for inoculum, same prior for all)
	K[orf,i] ~ dgamma(((k1[orf]^2)/(k2[orf]^2)),(k1[orf]/(k2[orf]^2))) 
	r[orf,i] ~ dgamma(((r1[orf]^2)/(r2[orf]^2)),(r1[orf]/(r2[orf]^2)))
	
	} #i
	
	## Priors for priors for parameters for each repeat
	# Specific to an experiment
	# Mean of K[orf,i]
	k1[orf] ~ dgamma(((k11^2)/(k12^2)),(k11/(k12^2)))
	# sd of K[orf,i]	
	k2[orf] ~ dgamma(((k21^2)/(k22^2)),(k21/(k22^2)))
	# Mean r[orf,i]
	r1[orf] ~ dgamma(((r11^2)/(r12^2)),(r11/(r12^2)))
	# sd r[orf,i]
	r2[orf] ~ dgamma(((r21^2)/(r22^2)),(r21/(r22^2)))
	# Prior for sigma[orf]
	sigma[orf] ~ dlnorm(log(s1),(s2^(-2)))

	} #orf

	## Priors for the parameters at the experiment level
	# so each experiment can have differing priors for each ORF's priors
	
	# K
	# Mean of (mean ORF K)
	k11 ~ dgamma((Kloc^2)/(Ksd^2),(Kloc/(Ksd^2)))
	# sd of (mean ORF K)
	k12m<-0.3*Kloc
	k12 ~ dlnorm(log(k12m),1)
	# Mean of (sd ORF K)
	k21m<-0.2*Kloc
	k21 ~ dlnorm(log(k21m),1)
	# sd of (sd ORF K)
	k22m<-0.085*Kloc
	k22 ~ dlnorm(log(k22m),1)
	
	# r
	# Mean of (mean ORF r) 
	r11 ~ dgamma(1.87,0.35)
	# sd of (mean ORF r)
	r12 ~ dlnorm(log(0.75),1)
	# Mean of (sd ORF r)
	r21 ~ dlnorm(log(0.75),0.85)
	# Sd of (sd ORF r)
	r22 ~ dlnorm(log(0.4),1)

	# sigma
	# Median of (median ORF sigma)
	s1m <- 0.4*Kloc
	s1 ~ dgamma(1,1/s1m)
	# sd of (median ORF sigma)
	s2m<-0.1*Kloc
	s2 ~ dlnorm(log(s2m),1)
	
	# Inoculum
	lbg<-(1-grange/100)*gloc
	ubg<-(1+grange/100)*gloc
	g ~ dunif(lbg,ubg)

} #model",sep=""),file="model.txt")

# Make JAGS model with data and model #
#jagd$ORF<-NULL
#jags.model(data=jagd,file="model.txt",...)
} #qfa.model


############ Check autocorrelations and plot correlograms to check convergence ############
qfa.acf<-function(posterior,lag.max=100,plot=TRUE,tol=0.2,printout=FALSE){
# Get maximum lag to calculate ACF for
if (lag.max==0){possize<-length(pos[,1])} else {possize<-lag.max}
# Make a list with each element giving acf statistics for each variable #
pars<-names(posterior)
acfs<-list(); nacmat<-c(); pnac<-c()
for (par in pars){acflist<-list()
	pos<-posterior[[par]]
	# Find ACF for each param & plot
	parnames<-dimnames(pos)[[2]]
	if (plot==TRUE){frame(); text(0.5,0.5,par,cex=5)}
	if (!is.null(parnames)){
	for (param in parnames){
		parn<-match(param,parnames)
		acflist[[parn]]<-acf(pos[,parn],lag.max=possize,main=param,plot=plot)
		if (plot==TRUE){abline(h=c(tol,-tol),v=10,col="red",lwd=1.5)}} #param
	names(acflist)<-parnames} else {acflist[[1]]<-acf(pos,lag.max=possize,main=par,plot=plot)
	if (plot==TRUE){abline(h=c(tol,-tol),v=10,col="red",lwd=1.5)}} #if.parnames
	# Number of autocorrelations with magnitude greater than tol #
	nacs<-sapply(acflist,nac,tol)
	# expressed as a percentage of the lags inspected
	pnacs<-100*nacs/(possize-10)
	if (printout==TRUE){print(paste(par,":",sep=""))
		print(paste("Number of autocorrelations (lag 10 ono.) with magnitude >",tol))
		print(nacs)
		print(paste("% of autocorrelations (lag 10 ono.) with magnitude >",tol))
		print(round(pnacs,3))} #print
	acfs[[match(par,pars)]]<-acflist; nacmat<-rbind(nacmat,nacs); pnac<-rbind(pnac,pnacs)} #par
names(acfs)<-pars; dimnames(nacmat)[[1]]<-pars; dimnames(pnac)[[1]]<-pars
list(ACF=acfs,nAC=nacmat,pAC=pnac)}


### Counts the number of magnitude>tol autocorrelations for an element of acflist ###
nac<-function(acs,tol){acs<-as.vector(acs[[1]])
acs<-acs[11:length(acs)]
length(acs[abs(acs)>tol])}

################# Makes posterior for each orf given list of orfs used ########################
########### and makes posterior organized by model variables. Suitable for analysis. ##########
qfa.simulate<-function(model,initial.update=5*10^4,samples=5*10^4,thins=10,orfs=c(),
params=c("k1","r1","g","k11","k12","k21","k22","r11","r12","r21","r22","s1","s2"),...){
# Do initial updating
print("Performing burn in and initial updating of model")
require(rjags)
update(model,initial.update)
# Sample from posterior using CODA
### Sample the posterior ###
print("Sampling the posterior")
x<-coda.samples(model,params,n.iter=samples,thin=thins,...)
# Get the matrix like part of the coda posterior
pm<-x[[1]]
### Build posterior object ###
varpos<-list()
# Get variable names from pm
posnames<-dimnames(pm)[[2]]
varnames<-strsplit(posnames,"\\[")
# Use this to extract variables' posteriors
for (param in params){paramposs<-sapply(varnames,varid,param)
varpos[[match(param,params)]]<-pm[,paramposs]}
names(varpos)<-params
# If we know what ORFs are in the model, add the names to relevant places
if (length(orfs)>0){
# Find out which parameters aren't there for each ORF
dims<-lapply(lapply(varpos,dim),colcheck)
specpars<-c(); norfs<-length(orfs)
for (k in 1:length(dims)){
	if (is.null(dims[[k]])){specpars<-append(specpars,params[k])}}
# Give the variables the names of their ORFs
orfpars<-params[!params%in%specpars]
for (par in orfpars){dimnames(varpos[[par]])[2][[1]]<-orfs}}
varpos}


# Add the posteriors for these to the orf organized posteriors
#for (spar in specpars){parn<-match(spar,specpars)
#orfpos$par[[parn]]<-varpos[[match(spar,names(varpos))]]}
#names(orfpos$par)[1:length(specpars)]<-specpars


#### Checks the identity of a split variable name ####
varid<-function(splitname,varname){splitv<-splitname[1]
if (splitv==varname){z=TRUE} else {z=FALSE}
z}


# checks number of columns in posterior matrix #
colcheck<-function(dim){dim[2]}  

# Removes null elements of orf posterior #
nullremove<-function(orfp,specpars){
for (par in specpars){orfp[[par]]<-NULL}
orfp}

############################### Likelihood Functions ################################
#### Upper and lower bounds for variables in optimization ####
lowerg<-function(inocguess){0.9*inocguess}; upperg<-function(inocguess){1.1*inocguess}
lowerr<-0; upperr<-function(rguess){2.5*rguess}
lowerK<-function(inocguess){0.9*inocguess}; upperK<-function(xdim,ydim){0.75*255*xdim*ydim}
# The lower bound on s seems to be causing problems with inf arising during optimisation
#lowers<-function(xdim,ydim){0.025*xdim*ydim}; uppers<-function(xdim,ydim){175*xdim*ydim}
lowers<-function(xdim,ydim){0.00025*xdim*ydim}; uppers<-function(xdim,ydim){50*xdim*ydim}

##### Does max. lik. fit for all colonies, given rod.read input #####
qfa.fit<-function(d,inocguess,ORF2gene="ORF2GENE.txt",fmt="%Y-%m-%d_%H-%M-%S",...){
# Remove edge colonies if edgestrip true
# Get xdim & ydim
xdim<-d[1,'Tile.Dimensions.X']; ydim<-d[1,'Tile.Dimensions.Y']
# Define optimization bounds reliant on xdim & ydim and inocguess #
lows<-lowers(xdim,ydim); ups<-uppers(xdim,ydim)
lowK<-lowerK(inocguess); upK<-upperK(xdim,ydim) 
lowg<-lowerg(inocguess); upg<-upperg(inocguess)
xybounds<-list(s=c(lows,ups),K=c(lowK,upK),g=c(lowg,upg))
# Create orf2gene dictionary
if (ORF2gene!=FALSE){gdict<-orf2gdict(ORF2gene)}
# Vector of barcodes
barcodes<-unique(d$Barcode); nbc<-length(barcodes)
# Get big data frame ready for results
results<-data.frame()
# For each barcode, optimize
bcount<-0
for (bcode in barcodes){bcount<-bcount+1
	# Say which plate you're optimizing
	print(paste("Optimizing Plate",bcount,"/",nbc,":",bcode))
	# Restrict data to just this barcode
	dbc<-d[d$Barcode==bcode,]
	# Get inoculation time for this plate
	inoctime<-as.POSIXlt(as.character(dbc[1,'Inoc.Time']),format=fmt)
	# Get list of unique colony positions
	positions<-lapply(1:length(dbc[,1]),index2pos,dbc)
	positions<-unique(positions)
	# Fit logistic model to each colony
	bcfit<-t(sapply(positions,colony.fit,dbc,inoctime,xybounds,fmt=fmt,...))
	info<-t(sapply(positions,colony.info,dbc))
	rows<-sapply(positions,rcget,"row")
	cols<-sapply(positions,rcget,"col")
	# Bind Data frame of barcode results to overall results
	results<-rbind(results,data.frame(Barcode=info[,1],Row=rows,Col=cols,
	Background=info[,9],Treatment=info[,2],Medium=info[,3],ORF=info[,4],
	K=bcfit[,1],r=bcfit[,2],g=bcfit[,3],s=bcfit[,4],Screen.Name=info[,5],
	Library.Name=info[,6],MasterPlate.Number=info[,7],Timeseries.order=info[,8],
	Inoc.Time=inoctime))
				} #bcode
if (ORF2gene!=FALSE){results$Gene<-sapply(as.character(results$ORF),orf2g,gdict)}
return(results)}

# Extrace row & col from list of position vectors
rcget<-function(posvec,rc){posvec[match(rc,c("row","col"))]}

### Function that does the optimization for one colony ###
colony.fit<-function(position,bcdata,inoctime,xybounds,fmt,...){
# Retrieve inoculation guess
inocguess<-(10/9)*xybounds$g[1]
# Get row & column to restrict data
row<-position[1]; col<-position[2]
d<-bcdata[(bcdata$Row==row)&(bcdata$Col==col),]
# Get growth data
growth<-sapply(d$Growth,nozero)
# Get time data
time<-sapply(d$Date.Time,tconv,inoctime,fmt)
assign("last.colony.growth",growth,envir=.GlobalEnv)
assign("last.colony.time",time,envir=.GlobalEnv)
# Check if worth optimizing
growththresh<-(1.5/255)*xybounds$K[2]
#if ((cor(time,log(growth))>0.1)&
if (max(growth)>growththresh){
	# Get initial guess for parameters
	init<-guess(time,growth,inocguess,xybounds)
	# Function to be optimized
	objf<-function(modpars){
	K<-modpars[1]; r<-modpars[2]
	g<-modpars[3]; s<-modpars[4]
	loglik(K,r,g,s,growth,time)}
	# Gradient of objf
	objfgr<-function(pars){
	K<-pars[1]; r<-pars[2]
	g<-pars[3]; s<-pars[4]
	loglikgr(K,r,g,s,growth,time)}
	# Perform optimization
	optsol<-optim(par=init,fn=objf,gr=objfgr,method="L-BFGS-B",
	lower=c(xybounds$K[1],lowerr,xybounds$g[1],xybounds$s[1]),
	upper=c(xybounds$K[2],upperr(init[2]),xybounds$g[2],xybounds$s[2]),...)
	pars<-optsol$par} else {pars<-c(inocguess,0,inocguess,Inf)}
pars}

# Logistic growth curve function #
logist<-function(K,r,g,t){
(K*g*exp(r*t))/(K+g*(exp(r*t)-1))}

# Log-likelihood (times(-1)) #
loglik<-function(K,r,g,s,growth,time){
sum(log(s)+(1/(2*s*s))*((growth-logist(K,r,g,time))^2))}

# partial derivative of -loglik wrt r #
dldr<-function(K,r,g,s,growth,time){
-(1/(s*s))*sum(
(growth-logist(K,r,g,time))*time*logist(K,r,g,time)*(1-logist(K,r,g,time)/K)
)}

# partial derivative of -loglik wrt g #
dldg<-function(K,r,g,s,growth,time){
-(1/(s*s))*sum(
(growth-logist(K,r,g,time))*((K*logist(K,r,g,time))/(g*(K+g*(exp(r*time)-1))))
)}

# partial derivative of -loglik wrt s #
dlds<-function(K,r,g,s,growth,time){
-sum(
(1/(s^3))*(growth-logist(K,r,g,time))-1/s
)}

# partial derivative of -loglik wrt K #
dldK<-function(K,r,g,s,growth,time){
-(1/(s*s))*sum(
(growth-logist(K,r,g,time))*((logist(K,r,g,time)*g*(exp(r*time)-1))/((K+g*exp(r*time))*K))
)}

# Gradient of log-likelihood
loglikgr<-function(K,r,g,s,growth,time){
c(dldK(K,r,g,s,growth,time),dldr(K,r,g,s,growth,time),
dldg(K,r,g,s,growth,time),dlds(K,r,g,s,growth,time))}
	
# Prevent zero growth from occurring
nozero<-function(growth){if (growth==0){growth<-1} else {growth<-growth}}

### Provide initial guess for logistic model parameters ###
guess<-function(time,growth,inocguess,xybounds){
n=length(time)
# g
G0g<-inocguess
# K
Kg<-max(growth)
#r
targ<-G0g+(Kg-G0g)/2.0
tmrate<-20.0
x=1
ct=0
while (x==1){
ct=ct+1
tmrate<-time[ct]+((targ-growth[ct])/(growth[ct+1]-growth[ct]))*(time[ct+1]-time[ct])
if (growth[ct+1]>=targ & growth[ct]<targ){x=0} else if (ct==(n-1)) {x=0
tmrate=20.0} else {x=1}
}#x
rg<-log(max(0.000000000001,(Kg-G0g)/G0g))/tmrate
rg<-max(lowerr,rg)
#s
s<-sqrt((1/n)*sum((growth-logist(Kg,rg,G0g,time))^2))
s<-max(xybounds$s[1],s); s<-min(s,xybounds$s[2])
# Out 
c(Kg,rg,G0g,s)					
}#guess


	
#### Get colony information ####
colony.info<-function(position,bcdata){
# Get row & column to restrict data
row<-position[1]; col<-position[2]
d<-bcdata[(bcdata$Row==row)&(bcdata$Col==col),]
d<-d[1,]
# Character vector of colony info
c(as.character(d$Barcode),
as.character(d$Treatments),as.character(d$Medium),as.character(d$ORF),
as.character(d$Screen.Name),as.character(d$Library.Name),as.character(d$MasterPlate.Number),
as.character(d$Timeseries.order),as.character(d$Background))}


##### Make PDFs #####
qfa.plot<-function(file,results,d,fmt="%Y-%m-%d_%H-%M-%S",barcodes=c(),master.plates=c(),
treatments=c(),screen.names=c(),backgrounds=c()){
# Get character vectors of requested barcodes, treatements,etc.; all if none specified
if (length(barcodes)==0){barcodes<-unique(results$Barcode)}
results<-results[results$Barcode%in%barcodes,]
d<-d[d$Barcode%in%barcodes,]
if (length(master.plates)==0){master.plates<-unique(results$MasterPlate.Number)} 
results<-results[results$MasterPlate.Number%in%master.plates,]
d<-d[d$MasterPlate.Number%in%master.plates,]
if (length(treatments)==0){treatments<-unique(results$Treatment)} 
results<-results[results$Treatment%in%treatments,]
d<-d[d$Treatment%in%treatments,]
if (length(screen.names)==0){screen.names<-unique(results$Screen.Name)} 
results<-results[results$Screen.Name%in%screen.names,]
d<-d[d$Screen.Name%in%screen.names,]
if (length(backgrounds)==0){backgrounds<-unique(results$Background)} 
results<-results[results$Background%in%backgrounds,]
d<-d[d$Background%in%backgrounds,]
print(paste("Plotting for",length(results[,1]),"colonies"))
# Produce PDF
colmax<-max(results$Col); rowmax<-max(results$Row)
pdf(file,4*colmax,4*rowmax)
cexfctr<-(rowmax*colmax)/384
bcount<-0; nbc<-length(barcodes)
for (bcode in barcodes){bcode<-as.character(bcode)
	# Say what you're doing
	bcount<-bcount+1
	print(paste("Plotting for Plate",bcount,"/",nbc,":",bcode))
	# Restricts results and data to this plate
	rbc<-results[as.character(results$Barcode)==bcode,]
	dbc<-d[as.character(d$Barcode)==bcode,]
	# Get starting time for this plate
	inoctime<-rbc$Inoc.Time[1]
	# Find number of rows and columns on this plate
	nrow<-max(rbc$Row)-min(rbc$Row)+1; ncol<-max(rbc$Col)-min(rbc$Col)+1
	# Find max growth to set y-axis for this plate
	maxg<-max(dbc$Growth)
	# Set graphics parameters for each plate
	par(mfrow=c(nrow,ncol),oma=c(13,15,22,1),
	mar=c(2,1,2,0.75),mgp=c(3,1,0),cex=cexfctr)
	## Plot for each row of results for that bcode ##
	z<-apply(rbc,1,rowplot,dbc,inoctime,maxg,fmt)
	# Title for the plate
	maintit<-paste(rbc$Barcode[1],"Treatment:",rbc$Treatment[1],
	"Medium:",rbc$Medium[1],"Plate:",rbc$MasterPlate.Number[1],sep=" ")
	cextit<-1008/length(strsplit(maintit,split="")[[1]])
	title(main=maintit,xlab="Time since inoculation (days)",line=7,
	ylab="Colony Size (Trimmed Greyscale)",cex.main=9,cex.lab=8,outer=TRUE)} #bcode
dev.off()}


#### Plot a colony's timecourse from a row of the results #####	
rowplot<-function(resrow,dbc,inoctime,maxg,fmt){
# Get logistic parameters, gene name and position
K<-as.numeric(resrow['K']); r<-as.numeric(resrow['r']); g<-as.numeric(resrow['g']);
row<-as.numeric(resrow['Row']); col<-as.numeric(resrow['Col'])
if ('Gene'%in%names(resrow)){gene<-resrow['Gene']} else {gene<-resrow['ORF']}
# Get data for that colony
dcol<-dbc[(dbc$Row==row)&(dbc$Col==col),]
growth<-sapply(dcol$Growth,nozero)
time<-sapply(dcol$Date.Time,tconv,inoctime,fmt)
# Draw the curves and data
logdraw(row,col,K,r,g,time,growth,gene,maxg,mdrmdp)}



### Converts row no. to position vector ###	
index2pos<-function(index,dbc){
c(dbc[index,'Row'],dbc[index,'Col'])}

### Do individual timecourse plot given parameters & data ###
logdraw<-function(row,col,K,r,g,time,growth,gene,maxg,fitfunct){
# Add logistic curve
curve((K*g*exp(r*x))/(K+g*(exp(r*x)-1)),n=100,xlim=c(0,ceiling(max(time))),ylim=c(0,1.2*maxg),
xlab="",ylab="",main=gene,frame.plot=0,cex.main=3,cex.axis=1,lwd=2.5)
# Add data points
points(time,growth,col="red",cex=2,pch=4,lwd=2)
# Add legend
# With fitness
F<-fitfunct(K,r,g)
legt1<-paste(c("K=","r=","g=","F="),c(round(K,0),round(r,3),round(g,0),round(F,2)),sep="")
legend(0,1.2*maxg,legt1,box.lty=0)
legend(0.9*max(time),1.2*maxg,c(row,col),box.lty=0,cex=0.5)
}

