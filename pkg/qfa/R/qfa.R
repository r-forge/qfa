#### Define Phenotype for fitness ####
mdrmdp<-function(K,r,g){MDR<-sapply(r/log((2*(K-g))/(K-2*g)),na2zero)
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

# Converts NAs to zeros # Should get rid of this - CONOR.
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
orf2g$ORF=toupper(orf2g$ORF)
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
qfa.epi<-function(double,control,qthresh,orfdict="ORF2GENE.txt",
GISthresh=0.0,plot=TRUE,modcheck=TRUE,fitfunct=mdrmdp,wctest=TRUE){
###### Get ORF median fitnesses for control & double #######
print("Calculating median (or mean) fitness for each ORF")
## LIK ##
# Get orfs in question
orfs<-as.character(double$ORF)
orfs<-orfs[orfs%in%as.character(control$ORF)]
orfs<-unique(orfs)
# Get lists with fitnesses for each repeat
cFstats<-lapply(orfs,orfstat,control,fitfunct)
names(cFstats)<-orfs
dFstats<-lapply(orfs,orfstat,double,fitfunct)
names(dFstats)<-orfs
# Get means or medians for each ORF
if(wctest){cFms<-sapply(cFstats,median)}else{cFms<-sapply(cFstats,mean)}
names(cFms)<-orfs
if(wctest){dFms<-sapply(dFstats,median)}else{dFms<-sapply(dFstats,mean)}
names(dFms)<-orfs
### Fit genetic independence model ###
m<-lm.epi(dFms,cFms,modcheck)
print(paste("Ratio of background mutant fitness to wildtype fitness =",round(m,4)))
###### Estimate probability of interaction #######
print("Calculating interaction probabilities")
pg<-sapply(orfs,pgis,m,cFstats,dFstats,cFms,dFms,wilcoxon=wctest)
pg<-as.data.frame(t(pg))
colnames(pg)=c("p","gis")
p<-pg$p

# Adjust for multiple comparisons
q<-p.adjust(p,"fdr")
# Get gene names if needed
# lik
genes<-as.character(double$Gene[match(orfs,double$ORF)]); 
#genes<-genes[genes%in%as.character(control$Gene)]
# If gene names missing, find them with orf2g
if (is.null(genes)){# Create orf-gene dictionary
	orfdict<-orf2gdict(orfdict)
	genes<-sapply(orfs,orf2g,orfdict)
}
	
# Get genetic interaction scores
meandiff<-mean(dFms-cFms)
#gis<-dFms/mean(dFms)-cFms/mean(cFms)
gis<-pg$gis/meandiff
# Put into data.frame
results<-data.frame(ORF=orfs,Gene=genes,P=p,Q=q,GIS=gis,Median.Double=dFms,Median.Control=cFms)
results$Type<-apply(results,1,typemake,m)
results<-results[order(results$GIS,results$Q,results$Type),]
# Get rid of duplicate entries in results
orflist<-unique(as.character(results$ORF))
results<-results[match(orflist,results$ORF),]
# Plot results
if (plot==TRUE){qfa.epiplot(results,qthresh,m)}
final<-list(Results=results,
Enhancers=gethits(results,qthresh,type="E",GISthresh=GISthresh),
Suppressors=gethits(results,qthresh,type="S",GISthresh=GISthresh))
return(final)
}

############### Epistasis Functions ##################
# Makes epistasis plot for a given fdr level #
qfa.epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",yylab="Query Fitness",mmain="Epistasis Plot",fmax=0){

enhancers<-gethits(results,qthresh,type="E")
suppressors<-gethits(results,qthresh,type="S")
others<-results[(!results$ORF%in%enhancers$ORF)&(!results$ORF%in%suppressors$ORF),]
if (fmax==0){
	# Get plot parameters
	ymax=1.1*max(results$Median.Double); ymin=0
	xmax=1.1*max(results$Median.Control); xmin=0
}else{
	ymin=0;ymax=fmax
	xmin=0;xmax=fmax
}
plot(NULL,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
xlab=xxlab,ylab=yylab,main=mmain,
col=8,pch=19,cex=0.5)
# Add line for genetic independence
if (fitratio!=FALSE){abline(0,fitratio,lwd=2,col=8)} else {
abline(0,lm.epi(results$Median.Double,results$Median.Control,modcheck=FALSE),lwd=2,col=8)}
# Add 1:1 fitness line
abline(0,1,lwd=2,lty=4,col=8)
# Add reference ORF fitnesses lines
if (ref.orf!=FALSE){reforf<-results[results$ORF==ref.orf,]
abline(v=reforf$Median.Control,col="lightblue",lwd=2)
abline(h=reforf$Median.Double,col="lightblue",lwd=2)}
# Add points for non-suppressors & non-enhancers
points(others$Median.Control,others$Median.Double,col="grey",cex=0.5,pch=19)
text(others$Median.Control,others$Median.Double,others$Gene,col="grey",pos=4,offset=0.1,cex=0.4)
# Add suppressors & enhancers
if(length(enhancers$ORF)>0){
	points(enhancers$Median.Control,enhancers$Median.Double,col='green',pch=19,cex=0.5)
	text(enhancers$Median.Control,enhancers$Median.Double,enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
}
if(length(suppressors$ORF)>0){
	points(suppressors$Median.Control,suppressors$Median.Double,col='red',pch=19,cex=0.5)
	text(suppressors$Median.Control,suppressors$Median.Double,suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
}
}

## Extract hits from epistasis results object ##
gethits<-function(results,qthresh,type="S",all.types=FALSE,GISthresh=0){
results<-results[results$Q<qthresh,]
if (all.types==FALSE){results<-results[results$Type==type,]}
results<-results[abs(results$GIS)>=GISthresh,]
results<-results[!is.na(results$ORF),]
return(results[order(results$GIS,results$Q,results$Type),])
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

# Estimates probability of no interaction and estimated strength of interaction for max lik method #
pgis<-function(orf,m,cFs,dFs,cFms,dFms,wilcoxon=TRUE){
	# If this orf is not present in both lists, return appropriate p,gis
	if((length(dFs[[orf]])==0)|(length(cFs[[orf]])==0)){return(c(1,0))}
	# If fitnesses are not unique in one condition (e.g. all dead) and only one repeat in another (e.g. after stripping)
	# This would cause wilcoxon test to fail due to ties, and t.test to fail due to insufficient y?
	ldFS=length(dFs[[orf]]); lcFS=length(cFs[[orf]])
	ludFS=length(unique(dFs[[orf]])); lucFS=length(unique(cFs[[orf]]))
	if (((ludFS==1)&(ldFS>1)&(lcFS==1))|((lucFS==1)&(lcFS>1)&(ldFS==1))){return(c(1,median(dFs[[orf]],na.rm=TRUE)-m*median(cFs[[orf]],na.rm=TRUE)))}
	# If both sets of cultures are dead (and fitnesses therefore equal) return appropriate p,gis
	if(sum(dFs[[orf]]==m*cFs[[orf]])==length(dFs[[orf]]==m*cFs[[orf]])){
		return(c(1,0))
	}else{
		if (wilcoxon){
			# Returns p-value for significance of difference, and estimate of difference between medians
			ctest<-wilcox.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
			p<-as.real(ctest$p.value)
			diff<-as.real(ctest$estimate)
			return(c(p,diff))
		}else{
			# t-test fails if only one element in either list, do one sample test
			if((ldFS<=1)|(lcFS<=1)){
				ctest<-t.test(dFs[[orf]]-m*cFs[[orf]])
				p<-as.real(ctest$p.value)
				diff<-as.real(ctest$estimate)
				return(c(p,diff))
			}else{
				# Returns p-value for significance of difference, and the difference between the means
				ctest<-t.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
				p<-as.real(ctest$p.value)
				diff<-as.real(ctest$estimate)
				diff<-diff[1]-diff[2]
				return(c(p,diff))
			}
		}
	}
}

# Function to make fitnesses
fitmake<-function(orf,poslist,fitfunct,g){
orfpos<-poslist[[orf]]
return(fitfunct(orfpos$k1,orfpos$r1,g))
}

# Get type of interaction
typemake<-function(row,m){md<-as.numeric(row['Median.Double']); c<-as.numeric(row['Median.Control'])
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

############################### Likelihood Functions ################################
#### Upper and lower bounds for variables in optimization ####
lowerg<-function(inocguess){0.9*inocguess}; upperg<-function(inocguess){1.1*inocguess}
lowerr<-0; upperr<-function(rguess){2.5*rguess}
lowerK<-function(inocguess){0.9*inocguess}; upperK<-function(xdim,ydim){0.75*255*xdim*ydim}
# The lower bound on s seems to be causing problems with inf arising during optimisation
#lowers<-function(xdim,ydim){0.025*xdim*ydim}; uppers<-function(xdim,ydim){175*xdim*ydim}
lowers<-function(xdim,ydim){0.00025*xdim*ydim}; uppers<-function(xdim,ydim){50*xdim*ydim}

##### Does max. lik. fit for all colonies, given colonyzer.read or rod.read input #####
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
	Inoc.Time=inoctime,TileX=info[,10],TileY=info[,11],XOffset=info[,12],YOffset=info[,13],
	Threshold=info[,14],EdgeLength=info[,15],EdgePixels=info[,16],RepQuad=info[,17]))
} #bcode
if (ORF2gene!=FALSE){results$Gene<-sapply(as.character(results$ORF),orf2g,gdict)}
return(results)
}

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
pars<-optsol$par
# Sanity check for fitted parameters (no negative growth)
if (pars[1]<pars[3]){
	pars[1]=pars[3]
	pars[2]=0
	n<-length(growth)
	pars[4]<-sqrt((1/n)*sum((growth-logist(pars[1],pars[2],pars[3],time))^2))
}
return(pars)
}

# Logistic growth curve function #
logist<-function(K,r,g,t){
(K*g*exp(r*t))/(K+g*(exp(r*t)-1))}

# Log-likelihood (times(-1)) #
loglik<-function(K,r,g,s,growth,time){
LL<-sum(log(s)+(1/(2*s*s))*((growth-logist(K,r,g,time))^2))
LL<-min(LL,999999999)
LL<-max(0.0000000001,LL)
if(is.na(LL)){LL<-0}
return(LL)
}

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
MaxRateFound=FALSE
ct=0
# Run through time points looking for maximum rate of change of density
while (MaxRateFound==FALSE){
ct=ct+1
tmrate<-time[ct]+((targ-growth[ct])/(growth[ct+1]-growth[ct]))*(time[ct+1]-time[ct])
if (growth[ct+1]>=targ & growth[ct]<targ){MaxRateFound=TRUE} 
# If we haven't found a maximum by the end, arbitrarily set tmrate to 20
else if (ct==(n-1)) {MaxRateFound=TRUE; tmrate=20;} 
else {MaxRateFound==FALSE}
}#MaxRateFound
rg<-log(max(0.000000000001,(Kg-G0g)/G0g))/tmrate
# Set initial guess to half of best estimate to bias against
# artificially large r
rg<-0.5*rg
# Sanity check for guessed parameter values
# If the data have low correlation, then set r=0 and K=g0
# If all elements of growth are equal (e.g. zero) then correlation function throws error...
if(length(unique(growth))==1){
	rg=0
	Kg=G0g
}else{
	if (cor(time,log(growth))<0.1){
		rg=0
		Kg=G0g
	}
}
rg<-max(lowerr,rg)
# Hardcode in an upper limit on r to remove floating point problems
rg<-min(50,rg)
#s
s<-sqrt((1/n)*sum((growth-logist(Kg,rg,G0g,time))^2))
s<-max(xybounds$s[1],s); s<-min(s,xybounds$s[2])
# Out 
guessparams=c(Kg,rg,G0g,s)					
return(guessparams)
}#guess
	
#### Get colony information ####
colony.info<-function(position,bcdata){
# Get row & column to restrict data
row<-position[1]; col<-position[2]
d<-bcdata[(bcdata$Row==row)&(bcdata$Col==col),]
# Want to use the LAST datapoint 
# Most relevant for edge length (possibly most reliable for position)
d<-d[order(d$Date.Time,decreasing=TRUE),]
d<-d[1,]
# Character vector of colony info
c(as.character(d$Barcode),
as.character(d$Treatments),as.character(d$Medium),as.character(d$ORF),
as.character(d$Screen.Name),as.character(d$Library.Name),as.character(d$MasterPlate.Number),
as.character(d$Timeseries.order),as.character(d$Background),as.numeric(d$Tile.Dimensions.X),as.numeric(d$Tile.Dimensions.Y),
as.numeric(d$X.Offset),as.numeric(d$Y.Offset),as.numeric(d$Threshold),as.numeric(d$Edge.length),as.numeric(d$Edge.Pixels),as.numeric(d$RepQuad))}


##### Make PDFs #####
qfa.plot<-function(file,results,d,fmt="%Y-%m-%d_%H-%M-%S",barcodes=c(),master.plates=c(),
treatments=c(),screen.names=c(),backgrounds=c(),maxg=0,maxt=0){
# Sort the data to be plotted sensibly, allowing easy comparison between repeats
results=results[order(results$MasterPlate.Number,results$Treatment,results$Screen.Name),]

# Get character vectors of requested barcodes, treatements,etc.; all if none specified
if (length(master.plates)==0){master.plates<-unique(results$MasterPlate.Number)} 
results<-results[results$MasterPlate.Number%in%master.plates,]
d<-d[d$MasterPlate.Number%in%master.plates,]
if (length(treatments)==0){treatments<-unique(results$Treatment)} 
results<-results[results$Treatment%in%treatments,]
d<-d[d$Treatments%in%treatments,]
if (length(screen.names)==0){screen.names<-unique(results$Screen.Name)} 
results<-results[results$Screen.Name%in%screen.names,]
d<-d[d$Screen.Name%in%screen.names,]
if (length(backgrounds)==0){backgrounds<-unique(results$Background)} 
results<-results[results$Background%in%backgrounds,]
d<-d[d$Background%in%backgrounds,]
if (length(barcodes)==0){barcodes<-unique(results$Barcode)}
results<-results[results$Barcode%in%barcodes,]
d<-d[d$Barcode%in%barcodes,]
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
	#as.POSIXlt(as.character(bigd$Date.Time),format=fmt)
	# Find number of rows and columns on this plate
	nrow<-max(rbc$Row)-min(rbc$Row)+1; ncol<-max(rbc$Col)-min(rbc$Col)+1
	# Find max growth to set y-axis for this plate
	if (maxg==0) maxg=max(dbc$Growth)
	# Set graphics parameters for each plate
	op<-par(mfrow=c(nrow,ncol),oma=c(13,15,22,1),
	mar=c(2,1,2,0.75),mgp=c(3,1,0),cex=cexfctr)
	## Plot for each row of results for that bcode ##
	z<-apply(rbc,1,rowplot,dbc,inoctime,maxg,fmt,maxt)
	# Title for the plate
	maintit<-paste(rbc$Barcode[1],"Treatment:",rbc$Treatment[1],
	"Medium:",rbc$Medium[1],"Plate:",rbc$MasterPlate.Number[1],sep=" ")
	cextit<-1008/length(strsplit(maintit,split="")[[1]])
	title(main=maintit,xlab="Time since inoculation (days)",line=7,
	ylab="Cell Density (AU)",cex.main=9,cex.lab=8,outer=TRUE)
      par(op)} #bcode
dev.off()}


#### Plot a colony's timecourse from a row of the results #####	
rowplot<-function(resrow,dbc,inoctime,maxg,fmt,maxt){
# Get logistic parameters, gene name and position
K<-as.numeric(resrow['K']); r<-as.numeric(resrow['r']); g<-as.numeric(resrow['g']);
row<-as.numeric(resrow['Row']); col<-as.numeric(resrow['Col'])
if ('Gene'%in%names(resrow)){gene<-resrow['Gene']} else {gene<-resrow['ORF']}
# Get data for that colony
dcol<-dbc[(dbc$Row==row)&(dbc$Col==col),]
growth<-sapply(dcol$Growth,nozero)
time<-sapply(dcol$Date.Time,tconv,inoctime,fmt)
# Draw the curves and data
logdraw(row,col,K,r,g,time,growth,gene,maxg,mdrmdp,maxt)}



### Converts row no. to position vector ###	
index2pos<-function(index,dbc){
c(dbc[index,'Row'],dbc[index,'Col'])}

### Do individual timecourse plot given parameters & data ###
logdraw<-function(row,col,K,r,g,time,growth,gene,maxg,fitfunct,maxt){
# Add logistic curve
if(maxt==0) maxt=ceiling(max(time))
curve((K*g*exp(r*x))/(K+g*(exp(r*x)-1)),n=31,xlim=c(0,maxt),ylim=c(0,1.2*maxg),
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

