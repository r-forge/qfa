#### Define Phenotype for fitness ####
#mdr<-function(K,r,g,v) sapply((r*v)/(log((K^v-g^v)/(K^v-(2^v)*(g^v)))+log(2)*v),na2zero)
mdr<-function(K,r,g,v) sapply((r*v)/log(1-(2^v-1)/((2^v)*(g/K)^v-1)),na2zero)
mdp<-function(K,r,g,v) sapply(log(K/g)/log(2),na2zero)
mdrmdp<-function(K,r,g,v) mdr(K,r,g,v)*mdp(K,r,g,v)

# Doubling time for generalised logistic function, as a function of time t
dtl<-function(K,r,g,v,t) sapply((-(r*t*v) + log((2**v*(1 - (K/g)**v)*((1 + (-1 + (K/g)**v)/exp(r*t*v))**(-1/v))**v)/(-1 + 2**v*((1 + (-1 + (K/g)**v)/exp(r*t*v))**(-1/v))**v)))/(r*v),na2zero)

# Logistic growth model #
logist<-function(K,r,g,t) (K*g*exp(r*t))/(K+g*(exp(r*t)-1))

# Generalised logistic growth model #
# Solution of dx/dt = r*x*(1-(x/K)^v)
Glogist<-function(K,r,g,v,t) K/(1+(-1+(K/g)^v)*exp(-r*v*t))^(1/v)

# Maximum slope of generalised logistic growth curve (on linear scale)
gen_logistic_maxslp<-function(K,r,g,v) r*v*K/(v+1)**((v+1)/v)

####### Normalise a column in a data frame by plate, grouping by the groupcol column (e.g. Treatment, Medium, or some combination)...
####### This function eradicates any plate effect within the group specified by groupcol
normalisePlates=function(d,column,groupcol="Treatment",fmin=0.000001){
	d[,groupcol]=as.character(d[,groupcol]); d[,column]=as.numeric(d[,column]); d$Barcode=as.character(d$Barcode)
	for(grouplabel in sort(unique(d[,groupcol]))){
		print(paste("Normalising plates by group",grouplabel))
		med=median(d[d[,groupcol]==grouplabel,column])
		for (b in unique(d$Barcode[d[,groupcol]==grouplabel])){
			datlst=as.numeric(d[(d$Barcode==b)&(d[,groupcol]==grouplabel),column])
			p.med=median(datlst)
			if(p.med>=fmin) datlst=datlst*med/p.med # Avoid division by zero median fitness...
			d[(d$Barcode==b)&(d[,groupcol]==grouplabel),column]=datlst
		}		
	}
	return(as.numeric(d[,column]))
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
na2zero<-function(num) {if (is.na(num)) num<-0; return(num)}

###### Create inoculation time environment #####
bctimes<-function(bctimef){
	startframe=read.delim(bctimef,colClasses=c("character"))
	starts<-new.env(hash=TRUE)
	for (k in 1:length(startframe$Barcode)){
		st=startframe$Start.Time[k]
		assign(startframe$Barcode[k],st,envir=starts)
	}
	return(starts)
} # bctimes

# Convert barcode to start time #
bc2st<-function(bc,inocenv) get(as.character(bc),envir=inocenv)

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

## Standard error of mean
sterr <- function(x) sqrt(var(x)/length(x))

summarise <- function(df,fdef="fit",summarystat=mean){
	summ=aggregate(formula(paste(fdef,"ORF",sep="~")),df,summarystat)
	return(setNames(summ[[fdef]],summ$ORF))
}

################################################## Epistasis Function ###########################################################
qfa.epi<-function(double,control,qthresh=0.05,orfdict="ORF2GENE.txt",
	GISthresh=0.0,plot=TRUE,modcheck=FALSE,wctest=TRUE,bootstrap=NULL,Nboot=5000,subSamp=Inf,reg="lmreg",fdef="fit"){
	###### Get ORF median fitnesses for control & double #######
	print("Calculating median (or mean) fitness for each ORF")
	## LIK ##
	# Get orfs in question
	dorfs<-unique(as.character(double$ORF))
	corfs<-unique(as.character(control$ORF))
	orfs<-sort(intersect(dorfs,corfs))
	control=control[control$ORF%in%orfs,]
	double=double[double$ORF%in%orfs,]
	
	# Get lists with fitnesses for each repeat
	cFstats<-split(control[[fdef]],control$ORF)
	dFstats<-split(double[[fdef]],double$ORF)
	# Get means or medians for each ORF
	if(!is.null(bootstrap)){
		cFms<-sapply(cFstats,bootstrap)
		dFms<-sapply(dFstats,bootstrap)	
	}else{
		if(wctest){
			cFms<-summarise(control,fdef,median)
			dFms<-summarise(double,fdef,median)
		}else{
			cFms<-summarise(control,fdef,mean)
			dFms<-summarise(double,fdef,mean)
		}
	}
	cSe<-summarise(control,fdef,sterr)
	dSe<-summarise(double,fdef,sterr)
	cCount<-summarise(control,fdef,length)
	dCount<-summarise(double,fdef,length)
	
	cFitDef<-findFit(control)
	dFitDef<-findFit(double)

	diffFrac=function(GIS,direction=`>`){
		return(length(GIS[direction(GIS,0)])/length(GIS))
	}
	
	### Fit genetic independence model ###
	m<-lm.epi(dFms,cFms,modcheck,reg=reg)
	print(paste("Ratio of background mutant fitness to wildtype fitness =",round(m,4)))
	###### Estimate probability of interaction #######
	print("Calculating interaction probabilities")
	if(!is.null(bootstrap)){
		GISlist<-sapply(orfs,pgis_bootstrap,cFstats,dFstats,sampSumm=bootstrap,Nreps=Nboot,subSamp=subSamp)
		g=apply(GISlist,2,bootstrap)
		greaterFrac=apply(GISlist,2,diffFrac,`>`)
		lessFrac=apply(GISlist,2,diffFrac,`<`)
		p=ifelse(g<0,1-lessFrac,1-greaterFrac)
		pg=rbind(p,g)
	}else{
		pg<-sapply(orfs,pgis,m,cFstats,dFstats,wilcoxon=wctest)
	}
	
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
	gis<-pg$gis
	# Put into data.frame
	nObs=length(p)
	if(wctest) {testType="wilcoxon"; sumType="median"}else{testType="t-test"; sumType="mean"}
	results<-data.frame(ORF=orfs,Gene=genes,P=p,Q=q,GIS=gis,
	QueryFitnessSummary=dFms,ControlFitnessSummary=cFms,QuerySE=dSe,ControlSE=cSe,QueryCount=dCount,ControlCount=cCount,QueryFit=dFitDef,ControlFit=cFitDef,
	TestType=rep(testType,nObs),SummaryType=rep(sumType,nObs),
	cTreat=rep(control$Treatment[1],nObs),cMed=rep(control$Medium[1],nObs),cScrID=rep(control$ScreenID[1],nObs),
	qTreat=rep(double$Treatment[1],nObs),qMed=rep(double$Medium[1],nObs),qScrID=rep(double$ScreenID[1],nObs),
	qGen=rep(double$Screen.Name[1],nObs),cGen=rep(control$Screen.Name[1],nObs),
	cLib=rep(paste(unique(control$Library.Name),collapse="_"),nObs),qLib=rep(paste(unique(double$Library.Name),collapse="_"),nObs))
	
	if("Client"%in%colnames(control)) results$cClient=rep(control$Client[1],nObs)
	if("Client"%in%colnames(double)) results$qClient=rep(double$Client[1],nObs)
	if("ExptDate"%in%colnames(control)) results$cDate=rep(gsub("'","",control$ExptDate[1]),nObs)
	if("ExptDate"%in%colnames(double)) results$qDate=rep(gsub("'","",double$ExptDate[1]),nObs)
	if("User"%in%colnames(control)) results$cUser=rep(control$User[1],nObs)
	if("User"%in%colnames(double)) results$qUser=rep(double$User[1],nObs)
	if("PI"%in%colnames(control)) results$cPI=rep(control$PI[1],nObs)
	if("PI"%in%colnames(double)) results$qPI=rep(double$PI[1],nObs)
	if("Condition"%in%colnames(control)) results$cCond=rep(control$Condition[1],nObs)
	if("Condition"%in%colnames(double)) results$qCond=rep(double$Condition[1],nObs)	
	if("Inoc"%in%colnames(control)) results$cInoc=rep(control$Inoc[1],nObs)
	if("Inoc"%in%colnames(double)) results$qInoc=rep(double$Inoc[1],nObs)

	#results$Type<-apply(results,1,typemake,m)
	results$Type=ifelse(results$GIS>0,"S","E") 
	results<-results[order(results$GIS,results$Q,results$Type),]
	# Get rid of duplicate entries in results
	orflist<-unique(as.character(results$ORF))
	results<-results[match(orflist,results$ORF),]
	# Plot results
	if (plot==TRUE){qfa.epiplot(results,qthresh,m,quantreg=quantreg)}
	final<-list(Results=results,
	Enhancers=gethits(results,qthresh,type="E",GISthresh=GISthresh),
	Suppressors=gethits(results,qthresh,type="S",GISthresh=GISthresh),
	slope=m,
	reg=reg)
	if(!is.null(bootstrap)) final$GISsummary=GISlist
	return(final)
}

############### Epistasis Functions ##################
# Makes epistasis plot for a given fdr level #
qfa.epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",yylab="Query Fitness",mmain="Epistasis Plot",fmax=0){
	enhancers<-gethits(results$Results,qthresh,type="E")
	suppressors<-gethits(results$Results,qthresh,type="S")
	others<-results$Results[(!results$Results$ORF%in%enhancers$ORF)&(!results$Results$ORF%in%suppressors$ORF),]
	if (fmax==0){
		# Get plot parameters
		ymax=1.1*max(results$Results$QueryFitnessSummary); ymin=0
		xmax=1.1*max(results$Results$ControlFitnessSummary); xmin=0
	}else{
		ymin=0;ymax=fmax
		xmin=0;xmax=fmax
	}
	plot(NULL,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
	xlab=xxlab,ylab=yylab,main=mmain,
	col=8,pch=19,cex=0.5)
	# Add line for genetic independence
	if (fitratio!=FALSE){abline(0,fitratio,lwd=2,col=8)} else {
		slope=results$slope
		abline(0,slope,lwd=2,col=8)}
	# Add 1:1 fitness line
	abline(0,1,lwd=2,lty=4,col=8)
	# Add reference ORF fitnesses lines
	if (ref.orf!=FALSE){
		reforf<-results$Results[results$Results$ORF==ref.orf,]
		abline(v=reforf$ControlFitnessSummary,col="lightblue",lwd=2)
		abline(h=reforf$QueryFitnessSummary,col="lightblue",lwd=2)}
	# Add points for non-suppressors & non-enhancers
	points(others$ControlFitnessSummary,others$QueryFitnessSummary,col="grey",cex=0.5,pch=19)
	text(others$ControlFitnessSummary,others$QueryFitnessSummary,others$Gene,col="grey",pos=4,offset=0.1,cex=0.4)
	# Add suppressors & enhancers
	if(length(enhancers$ORF)>0){
		points(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,col='green',pch=19,cex=0.5)
		text(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
	if(length(suppressors$ORF)>0){
		points(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,col='red',pch=19,cex=0.5)
		text(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
	Corr=signif(cor(results$Results$QueryFitnessSummary,results$Results$ControlFitnessSummary),3)
	Slope=signif(slope,3)
	legend("topleft",c(paste("Correlation: ",Corr),paste("Slope: ",Slope)),bty = "n")
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
	return(fitfunct(orfd$K,orfd$r,orfd$g,orfd$v))
}

# Regression through origin and point that splits dataset in two
splitRegression=function(x,y){
  # Find line that goes through one of the points (& origin) and has
  # about half of data above and about half of data below line
  fdat=data.frame(x=x,y=y)
  fdat=fdat[(fdat$x>0)&(fdat$y>0),]
  slopes=fdat$y/fdat$x
  pred=slopes%*%t(fdat$x)
  vert=t(pred)-fdat$y
  vert[abs(vert)<0.00000001]=0
  above=sign(vert)

  fracs=apply(above>0,1,sum)/(length(slopes)-1)
  best=which.min(abs(fracs-0.5))
  slopes[best]
}
# Regression minimising perpendicular distance between points and line
perpRegression=function(x,y){

  perpDist=function(x,y,m){
    sqrt(((x+m*y)/(m^2+1)-x)^2+(m*(x+m*y)/(m^2+1)-y)^2)
  }
  
  obf=function(m){
     return(sum(perpDist(x,y,m)))
  }

  perpRes=optim(50,obf,lower=0,upper=Inf)$par

  return(perpRes)
}

# Linear regression genetic ind. model fit
lm.epi<-function(doubles,controls,modcheck=FALSE,reg="splitreg"){
	if(reg=="lmreg"){
		indmod=lm(doubles~0+controls)
		m=as.numeric(indmod$coefficients)
	}
	if(reg=="quantreg"){
		indmod=rq(doubles~0+controls)
		m=as.numeric(indmod$coefficients)
	}
	if(reg=="splitreg"){
		m=splitRegression(controls,doubles)
	}
	if(reg=="perpreg"){
		m=perpRegression(controls,doubles)
	}
	
	# Check if linear regression OK
	if ((modcheck==TRUE)&(reg%in%c("quantreg","lmreg"))){
		pdf("ModelCheck.pdf")
		resids<-indmod$residuals
		hist(resids,xlab="Fitness Residuals",main="Histogram of Residuals")
		qqnorm(resids/sd(resids),main="Normal QQ Plot")
		abline(0,1,lwd=2,col="red")
		dev.off()
	}
	return(m)
}

# Estimates p-value and estimated strength of interaction
pgis<-function(orf,m,cFs,dFs,wilcoxon=TRUE){
	# If this orf is not present in both lists, return appropriate p,gis
	if((length(dFs[[orf]])==0)|(length(cFs[[orf]])==0)){return(c(1,0))}
	# If fitnesses are not unique in one condition (e.g. all dead) and only one repeat in another (e.g. after stripping)
	# This would cause wilcoxon test to fail due to ties, and t.test to fail due to insufficient y?
	ldFS=length(dFs[[orf]]); lcFS=length(cFs[[orf]])
	ludFS=length(unique(dFs[[orf]])); lucFS=length(unique(cFs[[orf]]))
	if (((ludFS==1)&(ldFS>1)&(lcFS==1))|((lucFS==1)&(lcFS>1)&(ldFS==1))){return(c(1,median(dFs[[orf]],na.rm=TRUE)-m*median(cFs[[orf]],na.rm=TRUE)))}
	if ((ludFS==1)&(lucFS==1)){return(c(1,median(dFs[[orf]],na.rm=TRUE)-m*median(cFs[[orf]],na.rm=TRUE)))}
	# If both sets of cultures are dead (and fitnesses therefore equal) return appropriate p,gis
	if(sum(dFs[[orf]]==m*cFs[[orf]])==length(dFs[[orf]]==m*cFs[[orf]])){
		return(c(1,0))
	}else{
		if (wilcoxon){
			valdiff=0
			p=1
			# Returns p-value for significance of difference, and estimate of difference between medians
			tryCatch({
				ctest<-wilcox.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
				valdiff<-as.numeric(ctest$estimate)
				p<-as.numeric(ctest$p.value)
			},error=function(e) {
				print(paste("Wilcox test failed for",orf,", using t-test instead for this gene"))
				ctest<-t.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
				valdiff<-as.numeric(ctest$estimate)
				valdiff<<-valdiff[1]-valdiff[2]
				p<<-as.numeric(ctest$p.value)
			})
			
			return(c(p,valdiff))
		}else{
			# t-test fails if only one element in either list, do one sample test
			if((ldFS<=1)|(lcFS<=1)){
				ctest<-t.test(dFs[[orf]]-m*cFs[[orf]])
				p<-as.numeric(ctest$p.value)
				valdiff<-as.numeric(ctest$estimate)
				return(c(p,valdiff))
			}else{
				# Returns p-value for significance of difference, and the difference between the means
				ctest<-t.test(dFs[[orf]],m*cFs[[orf]],alternative="two.sided",conf.int=TRUE)
				p<-as.numeric(ctest$p.value)
				valdiff<-as.numeric(ctest$estimate)
				valdiff<-valdiff[1]-valdiff[2]
				return(c(p,valdiff))
			}
		}
	}
}

# Estimates p-value and estimated strength of interaction
pgis_bootstrap<-function(orf,cFs,dFs,sampSumm=mean,wt="YOR202W",Nreps=10000,subSamp=Inf){
	# Need to come up with a good estimate for a minumum non-zero fitness
	# in order to avoid problems with division by zero later
	fitMin=median(c(unlist(cFs,use.names=FALSE),unlist(dFs,use.names=FALSE)))/200

	# If this orf is not present in both lists, return appropriate p,gis
	if((length(dFs[[orf]])==0)|(length(cFs[[orf]])==0)){return(c(1,0))}
	
	# Bootstrap estimates of fitness summary distribution for relevant genotypes
	obsDoubleMut=bsSamp(dFs[[orf]],Nreps,sampSumm,subSamp)
	arrayMut=bsSamp(cFs[[orf]],Nreps,sampSumm,subSamp)
	backMut=bsSamp(dFs[[wt]],Nreps,sampSumm,subSamp)
	wtMut=bsSamp(cFs[[wt]],Nreps,sampSumm,subSamp)
	
	# Uncertainty about summary of predicted double mutant fitness
	predDoubleMut=backMut*arrayMut/pmax(wtMut,fitMin)
	GIS=obsDoubleMut-predDoubleMut
	return(GIS)
}

bsSamp=function(A,Nrep=100000,sampSumm=mean,subSamp=Inf){
	bsreps=replicate(Nrep,sampSumm(sample(as.numeric(A),min(subSamp,length(A)),replace=TRUE)))
	return(bsreps)
}

# Get type of interaction (DEPRECATED)
typemake<-function(row,m){
	summd<-as.numeric(row['QueryFitnessSummary'])
	summc<-as.numeric(row['ControlFitnessSummary'])
	if (m<1){
		if (summd>m*summc){type<-"S"} else {type<-"E"}
	} 
	else {
		if (summd>m*summc){type<-"E"} else {type<-"S"}
	}
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

##### Does max. lik. fit for all colonies, given colonyzer.read or rod.read input #####
qfa.fit<-function(d,inocguess,ORF2gene="ORF2GENE.txt",fmt="%Y-%m-%d_%H-%M-%S",minK=0.025,detectThresh=0.0005,globalOpt=FALSE,logTransform=FALSE,fixG=TRUE,AUCLim=5,STP=20,nCores=1,glog=TRUE,modelFit=TRUE,checkSlow=TRUE,...){
	
	if(!"Column"%in%colnames(d)) d$Column=d$Col

	# ORF2gene argument is now deprecated, this link should be made with colony.read function instead
	if(!is.null(inocguess)){
		if(length(inocguess)==0) {
			print("ERROR: must specify an inoculum density guess (or NULL)")
			return()
		}
	}

	if(nCores>1){
		cl=makeCluster(nCores)
		clusterCall(cl,function() library(qfa))
	}else{cl=NULL}
	
	# Rename columns if necessary
	if(!"Tile.Dimensions.X"%in%colnames(d)) d[,"Tile.Dimensions.X"]=d$Diameter
	if(!"Tile.Dimensions.Y"%in%colnames(d)) d[,"Tile.Dimensions.Y"]=d$Diameter
	if(!"X.Offset"%in%colnames(d)) d[,"X.Offset"]=round(d$x-d$Diameter/2.0)
	if(!"Y.Offset"%in%colnames(d)) d[,"Y.Offset"]=round(d$y-d$Diameter/2.0)
	if(!"Edge.length"%in%colnames(d)) d[,"Edge.length"]=d$Perimeter
	if(!"Edge.Pixels"%in%colnames(d)) d[,"Edge.Pixels"]=d$Perimeter
	
	# Vector of barcodes
	barcodes<-unique(d$Barcode); nbc<-length(barcodes)
	# Get big data frame ready for results
	results<-data.frame()
	# For each barcode, optimize
	bcount<-0
	for (bcode in barcodes){
		bcount<-bcount+1
		print(paste("Optimizing Plate",bcount,"/",nbc,":",bcode))
		# Restrict data to just this barcode
		dbc<-d[d$Barcode==bcode,]
		inoctime=dbc$Inoc.Time[1]
		# Get list of unique colony positions
		positions<-lapply(1:length(dbc[,1]),index2pos,dbc)
		positions<-unique(positions)
		# Fit logistic model to each colony
		if(nCores>1){
			print(as.vector(lsf.str(envir=.GlobalEnv)))
			ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv))
			clusterExport(cl, ex)
			clusterExport(cl, as.vector(lsf.str(envir=.GlobalEnv)))
			bcfit<-t(parSapply(cl,positions,colony.fit,dbc,inocguess,fixG,globalOpt,detectThresh,minK,logTransform,AUCLim,STP,glog,modelFit,checkSlow,...))
		}else{
			bcfit<-t(sapply(positions,colony.fit,dbc,inocguess,fixG,globalOpt,detectThresh,minK,logTransform,AUCLim,STP,glog,modelFit,checkSlow,...))
		}
		info<-data.frame(t(sapply(positions,colony.info,dbc)))
		rows<-sapply(positions,rcget,"row")
		cols<-sapply(positions,rcget,"col")
		# Bind Data frame of barcode results to overall results
		barcMetadata<-data.frame(
		Barcode=as.character(info$Barcode),
		Row=rows,
		Column=cols,
		Col=cols,
		ScreenID=as.character(info$ScreenID),
		Treatment=as.character(info$Treatments),
		Medium=as.character(info$Medium),
		ORF=as.character(info$ORF),
		Screen.Name=as.character(info$Screen.Name),
		Library.Name=as.character(info$Library.Name),
		MasterPlate.Number=as.numeric(info$MasterPlate.Number),
		Timeseries.order=as.numeric(info$Timeseries.order),
		Inoc.Time=inoctime,
		TileX=as.numeric(info$Tile.Dimensions.X),
		TileY=as.numeric(info$Tile.Dimensions.Y),
		XOffset=as.numeric(info$X.Offset),
		YOffset=as.numeric(info$Y.Offset),
		Threshold=as.numeric(info$Threshold),
		EdgeLength=as.numeric(info$Edge.length),
		EdgePixels=as.numeric(info$Edge.Pixels),
		RepQuad=as.numeric(info$RepQuad),
		stringsAsFactors=FALSE)
		barcFitness<-data.frame(bcfit)
		barcResults=cbind(barcMetadata,barcFitness)
		if("Client"%in%colnames(info)){barcResults$Client=as.character(info$Client)}
		if("ExptDate"%in%colnames(info)){barcResults$ExptDate=as.character(info$ExptDate)}
		if("User"%in%colnames(info)){barcResults$User=as.character(info$User)}
		if("PI"%in%colnames(info)){barcResults$PI=as.character(info$PI)}
		if("Condition"%in%colnames(info)){barcResults$Condition=as.character(info$Condition);barcResults$Condition[is.na(barcResults$Condition)]=""}
		if("Inoc"%in%colnames(info)){barcResults$Inoc=as.character(info$Inoc)}
		if("Gene"%in%colnames(info)){barcResults$Gene=as.character(info$Gene)}
		results=rbind(results,barcResults)
		} #bcode
		if(nCores>1){
			stopCluster(cl)
		}
		return(results)
}

makeFitness<-function(results,AUCLim=5,dtmax=25){
	# Fitness definitions from Addinall et al. 2011
	# Update the inoculum density parameter in the case where tshift!=0
	#if("tshift"%in%names(results))	results$g=Glogist(results$K,results$r,results$g,results$v,0-results$tshift)
	results$MDP=mdp(results$K,results$r,results$g,results$v)
	results$MDR=mdr(results$K,results$r,results$g,results$v)
	results$MDRMDP=mdrmdp(results$K,results$r,results$g,results$v)
	results$glog_maxslp=gen_logistic_maxslp(results$K,results$r,results$g,results$v)
	# Set spurious fitnesses to zero
	results$MDR[(results$r>7)&(results$K<0.0275)]=0
	results$MDRMDP[(results$r>7)&(results$K<0.0275)]=0
	# Doubling time (doublings per hour)
	if("t0"%in%rownames(results)){
		results$DT=dtl(results$K,results$r,results$g,results$v,results$t0)*24
	}else{
		results$DT=(1/results$MDR)*24
	}
	# If doubling time is Inf, set to dtmax
	results$DT[abs(results$DT)>dtmax]=dtmax
	# Area under curve
	Glog=function(K,r,g,v,t) {
		# Eliminate problems with division by zero
		g=max(g,1E-9) 
		K=max(K,1E-9)
		ifelse(is.na(Glogist(K,r,g,v,t)),0,Glogist(K,r,g,v,t)) # Temporarily replace NA with zero here to allow integrate function below to handle modelFit=FALSE
	}
	AUC<-function(dno,tstar,dat) max(0,integrate(Glog,lower=0,upper=tstar,K=dat$K[dno],r=dat$r[dno],g=dat$g[dno],v=dat$v[dno],subdivisions=1000)$value - tstar*dat$g[dno])
	results$AUC=sapply(1:length(results[,1]),AUC,tstar=AUCLim,dat=results)
	return(results)
}

# Extract row & col from list of position vectors
rcget<-function(posvec,rc) posvec[match(rc,c("row","col"))] 

# Closure returning a loess-based approximating function for a timeseries
loapproxfun=function(t,g,span=0.2){
	lo=loess(g~t,span=span)
	# If loess smoothing with specified span does not work (e.g. too few points), revert to linear interpolation
	if(is.na(lo$fitted[1])){
		loex=approxfun(t,g,method="linear",rule=2)
	}else{
		loex=function(t){
			cdens=predict(lo,t)
			cdens[t<min(lo$x)]=predict(lo,min(lo$x))
			cdens[t>max(lo$x)]=predict(lo,max(lo$x))
			return(cdens)
		}
	}
	return(loex)
}

numericalfitness<-function(obsdat,AUCLim,STP){
	# Generate numerical AUC
	if(length(obsdat$Growth)>1){
			loapproxfree=loapproxfun(obsdat$Expt.Time,obsdat$Growth,span=0.5)
			loapprox=function(x) pmax(0,loapproxfree(x))
			nAUCfn=function(AUCLim) as.numeric(integrate(loapprox,0,AUCLim)$value)
			nSTPfn=function(STP) as.numeric(loapprox(STP))
			nAUC=sapply(AUCLim,nAUCfn)
			nSTP=sapply(STP,nSTPfn)
	}else{
			# If there's only one photograph (e.g. single time point 1536 assay)
			nAUC=rep(NA,length(AUCLim))
			nSTP=rep(obsdat$Growth[1],length(STP))
	}
	res=c(nAUC,nSTP)
	if(length(AUCLim)==1) {nAUCnames=c("nAUC")}else{nAUCnames=paste("nAUC",sprintf("%04d",round(AUCLim*24*60)),sep="")}
	if(length(STP)==1) {nSTPnames=c("nSTP")}else{nSTPnames=paste("nSTP",sprintf("%04d",round(STP*24*60)),sep="")}
	names(res)=c(nAUCnames,nSTPnames)
	if(length(AUCLim)>1) res["nAUC"]=res[nAUCnames[length(nAUCnames)]]
	if(length(STP)>1) res["nSTP"]=res[nSTPnames[length(nSTPnames)]]
	numr_lst=numerical_r(obsdat)
	res["nr"]=numr_lst$nr
	res["nr_t"]=numr_lst$nr_t
	res["maxslp"]=numr_lst$mslp
	res["maxslp_t"]=numr_lst$mslp_t	
	return(res)
}

numerical_r=function(obsdat,mkPlots=FALSE,span=0.3,nBrute=1000,cDiffDelta=0.0001,mlab=""){
	# Generate numerical (model-free) estimate for nr_t intrinsic growth rate
	tims=obsdat$Expt.Time
	gdat=obsdat$Growth
	tmax=max(tims)

	# Smooth data, find slope as function of time and nr_t slope
	lgdat=log(gdat)
	ltims=tims[!is.na(lgdat)]
	lgdat=lgdat[!is.na(lgdat)]
	la=NA
	try(a<-loapproxfun(tims,gdat,span=span),silent=TRUE)
	try(la<-loapproxfun(ltims,lgdat,span=span),silent=TRUE)
	if(!is.function(la)|!is.function(a)) return(list(nr=0,nr_t=NA,mslp=0,mslp_t=NA))
	centralDiff=function(f,delta) return(function(x) (f(x+delta/2.0)-f(x-delta/2.0))/delta)
	lslp=centralDiff(la,cDiffDelta)
	slp=centralDiff(a,cDiffDelta)
	# Brute force optimization
	stimes=seq(min(ltims),max(ltims),length.out=nBrute)
	vals=a(stimes)
	slps=slp(stimes)
	lvals=la(stimes)
	lslps=lslp(stimes)
	# Discard points too close to t=0 to avoid artificially high slopes
	opt=which.max(slps)
	lopt=which.max(lslps)
	res=list(nr=lslps[lopt],nr_t=stimes[lopt],mslp=slps[opt],mslp_t=stimes[opt])
	#maxlslp=optimize(lslp,lower=min(ltims),upper=max(ltims),nr_t=TRUE)
	#res$nr=maxlslp$objective
	#res$nr_t=maxlslp$maximum
	maxslope=res$nr
	
	if(mkPlots){
		mainlab=paste("Max. intrinsic growth rate =",formatC(res$nr,3),"(estimated on log scale)",mlab)
		# Plot synthetic data, Loess approximation and estimated slope
		op=par(mar=c(5,4,4,5)+.1,mfrow=c(1,2))
		plot(NULL,xlab="Time (d)",ylab="",xlim=c(-0.1*tmax,tmax),ylim=range(lgdat),main=mainlab)
		abline(v=res$nr_t,lwd=3,col="green",lty=2)
		slope=res$nr
		intercept=la(res$nr_t)-slope*res$nr_t
		abline(a=intercept,b=slope,lwd=3,col="green")
		points(ltims,lgdat)
		mtext("Log Cell density (AU)",side=2,line=3,col="red")
		curve(la,from=-0.1*tmax,to=tmax,add=TRUE,col="red",lwd=2,xlim=c(-0.1*tmax,tmax))
		par(new=TRUE)
		
		curve(lslp(x),from=-0.1*tmax,to=tmax,col="blue",xlab="",ylab="",xaxt="n",yaxt="n",lwd=2,xlim=c(-0.1*tmax,tmax))
		axis(4)
		mtext("Slope of Log Cell Density",side=4,line=3,col="blue")

		mainlab=paste("Max. slope =",formatC(res$mslp,3),"(estimated on linear scale)",mlab)
		# Plot synthetic data, Loess approximation and estimated slope
		op=par(mar=c(5,4,4,5)+.1)
		plot(NULL,xlab="Time (d)",ylab="",xlim=c(-0.1*tmax,tmax),ylim=range(obsdat$Growth),main=mainlab)
		abline(v=res$mslp_t,lwd=3,col="green",lty=2)
		slope=slp(res$mslp_t)
		intercept=a(res$mslp_t)-slope*res$mslp_t
		abline(a=intercept,b=slope,lwd=3,col="green",untf=FALSE)
		points(obsdat$Expt.Time,obsdat$Growth)
		mtext("Cell density (AU)",side=2,line=3,col="red")
		curve(a(x),from=-0.1*tmax,to=tmax,add=TRUE,col="red",lwd=2,xlim=c(-0.1*tmax,tmax))
		par(new=TRUE)
		
		curve(slp(x),from=-0.1*tmax,to=tmax,col="blue",xlab="",ylab="",xaxt="n",yaxt="n",lwd=2,xlim=c(-0.1*tmax,tmax))
		axis(4)
		mtext("Slope of Cell Density",side=4,line=3,col="blue")

		par(op)
	}
	return(res)
}

### Growth model fitting for one colony ###
colony.fit<-function(position,bcdata,inocguess,fixG=TRUE,globalOpt=FALSE,detectThresh=0,minK=0,logTransform=FALSE,AUCLim=5,STP=10,glog=TRUE,modelFit=TRUE,checkSlow=TRUE,...){
	# Get row & column to restrict data
	row<-position[1]; col<-position[2]
	#print(paste(bcdata$Barcode[1],"Row:",row,"Col:",col))
	do<-bcdata[(bcdata$Row==row)&(bcdata$Column==col),]
	obsdat=data.frame(Expt.Time=as.numeric(do$Expt.Time),Growth=as.numeric(do$Growth))
	pars=makefits(obsdat,inocguess,fixG,globalOpt,detectThresh,minK,logTransform,AUCLim,STP,glog,modelFit,checkSlow)
	return(pars)
}

makefits<-function(obsdat,inocguess,fixG=TRUE,globalOpt=FALSE,detectThresh=0,minK=0,logTransform=FALSE,AUCLim=5,STP=10,glog=TRUE,modelFit=TRUE,checkSlow=TRUE,...){
	numfit=numericalfitness(obsdat,AUCLim,STP)
	
	if(modelFit){
		pars=growthcurve(obsdat,inocguess,fixG,globalOpt,detectThresh,minK,logTransform,glog,checkSlow)
	}else{
		pars=c(max(obsdat$Growth),NA,NA,NA,NA,NA)
		names(pars)=c("K","r","g","v","objval","tshift")
	}
	# Add on time of first valid obs. and numerical AUC, STP
	valid=obsdat[obsdat$Growth>=detectThresh,]
	if(dim(valid)[1]>0) {
		t0=min(valid$Expt.Time)
		# Calculate mean here in case multiple data for same timepoint...
		d0=mean(valid$Growth[valid$Expt.Time==t0]) 
	}else{
		t0=Inf
		d0=NA
	}
	parsadd=c(t0,d0)
	names(parsadd)=c("t0","d0")
	pars=c(pars,parsadd,numfit)
	#if(length(pars)!=9)	{dput(pars); dput(obsdat)}
	return(pars)
}

makeBoundsQFA<-function(inocguess,d,minK=0,fixG=FALSE,globalOpt=FALSE,glog=TRUE){
	# Define optimization bounds based on inocguess #
	#lowr<-0; upr<-25
	lowr<-0; upr<-500
	if(glog){
		#lowv<-0.1; upv<-10.0
		lowv<-0.25; upv<-4.0
		lowv<-0.000001
	}else{
		lowv<-1; upv<-1
	}
	lowK<-max(0.9*inocguess,minK); upK<-1.0
		
	# We often fix inoculation density, but users might prefer to infer it from growth curves
	if(fixG) {lowg<-0.9*inocguess; upg<-1.1*inocguess}else{lowg<-0.01*inocguess; upg<-100.0*inocguess}
	if(globalOpt) {lowg<-1e-3*inocguess; upg<-1e3*inocguess}
	
	lowg<-max(0,lowg); upg<-min(upK,upg) 
	
	return(list(K=c(lowK,upK),r=c(lowr,upr),g=c(lowg,upg),v=c(lowv,upv),inocguess=inocguess))
}

growthcurve<-function(obsdat,iguess,fixG=TRUE,globalOpt=FALSE,detectThresh=0,minK=0,logTransform=FALSE,glog=TRUE,checkSlow=TRUE){
	len1=length(obsdat$Growth)
	# Throw away observations below the detectable threshold
	d=obsdat[obsdat$Growth>=detectThresh,]
	# Throw away observations occurring before inoculation date (negative times...)
	d=d[d$Expt.Time>=0,]
	d=d[!is.na(d$Growth),]
	len2=length(d$Growth)
	# If this has left us with too few points, return "dead colony"
	if ((len2/len1<0.25)|(len2<3)) {
		if(!is.null(iguess)) {pars=c(iguess,0,iguess,1,Inf,0)}else{pars=c(minK,0,minK,1,Inf,0)}
		names(pars)=c("K","r","g","v","objval","tshift")
		return(pars)
	}
	
	# Deal with situation where no inoculum guess specified
	tshift=0
	if(is.null(iguess)){
		# Without a sensible independent estimate for inoculum density, the best we can do is to estimate it based on observed data.
		# This strategy will only work well if the inoculum density is high enough to be measurable (e.g. pinned cultures or 
		# conc. spotted) and is clearly observed.  Clearly observed means: no condensation on plates immediately after they are 
		# placed in incubator for example.
		# If we are making an independent estimate of inoculum density, then we should also reset the time at which the experiment "begins".
		# This experiment start time should be the time at which the inoculum density is observed.
		if(length(d$Growth)>0) {
			candidate=min(d$Growth)
			tshift=d$Expt.Time[match(candidate,d$Growth)]
		}else{candidate=0}
		iguess=max(0.0001,candidate)
	}
	d$Expt.Time=d$Expt.Time-tshift
	# In the unlikely event that the smallest cell density did not occur at the earliest timepoint, need to delete earlier timepoints
	d=d[d$Expt.Time>=0,]
	
	bounds=makeBoundsQFA(iguess,d,minK,fixG,globalOpt,glog)
	inocguess=bounds$inocguess
	bounds$inocguess=NULL
	xybounds=bounds	

	# First *detectable* observation at time t0
	t0=min(d$Expt.Time)
	#xybounds$K=c(0.9*max(d$Growth),1.1*max(d$Growth))
	if(globalOpt) {
		pars=de.fit(d$Expt.Time,d$Growth,inocguess,xybounds,initPop=TRUE,logTransform=logTransform)
		# Check for high fraction of K at end of modelled experiment
		GEnd=Glogist(pars[1],pars[2],pars[3],pars[4],max(d$Expt.Time))
		if(GEnd/pars[1]<0.75){ # If experiment not quite finished...
			#print("Modelled growth at end of experiment is less than 0.75 of K estimate...")
			# Put tight bounds on K and optimise again
			Kmin=max(0.95*1.5*GEnd,minK); Kmax=max(1.05*1.5*GEnd,minK)
			xybounds$K=c(Kmin,Kmax)
			pars=de.fit(d$Expt.Time,d$Growth,inocguess,xybounds,initPop=TRUE,logTransform=logTransform)
		}
	}else{
		pars=data.fit(d$Expt.Time,d$Growth,inocguess,xybounds,logTransform=logTransform)
		# Check for high fraction of K at end of modelled experiment
		GEnd=Glogist(pars[1],pars[2],pars[3],pars[4],max(d$Expt.Time))
		if(GEnd/pars[1]<0.75){ # If experiment not quite finished...
			#print("Modelled growth at end of experiment is less than 0.75 of K estimate...")
			# Put tight bounds on K and optimise again
			Kmin=max(0.95*1.5*GEnd,minK); Kmax=max(1.05*1.5*GEnd,minK)
			xybounds$K=c(Kmin,Kmax)
			pars=data.fit(d$Expt.Time,d$Growth,inocguess,xybounds,logTransform=logTransform)
		}
	}
	# Check for spurious combination of relatively high r, low K, spend more time optimising...
	if(checkSlow&(mdr(pars[1],pars[2],pars[3],pars[4])>1.0)&((pars[1]<0.05)|(max(d$Growth)<0.05)|(tail(d$Growth,1)<0.05))){ # Try optimising with sick colony as guess
		 #print("Attempting to do global fit alternative sick colony growth curve")
		 Kmin=max(0.9*inocguess,minK); Kmax=1.5*max(pars[1],minK); xybounds$K=c(Kmin,Kmax); 
		 xybounds$r=c(0,3) # Slow growth
		 if(glog){
			xybounds$v=c(0.75,1.5) # More logistic growth
	     }else{
		    xybounds$v=c(1,1)
		 }
		 inits=list(K=pars[1],r=0.6,g=inocguess,v=1)
		 newpars=de.fit(d$Expt.Time,d$Growth,inocguess,xybounds,inits=inits,initPop=TRUE,widenr=FALSE,logTransform=logTransform)			
		 if(newpars[5]<=pars[5]) pars=newpars				
	 }
	opt=sumsq(pars[1],pars[2],pars[3],pars[4],obsdat$Growth,obsdat$Expt.Time,logTransform=logTransform) # Use all data (not just data below detection thresh)
	dead=sumsq(inocguess,0,inocguess,1,obsdat$Growth,obsdat$Expt.Time,logTransform=logTransform)
	if(dead<=opt) {
		pars=c(inocguess,0,inocguess,1,dead)
		names(pars)=c("K","r","g","v","objval") # Try dead colony
	}
	if(!is.finite(pars[["K"]])) pars[["K"]]=minK
	if(!is.finite(pars[["r"]])) pars[["r"]]=0
	if(!is.finite(pars[["g"]])) pars[["g"]]=0
	if("v"%in%names(pars)) if(!is.finite(pars[["v"]])) pars[["v"]]=1
	pars[["tshift"]]=tshift
	#if(length(pars)!=5) dput(pars)
	return(pars)
}

### Function that fits model to a timecourse with genetic optimisation algorithm
de.fit<-function(tim,growth,inocguess,xybounds,inits=list(),initPop=FALSE,widenr=TRUE,logTransform=FALSE,mxit=2000){
	# Fit to growth curve with differential evolution
	# Get initial guess for parameters
	if(length(inits)==0){
		# Get initial guess for parameters
		init<-guess(tim,growth,inocguess,xybounds)
	}else{init=inits}
	if (widenr) xybounds$r=c(0.01*init$r,10.0*init$r)

	# Function to be optimized
	objf<-function(modpars){
		K<-modpars[1]; r<-modpars[2]
		g<-modpars[3]; v<-abs(modpars[4])
		res=sumsq(K,r,g,v,growth,tim,logTransform=logTransform)
		return(res)
	}

	# Make results repeatable
	#set.seed(1234)
	# Format bounds
	low=c(xybounds$K[1],xybounds$r[1],xybounds$g[1],xybounds$v[1])
	up=c(xybounds$K[2],xybounds$r[2],xybounds$g[2],xybounds$v[2])

	NumParticles=11*length(low)

	if(initPop){
	# Randomly sample from within bounds for initial population
	Klist=runif(NumParticles,min=low[1],max=up[1])
	rlist=runif(NumParticles,min=low[2],max=up[2])
	glist=runif(NumParticles,min=low[3],max=up[3])
	vlist=runif(NumParticles,min=low[4],max=up[4])
	pop=data.frame(K=Klist,r=rlist,g=glist,v=vlist)
	pop=as.matrix(pop)
	# Set first particle equal to initial guess
	pop[1,]=as.numeric(t(init))
	# Set second particle equal to a dead colony
	pop[2,]=c(0.025,0,low[3],1)
	}else{pop=NULL}

	optsol=DEoptim(objf,lower=low,upper=up,
	DEoptim.control(trace = FALSE,NP=NumParticles,initialpop=pop,itermax=mxit)
	)
	pars=abs(as.numeric(optsol$optim$bestmem))
	objval=objf(pars)
	# Sanity check for fitted parameters (no negative growth)
	if (pars[1]<pars[3]){
		pars[1]=pars[3]; pars[2]=0; pars[4]=1
		objval=objf(pars)}
	pars=c(pars,objval)
	names(pars)=c("K","r","g","v","objval")
	return(pars)
}

### Function that fits to for a timecourse by maximum likelihood (least squares)
data.fit<-function(tim,growth,inocguess,xybounds,inits=list(),logTransform=FALSE,verbose=FALSE){
	if(length(inits)==0){
		# Get initial guess for parameters
		init<-guess(tim,growth,inocguess,xybounds)
		if(xybounds$v[1]==xybounds$v[2]) init$v=NULL
		if(xybounds$g[1]==xybounds$g[2]) init$g=NULL
	}else{init=inits}
	#xybounds$r=c(0.75*init$r,1.25*init$r)
	#xybounds$r=c(0.1*init$r,10.0*init$r)
	# Set initial guess to twice best estimate to bias against small r
	#init$r=2.0*init$r
	# Try to stop L-BFGS-B errors by moving away from minK
	xybounds$K[1]=1.01*xybounds$K[1]
	# Function to be optimized
	# Logistic only will have same upper and lower bounds for v...
	paramscale=c(0.2,10,inocguess,1)
	names(paramscale)=c("K","r","g","v")
	if(xybounds$v[1]==xybounds$v[2]){
		if(xybounds$g[1]==xybounds$g[2]){
			objf<-function(modpars){
				return(sumsq(modpars[["K"]],modpars[["r"]],xybounds$g[1],xybounds$v[1],growth,tim,logTransform=logTransform))
			}
			lbounds=c(xybounds$K[1],xybounds$r[1])
			ubounds=c(xybounds$K[2],xybounds$r[2])
			pscl=as.numeric(paramscale[c("K","r")])	
		}else{
			objf<-function(modpars){
				return(sumsq(modpars[["K"]],modpars[["r"]],modpars[["g"]],xybounds$v[1],growth,tim,logTransform=logTransform))
			}
			lbounds=c(xybounds$K[1],xybounds$r[1],xybounds$g[1])
			ubounds=c(xybounds$K[2],xybounds$r[2],xybounds$g[2])
			pscl=as.numeric(paramscale[c("K","r","g")])	
		}
	}else{
		if(xybounds$g[1]==xybounds$g[2]){
			objf<-function(modpars){
				return(sumsq(modpars[["K"]],modpars[["r"]],xybounds$g[1],modpars[["v"]],growth,tim,logTransform=logTransform))
			}
			lbounds=c(xybounds$K[1],xybounds$r[1],xybounds$v[1])
			ubounds=c(xybounds$K[2],xybounds$r[2],xybounds$v[2])
			pscl=as.numeric(paramscale[c("K","r","v")])				
		}else{
			objf<-function(modpars){
				return(sumsq(modpars[["K"]],modpars[["r"]],modpars[["g"]],modpars[["v"]],growth,tim,logTransform=logTransform))
			}
			lbounds=c(xybounds$K[1],xybounds$r[1],xybounds$g[1],xybounds$v[1])
			ubounds=c(xybounds$K[2],xybounds$r[2],xybounds$g[2],xybounds$v[2])
			pscl=as.numeric(paramscale[c("K","r","g","v")])	
		}
	}
	# Perform optimization
	optsol<-optim(par=unlist(init),fn=objf,gr=NULL,method="L-BFGS-B",
	lower=lbounds,
	upper=ubounds,
	control=list(maxit=1000,factr=1e7,trace=0,parscale=pscl))
	pars=abs(optsol$par)
	objval=objf(pars)
	if(!"g"%in%names(pars)) pars[["g"]]=xybounds$g[1]
	if(!"v"%in%names(pars)) pars[["v"]]=xybounds$v[1]
	# Sanity check for fitted parameters (no negative growth)
	if (pars[["K"]]<pars[["g"]]){
		pars[["K"]]=pars[["g"]]; pars[["r"]]=0; pars[["v"]]=1
		objval=objf(pars)
	}
	pars[["objval"]]=objval
	pars=pars[c("K","r","g","v","objval")]
	if(verbose) {if(optsol$message!="CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH") print(optsol$message)}
	return(pars)
}

### Provide initial guess for logistic model parameters ###
guess<-function(tim,growth,inocguess,xybounds,minK=0.025){
	# Sort time and growth
	growth=growth[order(tim)]
	tim=tim[order(tim)]
	n=sum((!is.na(tim))&(!is.na(growth)))

	# Enforce positivity and monotonic increasing behaviour in growth
	#growth[1]=max(c(growth[1],0.000000001))
	#for (x in 2:length(growth)) growth[x]=max(c(max(growth[1:(x-1)]),growth[x],0.00000001))
	if(is.null(inocguess)){G0g<-max(min(growth,na.rm=TRUE),0.0000001)}else{G0g<-inocguess}
	Kg<-max(max(growth,na.rm=TRUE),minK)
	vg=1 # Assume logistic model is adequate
	rg=0
	if(n>3){
		# Guess for r is a bit more complicated
		targ<-min(growth)+(max(growth)-min(growth))/2.0
		approxspline=smooth.spline(tim,growth)
		approxcurve=approxfun(approxspline$x,approxspline$y,rule=2)
		solvefn<-function(x) approxcurve(x)-targ
		if((approxcurve(max(tim))>=approxcurve(min(tim)))&(sign(solvefn(min(tim)))!=sign(solvefn(tim[which.max(growth)])))){
			sol=uniroot(solvefn,c(min(tim),tim[which.max(growth)]))
			tmrate=sol$root
		}else{ # Too few points to fit spline curve, still do something sensible...
			tmrate=(min(tim)+max(tim))/2.0;
		}
		if(Kg>G0g) rg=log((Kg-G0g)/G0g)/tmrate
	}
	# Sanity check for guessed parameter values
	# If the data have low correlation, then set r=0
	# If all elements of growth are equal (e.g. zero) then correlation function throws error...
	if((length(unique(growth))==1)|(cor(tim,growth)<0.1)){ rg=0; Kg=G0g }
	return(list(K=Kg,r=rg,g=G0g,v=vg))
}#guess

# Sum of squared error
sumsq<-function(K,r,g,v,growth,tim,logTransform=FALSE){
	# For generalised logistic function, K/g must be positive to avoid complex cell density estimates
	if((K/g)<0) return(Inf)
	if(logTransform){
		glog=growth[growth>0]; tlog=tim[growth>0]
		ss=sqrt(sum((log(glog)-log(Glogist(K,r,g,v,tlog)))^2))/length(glog)
	}else{
		ss=sqrt(sum((growth-Glogist(K,r,g,v,tim))^2))/length(growth)
	}
	if(is.na(ss)){
		return(Inf)
	}else{return(ss)}
}

# Prevent zero or negative growth from occurring
nozero<-function(growth){if (growth<=0){growth<-0.0000000000001} else {growth<-growth}}

#### Get colony information ####
colony.info<-function(position,bcdata){
	# Get row & column to restrict data
	row<-position[1]; col<-position[2]
	d<-bcdata[(bcdata$Row==row)&(bcdata$Column==col),]
	# Want to use the LAST datapoint 
	# Most relevant for edge length (possibly most reliable for position)
	d<-d[order(d$Date.Time,decreasing=TRUE),]
	d<-d[1,]
	# Character vector of colony info
	cvec=list(
	Barcode=as.character(d$Barcode),
	Treatments=as.character(d$Treatments),
	Medium=as.character(d$Medium),
	ORF=as.character(d$ORF),
	Screen.Name=as.character(d$Screen.Name),
	Library.Name=as.character(d$Library.Name),
	MasterPlate.Number=as.numeric(d$MasterPlate.Number),
	Timeseries.order=as.numeric(d$Timeseries.order),
	ScreenID=as.character(d$ScreenID),
	Tile.Dimensions.X=as.numeric(d$Tile.Dimensions.X),
	Tile.Dimensions.Y=as.numeric(d$Tile.Dimensions.Y),
	X.Offset=as.numeric(d$X.Offset),
	Y.Offset=as.numeric(d$Y.Offset),
	Threshold=as.numeric(d$Threshold),
	Edge.length=as.numeric(d$Edge.length),
	Edge.Pixels=as.numeric(d$Edge.Pixels),
	RepQuad=as.numeric(d$RepQuad)
	)
	if("Client"%in%colnames(d)){cvec$Client=as.character(d$Client)}
	if("ExptDate"%in%colnames(d)){cvec$ExptDate=as.character(d$ExptDate)}
	if("User"%in%colnames(d)){cvec$User=as.character(d$User)}
	if("PI"%in%colnames(d)){cvec$PI=as.character(d$PI)}
	if("Condition"%in%colnames(d)){cvec$Condition=as.character(d$Condition)}
	if("Inoc"%in%colnames(d)){cvec$Inoc=as.character(d$Inoc)}
	if("Gene"%in%colnames(d)){cvec$Gene=as.character(d$Gene)}
	return(cvec)
}

##### Make PDFs #####
qfa.plot<-function(file,results,d,fmt="%Y-%m-%d_%H-%M-%S",barcodes=c(),master.plates=c(),treatments=c(),screen.names=c(),screenIDs=c(),maxg=0,maxt=0,logify=FALSE,densityCol="Growth",curves=TRUE,ylabel="Cell density (AU)",ptype="p",ming=0.00001){
	if(!"Column"%in%colnames(d)) d$Column=d$Col

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
	if (length(screenIDs)==0){screenIDs<-unique(results$ScreenID)} 
	results<-results[results$ScreenID%in%screenIDs,]
	d<-d[d$ScreenID%in%screenIDs,]
	if (length(barcodes)==0){barcodes<-unique(results$Barcode)}
	results<-results[results$Barcode%in%barcodes,]
	d<-d[d$Barcode%in%barcodes,]

	print(paste("Plotting for",length(results[,1]),"colonies"))
	# Produce PDF
	colmin<-min(results$Col); rowmin<-min(results$Col)
	colmax<-max(results$Col); rowmax<-max(results$Row)
	pdf(file,4*(colmax-colmin+1),4*(rowmax-rowmin+1))
	cexfctr<-(rowmax*colmax)/384; marge=c(2,1,2,0.75)
	if(rowmax*colmax==96) {cexfctr=1.0;marge=c(2,1,2,0.75)}  
	if(rowmax*colmax>384) {cexfctr<-(rowmax*colmax)/1536; marge=c(2.0,1.5,2.0,1.5)}
	#cexfctr=1.0
	bcount<-0; nbc<-length(barcodes)
	for (bcode in barcodes){
		bcode<-as.character(bcode)
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
		#nrow<-max(rbc$Row)-min(rbc$Row)+1; ncol<-max(rbc$Column)-min(rbc$Column)+1
		nrow<-rowmax-rowmin+1; ncol=colmax-colmin+1
		#nrow=length(unique(rbc$Row)); ncol=length(unique(rbc$Col))
		# Find max growth to set y-axis for this plate
		if (maxg==0) maxg=max(dbc$Growth)
		# Set graphics parameters for each plate
		op<-par(mfrow=c(nrow,ncol),oma=c(13,15,22,1),
		mar=marge,mgp=c(3,1,0),cex=cexfctr)
		## Plot for each row of results for that bcode ##
		for(rno in rowmin:rowmax){
		  for(cno in colmin:colmax){
				params=rbc[(rbc$Row==rno)&(rbc$Col==cno),]
				if(nrow(params)>0){
					rowplot(params,dbc,inoctime,maxg,fmt,maxt,logify,densityCol=densityCol,curves=curves,ptype=ptype,ming=ming)
				}else{
					frame()
				}
			}
		}
		#z<-apply(rbc,1,rowplot,dbc,inoctime,maxg,fmt,maxt,logify,densityCol=densityCol,curves=curves,ptype=ptype)
		# Title for the plate
		maintit<-paste(rbc$Barcode[1],"Treatment:",rbc$Treatment[1],
		"Medium:",rbc$Medium[1],"Plate:",rbc$MasterPlate.Number[1],sep=" ")
		#cextit<-1008/length(strsplit(maintit,split="")[[1]])
		cextit<-500/nchar(maintit)
		title(main=maintit,xlab="Time since inoculation (days)",line=7,
		ylab=ylabel,cex.main=cextit,cex.lab=8,outer=TRUE)
		  par(op)} #bcode
	dev.off()
}

#### Plot a colony's timecourse from a row of the results #####	
rowplot<-function(resrow,dbc,inoctime,maxg,fmt,maxt,logify,densityCol="Growth",curves=TRUE,ptype="p",ming=0.00001){
	row<-as.numeric(resrow['Row']); col<-as.numeric(resrow['Col'])
	if ('Gene'%in%names(resrow)){gene<-resrow['Gene']} else {gene<-resrow['ORF']}
	# Get data for that colony
	dcol<-dbc[(dbc$Row==row)&(dbc$Column==col),]
	growth<-sapply(dcol[[densityCol]],nozero)
	tim<-dcol$Expt.Time
	if("tshift"%in%names(resrow)){tshift<-unique(as.numeric(resrow[["tshift"]]))[1]}else{tshift=0}
	# Draw the curves and data
	logdraw(row,col,resrow,tim,growth,gene,maxg,maxt=maxt,logify=logify,densityCol=densityCol,curves=curves,ptype=ptype,tshift=tshift,logmin=ming)
}

### Converts row no. to position vector ###	
index2pos<-function(index,dbc) c(dbc[[index,'Row']],dbc[[index,'Column']])

### Do individual timecourse plot given parameters & data ###
logdraw<-function(row,col,resrow,tim,growth,gene,maxg,fitfunct,maxt=0,scaleT=1.0,logify=FALSE,densityCol="Growth",curves=TRUE,ptype="p",tshift=0,logmin=0.00001){
	if(logify) {
	  ylog="y"
	  ylim=c(logmin,1.2*maxg)
	}else{
	  ylog=""
	  ylim=c(0,1.2*maxg)
	}
	plot(NULL,type="n",xlim=c(0,maxt),ylim=ylim,log=ylog,xlab="",ylab="",main=gene,frame.plot=0,cex.main=3*scaleT,cex.axis=1*scaleT)
	# Add data points
	points(tim,growth,col="red",cex=1.25*scaleT,pch=4,lwd=2,type=ptype,xlim=c(0,maxt),ylim=c(logmin,1.2*maxg))
	if(curves){
		if("nr_t"%in%names(resrow)) abline(v=resrow['nr_t'],col="blue")
		if("maxslp_t"%in%names(resrow)) abline(v=resrow['maxslp_t'],col="blue",lty=2)
		# Get logistic parameters, gene name and position
		K<-as.numeric(resrow['K']); r<-as.numeric(resrow['r']); g<-as.numeric(resrow['g']); v<-as.numeric(resrow['v']);
		MDR<-as.numeric(resrow['MDR']); MDP<-as.numeric(resrow['MDP']); AUC<-as.numeric(resrow['AUC']); DT<-as.numeric(resrow['DT']);
		if(logify) {ylog="y"}else{ylog=""}
		# Add logistic curve
		if(maxt==0) maxt=ceiling(max(tim))
		x=0
		curve(Glogist(K,r,g,v,x-tshift),n=31,lwd=2.5,add=TRUE,from=tshift,to=maxt,xlim=c(0,maxt),ylim=c(0.00001,1.2*maxg))
	}

	if(curves){
		# Add legend
		legt1<-paste(c("K=","r=","g=","v=","MDR=","MDP=","AUC=","DT="),c(signif(K,3),signif(r,3),signif(g,3),signif(v,3),signif(MDR,3),signif(MDP,3),signif(AUC,3),signif(DT,3)),sep="")
		if(logify){legend("bottomright",legt1,box.lty=0,cex=scaleT)}else{legend("topleft",legt1,box.lty=0,cex=scaleT)}
	}
	legend("topright",sprintf("R%02dC%02d",row,col),box.lty=0,cex=0.5*scaleT)
}

### Get the mode of a vector (why isn't there such a function in base?)
getMode <- function(x) {
  ux = unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}



