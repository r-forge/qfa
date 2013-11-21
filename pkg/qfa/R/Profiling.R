makeProfiles=function(flist,genChar){
	# Read fitnessReport or GIS files from flist (defined elsewhere)
	# Check if file valid output from qfa R packages
	# Check how many rows to skip before reading data (lskip)
	# Assume same filetype for all files in list
	con=file(flist[[1]],open="r")
	lskip=1
	while(!grepl("##########",readLines(con,n=1,warn=FALSE))) lskip=lskip+1
	close(con)

	dfs=list()
	for(f in names(flist)){
		dfs[[f]]=read.delim(flist[[f]],skip=lskip)
	}

	# Merge all datasets by ORF
	fits=Reduce(function(...) {merge(...,by=c("ORF","Gene"))}, dfs)
	colroots=colnames(dfs[[1]])[3:length(colnames(dfs[[1]]))] 
	cnames=c("ORF","Gene",c(t(sapply(colroots,paste,names(flist),sep="."))))
	colnames(fits)=cnames

	# First attempt, ignoring uncertainty on means
	df=list()
	for(f in names(flist)){
		df[[f]]=fits[[paste(genChar,f,sep=".")]]
	}
	df=data.frame(df)
	rownames(df)=fits$Gene
	# Cannot calculate correlation between rows where one row is filled with identical values 
	# (e.g. all strains have zero fitness)
	df$SD=apply(df,1,sd)
	df$ORF=fits$ORF
	print(df[df$SD==0,])
	df=df[df$SD>0,]
	return(df)
}

# http://stats.stackexchange.com/questions/65705/pairwise-mahalanobis-distance
fastPwMahal = function(x1_raw,invCovMat) {
	x1=as.matrix(x1_raw)
	SQRT = with(svd(invCovMat), u %*% diag(d^0.5) %*% t(v))
	as.matrix(dist(x1 %*% SQRT))
}

# Generate correlation, distance and Mahalanobis matrices
similarities=function(prof){
	invCovMat = solve(cov(prof[,1:length(flist)]))
	res=list(
	corMat=cor(t(as.matrix(prof[,1:length(flist)])),use="complete.obs"),
	distMat=as.matrix(dist(prof[,1:length(flist)])),
	mahMat=fastPwMahal(prof[,1:length(flist)],invCovMat)
	)
	return(res)
}

# Find genes most similar to a target
findSimilar=function(target,sims){
	hitCorr=sims$corMat[,target][order(sims$corMat[,target],decreasing=TRUE)]
	hitDist=sims$distMat[,target][order(sims$distMat[,target],decreasing=FALSE)]
	hitMah=sims$mahMat[,target][order(sims$mahMat[,target],decreasing=FALSE)]
	res=list()
	res[["Corr"]]=hitCorr
	res[["Dist"]]=hitDist
	res[["Mah"]]=hitMah
	return(res)
}

# Plot fitness profiles for most similar and most highly correlated genes
plotSimilarFit=function(targ,prof,sims,cutoff=50,ylim=c()){
	# Number of experiments
	nExp=dim(prof)[2]-2

	if(length(ylim)==0) ylim=quantile(as.numeric(unlist(prof[,1:length(flist)])),c(0,1))
	sim=findSimilar(targ,sims)
	corrs=head(sim$Corr,cutoff)
	dists=head(sim$Dist,cutoff)
	mahs=head(sim$Mah,cutoff)

	# Visualise similarity of similar genes
	op=par(mfrow=c(1,3),mai=c(1,0.5,0.5,0.25))
	cprof=prof[rownames(prof)%in%names(corrs),1:nExp]
	dprof=prof[rownames(prof)%in%names(dists),1:nExp]
	mprof=prof[rownames(prof)%in%names(mahs),1:nExp]

	# Correlation
	boxplot(prof[,1:nExp],col="lightblue",xlab="",ylab="Fitness",main=paste("Correlates with",targ),notch=FALSE,outline=FALSE,border="darkgrey",ylim=ylim,las=2)
	for(i in 1:dim(cprof)[1]) points(1:dim(cprof)[2],cprof[i,],type="l")
	points(1:dim(cprof)[2],cprof[targ,],col="red",type="l",lwd=3)
	text(dim(cprof)[2],cprof[,dim(cprof)[2]],tolower(rownames(cprof)),cex=0.5,pos=4)
	text(1,cprof[,1],tolower(rownames(cprof)),cex=0.5,pos=2)

	# Distance
	boxplot(prof[,1:nExp],col="lightblue",xlab="",ylab="Fitness",main=paste("Near to",targ),notch=FALSE,outline=FALSE,border="darkgrey",ylim=ylim,las=2)
	for(i in 1:dim(dprof)[1]) points(1:dim(dprof)[2],dprof[i,],type="l")
	points(1:dim(dprof)[2],dprof[targ,],col="red",type="l",lwd=3)
	text(dim(dprof)[2],dprof[,dim(dprof)[2]],tolower(rownames(dprof)),cex=0.5,pos=4)
	text(1,dprof[,1],tolower(rownames(dprof)),cex=0.5,pos=2)

	# Mahalanobis
	boxplot(prof[,1:nExp],col="lightblue",xlab="",ylab="Fitness",main=paste("Mahalanobis with",targ),notch=FALSE,outline=FALSE,border="darkgrey",ylim=ylim,las=2)
	for(i in 1:dim(mprof)[1]) points(1:dim(mprof)[2],mprof[i,],type="l")
	points(1:dim(mprof)[2],mprof[targ,],col="red",type="l",lwd=3)
	text(dim(mprof)[2],mprof[,dim(mprof)[2]],tolower(rownames(mprof)),cex=0.5,pos=4)
	text(1,mprof[,1],tolower(rownames(mprof)),cex=0.5,pos=2)
	
	par(op)

	return(list(correlation=corrs,distance=dists,mahalanobis=mahs))
}

#prof=makeProfiles(flist,"MedianFit")
#sims=similarities(prof)
#plotSimilarFit("SPE1",prof,sims,5)

# Just plot fitness profiles for a list of genes
plotGroupFits=function(glists,prof,clist=c(),mtitle="",ylim=c()){
	# Number of experiments
	nExp=dim(prof)[2]-2
	if(length(ylim)==0) ylim=quantile(as.numeric(unlist(prof[,1:nExp])),c(0,1))
	if(length(clist)<length(glists)) clist=rainbow(length(glists))
	#if(length(ltypes)<length(glists)) ltypes=rep(1,length(glists))
	op=par(mai=c(1.25,0.85,0.5,0.25))
	boxplot(prof[,1:nExp],col="lightgrey",xlab="",ylab="Median Fitness",main=mtitle,notch=FALSE,outline=FALSE,border="black",ylim=ylim,lwd=1.5,las=2)
	
	for(k in 1:length(names(glists))){	
		glist=glists[[names(glists)[k]]]
		cprof=prof[rownames(prof)%in%glist,]
		for(i in 1:dim(cprof)[1]) points(1:length(flist),cprof[i,1:length(flist)],type="l",col=clist[k],lwd=3,lty=2)
		text(length(flist),cprof[,length(flist)],tolower(rownames(cprof)),cex=0.5,pos=4,col=clist[k])
		text(1,cprof[,1],tolower(rownames(cprof)),cex=0.5,pos=2,col=clist[k])
	}
	legend("topright",names(glists),text.col=clist,col=clist,lty=2,lwd=2)
	par(op)

}

#plotGroupFits(list(spe=c("SPE1","SPE2","SPE3"),mrx=c("MRE11","XRS2","RAD50")),prof,clist=c("cornflowerblue","chartreuse4","darkorchid","orangered"))

# Just plot fitness profiles for a list of genes
plotFitsPoints=function(targlist,prof,pindex,colour="black",textleft=FALSE,alphascl=0.5,lty=1,mtitle="",useORFs=FALSE,ylim=c()){
	# Number of experiments
	nExp=dim(prof)[2]-2
	if(length(ylim)==0) ylim=quantile(as.numeric(unlist(prof[,1:nExp])),c(0,1))
	lightcolour=sapply(colour,col2rgb,alpha=TRUE)/255
	lightcolour[4,]=alphascl
	lightcolour=rgb(lightcolour[1,],lightcolour[2,],lightcolour[3,],lightcolour[4,])

	if(useORFs){
		cprof=prof[prof$ORF%in%targlist,]
	}else{
		cprof=prof[rownames(prof)%in%targlist,]
	}
	#if(mtitle=="") mtitle=paste(targlist,collapse=" ")
	boxplot(prof[,1:nExp],col="lightgrey",xlab="",ylab="Median Fitness",main=mtitle,notch=FALSE,outline=FALSE,border="darkgrey",ylim=ylim,cex.lab=1.5,lwd=1.5,las=2)
	if(dim(cprof)[1]>0) {
		pchlist=c()
		for(i in 1:dim(cprof)[1]){
				points(1:nExp,cprof[i,1:nExp],type="b",lty=lty,pch="",lwd=3,cex=2,col=lightcolour[(i-1)%%(length(colour))+1])
				points(1:nExp,cprof[i,1:nExp],type="p",lty=lty,pch=pindex,lwd=2,cex=2,col=colour[(i-1)%%(length(colour))+1])
			pchlist=c(pchlist,pindex)
			pindex=pindex+1
		}
		text(1,cprof[,1],tolower(rownames(cprof)),cex=1,pos=2,col=colour[(1:dim(cprof)[1]-1)%%(length(colour))+1])
		text(length(flist),cprof[,length(flist)],tolower(rownames(cprof)),cex=1,pos=4,col=colour[(1:dim(cprof)[1]-1)%%(length(colour))+1])
	}
	if(!textleft) legend("topright",tolower(rownames(cprof)),pch=pchlist,col=colour[(1:dim(cprof)[1]-1)%%(length(colour))+1],cex=1,pt.lwd=2,pt.cex=2)
	return(pindex)
}

#plotFitsPoints(c("SPE1","SPE2","SPE3"),prof,15,clist,alphascl=0.5,lty=6,mtitle="")

