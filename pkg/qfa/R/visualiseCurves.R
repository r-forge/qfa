plotRepsDF=function(gdf,target,type="l",lwd=2,mlab="",dmax=0.3,tmax=3){
	gdf$ID=paste(gdf$Barcode,sprintf("%02d",gdf$MasterPlate.Number),sprintf("%02d",gdf$Row),sprintf("%02d",gdf$Col))
	gdf=gdf[order(gdf$ORF,gdf$ID,gdf$Expt.Time),]
	clist=rainbow(length(unique(gdf$ID)))
	names(clist)=unique(gdf$ID)
	plot(gdf$Expt.Time,gdf$Growth,type="n",xlab="time(d)",ylab="Growth (AU)",main=paste(target,mlab),ylim=c(0,dmax),xlim=c(0,min(tmax,max(gdf$Expt.Time))))
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
	for(id in unique(gdf$ID))points(gdf$Expt.Time[gdf$ID==id],gdf$Growth[gdf$ID==id],type=type,col=clist[[id]],lwd=lwd)
}

plotAllReps=function(df,target,mlab="",type="l",lwd=2,returnDat=FALSE,dmax=0.3,tmax=3){
	# Plot growth curves for all available replicates of gene/ORF target in colonyzer.read dataframe df
	if(target%in%df$Gene){
		gdf=df[df$Gene==target,]
	}else{
		gdf=df[df$ORF==target,]
	}
	
	gdf$TrtMed=paste(gdf$Treatments,gdf$Medium)
	TrtMedList=unique(gdf$TrtMed)
	TrtMedList=sort(TrtMedList)
	print("Possible combinations of treatment and medium:")
	print(TrtMedList)
	for(tm in TrtMedList){
		plotRepsDF(gdf[gdf$TrtMed==tm,],target,type=type,lwd=lwd,mlab=paste(mlab,tm),dmax=dmax)
	}

	#gname="SPE1"
	#op=par(mfrow=c(1,2))
	#ares=plotAllReps(a,gname,"Control")
	#bres=plotAllReps(b,gname,"Query")
	#par(op)
}

findImages=function(path,barcRange=c(1,15),tRange=c(17,35),nchar=39,fmt="%Y-%m-%d_%H-%M-%S"){
	# Search for QFA images in file structure, extract as much information as possible from filenames
	flist=list.files(path=path,pattern="\\.(jpg|JPG)$",recursive=TRUE,full.name=TRUE)
	basename(flist)
	dat=data.frame(flist=flist,file=basename(flist),stringsAsFactors=FALSE)
	dat=dat[nchar(dat$file)==nchar,]
	dat$barcode=sapply(dat$file,substr,barcRange[1],barcRange[2])
	dat$time=as.POSIXlt(as.character(sapply(dat$file,substr,tRange[1],tRange[2])),format=fmt)
	dat=dat[order(dat$file),]
	starts=tapply(dat$time,dat$barcode,min,simplify=FALSE)
	dat$starts=do.call(c,starts[dat$barcode])
	dat$exptime=as.numeric(dat$time-dat$starts)/(60*60)
	dat$image=sapply(dat$file,function(x) strsplit(x,"\\.")[[1]][1])
	return(dat)
}