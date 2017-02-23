plotRepsDF=function(gdf,target,type="l",lwd=2,mlab="",dmax=0.3,tmax=3){
	gdf$ID=paste(gdf$Barcode,sprintf("%02d",gdf$MasterPlate.Number),sprintf("%02d",gdf$Row),sprintf("%02d",gdf$Column),sep="_")
	gdf=gdf[order(gdf$ORF,gdf$ID,gdf$Expt.Time),]
	clist=rainbow(length(unique(gdf$ID)))
	names(clist)=unique(gdf$ID)
	plot(gdf$Expt.Time,gdf$Growth,type="n",xlab="time(d)",ylab="Growth (AU)",main=paste(target,mlab),ylim=c(0,dmax),xlim=c(0,min(tmax,max(gdf$Expt.Time))))
	rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey")
	for(id in unique(gdf$ID))points(gdf$Expt.Time[gdf$ID==id],gdf$Growth[gdf$ID==id],type=type,col=clist[[id]],lwd=lwd)
	legend("topleft",lwd=1,legend=unique(gdf$ID),col=clist)
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

plotRow=function(i,gdf,x0=0,y0=0,dx=0.5,dy=0.5,byTime=TRUE){
	im=readJPEG(gdf$imFile[i])
	if(byTime){
		xstart=gdf$Expt.Time[i]
		xfin=gdf$Expt.Time[i]+dx
	}else{
		xstart=x0+(i-1)*dx
		xfin=x0+i*dx
	}
	rasterImage(im[median(gdf$Y.Offset):(median(gdf$Y.Offset)+median(gdf$Tile.Dimensions.Y)),median(gdf$X.Offset):(median(gdf$X.Offset)+median(gdf$Tile.Dimensions.X)),], xstart, y0, xfin, y0+dy)
}

prepareTimecourse=function(df,target,tmax=5.0,freqQuant=0.1,repMax=10,byTime=TRUE){

	if(target%in%df$Gene){
		gdf=df[df$Gene==target,]
	}else{
		gdf=df[df$ORF==target,]
	}
	gdf=gdf[gdf$Expt.Time<=tmax,]

	gdf$PINID=paste(gdf$Barcode,paste("R",sprintf("%02d",gdf$Row),sep=""),paste("C",sprintf("%02d",gdf$Column),sep=""),sep="_")
	listGrow=split(gdf,gdf$PINID)
	# Need to randomly sample replicates if there are too many
	if(length(listGrow)>repMax) listGrow=listGrow[sample(1:length(listGrow),repMax,replace=FALSE)]
	gdf=gdf[gdf$PINID%in%names(listGrow),]

	maxPhotos=max(by(gdf,gdf$PINID,function(x) length(x$ORF)))
	photoFreqs=unlist(by(gdf,gdf$PINID,function(x) diff(x$Expt.Time)))
	photoFreq=quantile(photoFreqs,freqQuant)


	reps=names(listGrow)
	delta=photoFreq
	print(delta)
	if(byTime){
		xmax=tmax
	}else{
		xmax=delta*maxPhotos
	}
	ymax=delta*length(reps)
	print(xmax)
	print(ymax)
	return(list(listGrow=listGrow,xmax=as.numeric(xmax),ymax=as.numeric(ymax),delta=as.numeric(delta),byTime=byTime,asp=as.numeric(ymax/xmax)))
	
}

timecourseImages=function(listGrow,target,xmax,ymax,delta,byTime,cexval=1.0){
	op=par(oma=c(0,0,0,0),mar=c(2,0,0,0))
	if(!byTime) {xtype="n"}else{xtype="s"}
	ps=12
	plot(NULL,xlim=c(0,xmax),ylim=c(0,ymax), type="n",xaxt=xtype,cex=2*cexval,ps=ps)
	text(-0.025*xmax,0.0,paste(listGrow[[1]]$Screen.Name[1],listGrow[[1]]$Condition[1],target),srt=90,cex=1.0*cexval,ps=ps,pos=4,col="darkgrey")
	for(i in 1:length(listGrow)){
		#listGrow[[i]]=listGrow[[i]][listGrow[[i]]$Expt.Time<=xmax,]
		nphotos=length(listGrow[[i]]$ORF)
		if(nphotos>0) {
			sapply(1:nphotos,plotRow,listGrow[[i]],dx=delta,dy=delta,x0=0,y0=(length(listGrow)-i)*delta,byTime=byTime)
			spid=listGrow[[i]]$PINID
			lab=paste(substr(spid,1,15),substr(spid,17,24),sep="\n")
			text(0.05*xmax,(length(listGrow)-i)*delta,lab,pos=3,cex=0.75*cexval,ps=ps)
		}
	}
	par(op)
}

makePDF=function(t1,raw,tmax=5.0,lab="",byTime=FALSE,cexval=1.0){
	p1=prepareTimecourse(raw,t1,tmax=tmax,byTime=byTime)
	pdf(paste(lab,"_",t1,"GrowthCurves.pdf",sep=""))
		plotAllReps(raw,t1,lwd=1)			
	dev.off()
	if(byTime){tlab="byTime"}else{tlab="byPhoto"}
	pdf(paste(lab,"_",t1,tlab,".pdf",sep=""),width=30,height=30*p1$asp)
	timecourseImages(p1$listGrow,t1,p1$xmax,p1$ymax,p1$delta,p1$byTime,cexval=cexval)
	dev.off()
}

findImages=function(path,barcRange=c(1,15),tRange=c(17,35),nchar=39,fmt="%Y-%m-%d_%H-%M-%S"){
	# Search for QFA images in file structure, extract as much information as possible from filenames
	flist=list.files(path=path,pattern="\\.(jpg|JPG)$",recursive=TRUE,full.names=TRUE)
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
