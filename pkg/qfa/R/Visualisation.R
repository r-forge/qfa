require(sp)

# Plotting function
epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",yylab="Query Fitness",mmain="Epistasis Plot",ymin=0,ymax=0,xmin=0,xmax=0,ecol="green",scol="red"){
	enhancers<-gethits(results,qthresh,type="E")
	suppressors<-gethits(results,qthresh,type="S")
	others<-results[(!results$ORF%in%enhancers$ORF)&(!results$ORF%in%suppressors$ORF),]
	plot(NULL,type="n",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
	xlab=xxlab,ylab=yylab,main=mmain,
	col=8,pch=19,cex=0.5)
	# Add line for genetic independence
	if (fitratio!=FALSE){abline(0,fitratio,lwd=2,col=8)} else {
		abline(0,lm.epi(results$QueryFitnessSummary,results$ControlFitnessSummary,modcheck=FALSE),lwd=2,col=8)}
	# Add 1:1 fitness line
	abline(0,1,lwd=2,lty=4,col=8)
	# Add reference ORF fitnesses lines
	if (ref.orf!=FALSE){
		reforf<-results[results$ORF==ref.orf,]
		abline(v=reforf$ControlFitnessSummary,col="lightblue",lwd=2)
		abline(h=reforf$QueryFitnessSummary,col="lightblue",lwd=2)}
	# Add points for non-suppressors & non-enhancers
	points(others$ControlFitnessSummary,others$QueryFitnessSummary,col="grey",cex=0.5,pch=19)
	#text(others$ControlFitnessSummary,others$QueryFitnessSummary,others$Gene,col="grey",pos=4,offset=0.1,cex=0.4)
	# Add suppressors & enhancers
	if(length(enhancers$ORF)>0){
		points(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,col=ecol,pch=19,cex=0.5)
		#text(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
	if(length(suppressors$ORF)>0){
		points(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,col=scol,pch=19,cex=0.5)
		#text(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
}

targFun=function(x,y,datobj){
	if(ratioPlot){
		d=((x-datobj[[ind]])/max(datobj[[ind]]))^2+(y-datobj$ratio)^2
	}else{
		d=(x-datobj$ControlFitnessSummary)^2+(y-datobj$QueryFitnessSummary)^2
	}
	best=which.min(d)
	return(best)
}

makePlot=function(datno,focusPlot=TRUE,...){
	if(compno>0){
		compnm=paste(compno,GROUPS$GroupName[compno],"\t",GROUPS$GroupID[compno])
	}else{
		compnm=""
	}
	maintitle=paste(datlist[[datno]]$datname,compnm,sep="\n")
	if((nchar(datlist[[datno]]$cMed)<22)&(nchar(datlist[[datno]]$qMed)<22)){
		xlab=paste(datlist[[datno]]$cScrID,datlist[[datno]]$cCli,datlist[[datno]]$cScrNm,datlist[[datno]]$cLibs,datlist[[datno]]$cUse,datlist[[datno]]$cDate,datlist[[datno]]$cTreat,datlist[[datno]]$cMed)
		ylab=paste(datlist[[datno]]$qScrID,datlist[[datno]]$qCli,datlist[[datno]]$qScrNm,datlist[[datno]]$qLibs,datlist[[datno]]$qUse,datlist[[datno]]$qDate,datlist[[datno]]$qTreat,datlist[[datno]]$qMed)
	}else{
		xlab=paste(datlist[[datno]]$cScrID,datlist[[datno]]$cCli,datlist[[datno]]$cScrNm,datlist[[datno]]$cLibs,datlist[[datno]]$cUse,datlist[[datno]]$cDate)
		ylab=paste(datlist[[datno]]$qScrID,datlist[[datno]]$qCli,datlist[[datno]]$qScrNm,datlist[[datno]]$qLibs,datlist[[datno]]$qUse,datlist[[datno]]$qDate)
	}
	orftargs=dat$ORF[targs]
	dat<<-datlist[[datno]]$res
	targs<<-match(orftargs,dat$ORF)
	epiplot(dat,0.05,mmain=maintitle,xxlab=xlab,yylab=ylab,ymin=fymin[datno],ymax=fymax[datno],xmin=fxmin[datno],xmax=fxmax[datno],...)
	if(compno>0){
		lst=strsplit(GROUPS$GroupORFs[compno]," ")[[1]]
		dcomp=dat[dat$ORF%in%lst,]
		if(length(dcomp$ORF)>0){
			points(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,col="purple",pch=16,cex=1)
			text(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,dcomp$Gene,pos=1,cex=0.75)
		}	
	}
	if(length(targs)>0){
		points(dat$ControlFitnessSummary[targs],dat$QueryFitnessSummary[targs],col="blue",pch=16,cex=0.2)
		text(dat$ControlFitnessSummary[targs],dat$QueryFitnessSummary[targs],dat$Gene[targs],pos=posits,cex=0.75)
	}
	if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
}

ratPlot=function(datno,focusPlot=TRUE,qthresh=0.05,ecol="green",scol="red"){
	if(compno>0){
		compnm=paste(compno,GROUPS$GroupName[compno],"\t",GROUPS$GroupID[compno])
	}else{
		compnm=""
	}
	maintitle=paste(datlist[[datno]]$datname,compnm,sep="\n")
	xlab=paste(datlist[[datno]]$cScrID)
	ylab=paste(datlist[[datno]]$qScrID)
	axislab=paste("Fit ",ylab,"/Fit ",xlab,sep="") 
	
	orftargs=dat$ORF[targs]
	if(ind=="gindex") ratlab="Ranked by gene name (alphabetical order)"
	if(ind=="oindex") ratlab="Ranked by Y-number (chromosome coordinate)"
	if(ind=="cindex") ratlab="Ranked by control fitness"
	if(ind=="qindex") ratlab="Ranked by query fitness"
	if(ind=="rindex") ratlab="Ranked by query/control fitness ratio"
	dat<<-datlist[[datno]]$res
	targs<<-match(orftargs,dat$ORF)
	plot(dat[[ind]],dat$ratio,xlab=ratlab,log="y",ylab=axislab,main=maintitle,col="grey",cex=0.5,pch=19)
	if(compno>0){
		lst=strsplit(GROUPS$GroupORFs[compno]," ")[[1]]
		dcomp=dat[dat$ORF%in%lst,]
		if(length(dcomp$ORF)>0){
			points(dcomp[[ind]],dcomp$ratio,col="purple",pch=16,cex=1)
			text(dcomp[[ind]],dcomp$ratio,dcomp$Gene,pos=1,cex=0.75)
		}	
	}
	if(length(targs)>0){
		points(dat[[ind]][targs],dat$ratio[targs],col="blue",pch=16,cex=0.2)
		text(dat[[ind]][targs],dat$ratio[targs],dat$Gene[targs],pos=posits,cex=0.75)
	}
	if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
}

drawSel=function(selx,sely){
	if(length(selx)==1) {
		points(selx,sely,col="black",pch=16,cex=0.4)
	}else{
		points(selx,sely,col="black",lty=2,type="l")
	}
}	

mouse=function(buttons,x,y){
	if(0%in%buttons) {
		# Problem with Cairo coordinate system here...
		xf=grconvertX(x,from="ndc",to="user")
		yf=grconvertY(y,from="ndc",to="user")
		if(zooming==TRUE){
			if(identical(zoomTL,c(-99,-99))) {
				zoomTL<<-c(xf,yf)
			}else{
				zoomBR<<-c(xf,yf)
				# Zoom in!
				fxmin[datno]<<-zoomTL[1]
				fxmax[datno]<<-zoomBR[1]
				fymin[datno]<<-zoomBR[2]
				fymax[datno]<<-zoomTL[2]
				if(ratioPlot){
					ratPlot(datno,ecol=Ecol,scol=Scol)
				}else{
					makePlot(datno,ecol=Ecol,scol=Scol)
				}
				# Reset zooming parameters
				zooming<<-FALSE
				zoomTL<<-c(-99,-99)
				zoomBR<<-c(-99,-99)
			}
		}else{
			if(sel==0){
				# Find nearest gene, highlight location, draw plot
				candidate=targFun(xf,yf,dat)
				if(!candidate%in%targs) {
					targs<<-c(targs,candidate)
					posits<<-c(posits,1)
				}else{
					targ=which(targs==candidate)
					pnew=posits[targ]
					posits[targ]<<-(pnew%%4)+1
				}
				if(ratioPlot){
					ratPlot(datno,ecol=Ecol,scol=Scol)
				}else{
					makePlot(datno,ecol=Ecol,scol=Scol)
				}
				if(length(selx)>0) drawSel(selx,sely)
			}else{
				selx<<-c(selx,xf)
				sely<<-c(sely,yf)
				if(ratioPlot){
					ratPlot(datno,ecol=Ecol,scol=Scol)
				}else{
					makePlot(datno,ecol=Ecol,scol=Scol)
				}
				drawSel(selx,sely)
			}
		}
	}
	if(1%in%buttons) {
		# Delete the last gene which was highlighted
		targs<<-targs[-length(targs)]
		posits<<-posits[-length(posits)]
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(2%in%buttons) {
		# Open gene webpage on SGD
		xf=grconvertX(x,from="ndc",to="user")
		yf=grconvertY(y,from="ndc",to="user")
		candidate=targFun(xf,yf,dat)
		browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",dat$ORF[candidate],sep=""))
	}
}

keybd=function(key){
	if(key=="q") return(invisible(1))
	if(key=="c") {
		# Clear highlighted/selected points
		targs<<-c()
		posits<<-c()
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="d"){
		# Unhighlight the last gene which was highlighted
		targs<<-targs[-length(targs)]
		posits<<-posits[-length(posits)]
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="w"){
		targ<<-targs[length(targs)]
		browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",dat$ORF[targ],sep=""))
	}
	if(key=="Right") {
		datno<<-(datno%%length(datlist))+1
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
		orftargs=dat$ORF[targs]
		dat<<-datlist[[datno]]$res
		targs<<-match(orftargs,dat$ORF)
	}
	if(key=="Left") {
		datno<<-((datno-2)%%length(datlist))+1
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
		orftargs=dat$ORF[targs]
		dat<<-datlist[[datno]]$res
		targs<<-match(orftargs,dat$ORF)
	}
	if(key=="Up") {
		compno<<-compno+1
		if(compno>length(GROUPS$GroupORFs)) compno<<-1
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="Down") {
		compno<<-compno-1
		if(compno<=0) compno<<-length(GROUPS$GroupORFs)
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="z") {
		if(sel==0){
			selx<<-c()
			sely<<-c()
		}
		sel<<-(sel+1)%%3
		if(sel==2){
			# Close polygon
			selx<<-c(selx,selx[1])
			sely<<-c(sely,sely[1])
			if(ratioPlot){
				ratPlot(datno,ecol=Ecol,scol=Scol)
			}else{
				makePlot(datno,ecol=Ecol,scol=Scol)
			}
			drawSel(selx,sely)
			sel<<-0
		}
	}
	if((key=="s")&(length(selx)>0)){
		if(ratioPlot){
			ind=point.in.polygon(dat[[ind]],dat$ratio,selx,sely)
		}else{
			ind=point.in.polygon(dat$ControlFitnessSummary,dat$QueryFitnessSummary,selx,sely)
		}
		tsel=which(ind==1)
		targs<<-c(targs,tsel)
		posits<<-c(posits,rep(1,length(tsel)))
		selx<<-c()
		sely<<-c()
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="p") {
		psname=sprintf("QFAVisualisation%04d.ps",plotno)
		cat("Printing plot to file:",file.path(getwd(),psname),"\n")
		cairo_ps(psname)
			if(ratioPlot){
				ratPlot(datno,ecol=Ecol,scol=Scol,focusPlot=FALSE)
			}else{
				makePlot(datno,ecol=Ecol,scol=Scol,focusPlot=FALSE)
			}
		dev.off()
		plotno<<-plotno+1
		if ((Sys.info()['sysname']=="Windows")) bringToTop()
	}
	if(key=="n") {
		pngname=sprintf("QFAVisualisation%04d.png",plotno)
		cat("Printing plot to file:",file.path(getwd(),pngname),"\n")
		png(pngname,width=480*2,height=480*2,pointsize=12*2)
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol,focusPlot=FALSE)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol,focusPlot=FALSE)
		}
		dev.off()
		plotno<<-plotno+1
		if ((Sys.info()['sysname']=="Windows")) bringToTop()
	}
	if(key=="m") {
		pdfname=sprintf("QFAVisualisation%04d.pdf",plotno)
		cat("Printing plot to file:",file.path(getwd(),pdfname),"\n")
		pdf(pdfname)
			if(ratioPlot){
				ratPlot(datno,ecol=Ecol,scol=Scol,focusPlot=FALSE)
			}else{
				makePlot(datno,ecol=Ecol,scol=Scol,focusPlot=FALSE)
			}
		dev.off()
		plotno<<-plotno+1
		if ((Sys.info()['sysname']=="Windows")) bringToTop()
	}
	if(key=="t") {
		ENH=Ecol
		SUP=Scol
		Scol<<-ENH
		Ecol<<-SUP
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="r") {
	# Zoom mode (click on TL and BR)
		zooming<<-TRUE
		zoomTL<<-c(-99,-99)
		zoomBR<<-c(-99,-99)
	}
	if(key=="u") {
	# Get user input for a new list of genes
		if (Sys.info()['sysname']=="Windows") bringToTop(which=-1)
		newgrp=getText(ORFGENE)
		newdf=data.frame(GroupName=newgrp[["Label"]],GroupID="R VisTool")
		newdf$GroupORFs=paste(newgrp[["ORFs"]],collapse=" ")
		GROUPS<<-rbind(newdf,GROUPS)
		compno<<-1
		if(ratioPlot){
			ratPlot(datno,ecol=Ecol,scol=Scol)
		}else{
			makePlot(datno,ecol=Ecol,scol=Scol)
		}
	}
	if(key=="l") {
	# Log ratio plot
	ratioPlot<<-!ratioPlot
	if(ratioPlot){
		ratPlot(datno,ecol=Ecol,scol=Scol)
	}else{
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	}
	if(key=="i") {
	# Change index for ratio plot
	indNo<<-indNo+1
	ind<<-indPoss[indNo%%length(indPoss)+1]
	if(ratioPlot){ratPlot(datno,ecol=Ecol,scol=Scol)}
	}
	return(NULL)
}

NAtoBlank=function(x){
	if(is.na(x)){
		return("")
	}else{
		return(x)
	}
}

getResults<-function(filename){
	res=read.delim(filename,skip=20,header=TRUE,stringsAsFactors=FALSE)
	hdr=read.delim(filename,nrows=20,header=FALSE,stringsAsFactors=FALSE)
	qfaVersion=NAtoBlank(strsplit(hdr[1,1],": ")[[1]][2])
	summType=NAtoBlank(strsplit(hdr[2,1],": ")[[1]][2])
	testType=NAtoBlank(strsplit(hdr[3,1],": ")[[1]][2])
	cTreat=NAtoBlank(strsplit(hdr[4,1],": ")[[1]][2])
	cMed=NAtoBlank(strsplit(hdr[5,1],": ")[[1]][2])
	cScrID=NAtoBlank(strsplit(hdr[6,1],": ")[[1]][2])
	cScrNm=NAtoBlank(strsplit(hdr[7,1],": ")[[1]][2])
	cLibs=NAtoBlank(strsplit(hdr[8,1],": ")[[1]][2])
	cCli=NAtoBlank(strsplit(hdr[9,1],": ")[[1]][2])
	cUse=NAtoBlank(strsplit(hdr[10,1],": ")[[1]][2])
	cDate=NAtoBlank(strsplit(hdr[11,1],": ")[[1]][2])
	qTreat=NAtoBlank(strsplit(hdr[12,1],": ")[[1]][2])
	qMed=NAtoBlank(strsplit(hdr[13,1],": ")[[1]][2])
	qScrID=NAtoBlank(strsplit(hdr[14,1],": ")[[1]][2])
	qScrNm=NAtoBlank(strsplit(hdr[15,1],": ")[[1]][2])
	qLibs=NAtoBlank(strsplit(hdr[16,1],": ")[[1]][2])
	qCli=NAtoBlank(strsplit(hdr[17,1],": ")[[1]][2])
	qUse=NAtoBlank(strsplit(hdr[18,1],": ")[[1]][2])
	qDate=NAtoBlank(strsplit(hdr[19,1],": ")[[1]][2])
	fMax=max(c(res$QueryFitnessSummary,res$ControlFitnessSummary))
	return(
	list(res=res,qfaVersion=qfaVersion,summType=summType,testType=testType,
	cTreat=cTreat,cMed=cMed,cScrID=cScrID,cScrNm=cScrNm,cLibs=cLibs,cCli=cCli,cUse=cUse,cDate=cDate,
	qTreat=qTreat,qMed=qMed,qScrID=qScrID,qScrNm=qScrNm,qLibs=qLibs,qCli=qCli,qUse=qUse,qDate=qDate,
	fMax=fMax)
	)
}

benschopFromSource<-function(){
	# Read in functionally related complexes
	# Largely from Benschopp Mol Cell 2010
	# Can add some complexes manually
	#compfile=system.file("/FunctionalComplexes.txt", package = "qfa")
	compfile=paste(system.file(package = "qfa"),"/FunctionalComplexes.txt",sep="")
	Benschopp=read.delim(compfile,stringsAsFactors=FALSE,sep="\t")
	Benschopp$CompList=gsub(";"," ",gsub("\\s","",Benschopp$Complex.members..systematic.name))
	res=data.frame(GroupName=Benschopp$X..Complex.name,GroupID="Func.")
	res$GroupORFs=Benschopp$CompList
	return(res)
}

buildBenschop<-function(){
	Benschop=read.delim(file.path(system.file(package = "qfa"),"Benschop.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	return(Benschop)
}

buildGO<-function(){
	#load(file=file.path(system.file(package = "qfa"),"GOAnnotation.Rda"))
	GO=read.delim(file.path(system.file(package = "qfa"),"GOAnnotation.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
	return(GO)
}

visToolDemo<-function(groupFun=buildBenschop){
	groups=groupFun()
	
	orfile=file.path(system.file(package = "qfa"),"ORF2GENE.txt")
	ORFGENE=read.delim(orfile,stringsAsFactors=FALSE,sep="\t",header=FALSE)
	colnames(ORFGENE)=c("ORF","Gene")
	ORFGENE=ORFGENE[!duplicated(ORFGENE$ORF),]

	# Read in GIS files
	filenames=list.files(system.file(package = "qfa"),pattern="*GIS.txt",full.names=TRUE)
	visTool(groups,ORFGENE,filenames)
}

visTool<-function(groups,orf2gene,GISfiles){
	# Function for generating interactive plot
	GROUPS<<-groups
	ORFGENE<<-orf2gene
	datlist<<-list()
	for (f in 1:length(GISfiles)){	
		report=getResults(GISfiles[f])
		report$res$ratio=pmax(0.005,report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
		report$res$gindex=rank(as.character(report$res$Gene))
		report$res$oindex=rank(as.character(report$res$ORF))
		report$res$cindex=rank(report$res$ControlFitnessSummary)
		report$res$qindex=rank(report$res$QueryFitnessSummary)
		report$res$rindex=rank(report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
		report$datname="Mean Fitnesses (MDR * MDP), t-test"
		datlist[[f]]<<-report
	}

	fxmax<<-c();fymax<<-c()
	fxmin<<-c();fymin<<-c()

	for(d in datlist){
		fxmax<<-c(fxmax,d$fMax)
		fymax<<-c(fymax,d$fMax)
		fxmin<<-c(0,fxmin)
		fymin<<-c(0,fymin)
	}
		
	n<<-0
	targs<<-c()
	posits<<-c()
	sel<<-0
	selx<<-c()
	sely<<-c()
	indPoss<<-c("gindex","oindex","rindex","cindex","qindex")
	indNo<<-0
	ind<<-indPoss[indNo%%length(indPoss)+1]
	
	# Toggling between colour of suppressors and enhancers
	Ecol<<-"green"
	Scol<<-"red"

	cat("Windows mouse\n")
	cat("~~~~~~~~~~~~~~~\n")
	cat("Left click: Highlight gene/Rotate text position\n")
	cat("Right click: Open SGD page for gene (alternatively, press 'w' on keyboard)\n")
	cat("Middle click: Remove last gene (alternatively, press 'd' on keyboard)\n\n")
	cat("Mac mouse\n")
	cat("~~~~~~~~~~~~~~~\n")
	cat("Click: Highlight gene/Rotate text position\n\n")
	cat("Keyboard\n")
	cat("~~~~~~~~~~~~~~~\n")
	cat("Left/Right arrow: Switch to next/previous fitness plot (comparison between different pair of experiments).\n")
	cat("Up/Down arrow: Change group of genes currently highlighted (in purple).\n")
	cat("u: Add new group of genes to list of highlightable groups.\n")
	cat("z: Toggle select tool on/off. Select genes for highlighting by drawing a polygon on plot.  Press 's' to add selection when finished.\n") 
	cat("s: Highlight genes encircled using select tool ('z').\n") 
	cat("c: Clear highlighting from currently selected genes.\n") 
	cat("w: Open SGD web-page for last gene highlighted.\n")
	cat("d: Remove highlighting from last gene highlighted.\n")
	cat("t: Toggle colours (red/green) indicating positive and negative interaction.\n")
	cat("l: Toggle plot style between fitness plot and log ratio plot.\n")
	cat("r: Enter zoom mode.  Click on top left and bottom right of area to inspect.\n")
	cat("p: Save current plot as vector graphic to QFAVisualisation.ps.  Other filetypes can also be generated - 'n': .png and 'm': .pdf\n")
	cat("q: Quit tool and print gene names currently selected to console window.\n")

	datno<<-1
	plotno<<-1
	compno<<-0
	pos<<-1
	zooming<<-FALSE
	zoomTL<<-c(-99,-99)
	zoomBR<<-c(-99,-99)
	ratioPlot<<-FALSE

	dat<<-datlist[[datno]]$res

	# Check if running under windows, if not, force X11
	# Note that event handling doesn't work under Linux :(
	sysinf=Sys.info()
	#if( sysinf["sysname"]!="Windows") x11()
	#if( sysinf["sysname"]=="Linux") Cairo()
	x11()
	if(ratioPlot){
		ratPlot(datno,ecol=Ecol,scol=Scol)
	}else{
		makePlot(datno,ecol=Ecol,scol=Scol)
	}

	getGraphicsEvent(prompt="L click: Highlight/Rotate, R click: SGD, M click: Remove, Left/Right: Change plot, z: select tool, s: add selection, c: clear, q: quit", onMouseDown=mouse, onKeybd=keybd)
	cat("~~~~~~~~~~~~~~~\n")
	cat("\n")
	cat("List of genes currently selected:\n")
	cat(dat$Gene[targs])
	cat("\n")
	cat("~~~~~~~~~~~~~~~\n")
	cat("\n")
	cat("Systematic names for genes currently selected:\n")
	cat(dat$ORF[targs])
	cat("\n")
	cat("\n")
	cat("\n")
	#if( sysinf["sysname"]!="Windows") dev.off()
	dev.off()
}

getText=function(ORFGENE){
	x=list()
	cat("Input a list of gene names, separated by spaces & press enter:\n")
	usergenes=scan(what=character(),nlines=1)
	usergenes=toupper(usergenes)
	cat("Please input a label for your list of genes (no spaces) & press enter:\n")
	lab=scan(what=character(),nlines=1,n=1)
	orfs=c()
	genes=c()
	for(gene in usergenes){
		if(gene%in%ORFGENE$Gene){
			orfs=c(orfs,ORFGENE$ORF[ORFGENE$Gene==gene])
			genes=c(genes,gene)
		} else if (gene%in%ORFGENE$ORF) {
			orfs=c(orfs,gene)
			genes=c(genes,ORFGENE$Gene[ORFGENE$ORF==gene])			
		}
	}
	x[["Genes"]]=genes
	x[["ORFs"]]=orfs
	x[["Label"]]=lab
	print(x)
	return(x)
}