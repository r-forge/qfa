#require(sp)
#globs<-new.env()
.onLoad = function(libname,pkgname){
	#assign("globs", new.env(), envir=parent.env(environment()))
	assign("globs", new.env(), envir="namespace:qfa")
	globs$test=1
}

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
	if(qfa:::globs$ratioPlot){
		d=((x-datobj[[qfa:::globs$ind]])/max(datobj[[qfa:::globs$ind]]))^2+(y-datobj$ratio)^2
	}else{
		d=(x-datobj$ControlFitnessSummary)^2+(y-datobj$QueryFitnessSummary)^2
	}
	best=which.min(d)
	return(best)
}

makePlot=function(datno,focusPlot=TRUE,...){
	if(qfa:::globs$compno>0){
		compnm=paste(qfa:::globs$compno,qfa:::globs$GROUPS$GroupName[qfa:::globs$compno],"\t",qfa:::globs$GROUPS$GroupID[qfa:::globs$compno])
	}else{
		compnm=""
	}
	maintitle=paste(qfa:::globs$datlist[[datno]]$datname,compnm,sep="\n")
	if((nchar(qfa:::globs$datlist[[datno]]$cMed)<22)&(nchar(qfa:::globs$datlist[[datno]]$qMed)<22)){
		xlab=paste(qfa:::globs$datlist[[datno]]$cScrID,qfa:::globs$datlist[[datno]]$cCli,qfa:::globs$datlist[[datno]]$cScrNm,qfa:::globs$datlist[[datno]]$cLibs,qfa:::globs$datlist[[datno]]$cUse,qfa:::globs$datlist[[datno]]$cDate,qfa:::globs$datlist[[datno]]$cTreat,qfa:::globs$datlist[[datno]]$cMed)
		ylab=paste(qfa:::globs$datlist[[datno]]$qScrID,qfa:::globs$datlist[[datno]]$qCli,qfa:::globs$datlist[[datno]]$qScrNm,qfa:::globs$datlist[[datno]]$qLibs,qfa:::globs$datlist[[datno]]$qUse,qfa:::globs$datlist[[datno]]$qDate,qfa:::globs$datlist[[datno]]$qTreat,qfa:::globs$datlist[[datno]]$qMed)
	}else{
		xlab=paste(qfa:::globs$datlist[[datno]]$cScrID,qfa:::globs$datlist[[datno]]$cCli,qfa:::globs$datlist[[datno]]$cScrNm,qfa:::globs$datlist[[datno]]$cLibs,qfa:::globs$datlist[[datno]]$cUse,qfa:::globs$datlist[[datno]]$cDate)
		ylab=paste(qfa:::globs$datlist[[datno]]$qScrID,qfa:::globs$datlist[[datno]]$qCli,qfa:::globs$datlist[[datno]]$qScrNm,qfa:::globs$datlist[[datno]]$qLibs,qfa:::globs$datlist[[datno]]$qUse,qfa:::globs$datlist[[datno]]$qDate)
	}
	qfa:::globs$orftargs=qfa:::globs$dat$ORF[qfa:::globs$targs]
	qfa:::globs$dat=qfa:::globs$datlist[[datno]]$res
	qfa:::globs$targs=match(qfa:::globs$orftargs,qfa:::globs$dat$ORF)
	epiplot(qfa:::globs$dat,0.05,mmain=maintitle,xxlab=xlab,yylab=ylab,ymin=qfa:::globs$fymin[datno],ymax=qfa:::globs$fymax[datno],xmin=qfa:::globs$fxmin[datno],xmax=qfa:::globs$fxmax[datno],...)
	if(qfa:::globs$compno>0){
		lst=strsplit(qfa:::globs$GROUPS$GroupORFs[qfa:::globs$compno]," ")[[1]]
		dcomp=qfa:::globs$dat[qfa:::globs$dat$ORF%in%lst,]
		if(length(dcomp$ORF)>0){
			points(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,col="purple",pch=16,cex=1)
			text(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,dcomp$Gene,pos=1,cex=0.75)
		}	
	}
	if(length(qfa:::globs$targs)>0){
		points(qfa:::globs$dat$ControlFitnessSummary[qfa:::globs$targs],qfa:::globs$dat$QueryFitnessSummary[qfa:::globs$targs],col="blue",pch=16,cex=0.2)
		text(qfa:::globs$dat$ControlFitnessSummary[qfa:::globs$targs],qfa:::globs$dat$QueryFitnessSummary[qfa:::globs$targs],qfa:::globs$dat$Gene[qfa:::globs$targs],pos=qfa:::globs$posits,cex=0.75)
	}
	if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
}

ratPlot=function(datno,focusPlot=TRUE,qthresh=0.05,ecol="green",scol="red"){
	if(qfa:::globs$compno>0){
		compnm=paste(qfa:::globs$compno,qfa:::globs$GROUPS$GroupName[qfa:::globs$compno],"\t",qfa:::globs$GROUPS$GroupID[qfa:::globs$compno])
	}else{
		compnm=""
	}
	maintitle=paste(qfa:::globs$datlist[[datno]]$datname,compnm,sep="\n")
	xlab=paste(qfa:::globs$datlist[[datno]]$cScrID)
	ylab=paste(qfa:::globs$datlist[[datno]]$qScrID)
	axislab=paste("Fit ",ylab,"/Fit ",xlab,sep="") 
	
	qfa:::globs$orftargs=qfa:::globs$dat$ORF[qfa:::globs$targs]
	if(qfa:::globs$ind=="gindex") ratlab="Ranked by gene name (alphabetical order)"
	if(qfa:::globs$ind=="oindex") ratlab="Ranked by Y-number (chromosome coordinate)"
	if(qfa:::globs$ind=="cindex") ratlab="Ranked by control fitness"
	if(qfa:::globs$ind=="qindex") ratlab="Ranked by query fitness"
	if(qfa:::globs$ind=="rindex") ratlab="Ranked by query/control fitness ratio"
	qfa:::globs$dat=qfa:::globs$datlist[[datno]]$res
	qfa:::globs$targs=match(qfa:::globs$orftargs,qfa:::globs$dat$ORF)
	plot(qfa:::globs$dat[[qfa:::globs$ind]],qfa:::globs$dat$ratio,xlab=ratlab,log="y",ylab=axislab,main=maintitle,col="grey",cex=0.5,pch=19)
	if(qfa:::globs$compno>0){
		lst=strsplit(qfa:::globs$GROUPS$GroupORFs[qfa:::globs$compno]," ")[[1]]
		qfa:::globs$dcomp=qfa:::globs$dat[qfa:::globs$dat$ORF%in%lst,]
		if(length(qfa:::globs$dcomp$ORF)>0){
			points(qfa:::globs$dcomp[[qfa:::globs$ind]],qfa:::globs$dcomp$ratio,col="purple",pch=16,cex=1)
			text(qfa:::globs$dcomp[[qfa:::globs$ind]],qfa:::globs$dcomp$ratio,qfa:::globs$dcomp$Gene,pos=1,cex=0.75)
		}	
	}
	if(length(qfa:::globs$targs)>0){
		points(qfa:::globs$dat[[qfa:::globs$ind]][qfa:::globs$targs],qfa:::globs$dat$ratio[qfa:::globs$targs],col="blue",pch=16,cex=0.2)
		text(qfa:::globs$dat[[qfa:::globs$ind]][qfa:::globs$targs],qfa:::globs$dat$ratio[qfa:::globs$targs],qfa:::globs$dat$Gene[qfa:::globs$targs],pos=qfa:::globs$posits,cex=0.75)
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
		qfa:::globs$xf=grconvertX(x,from="ndc",to="user")
		qfa:::globs$yf=grconvertY(y,from="ndc",to="user")
		if(qfa:::globs$zooming==TRUE){
			if(identical(qfa:::globs$zoomTL,c(-99,-99))) {
				qfa:::globs$zoomTL=c(qfa:::globs$xf,qfa:::globs$yf)
			}else{
				qfa:::globs$zoomBR=c(qfa:::globs$xf,qfa:::globs$yf)
				# Zoom in!
				qfa:::globs$fxmin[qfa:::globs$datno]=qfa:::globs$zoomTL[1]
				qfa:::globs$fxmax[qfa:::globs$datno]=qfa:::globs$zoomBR[1]
				qfa:::globs$fymin[qfa:::globs$datno]=qfa:::globs$zoomBR[2]
				qfa:::globs$fymax[qfa:::globs$datno]=qfa:::globs$zoomTL[2]
				if(qfa:::globs$ratioPlot){
					ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
				}else{
					makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
				}
				# Reset zooming parameters
				qfa:::globs$zooming=FALSE
				qfa:::globs$zoomTL=c(-99,-99)
				qfa:::globs$zoomBR=c(-99,-99)
			}
		}else{
			if(qfa:::globs$sel==0){
				# Find nearest gene, highlight location, draw plot
				candidate=targFun(qfa:::globs$xf,qfa:::globs$yf,qfa:::globs$dat)
				if(!candidate%in%qfa:::globs$targs) {
					qfa:::globs$targs=c(qfa:::globs$targs,candidate)
					qfa:::globs$posits=c(qfa:::globs$posits,1)
				}else{
					qfa:::globs$targ=which(qfa:::globs$targs==candidate)
					qfa:::globs$pnew=qfa:::globs$posits[qfa:::globs$targ]
					qfa:::globs$posits[qfa:::globs$targ]=(qfa:::globs$pnew%%4)+1
				}
				if(qfa:::globs$ratioPlot){
					ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
				}else{
					makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
				}
				if(length(qfa:::globs$selx)>0) drawSel(qfa:::globs$selx,qfa:::globs$sely)
			}else{
				qfa:::globs$selx=c(qfa:::globs$selx,qfa:::globs$xf)
				qfa:::globs$sely=c(qfa:::globs$sely,qfa:::globs$yf)
				if(qfa:::globs$ratioPlot){
					ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
				}else{
					makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
				}
				drawSel(qfa:::globs$selx,qfa:::globs$sely)
			}
		}
	}
	if(1%in%buttons) {
		# Delete the last gene which was highlighted
		qfa:::globs$targs=qfa:::globs$targs[-length(qfa:::globs$targs)]
		qfa:::globs$posits=qfa:::globs$posits[-length(qfa:::globs$posits)]
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(2%in%buttons) {
		# Open gene webpage on SGD
		qfa:::globs$xf=grconvertX(x,from="ndc",to="user")
		qfa:::globs$yf=grconvertY(y,from="ndc",to="user")
		candidate=targFun(qfa:::globs$xf,qfa:::globs$yf,qfa:::globs$dat)
		browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",qfa:::globs$dat$ORF[candidate],sep=""))
	}
}

keybd=function(key){
	if(key=="q") return(invisible(1))
	if(key=="c") {
		# Clear highlighted/selected points
		qfa:::globs$targs=c()
		qfa:::globs$posits=c()
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="d"){
		# Unhighlight the last gene which was highlighted
		qfa:::globs$targs=qfa:::globs$targs[-length(qfa:::globs$targs)]
		qfa:::globs$posits=qfa:::globs$posits[-length(qfa:::globs$posits)]
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="w"){
		qfa:::globs$targ=qfa:::globs$targs[length(qfa:::globs$targs)]
		browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",qfa:::globs$dat$ORF[qfa:::globs$targ],sep=""))
	}
	if(key=="Right") {
		qfa:::globs$datno<<-(qfa:::globs$datno%%length(qfa:::globs$datlist))+1
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
		qfa:::globs$orftargs=qfa:::globs$dat$ORF[qfa:::globs$targs]
		qfa:::globs$dat=qfa:::globs$datlist[[qfa:::globs$datno]]$res
		qfa:::globs$targs=match(qfa:::globs$orftargs,qfa:::globs$dat$ORF)
	}
	if(key=="Left") {
		qfa:::globs$datno=((qfa:::globs$datno-2)%%length(qfa:::globs$datlist))+1
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
		orftargs=qfa:::globs$dat$ORF[qfa:::globs$targs]
		qfa:::globs$dat=qfa:::globs$datlist[[qfa:::globs$datno]]$res
		qfa:::globs$targs=match(qfa:::globs$orftargs,qfa:::globs$dat$ORF)
	}
	if(key=="Up") {
		qfa:::globs$compno=qfa:::globs$compno+1
		if(qfa:::globs$compno>length(qfa:::globs$GROUPS$GroupORFs)) qfa:::globs$compno=1
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="Down") {
		qfa:::globs$compno=qfa:::globs$compno-1
		if(qfa:::globs$compno<=0) qfa:::globs$compno=length(qfa:::globs$GROUPS$GroupORFs)
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="z") {
		if(qfa:::globs$sel==0){
			qfa:::globs$selx=c()
			qfa:::globs$sely=c()
		}
		qfa:::globs$sel=(qfa:::globs$sel+1)%%3
		if(qfa:::globs$sel==2){
			# Close polygon
			qfa:::globs$selx=c(qfa:::globs$selx,qfa:::globs$selx[1])
			qfa:::globs$sely=c(qfa:::globs$sely,qfa:::globs$sely[1])
			if(qfa:::globs$ratioPlot){
				ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
			}else{
				makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
			}
			drawSel(qfa:::globs$selx,qfa:::globs$sely)
			qfa:::globs$sel=0
		}
	}
	if((key=="s")&(length(qfa:::globs$selx)>0)){
		if(qfa:::globs$ratioPlot){
			ptsind=point.in.polygon(qfa:::globs$dat[[qfa:::globs$ind]],qfa:::globs$dat$ratio,qfa:::globs$selx,qfa:::globs$sely)
		}else{
			ptsind=point.in.polygon(qfa:::globs$dat$ControlFitnessSummary,qfa:::globs$dat$QueryFitnessSummary,qfa:::globs$selx,qfa:::globs$sely)
		}
		qfa:::globs$tsel=which(ptsind==1)
		qfa:::globs$targs=c(qfa:::globs$targs,qfa:::globs$tsel)
		qfa:::globs$posits=c(qfa:::globs$posits,rep(1,length(qfa:::globs$tsel)))
		qfa:::globs$selx=c()
		qfa:::globs$sely=c()
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="p") {
		psname=sprintf("QFAVisualisation%04d.ps",qfa:::globs$plotno)
		cat("Printing plot to file:",file.path(getwd(),psname),"\n")
		cairo_ps(psname)
			if(qfa:::globs$ratioPlot){
				ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol,focusPlot=FALSE)
			}else{
				makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol,focusPlot=FALSE)
			}
		dev.off()
		qfa:::globs$plotno=qfa:::globs$plotno+1
		if ((Sys.info()['sysname']=="Windows")) bringToTop()
	}
	if(key=="n") {
		pngname=sprintf("QFAVisualisation%04d.png",qfa:::globs$plotno)
		cat("Printing plot to file:",file.path(getwd(),pngname),"\n")
		png(pngname,width=480*2,height=480*2,pointsize=12*2)
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol,focusPlot=FALSE)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol,focusPlot=FALSE)
		}
		dev.off()
		qfa:::globs$plotno=qfa:::globs$plotno+1
		if ((Sys.info()['sysname']=="Windows")) bringToTop()
	}
	if(key=="m") {
		pdfname=sprintf("QFAVisualisation%04d.pdf",qfa:::globs$plotno)
		cat("Printing plot to file:",file.path(getwd(),pdfname),"\n")
		pdf(pdfname)
			if(qfa:::globs$ratioPlot){
				ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol,focusPlot=FALSE)
			}else{
				makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol,focusPlot=FALSE)
			}
		dev.off()
		qfa:::globs$plotno=qfa:::globs$plotno+1
		if ((Sys.info()['sysname']=="Windows")) bringToTop()
	}
	if(key=="t") {
		ENH=qfa:::globs$Ecol
		SUP=qfa:::globs$Scol
		qfa:::globs$Scol=ENH
		qfa:::globs$Ecol=SUP
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="r") {
	# Zoom mode (click on TL and BR)
		qfa:::globs$zooming=TRUE
		qfa:::globs$zoomTL=c(-99,-99)
		qfa:::globs$zoomBR=c(-99,-99)
	}
	if(key=="u") {
	# Get user input for a new list of genes
		if (Sys.info()['sysname']=="Windows") bringToTop(which=-1)
		newgrp=getText(qfa:::globs$ORFGENE)
		newdf=data.frame(GroupName=newgrp[["Label"]],GroupID="R VisTool")
		newdf$GroupORFs=paste(newgrp[["ORFs"]],collapse=" ")
		qfa:::globs$GROUPS=rbind(newdf,qfa:::globs$GROUPS)
		qfa:::globs$compno=1
		if(qfa:::globs$ratioPlot){
			ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}else{
			makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
		}
	}
	if(key=="l") {
	# Log ratio plot
	qfa:::globs$ratioPlot=!qfa:::globs$ratioPlot
	if(qfa:::globs$ratioPlot){
		ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
	}else{
		makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
	}
	}
	if(key=="i") {
	# Change index for ratio plot
	qfa:::globs$indNo=qfa:::globs$indNo+1
	qfa:::globs$ind=qfa:::globs$indPoss[qfa:::globs$indNo%%length(qfa:::globs$indPoss)+1]
	if(qfa:::globs$ratioPlot){ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)}
	}
	if(key=="b") {
	# Print data to console
	cat("*************\n")
	print(qfa:::globs$datlist[[qfa:::globs$datno]]$rept$V1)
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
	hdr=read.delim(filename,nrows=19,header=FALSE,stringsAsFactors=FALSE)
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
	fMax=fMax,rept=hdr)
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
	# Create an environment for storing "global" variables
	#globs<-new.env(parent=parent.env(environment()))
	globs=new.env()
	#environment(globs)=asNamespace("qfa")
	environment(globs)=as.environment("package:qfa")

	# Function for generating interactive plot
	qfa:::globs$GROUPS=groups
	qfa:::globs$ORFGENE=orf2gene
	qfa:::globs$datlist=list()
	
	for (f in 1:length(GISfiles)){	
		report=getResults(GISfiles[f])
		report$res$ratio=pmax(0.005,report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
		report$res$gindex=rank(as.character(report$res$Gene))
		report$res$oindex=rank(as.character(report$res$ORF))
		report$res$cindex=rank(report$res$ControlFitnessSummary)
		report$res$qindex=rank(report$res$QueryFitnessSummary)
		report$res$rindex=rank(report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
		report$datname="Mean Fitnesses (MDR * MDP), t-test"
		qfa:::globs$datlist[[f]]=report
	}

	qfa:::globs$fxmax=c();qfa:::globs$fymax=c()
	qfa:::globs$fxmin=c();qfa:::globs$fymin=c()

	for(d in qfa:::globs$datlist){
		qfa:::globs$fxmax=c(qfa:::globs$fxmax,d$fMax)
		qfa:::globs$fymax=c(qfa:::globs$fymax,d$fMax)
		qfa:::globs$fxmin=c(0,qfa:::globs$fxmin)
		qfa:::globs$fymin=c(0,qfa:::globs$fymin)
	}
	
	print(qfa:::globs$fxmit)
		
	qfa:::globs$n=0
	qfa:::globs$targs=c()
	qfa:::globs$posits=c()
	qfa:::globs$sel=0
	qfa:::globs$selx=c()
	qfa:::globs$sely=c()
	qfa:::globs$indPoss=c("gindex","oindex","rindex","cindex","qindex")
	qfa:::globs$indNo=0
	qfa:::globs$ind=qfa:::globs$indPoss[qfa:::globs$indNo%%length(qfa:::globs$indPoss)+1]
	
	# Toggling between colour of suppressors and enhancers
	qfa:::globs$Ecol="green"
	qfa:::globs$Scol="red"

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
	cat("i: Toggle between log ratio plot styles.\n")
	cat("r: Enter zoom mode.  Click on top left and bottom right of area to inspect.\n")
	cat("p: Save current plot as vector graphic to QFAVisualisation.ps.  Other filetypes can also be generated - 'n': .png and 'm': .pdf\n")
	cat("b: Print experimental metadata report to console window\n")
	cat("q: Quit tool and print gene names currently selected to console window.\n")

	qfa:::globs$datno=1
	qfa:::globs$plotno=1
	qfa:::globs$compno=0
	qfa:::globs$pos=1
	qfa:::globs$zooming=FALSE
	qfa:::globs$zoomTL=c(-99,-99)
	qfa:::globs$zoomBR=c(-99,-99)
	qfa:::globs$ratioPlot=FALSE

	qfa:::globs$dat=qfa:::globs$datlist[[qfa:::globs$datno]]$res

	# Check if running under windows, if not, force X11
	# Note that event handling doesn't work under Linux :(
	sysinf=Sys.info()
	#if( sysinf["sysname"]!="Windows") x11()
	#if( sysinf["sysname"]=="Linux") Cairo()
	x11()
	print(qfa:::globs$datno)
	if(qfa:::globs$ratioPlot){
		ratPlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
	}else{
		makePlot(qfa:::globs$datno,ecol=qfa:::globs$Ecol,scol=qfa:::globs$Scol)
	}

	getGraphicsEvent(prompt="L click: Highlight/Rotate, R click: SGD, M click: Remove, Left/Right: Change plot, z: select tool, s: add selection, c: clear, q: quit", onMouseDown=mouse, onKeybd=keybd)
	cat("~~~~~~~~~~~~~~~\n")
	cat("\n")
	cat("List of genes currently selected:\n")
	cat(qfa:::globs$dat$Gene[qfa:::globs$targs])
	cat("\n")
	cat("~~~~~~~~~~~~~~~\n")
	cat("\n")
	cat("Systematic names for genes currently selected:\n")
	cat(qfa:::globs$dat$ORF[qfa:::globs$targs])
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