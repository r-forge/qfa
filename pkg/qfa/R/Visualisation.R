# Huge function closure to keep visualisation tool environment out of global environment and keep CRAN happy

makeVisTool=function(){
	globs<-new.env()

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

	makePlot=function(datno,focusPlot=TRUE,...){
		if(globs$compno>0){
			lst=strsplit(globs$GROUPS$GroupORFs[globs$compno]," ")[[1]]
			Ngrp=paste("(N=",length(lst),")",sep="")
			compnm=paste("Group",paste(globs$compno,":",sep=""),globs$GROUPS$GroupName[globs$compno],Ngrp,globs$GROUPS$GroupID[globs$compno])
		}else{
			compnm=""
		}
		#plotdesc=paste(globs$datlist[[datno]]$datname,globs$datlist[[datno]]$cScrNm,globs$datlist[[datno]]$cTreat,":",globs$datlist[[datno]]$qScrNm,globs$datlist[[datno]]$qTreat)
		plotdesc=globs$datlist[[datno]]$datname
		maintitle=paste(plotdesc,compnm,sep="\n")
		if (globs$datlist[[datno]]$cCond=="") {ControlCondition=globs$datlist[[datno]]$cMed}else{ControlCondition=globs$datlist[[datno]]$cCond}
		if (globs$datlist[[datno]]$qCond=="") {QueryCondition=globs$datlist[[datno]]$qMed}else{QueryCondition=globs$datlist[[datno]]$qCond}
		xlab=paste(globs$datlist[[datno]]$cScrNm,globs$datlist[[datno]]$cTreat,ControlCondition,globs$datlist[[datno]]$cFit,globs$datlist[[datno]]$cInoc,globs$datlist[[datno]]$cScrID,globs$datlist[[datno]]$cCli,globs$datlist[[datno]]$cDate)
		ylab=paste(globs$datlist[[datno]]$qScrNm,globs$datlist[[datno]]$qTreat,QueryCondition,globs$datlist[[datno]]$qFit,globs$datlist[[datno]]$qInoc,globs$datlist[[datno]]$qScrID,globs$datlist[[datno]]$qCli,globs$datlist[[datno]]$qDate)

		globs$orftargs=globs$dat$ORF[globs$targs]
		globs$dat=globs$datlist[[datno]]$res
		globs$targs=match(globs$orftargs,globs$dat$ORF)
		epiplot(globs$dat,0.05,mmain=maintitle,xxlab=xlab,yylab=ylab,ymin=globs$fymin[datno],ymax=globs$fymax[datno],xmin=globs$fxmin[datno],xmax=globs$fxmax[datno],...)
		if(globs$compno>0){
			lst=strsplit(globs$GROUPS$GroupORFs[globs$compno]," ")[[1]]
			dcomp=globs$dat[globs$dat$ORF%in%lst,]
			if(length(dcomp$ORF)>0){
				points(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,col="purple",pch=16,cex=1)
				text(dcomp$ControlFitnessSummary,dcomp$QueryFitnessSummary,dcomp$Gene,pos=1,cex=0.75)
			}	
		}
		if(length(globs$targs)>0){
			points(globs$dat$ControlFitnessSummary[globs$targs],globs$dat$QueryFitnessSummary[globs$targs],col="blue",pch=16,cex=0.2)
			text(globs$dat$ControlFitnessSummary[globs$targs],globs$dat$QueryFitnessSummary[globs$targs],globs$dat$Gene[globs$targs],pos=globs$posits,cex=0.75)
		}
		if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
	}

	ratPlot=function(datno,focusPlot=TRUE,qthresh=0.05,ecol="green",scol="red"){
		if(globs$compno>0){
			compnm=paste(globs$compno,globs$GROUPS$GroupName[globs$compno],"\t",globs$GROUPS$GroupID[globs$compno])
		}else{
			compnm=""
		}
		maintitle=paste(globs$datlist[[datno]]$datname,compnm,sep="\n")
		xlab=paste(globs$datlist[[datno]]$cScrID)
		ylab=paste(globs$datlist[[datno]]$qScrID)
		axislab=paste("Fit ",ylab,"/Fit ",xlab,sep="") 
		
		globs$orftargs=globs$dat$ORF[globs$targs]
		if(globs$ind=="gindex") ratlab="Ranked by gene name (alphabetical order)"
		if(globs$ind=="oindex") ratlab="Ranked by Y-number (chromosome coordinate)"
		if(globs$ind=="cindex") ratlab="Ranked by control fitness"
		if(globs$ind=="qindex") ratlab="Ranked by query fitness"
		if(globs$ind=="rindex") ratlab="Ranked by query/control fitness ratio"
		globs$dat=globs$datlist[[datno]]$res
		globs$targs=match(globs$orftargs,globs$dat$ORF)
		plot(globs$dat[[globs$ind]],globs$dat$ratio,xlab=ratlab,log="y",ylab=axislab,main=maintitle,col="grey",cex=0.5,pch=19)
		if(globs$compno>0){
			lst=strsplit(globs$GROUPS$GroupORFs[globs$compno]," ")[[1]]
			globs$dcomp=globs$dat[globs$dat$ORF%in%lst,]
			if(length(globs$dcomp$ORF)>0){
				points(globs$dcomp[[globs$ind]],globs$dcomp$ratio,col="purple",pch=16,cex=1)
				text(globs$dcomp[[globs$ind]],globs$dcomp$ratio,globs$dcomp$Gene,pos=1,cex=0.75)
			}	
		}
		if(length(globs$targs)>0){
			points(globs$dat[[globs$ind]][globs$targs],globs$dat$ratio[globs$targs],col="blue",pch=16,cex=0.2)
			text(globs$dat[[globs$ind]][globs$targs],globs$dat$ratio[globs$targs],globs$dat$Gene[globs$targs],pos=globs$posits,cex=0.75)
		}
		if ((Sys.info()['sysname']=="Windows")&focusPlot) bringToTop()
	}

	targFun=function(x,y,datobj){
		if(globs$ratioPlot){
			d=((x-datobj[[globs$ind]])/max(datobj[[globs$ind]]))^2+(y-datobj$ratio)^2
		}else{
			d=(x-datobj$ControlFitnessSummary)^2+(y-datobj$QueryFitnessSummary)^2
		}
		best=which.min(d)
		return(best)
	}

	mouse=function(buttons,x,y){
		if(0%in%buttons) {
			# Problem with Cairo coordinate system here...
			globs$xf=grconvertX(x,from="ndc",to="user")
			globs$yf=grconvertY(y,from="ndc",to="user")
			if(globs$zooming==TRUE){
				if(identical(globs$zoomTL,c(-99,-99))) {
					globs$zoomTL=c(globs$xf,globs$yf)
				}else{
					globs$zoomBR=c(globs$xf,globs$yf)
					# Zoom in!
					globs$fxmin[globs$datno]=globs$zoomTL[1]
					globs$fxmax[globs$datno]=globs$zoomBR[1]
					globs$fymin[globs$datno]=globs$zoomBR[2]
					globs$fymax[globs$datno]=globs$zoomTL[2]
					if(globs$ratioPlot){
						ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}else{
						makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}
					# Reset zooming parameters
					globs$zooming=FALSE
					globs$zoomTL=c(-99,-99)
					globs$zoomBR=c(-99,-99)
				}
			}else{
				if(globs$sel==0){
					# Find nearest gene, highlight location, draw plot
					candidate=targFun(globs$xf,globs$yf,globs$dat)
					if(!candidate%in%globs$targs) {
						globs$targs=c(globs$targs,candidate)
						globs$posits=c(globs$posits,1)
					}else{
						globs$targ=which(globs$targs==candidate)
						globs$pnew=globs$posits[globs$targ]
						globs$posits[globs$targ]=(globs$pnew%%4)+1
					}
					if(globs$ratioPlot){
						ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}else{
						makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}
					if(length(globs$selx)>0) drawSel(globs$selx,globs$sely)
				}else{
					globs$selx=c(globs$selx,globs$xf)
					globs$sely=c(globs$sely,globs$yf)
					if(globs$ratioPlot){
						ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}else{
						makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
					}
					drawSel(globs$selx,globs$sely)
				}
			}
		}
		if(1%in%buttons) {
			# Delete the last gene which was highlighted
			globs$targs=globs$targs[-length(globs$targs)]
			globs$posits=globs$posits[-length(globs$posits)]
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(2%in%buttons) {
			# Open gene webpage on SGD
			globs$xf=grconvertX(x,from="ndc",to="user")
			globs$yf=grconvertY(y,from="ndc",to="user")
			candidate=targFun(globs$xf,globs$yf,globs$dat)
			directurl="http://www.this-page-intentionally-left-blank.org/"
			if(!is.na(pmatch("Y",globs$dat$ORF[candidate]))) directurl=paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",globs$dat$ORF[candidate],sep="")
			if(!is.na(pmatch("SP",globs$dat$ORF[candidate]))) directurl=paste("http://www.pombase.org/spombe/result/",globs$dat$ORF[candidate],sep="")
			browseURL(directurl)
		}
	}

	keybd=function(key){
		if(key=="q") return(invisible(1))
		if(key=="c") {
			# Clear highlighted/selected points
			globs$targs=c()
			globs$posits=c()
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="d"){
			# Unhighlight the last gene which was highlighted
			globs$targs=globs$targs[-length(globs$targs)]
			globs$posits=globs$posits[-length(globs$posits)]
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="w"){
			globs$targ=globs$targs[length(globs$targs)]
			#browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",globs$dat$ORF[globs$targ],sep=""))
			directurl="http://www.this-page-intentionally-left-blank.org/"
			if(!is.na(pmatch("Y",globs$dat$ORF[globs$targ]))) directurl=paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",globs$dat$ORF[globs$targ],sep="")
			if(!is.na(pmatch("SP",globs$dat$ORF[globs$targ]))) directurl=paste("http://www.pombase.org/spombe/result/",globs$dat$ORF[globs$targ],sep="")
			browseURL(directurl)
		}
		if(key=="Right") {
			globs$datno<<-(globs$datno%%length(globs$datlist))+1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
			globs$orftargs=globs$dat$ORF[globs$targs]
			globs$dat=globs$datlist[[globs$datno]]$res
			globs$targs=match(globs$orftargs,globs$dat$ORF)
		}
		if(key=="Left") {
			globs$datno=((globs$datno-2)%%length(globs$datlist))+1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
			orftargs=globs$dat$ORF[globs$targs]
			globs$dat=globs$datlist[[globs$datno]]$res
			globs$targs=match(globs$orftargs,globs$dat$ORF)
		}
		if(key=="Up") {
			globs$compno=(globs$compno+1)%%(length(globs$GROUPS$GroupORFs)+1)
			#if(globs$compno>length(globs$GROUPS$GroupORFs)) globs$compno=1
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="Down") {
			globs$compno=(globs$compno-1)%%(length(globs$GROUPS$GroupORFs)+1)
			#if(globs$compno<=0) globs$compno=length(globs$GROUPS$GroupORFs)
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="z") {
			if(globs$sel==0){
				globs$selx=c()
				globs$sely=c()
			}
			globs$sel=(globs$sel+1)%%3
			if(globs$sel==2){
				# Close polygon
				globs$selx=c(globs$selx,globs$selx[1])
				globs$sely=c(globs$sely,globs$sely[1])
				if(globs$ratioPlot){
					ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
				}else{
					makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
				}
				drawSel(globs$selx,globs$sely)
				globs$sel=0
			}
		}
		if((key=="s")&(length(globs$selx)>0)){
			if(globs$ratioPlot){
				ptsind=point.in.polygon(globs$dat[[globs$ind]],globs$dat$ratio,globs$selx,globs$sely)
			}else{
				ptsind=point.in.polygon(globs$dat$ControlFitnessSummary,globs$dat$QueryFitnessSummary,globs$selx,globs$sely)
			}
			globs$tsel=which(ptsind==1)
			globs$targs=c(globs$targs,globs$tsel)
			globs$posits=c(globs$posits,rep(1,length(globs$tsel)))
			globs$selx=c()
			globs$sely=c()
			globs$sel=0
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="p") {
			psname=sprintf("QFAVisualisation%04d.ps",globs$plotno)
			cat("Printing plot to file:",file.path(getwd(),psname),"\n")
			cairo_ps(psname)
				if(globs$ratioPlot){
					ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}else{
					makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}
			dev.off()
			globs$plotno=globs$plotno+1
			if ((Sys.info()['sysname']=="Windows")) bringToTop()
		}
		if(key=="n") {
			pngname=sprintf("QFAVisualisation%04d.png",globs$plotno)
			cat("Printing plot to file:",file.path(getwd(),pngname),"\n")
			png(pngname,width=480*2,height=480*2,pointsize=12*2)
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
			}
			dev.off()
			globs$plotno=globs$plotno+1
			if ((Sys.info()['sysname']=="Windows")) bringToTop()
		}
		if(key=="o"){
			printSelected(globs)
		}
		
		if(key=="m") {
			pdfname=sprintf("QFAVisualisation%04d.pdf",globs$plotno)
			cat("Printing plot to file:",file.path(getwd(),pdfname),"\n")
			pdf(pdfname)
				if(globs$ratioPlot){
					ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}else{
					makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol,focusPlot=FALSE)
				}
			dev.off()
			globs$plotno=globs$plotno+1
			if ((Sys.info()['sysname']=="Windows")) bringToTop()
		}
		if(key=="t") {
			ENH=globs$Ecol
			SUP=globs$Scol
			globs$Scol=ENH
			globs$Ecol=SUP
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="r") {
		# Zoom mode (click on TL and BR)
			globs$zooming=TRUE
			globs$zoomTL=c(-99,-99)
			globs$zoomBR=c(-99,-99)
		}
		if(key=="u") {
		# Get user input for a new list of genes
			if (Sys.info()['sysname']=="Windows") bringToTop(which=-1)
			newgrp=getText(globs$ORFGENE)
			newdf=data.frame(Notes="Generated interactively by user",GroupName=newgrp[["Label"]],GroupID="VisTool")
			newdf$GroupORFs=paste(newgrp[["ORFs"]],collapse=" ")
			globs$GROUPS=rbind(globs$GROUPS,newdf)
			globs$compno=length(globs$GROUPS$GroupORFs)
			if(globs$ratioPlot){
				ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}else{
				makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
			}
		}
		if(key=="l") {
		# Log ratio plot
		globs$ratioPlot=!globs$ratioPlot
		if(globs$ratioPlot){
			ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}else{
			makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}
		}
		if(key=="i") {
		# Change index for ratio plot
		globs$indNo=globs$indNo+1
		globs$ind=globs$indPoss[globs$indNo%%length(globs$indPoss)+1]
		if(globs$ratioPlot){ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)}
		}
		if(key=="b") {
		# Print data to console
		cat("\nExperimental metadata: plot",globs$datno,"\n~~~~~~~~~~~~~~~\n")
		# Split report string and patch in replicate numbers for Query and Control
		rept=globs$datlist[[globs$datno]]$rept
		fcont=grep("Control",rept)
		if (length(fcont)==0) fcont=grep("x-axis",rept)
		findx=tail(fcont,1)
		rept2=c(rept[1:findx],paste("x-axis number of replicates:",globs$datlist[[globs$datno]]$cRep),rept[(findx+1):length(rept)],paste("y-axis number of replicates:",globs$datlist[[globs$datno]]$qRep))
		cat(rept2,sep="\n")
		}
		return(NULL)
	}

	visTool<-function(groups,orf2gene,GISfiles,metaRep="MetaReport.txt"){
		# Function for generating interactive plot
		globs$GROUPS=groups()
		globs$ORFGENE=orf2gene
		globs$datlist=list()
		clients=c()
		dates=c()
		
		for (f in 1:length(GISfiles)){	
			report=getResults(GISfiles[f])
			report$res$ratio=pmax(0.005,report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
			report$res$gindex=rank(as.character(report$res$Gene))
			report$res$oindex=rank(as.character(report$res$ORF))
			report$res$cindex=rank(report$res$ControlFitnessSummary)
			report$res$qindex=rank(report$res$QueryFitnessSummary)
			report$res$rindex=rank(report$res$QueryFitnessSummary/report$res$ControlFitnessSummary)
			#report$datname=paste("Plot",paste(f,":",sep=""),report$summType,"fitness",paste("(",report$testType,")",sep=""))
			globs$datlist[[f]]=report
			clients=c(clients,report$qCli)
			dates=c(dates,report$qDate)
		}
		
		# Re-order plots by query client and then by query date
		# Takes 22s for 138 datasets...
		dlist=unlist(globs$datlist)
		dlist=dlist[order(clients,dates)]
		dlist=as.list(dlist)
		neworder=order(clients,dates)
		newDatList=list()
		for (i in 1:length(neworder)){
			newDatList[[i]]=globs$datlist[[neworder[i]]]
			newDatList[[i]]$datname=paste("Plot",paste(i,":",sep=""),newDatList[[i]]$summType,"fitness",paste("(",newDatList[[i]]$testType,")",sep=""))
		}
		globs$datlist=newDatList
		
		reportExpts(globs,fname=metaRep)

		globs$fxmax=c();globs$fymax=c()
		globs$fxmin=c();globs$fymin=c()

		for(d in globs$datlist){
			globs$fxmax=c(globs$fxmax,d$fMax)
			globs$fymax=c(globs$fymax,d$fMax)
			globs$fxmin=c(0,globs$fxmin)
			globs$fymin=c(0,globs$fymin)
		}
			
		globs$n=0
		globs$targs=c()
		globs$posits=c()
		globs$sel=0
		globs$selx=c()
		globs$sely=c()
		globs$indPoss=c("gindex","oindex","rindex","cindex","qindex")
		globs$indNo=0
		globs$ind=globs$indPoss[globs$indNo%%length(globs$indPoss)+1]
		
		# Toggling between colour of suppressors and enhancers
		globs$Ecol="green"
		globs$Scol="red"

		cat("\nControls: Windows mouse\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("Left click: Highlight gene/Rotate text position\n")
		cat("Right click: Open SGD page for gene (alternatively, press 'w' on keyboard)\n")
		cat("Middle click: Remove last gene (alternatively, press 'd' on keyboard)\n\n")
		cat("Controls: Mac mouse\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("Click: Highlight gene/Rotate text position\n\n")
		cat("Controls: Keyboard\n")
		cat("~~~~~~~~~~~~~~~\n")
		cat("Left/Right arrow: Switch to next/previous fitness plot (comparison between different pair of experiments).\n")
		cat("Up/Down arrow: Change group of genes currently highlighted (in purple).\n")
		cat("b: Print experimental metadata report to console window\n")
		cat("c: Clear highlighting from currently selected genes.\n") 
		cat("d: Remove highlighting from last gene highlighted.\n")
		cat("i: Toggle between log ratio plot styles.\n")
		cat("l: Toggle plot style between fitness plot and log ratio plot.\n")
		cat("o: Print list of manually selected genes currently highlighted to console.\n")
		cat("p: Save current plot as vector graphic to QFAVisualisation.ps.  Other filetypes can also be generated - 'n': .png and 'm': .pdf\n")
		cat("q: Quit tool and print gene names currently selected to console window.\n")
		cat("r: Enter zoom mode.  Click on top left and bottom right of area to inspect (experimental, may have to quit to zoom out...).\n")
		cat("s: Highlight genes encircled using select tool ('z').\n") 
		cat("t: Toggle colours (red/green) indicating positive and negative interaction.\n")
		cat("u: Add new group of genes to list of highlightable groups.\n")
		cat("w: Open SGD web-page for last gene highlighted.\n")
		cat("z: Toggle select tool on/off. Select genes for highlighting by drawing a polygon on plot.  Press 's' to add selection when finished.\n") 

		globs$datno=1
		globs$plotno=1
		globs$compno=0
		globs$pos=1
		globs$zooming=FALSE
		globs$zoomTL=c(-99,-99)
		globs$zoomBR=c(-99,-99)
		globs$ratioPlot=FALSE

		globs$dat=globs$datlist[[globs$datno]]$res

		# Check if running under windows, if not, force X11
		# Note that event handling doesn't work under Linux :(
		sysinf=Sys.info()
		#if( sysinf["sysname"]!="Windows") x11()
		#if( sysinf["sysname"]=="Linux") Cairo()
		x11()
		#dev.new()

		if(globs$ratioPlot){
			ratPlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}else{
			makePlot(globs$datno,ecol=globs$Ecol,scol=globs$Scol)
		}

		getGraphicsEvent(prompt="Left/Right: Change plot, b: Expt. metadta, L-click: Highlight, M-click: Remove, q: quit", onMouseDown=mouse, onKeybd=keybd)
		printSelected(globs)
		#if( sysinf["sysname"]!="Windows") dev.off()
		dev.off()
	}
	
	return(visTool)
}

complexesFromSource<-function(){
	# Read in functionally related complexes
	# Largely from Benschopp Mol Cell 2010
	# Can add some complexes manually
	compfile=file.path(system.file(package = "qfa"),"extdata","FunctionalComplexes.txt",sep="")
	Complexes=read.delim(compfile,stringsAsFactors=FALSE,sep="\t")
	Complexes$CompList=gsub(";"," ",gsub("\\s","",Complexes$Complex.members..systematic.name))
	res=data.frame(GroupName=Complexes$X..Complex.name,GroupID="Func.")
	res$GroupORFs=Complexes$CompList
	return(res)
}

buildComplexes<-function(){
	fname=file.path(system.file(package = "qfa"),"extdata","FunctionalComplexes.txt")
	cat("\nGroups of functionally related complexes are specified in this file which you can edit with any text editor:\n~~~~~~~~~~~~~~~\n")
	cat(paste(fname,"\n"))
	ComplexesData=read.delim(fname,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	return(ComplexesData)
}

buildBenschop<-function() buildComplexes()

visToolDemo<-function(groupFun=buildComplexes){
	orfile=file.path(file.path(system.file(package = "qfa"),"extdata","ORF2GENE.txt.gz"))
	ORFGENE=read.delim(orfile,stringsAsFactors=FALSE,sep="\t",header=FALSE)
	colnames(ORFGENE)=c("ORF","Gene")
	ORFGENE=ORFGENE[!duplicated(ORFGENE$ORF),]

	# Read in GIS files
	filenames=list.files(file.path(system.file(package = "qfa"),"extdata"),pattern="*GIS.txt.gz",full.names=TRUE)
	visTool=makeVisTool()
	visTool(groupFun,ORFGENE,filenames,"MetaReport.txt")
}

NAtoBlank=function(x){
	return(ifelse((is.na(x))|(x=="NA"),"",x))
}

getResults<-function(filename){
	# Need to detect at least two different kinds of header now...
	fopen=file(filename,"rt")
	i=1
	while(!grepl("############",readLines(fopen,1))) i=i+1
	close(fopen)
	
	cStart=-1 # At what line does control metadata start?
	qStart=-1 # At what line does query metadata start?

	res=read.delim(filename,skip=i,header=TRUE,stringsAsFactors=FALSE)
	hdr=read.delim(filename,nrows=i-1,header=FALSE,stringsAsFactors=FALSE)$V1
	
	for(hdrow in 1:length(hdr)){
		if(((grepl("x-axis ",hdr[hdrow]))||(grepl("Control ",hdr[hdrow])))&&(cStart<0)) cStart=hdrow
		if(((grepl("y-axis ",hdr[hdrow]))||(grepl("Query ",hdr[hdrow])))&&(qStart<0)) qStart=hdrow
	}
	
	qfaVersion=NAtoBlank(strsplit(hdr[1],": ")[[1]][2])
	summType=NAtoBlank(strsplit(hdr[2],": ")[[1]][2])
	testType=NAtoBlank(strsplit(hdr[3],": ")[[1]][2])
	
	cTreat=NAtoBlank(strsplit(hdr[cStart+0],": ")[[1]][2])
	cMed=NAtoBlank(strsplit(hdr[cStart+1],": ")[[1]][2])
	cScrID=NAtoBlank(strsplit(hdr[cStart+2],": ")[[1]][2])
	cScrNm=NAtoBlank(strsplit(hdr[cStart+3],": ")[[1]][2])
	cLibs=NAtoBlank(strsplit(hdr[cStart+4],": ")[[1]][2])
	cCli=NAtoBlank(strsplit(hdr[cStart+5],": ")[[1]][2])
	cUse=NAtoBlank(strsplit(hdr[cStart+6],": ")[[1]][2])
	cDate=NAtoBlank(strsplit(hdr[cStart+7],": ")[[1]][2])
	cDate=gsub("'","",cDate)
	cDate=gsub('"',"",cDate)
	cPI=cCond=cInoc=cFit=""
	if(cStart+8<qStart) cPI=NAtoBlank(strsplit(hdr[cStart+8],": ")[[1]][2])
	if(cStart+9<qStart) cCond=NAtoBlank(strsplit(hdr[cStart+9],": ")[[1]][2])
	if(cStart+10<qStart) cInoc=NAtoBlank(strsplit(hdr[cStart+10],": ")[[1]][2])
	if(cStart+11<qStart) cFit=NAtoBlank(strsplit(hdr[cStart+11],": ")[[1]][2])
	
	qTreat=NAtoBlank(strsplit(hdr[qStart+0],": ")[[1]][2])
	qMed=NAtoBlank(strsplit(hdr[qStart+1],": ")[[1]][2])
	qScrID=NAtoBlank(strsplit(hdr[qStart+2],": ")[[1]][2])
	qScrNm=NAtoBlank(strsplit(hdr[qStart+3],": ")[[1]][2])
	qLibs=NAtoBlank(strsplit(hdr[qStart+4],": ")[[1]][2])
	qCli=NAtoBlank(strsplit(hdr[qStart+5],": ")[[1]][2])
	qUse=NAtoBlank(strsplit(hdr[qStart+6],": ")[[1]][2])
	qDate=NAtoBlank(strsplit(hdr[qStart+7],": ")[[1]][2])
	qDate=gsub("'","",qDate)
	qDate=gsub('"',"",qDate)
	qPI=qCond=qInoc=qFit=""
	if(qStart+8<=length(hdr)) qPI=NAtoBlank(strsplit(hdr[qStart+8],": ")[[1]][2])
	if(qStart+9<=length(hdr)) qCond=NAtoBlank(strsplit(hdr[qStart+9],": ")[[1]][2])
	if(qStart+10<=length(hdr)) qInoc=NAtoBlank(strsplit(hdr[qStart+10],": ")[[1]][2])
	if(qStart+11<=length(hdr)) qFit=NAtoBlank(strsplit(hdr[qStart+11],": ")[[1]][2])

	fMax=max(c(res$QueryFitnessSummary,res$ControlFitnessSummary))
	qRep=getMode(res$QueryCount)
	cRep=getMode(res$ControlCount)
	
	return(
	list(res=res,qfaVersion=qfaVersion,summType=summType,testType=testType,
	cTreat=cTreat,cMed=cMed,cScrID=cScrID,cScrNm=cScrNm,cLibs=cLibs,cCli=cCli,cUse=cUse,cDate=cDate,cPI=cPI,cCond=cCond,cInoc=cInoc,cFit=cFit,
	qTreat=qTreat,qMed=qMed,qScrID=qScrID,qScrNm=qScrNm,qLibs=qLibs,qCli=qCli,qUse=qUse,qDate=qDate,qPI=qPI,qCond=qCond,qInoc=qInoc,qFit=qFit,
	fMax=fMax,rept=hdr,qRep=qRep,cRep=cRep)
	)
}

buildGO<-function(){
	# Note we do not recommend that users edit this particular set of functionally related genes
	fname=file.path(system.file(package = "qfa"),"extdata","GOAnnotation.txt.gz")
	cat("\nGrouping of genes by GO terms specified in this text file (note: must be uncompressed before viewing with text editor)\n~~~~~~~~~~~~~~~\n")
	cat(paste(fname,"\n"))
	GO=read.delim(fname,sep="\t",header=TRUE,stringsAsFactors=FALSE)
	return(GO)
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
		if (gene%in%ORFGENE$ORF) {
			orfs=c(orfs,gene)
			genes=c(genes,ORFGENE$Gene[ORFGENE$ORF==gene])			
		} else {
			if(gene%in%ORFGENE$Gene){
				# Standard gene names are ambiguous
				# In this scenario, they do not unambiguously specify species (e.g. RIF1 in cerevisiae and pombe)
				orfs=c(orfs,ORFGENE$ORF[ORFGENE$Gene==gene])
				genes=c(genes,gene)
			}
		}
	}
	x[["Genes"]]=genes
	x[["ORFs"]]=orfs
	x[["Label"]]=lab
	print(x)
	return(x)
}

drawSel=function(selx,sely){
	if(length(selx)==1) {
		points(selx,sely,col="black",pch=16,cex=0.4)
	}else{
		points(selx,sely,col="black",lty=2,type="l")
	}
}

printSelected=function(globs){
	cat("~~~~~~~~~~~~~~~\n")
	cat("\n")
	cat("List of genes currently selected:\n")
	cat(globs$dat$Gene[globs$targs])
	cat("\n")
	cat("~~~~~~~~~~~~~~~\n")
	cat("\n")
	cat("Systematic names for genes currently selected:\n")
	cat(globs$dat$ORF[globs$targs])
	cat("\n")
	cat("\n")
	cat("\n")
}

reportExpts=function(globs,fname="MetaReport.txt"){
	if(fname!="") {
		cat("\nDetailed plot metadata will be written to file:\n~~~~~~~~~~~~~~~\n")
		cat(file.path(getwd(),fname),"\n")
		sink(fname)
	}
	metaSumm=data.frame(PlotNo=numeric(),x.axis=character(),y.axis=character(),stringsAsFactors=FALSE)
	for (datno in 1:length(globs$datlist)){
		metaSumm[datno,]=c(datno,paste(globs$datlist[[datno]]$cScrNm,globs$datlist[[datno]]$cTreat,globs$datlist[[datno]]$cCond,globs$datlist[[datno]]$cFit),paste(globs$datlist[[datno]]$qScrNm,globs$datlist[[datno]]$qTreat,globs$datlist[[datno]]$qCond,globs$datlist[[datno]]$qFit))
		# Print data to console
		cat("\nExperimental metadata: plot",datno,"\n~~~~~~~~~~~~~~~\n")
		# Split report string and patch in replicate numbers for Query and Control
		rept=globs$datlist[[datno]]$rept
		fcont=grep("Control",rept)
		if (length(fcont)==0) fcont=grep("x-axis",rept)
		findx=tail(fcont,1)
		rept2=c(rept[1:findx],paste("x-axis number of replicates:",globs$datlist[[datno]]$cRep),rept[(findx+1):length(rept)],paste("y-axis number of replicates:",globs$datlist[[datno]]$qRep))
		strt=4
		r2=rept2[strt:length(rept2)]
		fin=length(r2)
		reptCont=r2[1:(fin/2)]
		reptCont=gsub("Control ","",reptCont)
		reptQuer=r2[((fin/2)+1):fin]
		reptQuer=gsub("Query ","",reptQuer)
		df=data.frame(Control=reptCont,Query=reptQuer)
		#print(df,row.names=FALSE)
		cat(rept2,sep="\n")
	}
	if(fname!="") sink()
	
	cat("\nExperimental metadata summary:\n~~~~~~~~~~~~~~~\n")
	flist=strsplit(fname,"\\.")[[1]]
	froot=head(flist,1)
	fext=tail(flist,1)
	froot2=paste(froot,"Summary",sep="")
	fname2=paste(froot2,fext,sep=".")
	print(metaSumm,row.names=FALSE)
	cat("\nExperimental metadata summary written to file:\n~~~~~~~~~~~~~~~\n")
	cat(fname2,"\n")
	write.table(metaSumm,file=fname2,sep="\t",row.names=FALSE,quote=FALSE)
}