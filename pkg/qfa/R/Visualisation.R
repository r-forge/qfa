require(sp)

# Plotting function
epiplot<-function(results,qthresh,fitratio=FALSE,ref.orf="YOR202W",xxlab="Control Fitness",yylab="Query Fitness",mmain="Epistasis Plot",ymin=0,ymax=0,xmin=0,xmax=0){
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
		points(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,col='green',pch=19,cex=0.5)
		#text(enhancers$ControlFitnessSummary,enhancers$QueryFitnessSummary,enhancers$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
	if(length(suppressors$ORF)>0){
		points(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,col='red',pch=19,cex=0.5)
		#text(suppressors$ControlFitnessSummary,suppressors$QueryFitnessSummary,suppressors$Gene,col=1,pos=4,offset=0.1,cex=0.4)
	}
}

targFun=function(x,y,dat){
	d=(x-dat$ControlFitnessSummary)^2+(y-dat$QueryFitnessSummary)^2
	return(which.min(d))
}

makePlot=function(datno,...){
	if(compno>0){
		compnm=paste(Benschopp$X..Complex.name[compno],"\t",Benschopp$comments[compno])
	}else{
		compnm=""
	}
	maintitle=paste(datlist[[datno]]$datname,compnm,sep="\n")
	if((nchar(datlist[[datno]]$cMed)<10)&(nchar(datlist[[datno]]$qMed)<10)){
		xlab=paste(datlist[[datno]]$cBack,datlist[[datno]]$cTreat,datlist[[datno]]$cMed)
		ylab=paste(datlist[[datno]]$qBack,datlist[[datno]]$qTreat,datlist[[datno]]$qMed)
	}else{
		xlab=paste(datlist[[datno]]$cBack,datlist[[datno]]$cTreat)
		ylab=paste(datlist[[datno]]$qBack,datlist[[datno]]$qTreat)	
	}
	dat=datlist[[datno]]$res
	epiplot(dat,0.05,mmain=maintitle,xxlab=xlab,yylab=ylab,ymin=fymin[datno],ymax=fymax[datno],xmin=fxmin[datno],xmax=fxmax[datno],...)
	if(compno>0){
		lst=Benschopp$CompList[compno][[1]]
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
				makePlot(datno)
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
				makePlot(datno)
				if(length(selx)>0) drawSel(selx,sely)
			}else{
				selx<<-c(selx,xf)
				sely<<-c(sely,yf)
				makePlot(datno)
				drawSel(selx,sely)
			}
		}
	}
	if(1%in%buttons) {
		# Delete the last gene which was highlighted
		targs<<-targs[-length(targs)]
		posits<<-posits[-length(posits)]
		makePlot(datno)
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
		makePlot(datno)
	}
	if(key=="d"){
		# Unhighlight the last gene which was highlighted
		targs<<-targs[-length(targs)]
		posits<<-posits[-length(posits)]
		makePlot(datno)
	}
	if(key=="w"){
		targ<<-targs[length(targs)]
		browseURL(paste("http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=",dat$ORF[targ],sep=""))
	}
	if(key=="Right") {
		orftargs=dat$ORF[targs]
		datno<<-(datno%%length(datlist))+1
		dat<<-datlist[[datno]]$res
		targs<<-match(orftargs,dat$ORF)
		makePlot(datno)
	}
	if(key=="Left") {
		orftargs=dat$ORF[targs]
		datno<<-((datno-2)%%length(datlist))+1
		dat<<-datlist[[datno]]$res
		targs<<-match(orftargs,dat$ORF)
		makePlot(datno)
	}
	if(key=="Up") {
		compno<<-compno+1
		if(compno>length(Benschopp$CompList)) compno<<-1
		makePlot(datno)
	}
	if(key=="Down") {
		compno<<-compno-1
		if(compno<=0) compno<<-length(Benschopp$CompList)
		makePlot(datno)
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
			makePlot(datno)
			drawSel(selx,sely)
			sel<<-0
		}
	}
	if((key=="s")&(length(selx)>0)){
		ind=point.in.polygon(dat$ControlFitnessSummary,dat$QueryFitnessSummary,selx,sely)
		tsel=which(ind==1)
		targs<<-c(targs,tsel)
		posits<<-c(posits,rep(1,length(tsel)))
		selx<<-c()
		sely<<-c()
		makePlot(datno)
	}
	if(key=="p") {
		pdf(sprintf("QFAVisualisation%04d.pdf",plotno))
			makePlot(datno)
			plotno<<-plotno+1
		dev.off()
	}
	if(key=="r") {
	# Zoom mode (click on TL and BR)
		zooming<<-TRUE
		zoomTL<<-c(-99,-99)
		zoomBR<<-c(-99,-99)
	}
	return(NULL)
}

getResults<-function(filename){
	res=read.delim(filename,skip=10,header=TRUE,stringsAsFactors=FALSE)
	hdr=read.delim(filename,nrows=10,header=FALSE,stringsAsFactors=FALSE)
	qfaVersion=strsplit(hdr[1,1],": ")[[1]][2]
	summType=strsplit(hdr[2,1],": ")[[1]][2]
	testType=strsplit(hdr[3,1],": ")[[1]][2]
	cTreat=strsplit(hdr[4,1],": ")[[1]][2]
	cMed=strsplit(hdr[5,1],": ")[[1]][2]
	cBack=strsplit(hdr[6,1],": ")[[1]][2]
	qTreat=strsplit(hdr[7,1],": ")[[1]][2]
	qMed=strsplit(hdr[8,1],": ")[[1]][2]
	qBack=strsplit(hdr[9,1],": ")[[1]][2]
	fMax=max(c(res$QueryFitnessSummary,res$ControlFitnessSummary))
	return(
	list(res=res,qfaVersion=qfaVersion,summType=summType,testType=testType,cTreat=cTreat,
	cMed=cMed,cBack=cBack,qTreat=qTreat,qMed=qMed,qBack=qBack,fMax=fMax)
	)
}

visTool<-function(){
	# Read in functionally related complex
	# Largely from Benschopp Mol Cell 2010
	# Can add some complexes manually
	#compfile=system.file("/FunctionalComplexes.txt", package = "qfa")
	compfile=paste(system.file(package = "qfa"),"/FunctionalComplexes.txt",sep="")
	Benschopp<<-read.delim(compfile,stringsAsFactors=FALSE,sep="\t")
	colnames(Benschopp)
	Benschopp$CompList<<-strsplit(Benschopp$Complex.members..systematic.name,"; ")

	# Read in GIS files
	datlist<<-list()
	filenames<<-c(
	"CDC13-1_20GIS.txt","CDC13-1_27GIS.txt","CDC13-1_36GIS.txt",
	"YKU70_23GIS.txt","YKU70_30GIS.txt","YKU70_37GIS.txt","YKU70_37.5GIS.txt")
	for (f in 1:length(filenames)){
		#sysfile=system.file(paste("/",filenames[f],sep=""), package = "qfa")
		sysfile=paste(system.file(package = "qfa"),"/",filenames[f],sep="")
		report=getResults(sysfile)
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

	print("Windows mouse")
	print("~~~~~~~~~~~~~~~")
	print("L click: Highlight gene/Rotate text position")
	print("R click: SGD (or press 'w' on keyboard)")
	print("M click: Remove last gene (or press 'm' on keyboard)")
	print("Mac mouse")
	print("~~~~~~~~~~~~~~~")
	print("Click: Highlight gene/Rotate text position")
	print("Keyboard")
	print("~~~~~~~~~~~~~~~")
	print("Left/Right arrow: Change plot")
	print("z: select tool (toggle on and off)") 
	print("s: add selection") 
	print("c: clear selection") 
	print("w: open last gene highlighted in SGD")
	print("d: unhighlight last gene highlighted")
	print("r: begin zoom (now click on top left and bottom right of new zoomed plot)")
	print("p: print current plot to QFAVisualisation.pdf")
	print("q: quit")

	datno<<-1
	plotno<<-1
	compno<<-0
	pos<<-1
	zooming<<-FALSE
	zoomTL<<-c(-99,-99)
	zoomBR<<-c(-99,-99)

	dat<<-datlist[[datno]]$res

	# Check if running under windows, if not, force X11
	sysinf=Sys.info()
	if( sysinf["sysname"]=="Mac") x11()
	if( sysinf["sysname"]=="Linux") Cairo()
		makePlot(datno)
		getGraphicsEvent(prompt="L click: Highlight/Rotate, R click: SGD, M click: Remove, Left/Right: Change plot, z: select tool, s: add selection, c: clear, q: quit", onMouseDown=mouse, onKeybd=keybd)
		print(dat$Gene[targs])
	if( sysinf["sysname"]%in%c("Linux","Mac")) dev.off()
}
