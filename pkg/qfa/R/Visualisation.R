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

targFun=function(x,y,dat){
	d=(x-dat$ControlFitnessSummary)^2+(y-dat$QueryFitnessSummary)^2
	return(which.min(d))
}

makePlot=function(datno,...){
	if(compno>0){
		compnm=paste(GROUPS$GroupName[compno],"\t",GROUPS$GroupID[compno])
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
		lst=GROUPS$GroupORFs[compno][[1]]
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
				makePlot(datno,ecol=Ecol,scol=Scol)
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
				makePlot(datno,ecol=Ecol,scol=Scol)
				if(length(selx)>0) drawSel(selx,sely)
			}else{
				selx<<-c(selx,xf)
				sely<<-c(sely,yf)
				makePlot(datno,ecol=Ecol,scol=Scol)
				drawSel(selx,sely)
			}
		}
	}
	if(1%in%buttons) {
		# Delete the last gene which was highlighted
		targs<<-targs[-length(targs)]
		posits<<-posits[-length(posits)]
		makePlot(datno,ecol=Ecol,scol=Scol)
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
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	if(key=="d"){
		# Unhighlight the last gene which was highlighted
		targs<<-targs[-length(targs)]
		posits<<-posits[-length(posits)]
		makePlot(datno,ecol=Ecol,scol=Scol)
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
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	if(key=="Left") {
		orftargs=dat$ORF[targs]
		datno<<-((datno-2)%%length(datlist))+1
		dat<<-datlist[[datno]]$res
		targs<<-match(orftargs,dat$ORF)
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	if(key=="Up") {
		compno<<-compno+1
		if(compno>length(GROUPS$GroupORFs)) compno<<-1
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	if(key=="Down") {
		compno<<-compno-1
		if(compno<=0) compno<<-length(GROUPS$GroupORFs)
		makePlot(datno,ecol=Ecol,scol=Scol)
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
			makePlot(datno,ecol=Ecol,scol=Scol)
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
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	if(key=="p") {
		pdf(sprintf("QFAVisualisation%04d.pdf",plotno))
			makePlot(datno,ecol=Ecol,scol=Scol)
			plotno<<-plotno+1
		dev.off()
	}
	if(key=="t") {
		ENH=Ecol
		SUP=Scol
		Scol<<-ENH
		Ecol<<-SUP
		makePlot(datno,ecol=Ecol,scol=Scol)
	}
	if(key=="r") {
	# Zoom mode (click on TL and BR)
		zooming<<-TRUE
		zoomTL<<-c(-99,-99)
		zoomBR<<-c(-99,-99)
	}
	if(key=="u") {
	# Get user input for a new list of genes
		newgrp=getText(ORFGENE)
		newdf=data.frame(GroupName=newgrp[["Label"]],GroupID="R VisTool")
		newdf$GroupORFs=list(newgrp[["ORFs"]])
		GROUPS<<-rbind(newdf,GROUPS)
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

buildBenschop<-function(){
	# Read in functionally related complexes
	# Largely from Benschopp Mol Cell 2010
	# Can add some complexes manually
	#compfile=system.file("/FunctionalComplexes.txt", package = "qfa")
	compfile=paste(system.file(package = "qfa"),"/FunctionalComplexes.txt",sep="")
	Benschopp=read.delim(compfile,stringsAsFactors=FALSE,sep="\t")
	Benschopp$CompList=strsplit(Benschopp$Complex.members..systematic.name,"; ")
	res=data.frame(GroupName=Benschopp$X..Complex.name,GroupID="Manual or Benschop Mol. Cell 2010")
	res$GroupORFs=Benschopp$CompList
	return(res)
}

buildGO<-function(){
	load(file=paste(system.file(package = "qfa"),"/GOAnnotation.Rda",sep=""))
	return(tooldf)
}

visToolDemo<-function(groupFun=buildBenschop){
	groups=groupFun()
	
	orfile=paste(system.file(package = "qfa"),"/ORF2GENE.txt",sep="")
	ORFGENE=read.delim(orfile,stringsAsFactors=FALSE,sep="\t",header=FALSE)
	colnames(ORFGENE)=c("ORF","Gene")
	ORFGENE=ORFGENE[!duplicated(ORFGENE$ORF),]

	# Read in GIS files
	flist=c(
	"CDC13-1_20GIS.txt","CDC13-1_27GIS.txt","CDC13-1_36GIS.txt",
	"YKU70_23GIS.txt","YKU70_30GIS.txt","YKU70_37GIS.txt","YKU70_37.5GIS.txt")
	filenames=c()
	for (f in flist){
		sysfile=paste(system.file(package = "qfa"),"/",f,sep="")
		filenames=c(filenames,sysfile)
	}
	visTool(groups,ORFGENE,filenames)
}

visTool<-function(groups,orf2gene,GISfiles){
	# Function for generating interactive plot
	GROUPS<<-groups
	ORFGENE<<-orf2gene
	datlist<<-list()
	for (f in 1:length(GISfiles)){	
		report=getResults(GISfiles[f])
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
	
	# Toggling between colour of suppressors and enhancers
	Ecol<<-"green"
	Scol<<-"red"

	print("Windows mouse")
	print("~~~~~~~~~~~~~~~")
	print("Left click: Highlight gene/Rotate text position")
	print("Right click: SGD (or press 'w' on keyboard)")
	print("Middle click: Remove last gene (or press 'd' on keyboard)")
	print("Mac mouse")
	print("~~~~~~~~~~~~~~~")
	print("Click: Highlight gene/Rotate text position")
	print("Keyboard")
	print("~~~~~~~~~~~~~~~")
	print("Left/Right arrow: change plot")
	print("Up/Down arrow: change group highlighted")
	print("u: add new group of genes to list of groups")
	print("z: select tool (toggle on and off)") 
	print("s: add selection") 
	print("c: clear selection") 
	print("w: open last gene highlighted in SGD")
	print("d: unhighlight last gene highlighted")
	print("t: toggle colours indicating positive and negative interaction")
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
	# Note that event handling doesn't work under Linux :(
	sysinf=Sys.info()
	#if( sysinf["sysname"]!="Windows") x11()
	#if( sysinf["sysname"]=="Linux") Cairo()
	x11()
	makePlot(datno,ecol=Ecol,scol=Scol)
	getGraphicsEvent(prompt="L click: Highlight/Rotate, R click: SGD, M click: Remove, Left/Right: Change plot, z: select tool, s: add selection, c: clear, q: quit", onMouseDown=mouse, onKeybd=keybd)
	print(dat$Gene[targs])
	#if( sysinf["sysname"]!="Windows") dev.off()
	dev.off()
}

getText=function(ORFGENE){
	x=list()
	cat("Input a list of gene names, separated by spaces & press enter:")
	usergenes=scan(what=character(),nlines=1)
	usergenes=toupper(usergenes)
	cat("Please input a label for your list of genes (no spaces) & press enter:")
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




