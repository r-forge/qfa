showDemo<-function(demoname="telomereCap"){
	# Function to just display text in demo
	file.show(system.file("demo", paste(demoname,".R",sep=""), package = "qfa"))
}

findFit<-function(df){
	# Find which column corresponds to the "fit" column in df
	if(!"fit"%in%colnames(df)) return("NoFitCol")
	cnames=colnames(df)[colnames(df)!="fit"]
	for(colname in cnames){
		compare<-df$fit==df[colname]
		if((!is.na(sum(compare)))&&(sum(compare)==length(df$fit))) return(colname)
	}
}
	
fitnessReport<-function(grp,outputfile,dataframe,groupcol="Treatment",stripORFs=c(),fdef="fit"){
	# Eliminate spurious precision to make smaller files
	#results[,3:9]=signif(results[,3:9],6)
	# Summarises mean and median fitnesses for all orfs in an .fit object
	print(outputfile)
	fitdf=dataframe[dataframe[[groupcol]]==grp,]
	if("ExptDate"%in%colnames(fitdf)){
		fitdf$ExptDate=gsub("'","",fitdf$ExptDate)
		fitdf$ExptDate=gsub('"',"",fitdf$ExptDate)
	}

	orflst=unique(as.character(fitdf$ORF))
	orflst=sort(orflst)
	if(!is.null(stripORFs)) fitdf=fitdf[!fitdf$ORF%in%stripORFs,]

	#meanRes=by(fitdf[,27:31],fitdf$ORF,median)
	medRes=tapply(fitdf[[fdef]],fitdf$ORF,median,na.rm=TRUE)
	meanRes=tapply(fitdf[[fdef]],fitdf$ORF,mean,na.rm=TRUE)
	varRes=tapply(fitdf[[fdef]],fitdf$ORF,var,na.rm=TRUE)
	numRes=tapply(fitdf[[fdef]],fitdf$ORF,length)
	seRes=tapply(fitdf[[fdef]],fitdf$ORF,sd,na.rm=TRUE)/sqrt(numRes)
	orflst=names(medRes)
	genelst=as.character(fitdf$Gene[match(orflst,as.character(fitdf$ORF))])

	results=data.frame(Gene=genelst,ORF=orflst,MedianFit=medRes,MeanFit=meanRes,VarianceFit=varRes,NumRepeats=numRes,SEFit=seRes)
	
	# Automatically extract qfa package version number
	sinfo=sessionInfo()
	vno=sinfo$otherPkgs$qfa$Version[1]
	QFAversion=paste("R package version:",vno)
	rTreat=paste("Treatment:",paste(unique(fitdf$Treatment),collapse=" "))
	rMed=paste("Medium:",paste(unique(fitdf$Medium),collapse=" "))
	rScrID=paste("Screen ID:",paste(unique(fitdf$ScreenID),collapse=" "))
	rGen=paste("Screen name:",paste(unique(fitdf$Screen.Name),collapse=" "))
	rLib=paste("Libraries:",paste(unique(fitdf$Library.Name),collapse=" "))
	rClient=paste("Client:",paste(unique(fitdf$Client),collapse=" "))
	rUser=paste("User:",paste(unique(fitdf$User),collapse=" "))
	rPI=paste("PI:",paste(unique(fitdf$PI),collapse=" "))
	rDate=paste("Date:",paste(unique(fitdf$ExptDate),collapse=" "))
	if(fdef=="fit") {rFit=paste("Fitness definition:",findFit(fitdf))}else{rFit=paste("Fitness definition:",fdef)}
	rCond=paste("Condition:",paste(unique(fitdf$Condition),collapse=" "))
	spacer="#########################################################"
	header=c(QFAversion,
	rTreat,rMed,rScrID,rGen,rLib,rClient,rUser,rPI,rDate,rFit,rCond,
	spacer)
	tmpfname=paste(outputfile,"tmp",sep="_")
	write.table(header,outputfile,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	write.table(results,tmpfname,sep="\t",quote=FALSE,row.names=FALSE)
	file.append(outputfile,tmpfname)
	file.remove(tmpfname)	
	return(results)
}

report.epi<-function(results,filename){
	# Eliminate spurious precision to make smaller files
	results[,3:9]=signif(results[,3:9],6)
	# Automatically extract qfa package version number
	sinfo=sessionInfo()
	vno=sinfo$otherPkgs$qfa$Version[1]

	QFAversion=paste("R package version:",vno)
	sumType=paste("Summary type:",paste(unique(results$SummaryType),collapse=" "))
	testType=paste("Test type:",paste(unique(results$TestType),collapse=" "))
	
	cTreat=paste("x-axis treatment:",paste(unique(results$cTreat),collapse=" "))
	cMed=paste("x-axis medium:",paste(unique(results$cMed),collapse=" "))
	cScrID=paste("x-axis screen ID:",paste(unique(results$cScrID),collapse=" "))
	cGen=paste("x-axis screen name:",paste(unique(results$cGen),collapse=" "))
	cLib=paste("x-axis libraries:",paste(unique(results$cLib),collapse=" "))
	cClient=paste("x-axis client:",paste(unique(results$cClient),collapse=" "))
	cUser=paste("x-axis user:",paste(unique(results$cUser),collapse=" "))
	cDate=paste("x-axis date:",paste(unique(results$cDate),collapse=" "))
	cPI=paste("x-axis PI:",paste(unique(results$cPI),collapse=" "))
	cCond=paste("x-axis condition:",paste(unique(results$cCond),collapse=" "))
	cInoc=paste("x-axis inoculation type:",paste(unique(results$cInoc),collapse=" "))
	cFit=paste("x-axis fitness definition:",paste(unique(results$ControlFit),collapse=" "))
	
	qTreat=paste("y-axis treatment:",paste(unique(results$qTreat),collapse=" "))
	qMed=paste("y-axis medium:",paste(unique(results$qMed),collapse=" "))
	qScrID=paste("y-axis screen ID:",paste(unique(results$qScrID),collapse=" "))
	qGen=paste("y-axis screen name:",paste(unique(results$qGen),collapse=" "))
	qLib=paste("y-axis libraries:",paste(unique(results$qLib),collapse=" "))
	qClient=paste("y-axis client:",paste(unique(results$qClient),collapse=" "))
	qUser=paste("y-axis user:",paste(unique(results$qUser),collapse=" "))
	qDate=paste("y-axis date:",paste(unique(results$qDate),collapse=" "))
	qPI=paste("y-axis PI:",paste(unique(results$qPI),collapse=" "))
	qCond=paste("y-axis condition:",paste(unique(results$qCond),collapse=" "))
	qInoc=paste("y-axis inoculation type:",paste(unique(results$qInoc),collapse=" "))
	qFit=paste("y-axis fitness definition:",paste(unique(results$QueryFit),collapse=" "))
	
	spacer="#########################################################"
	
	header=c(QFAversion,sumType,testType,
	cTreat,cMed,cScrID,cGen,cLib,cClient,cUser,cDate,cPI,cCond,cInoc,cFit,
	qTreat,qMed,qScrID,qGen,qLib,qClient,qUser,qDate,qPI,qCond,qInoc,qFit,
	spacer)
	write.table(header,filename,sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
	results$SummaryType=results$testType=NULL
	results$cTreat=results$cMed=results$cScrID=results$cGen=results$cLib=results$cClient=results$cUser=results$cDate=results$cPI=results$cCond=results$cInoc=results$ControlFit=NULL
	results$qTreat=results$qMed=results$qScrID=results$qGen=results$qLib=results$qClient=results$qUser=results$qDate=results$qPI=results$qCond=results$qInoc=results$QueryFit=NULL
	tmpfname=paste(filename,"tmp",sep="_")
	write.table(results,tmpfname,sep="\t",quote=FALSE,row.names=FALSE)
	file.append(filename,tmpfname)
	file.remove(tmpfname)	
}

correlationReport<-function(scrnms,dataframe,outputfile,aw=4,ah=4,fitmax=185,fdef="fit"){
	# Tests the correlation of all possible repeats of the same plate/treatment
	# Useful tool for searching for incorrect plate orientation, or misplaced/mislabelled plates
	# Can also give clues about plates with incorrect medium

	# This might fail for various reasons
	# But probably not a critical part of any workflow, so wrap in "try"
	try({
		dataframe=dataframe[dataframe$Screen.Name%in%scrnms,]
		# Generate correlation plot for every possible 2-way comparison
		dataframe$TrtMed=paste(dataframe$Treatment,dataframe$Medium)
		trts=unique(as.character(dataframe$TrtMed))
		plates=unique(as.numeric(as.character(dataframe$MasterPlate.Number)))
		plates=plates[order(plates)]
		pnum=length(plates)
		pdf(outputfile)
		corrs=c()
		op<-par(mfrow=c(aw,ah),cex.main=0.8,mar=c(1,1,1,1),mgp=c(0,0,0))

		for(t in trts){
			for(p in plates){
				reps=dataframe[(as.character(dataframe$TrtMed)==t)&(as.numeric(as.character(dataframe$MasterPlate.Number))==p),]
				reps=reps[order(reps$Screen.Name,reps$MasterPlate.Number,reps$Gene),]
				barcs=unique(reps$Barcode)
				for (bpair in data.frame(combn(barcs,2),stringsAsFactors=FALSE)){
					r1=reps[reps$Barcode==bpair[1],]
					r2=reps[reps$Barcode==bpair[2],]

					if((length(r1[[fdef]])==length(r2[[fdef]]))&(length(r1[[fdef]])>0)){
						cols=rainbow(max(as.numeric(dataframe$MasterPlate.Number),na.rm=TRUE))
						correlate=cor(r1[[fdef]],r2[[fdef]])
						corrs=rbind(corrs,c(bpair[1],bpair[2],t,p,correlate))
						ptitle=paste("Plate:",p,"Corr:",formatC(correlate,4))
						plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab=bpair[1],ylab=bpair[2],main=ptitle,axes=FALSE,cex.main=1.2)
						points(r1[[fdef]],r2[[fdef]],col=cols[r1$MasterPlate.Number],pch=16,cex=0.4)
					}else{
						plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab=bpair[1],ylab=bpair[2],main=paste("Trt:",t,"Plate:",p),axes=FALSE)
					}
				}
			}
		}
		par(op)

		corrs=as.data.frame(corrs,stringsAsFactors=FALSE)
		colnames(corrs)=c("rep1","rep2","trt","p","correlate")
		corrs$p=as.numeric(corrs$p)
		corrs$correlate=as.numeric(corrs$correlate)
		hist(corrs$correlate,xlim=c(0,1),xlab="Correlation Coefficient",ylab="Frequency",main=outputfile)
	})
	dev.off()
}

plateBoxplots<-function(dataframe,outputfile,fitmax=185,groupcol="Treatment",fdef="fit"){
	# Generates plate by plate boxplots of culture fitnesses
	# Useful for identifying plates with medium problems
	dataframe$ScreenRepQuad=paste(dataframe$Screen.Name,"RQ",sprintf("%02d",dataframe$RepQuad),sep="_")
	grps=unique(dataframe[[groupcol]])
	scrns=unique(dataframe$ScreenRepQuad)

	grps=grps[order(grps)]
	scrns=scrns[order(scrns)]
	pdf(outputfile,width=12,height=8)
	for (grp in grps){
		for (scr in scrns){
			dt=dataframe[(dataframe$ScreenRepQuad==scr)&(dataframe[[groupcol]]==grp),]
			dt=dt[order(as.numeric(as.character(dt$MasterPlate.Number))),]
			plateNums=unique(as.numeric(as.character(dt$MasterPlate.Number)))
			plateNums=as.character(plateNums)
			dt$MasterPlate.Number<-factor(dt$MasterPlate.Number,levels=plateNums)
			if(length(dt[[fdef]])>0) boxplot(dt[[fdef]]~dt$MasterPlate.Number,notch=TRUE,xlab="Plate Number",ylab="Fitness",col=rainbow(23),main=paste(grp,scr),ylim=c(0,fitmax),cex.axis=1)
		}
	}
	if(fdef=="fit"){fitdef=findFit(dataframe)}else{fitdef=fdef}
	by(dataframe,dataframe$Barcode,summaryPlots,fdef=fitdef,fdef2="MDRMDP",label="")
	dev.off()
}

summaryPlots=function(datFit,fdef="nAUC",fdef2="MDRMDP",label=""){
	datFit$Quad=(((datFit$Row-1)%%2)*2+(datFit$Column-1)%%2)+1
	bc=unique(datFit$Barcode)
	print(paste("Summary plot for",bc))
	op=par(mfrow=c(2,3))
	hist(datFit[[fdef]],breaks=100,xlab=fdef,main="",freq=TRUE, cex=1.75)
	boxplot(datFit[[fdef]]~datFit$Row,xlab="Row",ylab=fdef, cex=1.75)
	boxplot(datFit[[fdef]]~datFit$Column,xlab="Column",ylab=fdef, cex=1.75)
	# Only really useful for plates with several replicate of a few genotypes
	#boxplot(datFit[[fdef]]~datFit$Gene,las=2,xlab="",ylab=fdef,cex.axis=0.5, cex=1.75) 
	boxplot(datFit[["d0"]]~datFit$Quad,notch=TRUE,xlab="Quadrant",ylab="Inoculum Density", cex=1.75)
	boxplot(datFit[[fdef]]~datFit$Quad,notch=TRUE,xlab="Quadrant",ylab=fdef, cex=1.75)
	boxplot(datFit[[fdef2]]~datFit$Quad,notch=TRUE,xlab="Quadrant",ylab=fdef2, cex=1.75)
	title(paste(bc,label))
	par(op)
}