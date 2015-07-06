readSGD=function(fname="http://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"){
	# Read in SGD features tab file and give columns appropriate headers
	sgdNames=c("SGDID1","FType","FQual","FName","Gene","Alias","Parent","SGDID2","Chr","Start","Stop","Strand","Pos","CoordVersion","SeqVersion","Description")
	# Note I had to strip a stray " from line 5542 to get whole file to read in, or read with quote=""
	sgd=read.delim(fname,header=FALSE,stringsAsFactors=FALSE,col.names=sgdNames,sep="\t",quote="")
	sgd$Mid=(sgd$Start+sgd$Stop)/2
	return(sgd)
}

checkTarg=function(targ,sgd){
	targ=toupper(targ)
	# If targ is not a y-number, try standard gene name or alias
	if(!targ%in%sgd$FName){
		tfinal="MISSING"
		if(targ%in%sgd$Gene){
			tfinal=sgd$FName[sgd$Gene==targ]
		}
		if(targ%in%sgd$Alias){
			tfinal=sgd$FName[sgd$Alias==targ]
		}
		return(tfinal)
	}else{
		return(targ)
	}
}

getNeighbours=function(gvec,kb,sgd,geneCoord="Mid",geneDist="Mid",geneDir="All"){
	# Can specify geneCoord = "Start", "Mid" or "Stop" for the definition of the position of target genes
	# Can specity geneDist = "Start", "Stop", "Mid" or "Min" (smallest of all alternatives) for definition of distance between genes
	# Can specify geneDir = "All", "Up" or "Down"
	geneDist=paste("d",geneDist,sep="")
	gvec=sapply(gvec,checkTarg,sgd)
	res=head(sgd,0)
	res=data.frame()
	
	for(g in gvec){
		targ=sgd[sgd$FName==g,]
		sgdChr=sgd[(sgd$Chr==targ$Chr),]
		sgdChr$dStart=sgdChr$Start-targ[[geneCoord]]
		sgdChr$dStop=sgdChr$Stop-targ[[geneCoord]]
		sgdChr$dMid=sgdChr$Mid-targ[[geneCoord]]
		sgdChr$dMin=sign(sgdChr$dMid)*pmin(abs(sgdChr$dStart),abs(sgdChr$dStop))
		sgdChr$All=abs(sgdChr[[geneDist]])
		sgdChr$Up=if(targ$Stop>targ$Start){sgdChr[[geneDist]]}else{-sgdChr[[geneDist]]}
		sgdChr$Down=if(targ$Stop>targ$Start){-sgdChr[[geneDist]]}else{sgdChr[[geneDist]]}
		near=sgdChr[(sgdChr[[geneDir]]>=0)&(sgdChr[[geneDir]]<kb*1000)&(sgdChr$FType=="ORF"),]
		near$Target_ORF=g
		near$Target_Gene=targ$Gene
		res=rbind(res,near)
	}
	return(res)
}

## DEMO Stripping genes found nearby other genes
# Parse SGD features.tab
# Find chromosome
# Find all genes on chromosome within cutoff kb of target
# Return Y-numbers and gene names
#sgd=readSGD("CommonAuxiliary/SGD_features.tab")
#neighbs=getNeighbours(c("cdc13","lyp1","yku70","ura3"),5,sgd)
#print(neighbs$Gene)
#print(neighbs$FName)

readMer=function(fname="CommonAuxiliary/QFA_Database.mer"){
	mer=read.delim(fname,sep=",",stringsAsFactors=FALSE)
	# Need to tidy up SGA Number column
	mer$SGA.Number=as.numeric(gsub("SGA","",mer$SGA.Number))
	return(mer)
}

getMissingSGA=function(SGAs,threshfrac){
	SGAs=sprintf("SGA%04d",SGAs)
	findFirst=function(SGA,threshfrac){
		# Want YPD_G fitness report, but can't rely on client or temperature (also in file name)
		SGAfiles=list.files(file.path("SGA_EXPERIMENTS",SGA,"ANALYSISOUT"))
		# Want only fitness reports
		SGAfiles=SGAfiles[grepl("FitnessReport",SGAfiles)]
		if(length(SGAfiles)>1){
			# Sometimes YPD_G/YPD_GN, sort and choose lowest
			SGAfile=sort(SGAfiles[grepl("YPD_G",SGAfiles)])[1]
		}else{
			SGAfile=SGAfiles
		}
		print(paste("Identifying dead strains in",SGAfile))
		SGApath=file.path("SGA_EXPERIMENTS",SGA,"ANALYSISOUT",SGAfile)
		f=file(SGApath, "r", blocking = FALSE)
		flines=readLines(f)
		close(f)
		delimiter=which(grepl("#######",flines))[1]
		SGAfit=read.delim(SGApath,sep="\t",stringsAsFactors=FALSE,skip=delimiter)
		SGAfit$MedianScaled=SGAfit$MedianFit/median(SGAfit$MedianFit)
		SGAfit$MeanScaled=SGAfit$MeanFit/mean(SGAfit$MedianFit)
		SGAfit$SGA=SGA
		return(SGAfit[SGAfit$MedianScaled<=threshfrac,])
	}
	reptList=lapply(SGAs,findFirst,threshfrac)
	combined=do.call(rbind,reptList)
	res=combined[!duplicated(combined$ORF),]
	return(res)
}

getMissingGenotypes=function(QFAs,threshfrac,mer){
	QFAnums=as.numeric(substr(QFAs,nchar(QFAs)-3,nchar(QFAs)))
	SGAnums=mer$SGA.Number[mer$Screen.Number%in%QFAnums]
	SGAnums=SGAnums[!is.na(SGAnums)]
	if(length(SGAnums)==0) {
		print("Warning!  No SGAs found in database when searching for dead cultures")
		return(NULL)
	}
	print(paste("Identifying cultures missing in",sprintf("SGA%04d",SGAnums)))
	res=getMissingSGA(SGAnums,threshfrac)
	return(res)
}

## DEMO Stripping genes that were observed to be dead in starting libraries
# Parse .mer file
# Find SGA screen corresponding to QFA
# Find strains whose median fitness lie below a threshold (what about genotypes with many replicate strains?)
# Alternatively find pinIDs whose median fitnesses lie below threshold...
#mer=readMer("CommonAuxiliary/QFA_Database.mer")
#dead=getMissingGenotypes(c(136,60),0.25,mer)
#print(dead$ORF)
#print(dead$Gene)