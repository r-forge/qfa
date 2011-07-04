fitnessReport<-function(treatment,outputfile,dataframe){
	# Summarises mean and median fitnesses for all orfs in an .fit object
	print(outputfile)
	fitdf=dataframe[dataframe$Treatment==treatment,]

	genelst=unique(as.character(fitdf$Gene))
	genelst=sort(genelst)
	orflst=as.character(fitdf$ORF[match(genelst,as.character(fitdf$Gene))])

	# Get median fitness and variance for each gene
	calcStuff <- function(queryORF){
		res=fitdf[as.character(fitdf$ORF)==queryORF,]
		res=res[!is.na(res$fit),]
		return(c(median(res$fit),mean(res$fit),var(res$fit),length(res$fit),sd(res$fit)/sqrt(length(res$fit))))
	}

	results=sapply(orflst,calcStuff)
	results=data.frame(t(results))
	results=cbind(genelst,orflst,results,rep(treatment,length(genelst)))
	colnames(results)=c("Gene","ORF","MedianFit","MeanFit","VarianceFit","NumRepeats","SEFit","Treatment")

	write.table(results,file=outputfile,sep="\t",quote=FALSE,row.names=FALSE)
	return(results)
}

correlationReport<-function(screennames,dataframe,outputfile,aw=4,ah=4,fitmax=130){
	# Tests the correlation of all possible repeats of the same plate/treatment
	# Useful tool for searching for incorrect plate orientation, or misplaced/mislabelled plates
	# Can also give clues about plates with incorrect medium
	
	# This might fail for various reasons
	# But probably not a critical part of any workflow, so wrap in "try"
	try({
	# Generate correlation plot for every possible 2-way comparison
	poss=combn(scrnms,2)
	trts=unique(dataframe$Treatment)
	pnum=max(as.numeric(dataframe$MasterPlate.Number))
	# Adjust aw, ah to get all plates on one page (e.g. 4*4 gives 16 panels, good for 15 plate lib)
	pdummy=aw*ah-pnum

	corrs=c()
	pdf(outputfile)
	op<-par(mfrow=c(aw,ah),cex.main=0.8,mar=c(1,1,1,1),mgp=c(0,0,0))

	for(t in trts){
	for(comb in 1:length(poss[1,])){
	for(p in 1:pnum){
		rep1=poss[1,comb]; rep2=poss[2,comb]
		r1=dataframe[(dataframe$Treatment==t)&(dataframe$Screen.Name==rep1)&(dataframe$MasterPlate.Number==p),]
		r2=dataframe[(dataframe$Treatment==t)&(dataframe$Screen.Name==rep2)&(dataframe$MasterPlate.Number==p),]
		r1=r1[order(r1$Screen.Name,r1$MasterPlate.Number,r1$Gene),]
		r2=r2[order(r2$Screen.Name,r2$MasterPlate.Number,r2$Gene),]
		
		cols=rainbow(max(as.numeric(dataframe$MasterPlate.Number),na.rm=TRUE))
		correlate=cor(r1$fit,r2$fit)
		corrs=rbind(corrs,c(rep1,rep2,t,p,correlate))
		ptitle=paste("Trt:",t,"Plate:",p,"Corr:",formatC(correlate,4))
		#print(paste(rep1,rep2,ptitle))
		plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab=rep1,ylab=rep2,main=ptitle,axes=FALSE)
		abline(0,1,lwd=3,col="grey")
		#text(r1$fit,r2$fit,r1$Gene,col="black",pos=4,offset=0.1,cex=0.4)
		points(r1$fit,r2$fit,col=cols[r1$MasterPlate.Number],pch=16,cex=0.4)
		#print(c(sum(r1$Gene==r2$Gene),length(r1$Gene)))
	}
	for(p in 1:pdummy) plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab="",ylab="",main="",axes=FALSE)
	}
	}
	par(op)
	corrs=as.data.frame(corrs,stringsAsFactors=FALSE)
	colnames(corrs)=c("rep1","rep2","trt","p","correlate")
	corrs$p=as.numeric(corrs$p)
	corrs$correlate=as.numeric(corrs$correlate)
	hist(corrs$correlate,xlim=c(0,1),xlab="Correlation Coefficient",ylab="Frequency",main=cfolder)
	dev.off()
	})
}

plateBoxplots<-function(dataframe,outputfile,fitmax=130){
# Generates plate by plate boxplots of culture fitnesses
# Useful for identifying plates with medium problems
pdf(outputfile)
for (trt in unique(dataframe$Treatment)){
	for (scr in unique(dataframe$Screen.Name)){
		dt=dataframe[(dataframe$Screen.Name==scr)&(dataframe$Treatment==trt),]
		dt$MasterPlate.Number=as.numeric(dt$MasterPlate.Number)
		dt=dt[order(dt$MasterPlate.Number),]
		dt$MasterPlate.Number=as.factor(dt$MasterPlate.Number)
		boxplot(dt$fit~dt$MasterPlate.Number,notch=TRUE,xlab="Plate Number",ylab="Fitness",col=rainbow(23),main=paste(trt,scr),ylim=c(0,fitmax),cex.axis=0.5)
}
}
dev.off()
}

