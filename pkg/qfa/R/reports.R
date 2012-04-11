fitnessReport<-function(treatment,outputfile,dataframe){
	# Summarises mean and median fitnesses for all orfs in an .fit object
	print(outputfile)
	fitdf=dataframe[dataframe$Treatment==treatment,]

	orflst=unique(as.character(fitdf$ORF))
	orflst=sort(orflst)

	#meanRes=by(fitdf[,27:31],fitdf$ORF,median)
	medRes=tapply(fitdf$fit,fitdf$ORF,median,na.rm=TRUE)
	meanRes=tapply(fitdf$fit,fitdf$ORF,mean,na.rm=TRUE)
	varRes=tapply(fitdf$fit,fitdf$ORF,var,na.rm=TRUE)
	numRes=tapply(fitdf$fit,fitdf$ORF,length)
	seRes=tapply(fitdf$fit,fitdf$ORF,sd,na.rm=TRUE)/sqrt(numRes)
	orflst=names(medRes)
	genelst=as.character(fitdf$Gene[match(orflst,as.character(fitdf$ORF))])

	results=data.frame(Gene=genelst,ORF=orflst,MedianFit=medRes,MeanFit=meanRes,VarianceFit=varRes,NumRepeats=numRes,SEFit=seRes,Treatment=rep(treatment,length(orflst)))

	write.table(results,file=outputfile,sep="\t",quote=FALSE,row.names=FALSE)
	return(results)
}

correlationReport<-function(scrnms,dataframe,outputfile,aw=4,ah=4,fitmax=185){
        # Tests the correlation of all possible repeats of the same plate/treatment
        # Useful tool for searching for incorrect plate orientation, or misplaced/mislabelled plates
        # Can also give clues about plates with incorrect medium
        
        # This might fail for various reasons
        # But probably not a critical part of any workflow, so wrap in "try"
        try({
        # Generate correlation plot for every possible 2-way comparison
        poss=combn(scrnms,2)
        trts=unique(as.character(dataframe$Treatment))
        plates=unique(as.numeric(as.character(dataframe$MasterPlate.Number)))
	  plates=plates[order(plates)]
	  pnum=length(plates)
        # Adjust aw, ah to get all plates on one page (e.g. 4*4 gives 16 panels, good for 15 plate lib)
        pdummy=aw*ah-pnum
	  pdf(outputfile)
        corrs=c()
        op<-par(mfrow=c(aw,ah),cex.main=0.8,mar=c(1,1,1,1),mgp=c(0,0,0))

        for(t in trts){
        for(comb in 1:length(poss[1,])){
        for(p in plates){
                rep1=poss[1,comb]; rep2=poss[2,comb]
                r1=dataframe[(as.character(dataframe$Treatment)==t)&(as.character(dataframe$Screen.Name)==rep1)&(as.numeric(as.character(dataframe$MasterPlate.Number))==p),]
                r2=dataframe[(as.character(dataframe$Treatment)==t)&(as.character(dataframe$Screen.Name)==rep2)&(as.numeric(as.character(dataframe$MasterPlate.Number))==p),]
                r1=r1[order(r1$Screen.Name,r1$MasterPlate.Number,r1$Gene),]
                r2=r2[order(r2$Screen.Name,r2$MasterPlate.Number,r2$Gene),]
                
		    if((length(r1$fit)==length(r2$fit))&(length(r1$fit)>0)){
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
		    }else{plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab=rep1,ylab=rep2,main=paste("Trt:",t,"Plate:",p),axes=FALSE)}
        }
        if(pdummy>0) for(p in 1:pdummy) plot(NULL,xlim=c(0,fitmax),ylim=c(0,fitmax),xlab="",ylab="",main="",axes=FALSE)
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

plateBoxplots<-function(dataframe,outputfile,fitmax=185){
# Generates plate by plate boxplots of culture fitnesses
# Useful for identifying plates with medium problems
trts=unique(dataframe$Treatment)
scrns=unique(dataframe$Screen.Name)
trts=trts[order(trts)]
scrns=scrns[order(scrns)]
pdf(outputfile)
for (trt in trts){
	for (scr in scrns){
		dt=dataframe[(dataframe$Screen.Name==scr)&(dataframe$Treatment==trt),]
		dt=dt[order(as.numeric(as.character(dt$MasterPlate.Number))),]
		plateNums=unique(as.numeric(as.character(dt$MasterPlate.Number)))
		plateNums=as.character(plateNums)
		dt$MasterPlate.Number<-factor(dt$MasterPlate.Number,levels=plateNums)
		boxplot(dt$fit~dt$MasterPlate.Number,notch=TRUE,xlab="Plate Number",ylab="Fitness",col=rainbow(23),main=paste(trt,scr),ylim=c(0,fitmax),cex.axis=0.5)
}
}
dev.off()
}

