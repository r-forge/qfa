#########################################################################
# Demonstration of the computational component of Quantitative Fitness Analysis.
# Repeating a subset of the genetic epistasis analysis from Addinall et al. PLoS Genetics
# http://dx.doi.org/10.1371/journal.pgen.1001362
# Looking at gene deletions interacting with cdc13-1 at 27C
# Using 8 replicates of a plate containing telomere specific deletions (plate 15)
# combined with a neutral, wild-type deletion (control: ura3D) or a temperature sensitive, 
# telomere capping mutation (query: cdc13-1).
#########################################################################

library(qfa)

# Change current working directory to root of the qfa package so that we can read in the demo files
# This block is NOT required for a normal workflow
current=getwd()
packroot=file.path(system.file(package = "qfa"),"extdata")
setwd(packroot)

paste("To read the source code for this demo, open this file:",file.path(system.file(package = "qfa"),"demo","telomereCap.R"))

# Assign plate ids, ORFs and gene names to observations
control=colonyzer.read(files=c("URA3Control.dat"),experiment="URA3ExptDescription.txt",ORF2gene="ORF2GENE.txt",libraries="LibraryDescription.txt",screenID="URA3D")
query=colonyzer.read(files=c("cdc13-1Query.dat"),experiment="cdc13-1ExptDescription.txt",ORF2gene="ORF2GENE.txt",libraries="LibraryDescription.txt",screenID="cdc13-1")

# Only use one of the 8 available replicates to save time.  Comment these lines out to do full analysis (takes about 30 mins)
control=control[(control$Library.Name=="BooneSDLV3")&(control$RepQuad==1),]
query=query[(query$Library.Name=="BooneSDLV3")&(query$RepQuad==1),]

# Strip non-experimental edge cultures
control=control[(control$Row!=1)&(control$Col!=1)&(control$Row!=16)&(control$Col!=24),]
query=query[(query$Row!=1)&(query$Col!=1)&(query$Row!=16)&(query$Col!=24),]

# Define which measure of cell density to use
control$Growth=control$Intensity
query$Growth=query$Intensity

# Fit generalised logistic model to observed data
control.fit<-qfa.fit(control,inocguess=1.4e-05,ORF2gene="ORF2GENE.txt",fixG=TRUE,detectThresh=0.001,AUCLim=4,STP=4,glog=FALSE)
query.fit<-qfa.fit(query,inocguess=1.4e-05,ORF2gene="ORF2GENE.txt",fixG=TRUE,detectThresh=0.001,AUCLim=4,STP=4,glog=FALSE)

# Construct QFA fitness measures
control.fit=makeFitness(control.fit,AUCLim=1:4)
query.fit=makeFitness(query.fit,AUCLim=1:4)

# Revert to original working directory to prepare for outputting some files
# Also not necessary for a regular workflow
setwd(current)
paste("Output files will be generated in this directory:",current)

# Produce pdfs of fitted curves & data
qfa.plot("cdc13-1_GrowthCurves.pdf",query.fit,query,maxg=0.25,maxt=6)
qfa.plot("ura3D_GrowthCurves.pdf",control.fit,control,maxg=0.25,maxt=6)

# Write summarised fitnesses to file
qresults=fitnessReport("27","cdc13-1_Fitnesses.txt",query.fit)
cresults=fitnessReport("27","ura3D_Fitnesses.txt",control.fit)

# Epistasis analysis
fitness_definition="r"
epi<-qfa.epi(query.fit,control.fit,0.05,plot=FALSE,wctest=TRUE,reg="splitreg",fdef=fitness_definition)
pdf("FitnessPlotGIS.pdf")
qfa.epiplot(epi,0.05,xxlab="ura3D Fitness",yylab="cdc13-1 Fitness",mmain=paste("Median fitness,",fitness_definition,"(Wilcoxon test for significance) at 27C"),fmax=max(control.fit[[fitness_definition]]))
dev.off()
report.epi(epi$Results,file="EpistasisReport.txt")
