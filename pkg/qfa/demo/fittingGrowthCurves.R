# Simulations using parameter values representing a typical conc. spotted growth curve
times=seq(0,5,0.25)
modelfn=function(t) Glogist(K=0.25,r=7.5,g=9e-4,v=0.3,t)
# Add some measurement error to simulate experimental observations
simobsLogNormal=exp(log(modelfn(times))+rnorm(length(times),mean=0,sd=0.02)) # Lognormal measurement error
simobsNormal=pmax(0.00001,modelfn(times)+rnorm(length(times),mean=0,sd=0.015)) # Normal measurement error
simobs=simobsLogNormal

# Use QFA package to fit model to data
# Constrain parameter search space (put physically realistic lower and upper bounds on parameters)
xybounds=list(K=c(0,1),r=c(0,15),g=c(5E-7,0.1),v=c(0.01,10))
fastfit=data.fit(times,simobs,0.01,xybounds)
fastfitlog=data.fit(times,simobs,0.01,xybounds,logTransform=TRUE)

# Draw some plots comparing fits to data
fitfn=function(t,obj) Glogist(K=obj[["K"]],r=obj[["r"]],g=obj[["g"]],v=obj[["v"]],t)
mkplot=function(logscale=NULL){
	ylim=c(0,1.1*max(simobs))
	if(is.null(logscale)){mlab="Cell density on linear scale"}else{mlab="Cell density on log scale"}
	curve(modelfn,from=min(times),to=max(times),xlab="Time (d)",ylab="Cell density (AU)",col="pink",lwd=2,log=logscale,main=mlab)
	points(times,simobs,type="b")
	curve(fitfn(x,fastfit),from=min(times),to=max(times),add=TRUE,col="blue")
	curve(fitfn(x,fastfitlog),from=min(times),to=max(times),add=TRUE,col="green")
}
op=par(mfrow=c(1,2))
mkplot()
mkplot(logscale="y")
legend("bottomright",legend=c("True process", "Simulated obs.","Fit linear scale","Fit log scale"),col=c("pink","black","blue","green"),lwd=2)
par(op)

# It is easy to find examples of observations simulated with normal measurement error for which it is not possible to find a fit 
# that looks good by eye.  This happens when simulated observations deviate from the true process by a large amount when cell density is low
# (at the start of the experiment).  According to the lognormal measurement error model, large deviations when cell density is low are very 
# unlikely, and so the mean process (the growth curve) is forced through these points to compensate.