### Creates an object with information in the format for running the SHM ###
### User should strip unwanted strains (e.g. edge spots) and ensure that ###
### data only include plates from one experiment (e.g. unique tret & med)###
### User should also ensure that density observations are sorted by Expt.Time ###
### and that any observations with "negative times" are stripped or time set to zero ###
### The SHM_preprocess function does not sort or filter input data. ###
### This allows us to attach results to data, preserving metadata.  ###
### Data should be sorted by ORF then by spot ID.
SHM_preprocess<-function(a)
{
  a=a[order(a$ORF,a$ID),]
  ORFuni=unique(a$ORF) 
  # Was a bug here, need to sort (added line above) to give same ordering as that returned from lapply...
  # Result was that wrong number of reps were being analysed for many genotypes (Reps incorrect below).

  gene<-a$Gene[match(ORFuni,a$ORF)]

  N<-length(ORFuni)
  M<-length(unique(a$ID))

  Reps<-as.numeric(lapply(split(a,a$ORF),function(x) length(unique(x$ID))))
  Times<-as.numeric(lapply(split(a,a$ID),nrow))
  SumReps<-cumsum(c(0,Reps))

  dimr<-max(Reps)
  dimc<-max(Times)

  gmat=array(data=NA,dim=c(max(Reps),max(Times),length(ORFuni)))
  tmat=array(data=NA,dim=c(max(Reps),max(Times),length(ORFuni)))
  for (orfno in seq_along(ORFuni)){
	df=a[a$ORF==ORFuni[orfno],]
	ids=unique(df$ID)
	for(idno in seq_along(ids)){
		g=df$Growth[df$ID==ids[idno]]
		t=df$Expt.Time[df$ID==ids[idno]]
		gmat[idno,1:length(g),orfno]=g
		tmat[idno,1:length(t),orfno]=t
	}
  }

  QFA.I<-list("NoORF"=Reps,"NoTime"=Times,"NoSum"=SumReps,
    "N"=N,"M"=M,"gene"=gene)
  QFA.D<-list(y=gmat,x=tmat,ORFuni=ORFuni)
  tmat[is.na(tmat)]=-999
  gmat[is.na(gmat)]=-999
  xx<-aperm(tmat,c(2,1,3))
  yy<-aperm(gmat,c(2,1,3))

  QFA.x=as.double(xx)
  QFA.y=as.double(yy)
  QFA.NoORF=Reps
  QFA.NoTIME=Times
  QFA.NoSUM=SumReps
  QFA.I=as.integer(c(N,max(Reps),max(Times),dimr*dimc*N,length(Times)))
   
  #QFA.I=as.integer(c(L,M,N,maxy,maxNoTIME)
  SHM=list(y=QFA.D$y,x=QFA.D$x,QFA.I=QFA.I,QFA.y=QFA.y,QFA.x=QFA.x,
    QFA.NoORF=QFA.NoORF,QFA.NoTIME=QFA.NoTIME,QFA.NoSUM=QFA.NoSUM,gene=gene,orf=ORFuni)
}

### Takes qfa Raw data (i.e. Colonyzer output with metadata) and qfaBayes SHM posterior
### Returns fitness object similar to the output from qfa.fit
### Optionally draws growth curves and summary curves using mean of posterior parameter values
SHM_makeQFAfits=function(dat,post,priors=NULL,fname=NULL,fdefs=c("r","K","MDR","MDP","MDRMDP","nAUC","nSTP","maxslp","nr")){
	dat=dat[order(dat$ORF,dat$ID),]
	# Values that vary for a given spot correspond to final timepoint
	meta=subset(aggregate(dat,by=list(dat$ID),FUN=tail,1),select=c(-Group.1))
	meta=meta[order(meta$ORF,meta$ID),]
	metanames=unique(meta[,c("ORF","Gene")])
	lup=as.vector(metanames$Gene)
	names(lup)=metanames$ORF
	scrname=unique(meta$Screen.Name)
	gvec=strsplit(scrname,"_")[[1]]
	if(length(gvec)==2){
		scrname=paste(scrname,"CDC13+")
	}else{
		scrname=paste(scrname,"cdc13-1")
	}
	bckgrnd=paste(scrname,unique(meta$Treatment))
	
	num=by(dat,dat$ID,FUN=numericalfitness,0.85,0.85)
	num=num[meta$ID]
	numdf=data.frame(matrix(unlist(num),nrow=length(num),byrow=TRUE))
	names(numdf)=names(num[[1]])

	postmean=colMeans(post)
	
	meta$K=exp(postmean[sprintf("K_lm[%i]",0:(dim(meta)[1]-1))])
	meta$r=exp(postmean[sprintf("r_lm[%i]",0:(dim(meta)[1]-1))])
	meta$g=exp(postmean["P"])
	meta$v=1
	meta$SHM_index=0:(length(meta$K)-1)
	meta=makeFitness(meta)
	meta=cbind(meta,numdf)
	
	plotDists=function(paramSumm="K_o_l[%i]",param="K_lm[%i]",pname="K",prior_mu="K_mu",prior_eta="eta_K_p",dfids=c(),orfno=1,xmax=NULL,Nsamp=100000){
			smmeta=meta[meta$ID%in%dfids,]
			if(is.null(xmax)) xmax=1.25*max(meta[[pname]])
			summVals=exp(post[,sprintf(paramSumm,orfno-1)])
			Dens=density(summVals,from=0,to=xmax)
			plot(Dens,xlim=c(0,xmax),col="red",lwd=3,type="n",main=paste(bckgrnd,pname,unique(smmeta$Gene),sep="\n"))
			if(!is.null(priors)){
				mln=priors[prior_mu]; vln=1/priors[prior_eta]
				PDist=exp(rnorm(Nsamp,mean=mln,sd=sqrt(vln)))
				points(density(PDist,from=0,to=xmax),type="l",col="blue",lwd=2)
			}
			if(!is.null(param)){
				SHMinds=smmeta$SHM_index
				for(ind in SHMinds) points(density(exp(post[,sprintf(param,ind)]),from=0,to=xmax),type="l",col="grey")
			}
			points(Dens,type="l",lwd=4,col="red")
	}
	
	if(!is.null(fname)){
		pdf(fname, height=7, width=9.898)
		tmax=max(dat$Expt.Time)
		gmax=max(dat$Growth)
		uniORF=unique(meta$ORF)
		for(f in fdefs){
			# Biological replicates
			bymedian = with(meta, reorder(Gene,-eval(parse(text=f)),median))
			boxplot(eval(parse(text=f))~bymedian,data=meta,ylab=f,las=2,cex.axis=0.45,cex=0.25,main=paste(bckgrnd,"replicate strains"),col=c("pink","lightblue"))
		}
		# MCMC samples
		plotMCMC=function(paramSumm="r_o_l[%i]",fdef="r"){
			cols=sprintf(paramSumm,0:(length(uniORF)-1))
			samps=post[,cols]
			newcols=sprintf("val.%i",0:(length(uniORF)-1))
			colnames(samps)=newcols
			samps$sampno=1:dim(post)[1]
			rpost=reshape(samps,varying=newcols,direction="long",idvar="sampno",timevar="ORFno")
			rpost$ORF=uniORF[rpost$ORFno+1]
			rpost$Gene=lup[rpost$ORF]
			bymed = with(rpost, reorder(Gene,-val,median))
			boxplot(val~bymed,data=rpost,ylab=fdef,las=2,cex.axis=0.45,cex=0.25,main=paste(bckgrnd,"MCMC summary samples"),col=c("pink","lightblue"))
		}
		plotMCMC(paramSumm="r_o_l[%i]",fdef="r")
		plotMCMC(paramSumm="K_o_l[%i]",fdef="K")

		for(orfno in seq_along(uniORF)){
			df=dat[dat$ORF==uniORF[orfno],]
			uniID=unique(df$ID)
			plot(NULL,xlim=c(0,tmax),ylim=c(0,gmax),main=paste(unique(df$Gene),unique(df$ORF)),xlab="Time (d)",ylab="Cell density (AU)")
			for(ID in uniID){
				points(df$Expt.Time[df$ID==ID],df$Growth[df$ID==ID],type="b",cex=0.5)
				# Curve
				K=meta$K[meta$ID==ID]
				r=meta$r[meta$ID==ID]
				g=meta$g[meta$ID==ID]
				curve((K*g*exp(r*x))/(K+g*(exp(r*x)-1)), from=0, to=tmax,col="blue",add=TRUE)
			}
			K_summ=exp(postmean[sprintf("K_o_l[%i]",orfno-1)])
			r_summ=exp(postmean[sprintf("r_o_l[%i]",orfno-1)])
			g_summ=exp(postmean["P"])
			DT=dtl(K_summ,r_summ,g_summ,1,0)*24*60
			curve((K_summ*g_summ*exp(r_summ*x))/(K_summ+g_summ*(exp(r_summ*x)-1)), from=0, to=tmax,col="red",lwd=3,add=TRUE)
			# Add legend
			legt1<-paste(c("K=","r=","g=","DT="),c(signif(K_summ,3),signif(r_summ,3),signif(g_summ,3),signif(DT,3)),sep="")
			legend("topleft",legt1,box.lty=0,cex=1)
			op=par(mfrow=c(2,2))
			plotDists(paramSumm="K_o_l[%i]",param="K_lm[%i]",pname="K",prior_mu="K_mu",prior_eta="eta_K_p",dfids=df$ID,orfno=orfno)
			plotDists(paramSumm="r_o_l[%i]",param="r_lm[%i]",pname="r",prior_mu="r_mu",prior_eta="eta_r_p",dfids=df$ID,orfno=orfno)
			plotDists(paramSumm="P",param=NULL,pname="g",prior_mu="P_mu",prior_eta="eta_P",dfids=dat$ID[dat$ORF==uniORF[orfno]],orfno=orfno,xmax=0.001)
			par(op)
		}
		dev.off()
		
	}
	
	return(meta)
}

### Calls the C code for running the SHM MCMC ###
SHM_main <- function(burn,iters,thin,adaptive_phase,
  QFA.I,QFA.y,QFA.x,QFA.NoORF,QFA.NoTIME,PRIORS,TUNING){
  if(adaptive_phase>burn){
    stop()
  }
  L=length(QFA.NoORF)
  LM<-sum(QFA.NoORF[1:L])
  NCOL=LM+L+L+1+1+1+LM+L+L+1+1+L+1+1+1+1+1+1
  tmp <- .C("main_SHM", as.integer(burn),as.integer(iters),as.integer(thin),
    as.integer(adaptive_phase),
	OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
    QFAI=as.integer(QFA.I),QFAy=as.double(QFA.y),QFAx=as.double(QFA.x),
	QFANoORF=as.integer(QFA.NoORF),QFANoTIME=as.integer(QFA.NoTIME),
    PRIORS=as.double(PRIORS),TUNING=as.double(TUNING),PACKAGE="qfaBayes")
  mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
  mat=data.frame(mat)
  names(mat)=tmp$HEADER
  mat
}




### Displays repeat level logistic growth curves for each ORF###
plot_SHM_simple<-function(SHM_output,SHM,outfile=NULL){
  samp<-SHM_output
  y<-SHM$y
  x<-SHM$x
  N=L=SHM$QFA.I[1]
  NoSum<-SHM$QFA.NoSUM
  NoORF<-SHM$QFA.NoORF
  NoTime<-SHM$QFA.NoTIME
  gene=SHM$gene

  M=sum(NoORF[1:L])
 
  K_lm=tau_K_l=K_o_l=sigma_K_o=K_p=P=r_lm=tau_r_l=r_o_l=sigma_r_o=r_p=nu_cl=nu_p=sigma_nu=0
  aa<-samp
  #K_lm[%i] # CORRECT
  t=1
  for (i in 1:M){
    j=i
    K_lm[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #tau_K_l[%i] # CORRECT.  ITERATES THROUGH i BUT TAKES jTH ELEMENT FROM samp, THEN INCREMENTS j
  j=M+1
  for (i in (2*M+3*N+8):(2*M+4*N+7)){
    tau_K_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #"K_o_l[%i]  # CORRECT.  ITERATES THROUGH i BUT TAKES jTH ELEMENT FROM samp, THEN INCREMENTS j
  j=M+N+1
  for (i in (M+1):(M+N)){
    K_o_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_K_o "); # CORRECT
  i=2*M+3*N+5
  j=M+2*N+1
  sigma_K_o=mean(samp[,j])
  
  t=1
  #K_p ");  # CORRECT
  i=M+1+N
  j=M+2*N+2
  K_p=mean(samp[,j])

  t=1
  #"P # CORRECT
  i=(M+N+2)
  j=M+2*N+3
  P=mean(samp[,j])

  t=1
  #r_lm[%i] # CORRECT
  j=M+2*N+4
  for (i in (M+2*N+4):(2*M+2*N+3)){
    r_lm[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #tau_r_l[%i] ",l); # CORRECT.  ITERATES THROUGH i BUT TAKES jTH ELEMENT FROM samp, THEN INCREMENTS j
  j=2*M+2*N+4
  for (i in (2*M+4*N+8):(2*M+5*N+7)){
    tau_r_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #r_o_l[%i] ",l); # CORRECT.  ITERATES THROUGH i BUT TAKES jTH ELEMENT FROM samp, THEN INCREMENTS j
  j=2*M+3*N+4
  for (i in (2*M+2*N+4):(2*M+3*N+3)){
    r_o_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }
  
  t=1
  #sigma_r_o "); # CORRECT 
  i=2*M+3*N+7
  j=2*M+4*N+4
  sigma_r_o=mean(samp[,j]);

  t=1
  #r_p "); # CORRECT
  i=2*M+3*N+4
  j=2*M+4*N+5
  r_p=mean(samp[,j]);

  t=1
  #"nu_cl[%i] ",l);  # CORRECT?  nu_l[0]
  j=2*M+4*N+6
  for (i in (M+N+3):(M+2*N+2)){
    nu_cl[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_nu "); # CORRECT
  i=2*M+3*N+6
  j=2*M+5*N+6
  sigma_nu=mean(samp[,j]);

  t=1
  #nu_p "); # CORRECT
  i=M+2*N+3
  j=2*M+5*N+7
  nu_p=mean(samp[,j]);

  ###
  K<-exp(K_p)
  K_i<-exp(K_o_l)
  K_ij<-exp(K_lm)
  P<-exp(P)
  r<-exp(r_p)
  r_i<-exp(r_o_l)
  r_ij<-exp(r_lm)
  taui<-exp(nu_cl)
  tau<-exp(nu_p)
  K_i_tau<-exp(sigma_K_o)
  r_i_tau<-exp(sigma_r_o)
  K_ij_tau<-exp(tau_K_l)
  r_ij_tau<-exp(tau_r_l)
  
  if(!is.null(outfile)) pdf(outfile)
  ylimmax=max(y,na.rm=TRUE)
  xlimmax=max(x,na.rm=TRUE)
  for (i in 1:N){
    plot(-1,-1,main=paste(gene[i],"Repeat Curves"),xlab="Time (days)",
	  ylab="Culture Density (AU)",xlim=c(0,xlimmax),ylim=c(0,ylimmax))
    for (j in 1:NoORF[i]){
	  srt=order(x[j,,i])
      points(x[j,,i][srt],y[j,,i][srt],type="b",cex=0.5)
      KK=K_ij[(j+NoSum[i])]
      rr=r_ij[(j+NoSum[i])]
      curve((KK*P*exp(rr*x))/(KK+P*(exp(rr*x)-1)), 0, xlimmax,add=TRUE) 
    }
    K=exp(K_o_l[i])
    r=exp(r_o_l[i])
	DT=dtl(K,r,P,1,0)*24*60
    curve((K*P*exp(r*x))/(K+P*(exp(r*x)-1)),lwd=3,col="red",add=T)
	legt1<-paste(c("K=","r=","g=","DT="),c(signif(K,3),signif(r,3),signif(P,3),signif(DT,3)),sep="")
	legend("topleft",legt1,box.lty=0,cex=1)
  }
  
  if(!is.null(outfile)) dev.off()
  
}

### Creates univariate MDRxMDP fitnesses measures from SHM output (x2) for input to the IHM ###
IHM_MDRxMDP_postpro<-function(SHM_a,SHM_output_a,SHM_b,SHM_output_b){
  QFA.yA=colMeans(SHM_output_a)
  L=SHM_a$QFA.I[1]
  SHIFTmn=SHM_a$QFA.NoSUM[L+1]
  t=0;
  K_lm=r_lm=P=numeric()
  
  for (i in 1:(SHIFTmn)){
   	t=t+1
	K_lm[i]=exp(QFA.yA[t])
  }
	
  t=t+2*L+3
  P=exp(QFA.yA[t])
 
  for (i in 1:SHIFTmn){
	t=t+1
	r_lm[i]=exp(QFA.yA[t])
  }
        
  for (i in 1:SHIFTmn){
	if(K_lm[i]<=2*P){
		K_lm[i]=2*P
		r_lm[i]=0
    }
  }
		
  for (i in 1:SHIFTmn){
	QFA.yA[i]=(r_lm[i]/log(2*max(0,K_lm[i]-P)/max(0,K_lm[i]-2*P)))*(log(K_lm[i]/P)/log(2));
  }

  QFA.yB=colMeans(SHM_output_b)
  L=SHM_b$QFA.I[1]
  SHIFTmn=SHM_b$QFA.NoSUM[L+1]
  t=0
  K_lm=r_lm=P=numeric()
  
  for (i in 1:(SHIFTmn)){
	t=t+1
	K_lm[i]=exp(QFA.yB[t])
  }
	
  t=t+2*L+3
  P=exp(QFA.yB[t])
 
  for (i in 1:SHIFTmn){
	t=t+1
	r_lm[i]=exp(QFA.yB[t])
  }
        
  for (i in 1:SHIFTmn){
	if(K_lm[i]<=2*P){
		K_lm[i]=2*P
		r_lm[i]=0
    }
  }
		
  for (i in 1:SHIFTmn){
	QFA.yB[i]=(r_lm[i]/log(2*max(0,K_lm[i]-P)/max(0,K_lm[i]-2*P)))*(log(K_lm[i]/P)/log(2));
  }

  list(QFA.yA=QFA.yA,QFA.yB=QFA.yB)
}

### Calls the C code for running the IHM MCMC ###
IHM_main <- function(burn,iters,thin,adaptive_phase,
  QFA.IA,QFA.yA,QFA.NoORFA,QFA.IB,QFA.yB,QFA.NoORFB,PRIORS,TUNING){
  if(adaptive_phase>burn){
    stop()
  }
  aa<-QFA.NoORFA
  bb<-QFA.NoORFB
  if(!(length(aa)==length(bb))){
    stop()
  }
  L=length(aa)
  NCOL=L+1+1+2*L+1+1+1+L+L+1
  tmp <- .C("main_IHM", as.integer(burn),as.integer(iters),as.integer(thin),
    as.integer(adaptive_phase),
    OUT=as.double(1:(NCOL*iters)),HEADER=as.character(rep("NULLNULL",NCOL)),
    QFAIA=QFA.IA,QFAyA=as.double(QFA.yA),QFANoORFA=as.integer(QFA.NoORFA),
    QFAIB=as.integer(QFA.IB),QFAyB=as.double(QFA.yB),
	QFANoORFB=as.integer(QFA.NoORFB),PRIORS=as.double(PRIORS),
	TUNING=as.double(TUNING),PACKAGE="qfaBayes")
  mat=matrix((tmp$OUT),nrow=iters,byrow=T)
  mat=data.frame(mat)
  names(mat)=tmp$HEADER
  mat
}

### Displays a fitness plot of the control vs query, using MDRxMDP ###
plot_IHM_simple<-function(IHM_output,SHM){
  N=SHM$QFA.I[1]
  gene=SHM$gene
  samp=IHM_output
  if(nrow(samp)>1){
    vecsamp<-colMeans(samp)
  } 
  else{
    vecsamp<-as.numeric(samp)
  }
  namesamp<-names(vecsamp)
  Z_l<-exp(vecsamp[1:(N)])
  sigma_Z<-exp(vecsamp[N+1])
  Z<-exp(vecsamp[N+2])
  nu_cl<-exp(vecsamp[(N+3):(3*N+2)])
  sigma_nu<-exp(vecsamp[3*N+3])
  nu<-exp(vecsamp[(3*N+4)])
  A1<-exp(0)
  A2<-exp(vecsamp[3*N+5])
  delta<-vecsamp[(3*N+6):(4*N+5)]
  gamma<-vecsamp[(4*N+6):(5*N+5)]
  sigma_gamma<-exp(vecsamp[(5*N+6)])
  if(nrow(samp)>1){
    delta_gamma<-colMeans(samp[,(3*N+6):(4*N+5)]*samp[,(4*N+6):(5*N+5)])
  }
  else{
    delta_gamma<-colMeans(samp[,(3*N+6):(4*N+5)]*(samp[,(4*N+6):(5*N+5)]))
  }
  delta_gamma=exp(delta_gamma)

  sig<-sum(rep(1,N)[delta>0.5])
  order<-order(1-delta)
  vecorder<-order(1-delta)[1:sig]

  limmin<-0
  limmax<-max(A2*Z_l*delta_gamma)
  limmaxx<-max(A1*Z_l)
  i=1:N
  plot(1,type="n",main=expression(paste("Treatment",Treatment,degree,"C",
  " (delta=Posterior Expectations)")),ylim=c(limmin,limmax),
  xlim=c(limmin,limmaxx),xlab="Control Fitness (=exp(Z_l))",
  ylab="Query Fitness (=exp(alpha+Z_l+delta_l*gamma_l))",col=8,pch=19,cex=0.5)
  lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
  lines(A1*c(-1000,10000),A2*c(-1000,10000),col="grey",lwd=2)
  points(A1*Z_l[i], A2*(Z_l[i]*delta_gamma[i]),col=8,pch=19,cex=0.5)
  i=vecorder[log(delta_gamma)[vecorder]>0]
  points(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),col=2,pch=19,cex=0.5)
  i=vecorder[log(delta_gamma)[vecorder]<=0]  
  points(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),col=3,pch=19,cex=0.5)
  i=vecorder
  text(A1*Z_l[i],A2*(Z_l[i]*delta_gamma[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}

### Creates an object with information in the format for running the JHM ###
JHM_postpro<-function(a,Treatment_a,Screen_a,MPlate_a,remove_row_a,remove_col_a,b,Treatment_b,Screen_b,MPlate_b,remove_row_b,remove_col_b)
{
  print("Stripping extra data and sorting...")
  # Discard any superfluous data from a
  a<-funcREMOVE(a,Screen_a,Treatment_a,MPlate_a)
  if (length(remove_row_a)>=1) a=a[!a$Row%in%remove_row_a,]
  if (length(remove_col_a)>=1) a=a[!a$Col%in%remove_col_a,]

  # Sort a by plate, row, col, time
  aRow<-sprintf("%02d",a$Row)
  aCol<-sprintf("%02d",a$Col)
  aPlate<-sprintf("%02d",a$MasterPlate.Number)
  a$ID<-paste(a$Barcode,aPlate,aRow,aCol,sep="_")
  a<-a[order(a$ORF,a$ID,a$Expt.Time), ]

  # Discard any superfluous data from b
  b<-funcREMOVE(b,Screen_b,Treatment_b,MPlate_b)
  if (length(remove_row_b)>=1) b=b[!b$Row%in%remove_row_b,]
  if (length(remove_col_b)>=1) b=b[!b$Col%in%remove_col_b,]

  # Sort b by plate, row, col, time
  bRow<-sprintf("%02d",b$Row)
  bCol<-sprintf("%02d",b$Col)
  bPlate<-sprintf("%02d",b$MasterPlate.Number)
  b$ID<-paste(b$Barcode,bPlate,bRow,bCol,sep="_")
  b<-b[order(b$ORF,b$ID,b$Expt.Time), ]
  
  # Only analyse ORFs common to both screens
  ORFuni_a<-unique(a$ORF)
  ORFuni_b<-unique(b$ORF)
  ORFuni<-sort(intersect(ORFuni_a,ORFuni_b))
  N<-length(ORFuni)
  a=a[a$ORF%in%ORFuni,]
  b=b[b$ORF%in%ORFuni,]
  
  # Map from ORF to gene name
  gene<-a$Gene[match(ORFuni,a$ORF)]

  print("Calculating number of repeats, times etc. for each orf...")
  Ma=length(unique(a$ID))
  NoORF_a<-as.numeric(lapply(split(a,a$ORF),funcNoORF))#no of repeats each orf
  NoTime_a<-c(0,as.numeric(lapply(split(a,a$ID),nrow)))# 0+ no of time each repeat
  NoSum_a<-c(0,cumsum(NoORF_a))

  Mb=length(unique(b$ID))
  NoORF_b<-as.numeric(lapply(split(b,b$ORF),funcNoORF))#no of repeats each orf
  NoTime_b<-c(0,as.numeric(lapply(split(b,b$ID),nrow)))# 0+ no of time each repeat
  NoSum_b<-c(0,cumsum(NoORF_b))

  dimr<-max(NoORF_a,NoORF_b)
  dimc<-max(NoTime_a,NoTime_b)

  # This code results in stray data points being added to the start of each timeseries
  # I don't understand how it is supposed to work, so can't fix, so replacing it instead CONOR Aug 2014
  #y<-funcXY_J(a$Growth,b$Growth,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)
  #x<-funcXY_J(a$Expt.Time,b$Expt.Time,Ma,Mb,N,NoTime_a,NoSum_a,NoTime_b,NoSum_b,dimr,dimc)
  
  print("Filling data arrays...")
  # Initialise 4D arrays for observation times (x) and colony sizes (y) 
  x=array(NA,dim=c(dimr,dimc,N,2))
  y=array(NA,dim=c(dimr,dimc,N,2))
  pb <- txtProgressBar(min = 1, max = N, style = 3)
  for(orfno in seq_along(ORFuni)){
	orf=ORFuni[orfno]
	setTxtProgressBar(pb, orfno)
	orfa=a[a$ORF==orf,]; orfb=b[b$ORF==orf,]
	orfa=split(orfa,orfa$ID); orfb=split(orfb,orfb$ID)
	# For each available replicate of each orf in the control and query condition, add data to arrays
	repind=1
	for(repl in orfa){
		tim=repl$Expt.Time
		siz=repl$Growth
		len=length(tim)
		x[repind,1:len,orfno,1]=tim
		y[repind,1:len,orfno,1]=siz
		repind=repind+1
	}
	repind=1
	for(repl in orfb){
		tim=repl$Expt.Time
		siz=repl$Growth
		len=length(tim)
		x[repind,1:len,orfno,2]=tim
		y[repind,1:len,orfno,2]=siz
		repind=repind+1
	}
  }
  close(pb)

  QFA.D<-list(x=x,y=y)

  x[is.na(x)]=-999
  y[is.na(y)]=-999
  xx_a<-aperm(x[,,,1],c(2,1,3))
  yy_a<-aperm(y[,,,1],c(2,1,3))

  xx_b<-aperm(x[,,,2],c(2,1,3))
  yy_b<-aperm(y[,,,2],c(2,1,3))

  list(y=QFA.D$y,x=QFA.D$x,QFA.IA=c(N,max(NoORF_a),max(NoTime_a),length(y)/2,
    length(NoTime_a[-1])),QFA.yA=c(yy_a),QFA.xA=c(xx_a),QFA.NoORFA=c(NoORF_a),
	QFA.NoTIMEA=c(NoTime_a)[-1],QFA.NoSUMA=c(NoSum_a),QFA.IB=c(N,max(NoORF_b),
	max(NoTime_b),length(y)/2,length(NoTime_b[-1])),QFA.yB=c(yy_b),
	QFA.xB=c(xx_b),QFA.NoORFB=c(NoORF_b),QFA.NoTIMEB=c(NoTime_b)[-1],
	QFA.NoSUMB=c(NoSum_b),gene=gene,orf=ORFuni)
}

### Calls the C code for running the JHM MCMC ###
JHM_main<- function(burn,iters,thin,adaptive_phase,QFA.IA,QFA.yA,QFA.xA,QFA.NoORFA,QFA.NoTIMEA,QFA.IB,QFA.yB,QFA.xB,QFA.NoORFB,QFA.NoTIMEB,PRIORS,TUNING) {
  if(adaptive_phase>burn){
    stop()
  }
  aa<-QFA.NoORFA
  bb<-QFA.NoORFB
  if(!(length(aa)==length(bb))){
    stop()
  }
  L=length(aa)
  LMa<-sum(aa)
  LMb<-sum(bb)
  NCOL=LMa+LMb+2*L+L+1+1+1+LMa+LMb+2*L+L+1+1+2*L+1+1+1+1+L+L+1+L+1+2*2+2*2
  tmp <- .C("main_JHM", as.integer(burn),as.integer(iters),as.integer(thin),
    as.integer(adaptive_phase),OUT=as.double(1:(NCOL*iters)),
    HEADER=as.character(rep("NULLNULL",NCOL)),
    QFAIA=as.integer(QFA.IA),QFAy=as.double(QFA.yA),QFAxA=as.double(QFA.xA),
    QFANoORFA=as.integer(QFA.NoORFA),QFANoTIMEA=as.integer(QFA.NoTIMEA),
    QFAIB=as.integer(QFA.IB),QFAy=as.double(QFA.yB),QFAxB=as.double(QFA.xB),
    QFANoORFB=as.integer(QFA.NoORFB),QFANoTIMEB=as.integer(QFA.NoTIMEB),
    PRIORS=as.double(PRIORS),TUNING=as.double(TUNING),PACKAGE="qfaBayes")
  mat=matrix(c(tmp$OUT),nrow=iters,byrow=T)
  mat=data.frame(mat)
  names(mat)=tmp$HEADER
  mat
}

### Displays a fitness plot of the control vs query, using MDRxMDP, K and then r ###
plot_JHM_simple<-function(JHM_output,JHM){
  samp<-JHM_output
  gene=JHM$gene
  L<-JHM$QFA.IA[1]
  M=sum(c(JHM$QFA.NoORFA,c(JHM$QFA.NoORFB)))

  K_clm=tau_K_cl=K_o_l=sigma_K_o=K_p=P=r_clm=tau_r_cl=r_o_l=sigma_r_o=r_p=
    nu_cl=sigma_nu=nu_p=alpha_c=beta_c=delta_l=gamma_cl=sigma_gamma=omega_cl=sigma_omega=upsilon_c=sigma_upsilon=0
  ####
  t=1
  #K_clm
  for (i in 1:c(M)){
    j=i
    K_clm[t]=mean(samp[,j]);t=t+1
  }

  t=1
  #tau_K_cl
  j=M+1
  for (i in (2*M+9*L+15):(2*M+11*L+14)){
    tau_K_cl[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #K_o_l
  j=M+2*L+1
  for (i in (M+1):(M+L)){
    K_o_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_K_o
  i=2*M+9*L+9
  j=M+3*L+1
  sigma_K_o=mean(samp[,j])

  t=1
  #K_p
  i=M+L+1
  j=M+3*L+2
  K_p=mean(samp[,j])

  t=1
  #P
  i=M+L+2
  j=M+3*L+3
  P=mean(samp[,j])

  t=1
  #r_clm
  j=M+3*L+4
  for (i in (M+8*L+8):(2*M+8*L+7)){
    r_clm[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #tau_r_cl
  j=2*M+3*L+4
  for (i in (2*M+11*L+15):(2*M+13*L+14)){
    tau_r_cl[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #r_o_l
  j=2*M+5*L+4
  for (i in (2*M+8*L+8):(2*M+9*L+7)){
    r_o_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_r_o
  i=2*M+9*L+13
  j=2*M+6*L+4
  sigma_r_o=mean(samp[,j])

  t=1
  #r_p
  i=2*M+9*L+8
  j=2*M+6*L+5
  r_p=mean(samp[,j])


  t=1
  #nu_cl
  j=2*M+6*L+6
  for (i in (M+5*L+7):(M+7*L+6)){
    nu_cl[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_nu
  i=2*M+9*L+11
  j=2*M+8*L+6
  sigma_nu=mean(samp[,j])

  t=1
  #nu_p
  i=M+6*L+7
  j=2*M+8*L+7
  nu_p=mean(samp[,j])

  t=1
  #alpha_c
  i=M+L+4
  j=2*M+8*L+8
  alpha_c=mean(samp[,j])

  t=1
  #beta_c
  i=M+L+6
  j=2*M+8*L+9
  beta_c=mean(samp[,j])

  t=1
  #delta_l
  j=2*M+8*L+10
  for (i in (M+2*L+7):(M+3*L+6)){
    delta_l[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #gamma_cl
  j=2*M+9*L+10
  for (i in (M+4*L+7):(M+5*L+6)){
    gamma_cl[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_gamma
  i=2*M+9*L+10
  j=2*M+10*L+10
  sigma_gamma=mean(samp[,j])

  t=1
  #omega_cl
  j=2*M+10*L+11
  for (i in (M+7*L+8):(M+8*L+7)){
    omega_cl[t]=mean(samp[,j]);t=t+1
    j=j+1
  }

  t=1
  #sigma_omega
  i=2*M+9*L+12
  j=2*M+11*L+11
  sigma_omega=mean(samp[,j])

  K<-exp(K_p)
  K_i<-exp(K_o_l)
  K_ij<-exp(K_clm)
  PO<-exp(P)
  r<-exp(r_p)
  r_i<-exp(r_o_l)
  r_ij<-exp(r_clm)

  taui<-exp(nu_cl)
  tau<-exp(nu_p)
  gam<-gamma_cl
  omega<-omega_cl
  nuc<-exp(upsilon_c)

  gamdelt=0
  j=2*M+8*L+10
  jj=2*M+9*L+10
  ii=M+4*L+7
  t=1
  for (i in (M+2*L+7):(M+3*L+6)){
    gamdelt[t]=mean(samp[,j]*samp[,jj]);t=t+1
    j=j+1
    ii=i+1
    jj=jj+1
  }

  omegadelt=0
  j=2*M+8*L+10
  jj=2*M+10*L+11
  ii=M+7*L+8
  t=1
  for (i in (M+2*L+7):(M+3*L+6)){
    omegadelt[t]=mean(samp[,j]*samp[,jj]);t=t+1
    j=j+1
    ii=i+1
    jj=jj+1
  }

  delta<-delta_l

  A1<-1
  A2<-exp(alpha_c)
  B1<-1
  B2<-exp(beta_c)
  sig<-sum(rep(1,L)[delta>0.5])
  order<-order(1-delta)
  vecorder<-order(1-delta)[1:sig]

  K_ij=vecK=c(exp(K_o_l),A2*exp(K_o_l+gamdelt))
  r_ij=c(exp(r_o_l),B2*exp(r_o_l+omegadelt))

  K_ij[K_ij<2*PO]=2*PO+0.001

  vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
  vecMDPa<-log(K_ij/PO)/log(2) #MDP
  vecMDRPMDR<-vecMDRa*vecMDPa
  vecMDRPMDR[vecK<2*PO]=0
  mu_a=(vecMDRPMDR)[1:L]
  mu_b=(vecMDRPMDR)[(1+L):(2*L)]

  limmin<-0
  limmax<-max(na.omit(c(mu_b)))
  limmaxx<-max(na.omit(c(mu_a)))

  plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",
    ylab="",pch=19,col=8,cex=0.5)
  lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)

  vecMDRa<-r_ij/log(2*(K_ij-PO)/(K_ij-2*PO)) #MDR
  vecMDPa<-log(K_ij/PO)/log(2) #MDP
  vecMDPa*vecMDRa

  i=1:L
  points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,
    pch=19,cex=0.5)
  i=vecorder[omegadelt[order][1:sig]>0]#######
  points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),
    xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
  i=vecorder[omegadelt[order][1:sig]<=0]  #######
  points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),
    xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
  i=vecorder
  text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)

  mu_a=exp(K_o_l)#####
  mu_b=A2*exp(K_o_l+gamdelt)#####
  limmin<-0
  limmax<-max(na.omit(c(mu_b)))
  limmaxx<-max(na.omit(c(mu_a)))

  plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",
  ylab="",pch=19,col=8,cex=0.5)
  lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
  lines(c(-1000,1000),A2*c(-1000,1000),col="grey",lty=2)########
  i=1:L
  points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,
    pch=19,cex=0.5)
  i=vecorder[gamdelt[order][1:sig]>0]####
  points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),
    xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
  i=vecorder[gamdelt[order][1:sig]<=0]  #######
  points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),
    xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
  i=vecorder
  text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)

  mu_a=exp(r_o_l)#####
  mu_b=B2*exp(r_o_l+omegadelt)#####
  limmin<-0
  limmax<-max(na.omit(c(mu_b)))
  limmaxx<-max(na.omit(c(mu_a)))

  plot(1,type="n",ylim=c(limmin,limmax),xlim=c(limmin,limmaxx),main="",xlab="",
    ylab="",pch=19,col=8,cex=0.5)
  lines(c(-1000,10000),c(-1000,10000),lwd=2,col="grey",lty=4)
  lines(c(-1000,1000),B2*c(-1000,1000),col="grey",lty=2)#######
  i=1:L
  points(mu_a[i],mu_b[i],ylim=c(limmin,limmax),xlim=c(limmin,limmax),col=8,
    pch=19,cex=0.5)
  i=vecorder[omegadelt[order][1:sig]>0]#######
  points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),
    xlab="Single",ylab="Double",col=2,pch=19,cex=0.5)
  i=vecorder[omegadelt[order][1:sig]<=0]  #######
  points(mu_a[i],(mu_b[i]),ylim=c(limmin,limmax),xlim=c(limmin,limmax),
    xlab="Single",ylab="Double",col=3,pch=19,cex=0.5)
  i=vecorder
  text(mu_a[i],(mu_b[i]),gene[i],pos=4,offset=0.1,cex=0.4)
}

### Gives a trace plot for one of every parameter type to check convergence qfaBayes output ###
# Generalised so same function applies to all models - CONOR.
visual_convergence_simple_check<-function(qfaBayes_output,verbose=TRUE){
	zeroes=grep("\\[0\\]",colnames(qfaBayes_output))
	unindexed=grep("\\[",colnames(qfaBayes_output),invert=TRUE)
	namelist=colnames(qfaBayes_output)[c(zeroes,unindexed)]
	print(namelist)
	cat("Plotting ",length(namelist)," out of ",length(colnames(qfaBayes_output))," possible variables")
	for (pname in namelist){
		op=par(mfrow=c(2,1))
		particles=qfaBayes_output[,pname]
		heidel_welch=heidel.diag(particles)
		eff=effectiveSize(particles)
		hpval=formatC(heidel_welch[3],3)
		effval=formatC(eff,3)
		if(verbose){
			print(pname)
			print(heidel_welch)
		}
		plot(particles,type="l",ylab=pname,xlab="iter",main=paste("Heidel-Welch p-value:",hpval,"ESS:",effval))
		if(var(particles)>0){
			acf(particles,main=pname)
		}else{
			plot(1,type="n",main=pname)
		}
		par(op)
	}
}

### Gives a trace plot for one of every parameter type to check convergence of SHM output ###
visual_convergence_simple_check_SHM<-function(SHM_output){
  visual_convergence_simple_check(SHM_output)
}

### Gives a trace plot for one of every parameter type to check convergence IHM output ###
visual_convergence_simple_check_IHM<-function(IHM_output){
  visual_convergence_simple_check(IHM_output)
}

### Gives a trace plot for one of every parameter type to check convergence IHM output ###
visual_convergence_simple_check_JHM<-function(JHM_output){
  visual_convergence_simple_check(JHM_output)
}
#############
