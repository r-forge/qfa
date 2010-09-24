require(qfa)
cdc.d<-qfa.data(cdc)
cdc.m<-qfa.model(cdc.d,5*10^5,200)
cdc.pos<-qfa.simulate(cdc.m,initial.update=5*10^3,samples=10^4,thins=5)

pdf("CDC27ACF.pdf")
cdc.acf<-qfa.acf(cdc.pos)
dev.off()

