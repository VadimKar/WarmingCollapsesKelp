


library(deSolve); library(rootSolve);
EQfun1d=function(parmsi,fun,xlims=0:1){
     EQs=rev(uniroot.all(fun,xlims,parms=parmsi))
     EQs=EQs[abs(fun(EQs,parmsi))<1] #singularities get classified as roots
     if(length(EQs)==4) EQs=EQs[1:3]; if(length(EQs)==3) return(EQs);
     if(EQs[1]==0) EQs=c(NA,NA,EQs[1]) else EQs=c(EQs[1],NA,NA); EQs;
}
Bfn1d=function(fun,xlims,parmsi,parv,parlim,xrev=FALSE,N=2e3,transfun=function(x) x,...){
     varseq=seq(parlim[1],parlim[2],len=N)
     xm=t(sapply(varseq, function(x){ parmsi[parv]=x; EQfun1d(parmsi,fun,xlims=xlims); }))
     pltfun(varseq,transfun(xm),xrev,...)
}
EQfun2d=function(parmsi,fun=CoralODEpp,funDE=CoralODE,xlims=0:1,ylims=0:1,resol=1e3,Buff=0.01){
     #Find equilibria of a competition model. Buff is buffer region to exclude trivial 0,0 eq
     dx=0.01*diff(xlims); dy=0.01*diff(ylims); 
     x=matrix(seq(xlims[1]-dx,xlims[2]+dx, length=resol), byrow=F, resol,resol);
     y=matrix(seq(ylims[1]-dy,ylims[2]+dy, length=resol),byrow=T, resol, resol);
     npts = resol*resol; z=fun(x,y,parmsi);
     xz=matrix(z[1:npts],resol,resol); yz=matrix(z[(npts+1):(2*npts)],resol,resol);
     dyxz=xz^2+yz^2; eqf=function(xbc,ybc){ targ=min(dyxz[xbc & ybc]); mn=which(dyxz==targ,arr.ind=TRUE)[1,]; c(x[mn[1],1],y[1,mn[2]]); }
     EQs=ltrim(rbind(xon=eqf(x>Buff,y<Buff), int=eqf(x>Buff,y>Buff), yon=eqf(x<Buff,y>Buff)),0.001)
     xon.st=sum(apply(ode(funDE,y=EQs["xon",],times=c(0,1000),parms=parmsi)[,-1],2,diff)^2)<1e-2
     yon.st=sum(apply(ode(funDE,y=EQs["yon",],times=c(0,1000),parms=parmsi)[,-1],2,diff)^2)<1e-2
     if(!xon.st) EQs[c("int","xon"),]=NA; if(!yon.st) EQs[c("int","yon"),]=NA; EQs;
}
pltfun=function(x,Y,xrev=FALSE,nmx=22,Buff=0.025,axCol="darkgray",...){
     x=x*c(1,-1)[1+(xrev)]; N=nrow(Y);
     # plot(rep(x,3),Y,col=rep(c(4,2,1),each=N),cex=0.5,pch=16,axes=TRUE,...)
     plot(rep(x,3),Y,col=rep(c(4,2,1),each=N),cex=0.5,pch=16,axes=FALSE,col.lab=axCol,...)
     axis(2,las=1,col=8,col.axis=axCol); box(bty="l",col=axCol); if(xrev) axis(1,axTicks(1),abs(axTicks(1)),col=axCol,col.axis=axCol) else axis(1,col=axCol,col.axis=axCol);
     whichs=rep(FALSE,N); whichs[round(seq(1,N,len=20))]=TRUE; 
     if(max(x)>0) td=head(which(complete.cases(Y) & whichs),-1) else td=tail(which(complete.cases(Y) & whichs),-1);
     arrows(x[td],Y[td,1]-Buff,x[td],Y[td,2]+Buff,col=rgb(0,0,1,alpha=0.5),length=0.05)
     arrows(x[td],Y[td,1]-2*Buff,x[td],Y[td,1]-Buff,col=rgb(0,0,1,alpha=0.25),length=0.05)
}




# Mumby et al. (2007) model. Parameter values in Blackwood et al 2012:
# https://link.springer.com/article/10.1007/s12080-010-0102-0
CoralODE=funDE=function(t,state,parms) {
     C=state[1]; M=state[2];
     dC=parms["r"]*C*(1-C-M) - parms["d"]*C - parms["a"]*M*C
     dM=parms["a"]*M*C + parms["y"]*M*(1-C-M) - parms["g"]*M/(1-C)
     list(c(dC,dM))
}
CoralODEpp=fun=function(C,M,parms) {
     dC=parms["r"]*C*(1-C-M) - parms["d"]*C - parms["a"]*M*C
     dM=parms["a"]*M*C + parms["y"]*M*(1-C-M) - parms["g"]*M/(1-C)
     c(dC,dM)
}
Cparms=c(r=1,a=0.1,y=0.8,d=0.44,g=0.3)
dv=seq(0.1,0.7,len=80); Cparms["r"]=1; xd=t(sapply(dv, function(x){ Cparms["d"]=x; EQfun2d(Cparms)[,1]; }))
pltfun(dv,xd,ylab="Coral cover",xlab="Coral mortality")
rv=seq(0.5,1.75,len=80); Cparms["d"]=0.44; xr=t(sapply(rv, function(x){ Cparms["r"]=x; EQfun2d(Cparms)[,1]; }))
pltfun(rv,xr,TRUE,ylab="Coral cover",xlab="Coral growth")





#Kelp model of Karatatev et al 2021 and present work
KelpODE=function(N,parms){ 
     graze=parms["U"]*parms["delta"]/(1 + parms["fg"]*(N^2))
     N*(parms["r"]*(1-N)/(1+graze*parms["RecEat"]) - graze - parms["m"])
}
Kparms=c(r=57.36,Eff=14.96,delta=8.94,RecEat=1.25,fg=12.18,U=0.6,m=0) #mean values from fitted model
Bfn1d(KelpODE,0:1,Kparms,"m",c(0,12),ylab="Kelp cover",xlab="Kelp mortality")
Bfn1d(KelpODE,0:1,Kparms,"r",c(1,60),TRUE,ylab="Kelp cover",xlab="Kelp growth")



#Noy-Meir 1975 / May 1977 model of overgrazing / overharvest of a logistic population. 
grazeDE=function(x,parms)   x*(parms["r"]*(1-x/parms["K"]) - parms["H"]*x / (1+x^2) - parms["m"])
Gparms=c(r=1.25,K=10,H=2,m=0)
Bfn1d(grazeDE,c(0.01,10),Gparms,"r",c(0.5,1.25),TRUE,ylab="Grass biomass",xlab="Grass growth")
# Bfn1d(grazeDE,c(0.1,10),Gparms,"m",c(0,0.5),ylab="Grass biomass",xlab="Grass mortality")
# Bfn1d(grazeDE,c(0.1,15),Gparms,"K",c(4,7),ylab="Grass biomass",xlab="Grass growth")



#A model of fishing intensity dynamics where allee effects arise from slow institutions by Tekwa et al 2018
#This comes from assuming a "fast" fish stock S: dS/dt=S(r-S/K-F)=0 --> S(F)=K*(r-F)
fishnDE=function(Ft,parms)   (parms["r"] - 2*Ft)*(parms["H"]/(Ft*(parms["r"]-Ft)) - parms["K"])
Fparms=c(K=6,H=1,r=1)
# Bfn1d(fishnDE,c(0,8),Fparms,"r",c(0.5,1.25),TRUE,transfun=function(x) ltrim(Fparms["K"]*(Fparms["r"]-x),0),ylab="Stock biomass",xlab="Stock growth")
Bfn1d(fishnDE,c(0,8),Fparms,"K",c(2,8),TRUE,transfun=function(x) ltrim(Fparms["K"]*(Fparms["r"]-t(apply(x,1,rev))),0),ylab="Stock biomass",xlab="Stock growth")








#Pollinator community collapse model from Bhandary et al 2023. Altered k 0.1 to 0.3, Ak 10^4 to 14^4 and a0 0.35 to 0.15
# https://royalsocietypublishing.org/rsos/article/10/2/221363/91981/Rising-temperature-drives-tipping-points-in
polinDE=function(A,parms){
     T=parms["Ta"]+293; T0=293; dTs=-0.5*(T-293)^2; 
     at=parms["a0"]*exp(dTs/parms["siga"]^2); ht=parms["h0"]*exp(dTs/parms["sigh"]^2);
     kt=parms["k0"]*exp(parms["Ak"]*((1/293)-1/T))
     P=(at + A*parms["Y"]/(1+ht*parms["Y"]*A))/parms["beta"]
     A*(at - kt - A*parms["beta"] + parms["Y"]*P/(1+ht*parms["Y"]*P))
}
Pparms=c(beta=1,Y=2, a0=0.15,Ta=3,siga=5,h0=0.15,sigh=15,k0=0.3,Ak=14^4)
Bfn1d(polinDE,c(0,8),Pparms,"Ta",c(0,4),ylab="Pollinator abundance",xlab="Temperature increase")




#Forest collapse model of van Nes et al 2014 with climate-forest feedback = 0, simplified via dP/dt=0
# https://onlinelibrary.wiley.com/doi/abs/10.1111/gcb.12398
treeDE=function(T,parms){
     P=parms["P"] + parms["b"]*T/parms["K"]
     r=parms["rm"]*P/(parms["hP"]+P); p=parms["p"];
     T*(r*(1-T/parms["K"]) - (parms["ma"]*parms["ha"]/(T+parms["ha"])) - (parms["mf"]*parms["hf"]^p)/((T^p)+parms["hf"]^p))
}
Tparms=c(rm=0.3,hP=0.5,K=0.9,P=4, ma=0.15,ha=0.1, mf=0.11,hf=0.64,p=7, b=0)
Bfn1d(treeDE,0:1,Tparms,"P",c(0,5),TRUE,ylab="Tree cover",xlab="Precipitation")


#Food web collapse model of Karatayev et al 2023 with generalized feedbacks
webDE=function(C,parms){
     N=(parms["r"]-parms["delta"]*C)*parms["K"]/parms["r"]
     C*(parms["b"]*parms["delta"]*N*(1-parms["f"]*N/parms["K"]) - parms["m"])
}
Wparms=c(r=1,K=1.35,delta=1.1,b=1,f=0.87,m=0.3)
Bfn1d(webDE,0:1,Wparms,"m",c(0.1,0.5),ylab="Consumer biomass",xlab="Consumer mortality")
Bfn1d(webDE,0:1,Wparms,"b",c(0.5,1.75),TRUE,ylab="Consumer biomass",xlab="Consumer growth")




#All together:
#mgp controls placement of tick labels relative to axis labels.  First is x, second is y; bigger = closer.
par(mgp=c(1.5,0.4,0),mar=c(3,2.5,1,1),mfrow=c(3,4),tck=-0.035)
pltfun(dv,xd,ylab="Coral cover",xlab="Coral mortality",main="Tropical reefs")
pltfun(rv,xr,TRUE,ylab="Coral cover",xlab="Coral growth",main="Tropical reefs")
Bfn1d(KelpODE,0:1,Kparms,"m",c(0,12),ylab="Kelp cover",xlab="Kelp mortality",main="Temperate reefs")
Bfn1d(KelpODE,0:1,Kparms,"r",c(1,60),TRUE,ylab="Kelp cover",xlab="Kelp growth",main="Temperate reefs")
Bfn1d(grazeDE,c(0.01,10),Gparms,"r",c(0.5,1.25),TRUE,ylab="Grass biomass",xlab="Grass growth",main="Grazed grasslands")
Bfn1d(fishnDE,c(0,8),Fparms,"K",c(2,8),TRUE,transfun=function(x) ltrim(Fparms["K"]*(Fparms["r"]-t(apply(x,1,rev))),0),ylab="Fish biomass",xlab="Population growth",main="Fished populations")
Bfn1d(polinDE,c(0,8),Pparms,"Ta",c(0,4),ylab="Pollinator abundance",xlab="Temperature increase",main="Pollinator communities")
Bfn1d(treeDE,0:1,Tparms,"P",c(0,5),TRUE,ylab="Tree cover",xlab="Precipitation",main="Temperate forests")
Bfn1d(webDE,0:1,Wparms,"m",c(0.1,0.5),ylab="Consumer biomass",xlab="Consumer mortality",main="Aquatic food webs")
Bfn1d(webDE,0:1,Wparms,"b",c(0.5,1.75),TRUE,ylab="Consumer biomass",xlab="Consumer growth",main="Aquatic food webs")
plot(1:10,axes=FALSE,col=0,ylab="",xlab="")
legend("top",c("Natural state","Collapsed state","Tipping point","Resilience"),cex=1.2,box.col=0,col=c(4,1,2,"blue"),lwd=3,seg.len=2)












