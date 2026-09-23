


#Script 5: Run model fitting, analysis, and validation (Figures 2-4)


library(deSolve); library(nloptr); library(mgcv); library(GPEDM); library(abind); library(MASS); library(doParallel); library(foreach); library(domir); #library(home);



setwd("/Users/aqeco/Dropbox/Projects/KelpAllee/CodebaseFinal/Data")
gamdatFutAll=readRDS("ProcessedData\\gamdatFutAll.rds")
gamdatFutCWAll=readRDS("ProcessedData\\gamdatFutCWAll.rds")
alldat0=readRDS("ProcessedData\\alldat0.rds");
usf=import('LTM data\\KFM_HabitatSizeFrequencyFull.csv')
ARMs=as.matrix(import("LTM data\\KFM_ARMs_RawData_1992-2017.csv")[,c(1,6,9,11:13)])


#### CORE MODEL  ####
datFitMake2=function(parmsG,gamGive=FALSE,plotgive=0,temp=0){
     ntp1=gamdat$res; nt=gamdat$restm1; C=parmsG["a"];
     Et=(C+(ntp1+C)/(nt+C)) / (parmsG["b"]-nt)
     gamRes=gam(Et~s(no3s,k=5)+s(n2,k=8)+s(tempL1,k=6)+s(tempL2,k=8)+s(tempL3,k=8)+s(NPGO,k=4),data=gamdat); rsq(gamRes);
     if(plotgive){ par(mfrow=c(2,3),mar=c(2.1,2.1,2.1,2.1)); plot(gamRes,shade=TRUE,shade.col="gray",las=1,ylim=c(-0.4,0.4)); par(mfrow=c(1,1)); }

     if(gamGive) return(gam(Et~s(no3s,k=5)+s(tempL1,k=6)+s(tempL2,k=8)+s(tempL3,k=8),data=gamdat))
     gamdat$resEpred=round(nrm(predict.gam(gamRes,newdata=gamdat[,c("temp","tempL1","tempL2","tempL3","no3s","NPGO","no3sL1","no3","no3L1","n2","n5","MEI")],type="response"))^parmsG["EXP"],2)
     NSA=parmsG["NSA"]; KI=makelags(gamdat,"res","trans",NSA+1,tau=1,time="yr")[,NSA]; 
     Elag=makelags(gamdat,"resEpred","trans",NSA+1,tau=1,time="yr")[,0:(NSA-1)];
     gamdat$Tpast=rowMeans(t(t(gamdat[,c("temp","tempL1","tempL2","tempL3","tempL4","tempL5")[pmax(temp,1)]])))
     datFit=cbind(U=gamdat$lusm,K=gamdat$res,KI=KI,Y=gamdat$yr,E=gamdat$resEpred,Elag,trans=gamdat$trans); 
     du=complete.cases(datFit); datFit=datFit[du,];
     if(temp[1]>0) datFit=cbind(datFit,T=gamdat$Tpast[du],gamdat[du,c("nPixColL1","nPixRecL1","nPixForL1","nPixBarL1","no3s","temp","tempL1","tempL2","tempL3","tempL4")])
     as.matrix(datFit);
}
odeFull1=function(t,start,parms,Vari,PE=0){    #Vari=c(lusm,envt)
     MIN=1e-6; n=pmax(start,MIN); Vari=matrix(Vari,length(start),2);
     r=parms["r"]*pmax(1 + parms["Eff"]*Vari[,2],0.01)*(1+PE)
     graze=Vari[,1]*parms["delta"]/(1 + parms["fg"]*(n^2))
     return(list(n*(r*(1-n)/(1+graze*parms["RecEat"]) - graze)))
}
stepIterv=function(Var,parms,PE=0,ICi=Var[,3],NT=1,NE=4,NS=2){ 
     for(t in 1:parms["NSA"]) ICi=lsoda(ICi,c(0,NT),odeFull1,parms=parms,Vari=Var[,c(1,NE+t)],PE=PE)[2,-1]; round(ICi,NS); }
rlnormT=function(n,mu,sig) rlnorm(n, meanlog=log(mu/sqrt(1+(sig/mu)^2)), sdlog=sqrt(log(1+(sig/mu)^2)) )
LoglikFunLS=function(parmsTry,parmsFit=parmnames[pfit],parmsG,fitgive=0,parmsUse=NULL,datFit=NULL){
     if(fitgive==0) print(round(c(parmsTry,Bst=REPORT[which.min(REPORT[,ncol(REPORT)]),c(head(1:ncol(REPORT),-7),ncol(REPORT))]),3))
     if(!is.null(parmsUse)) parms=parmsUse; parms[parmsFit]=parmsTry; fg0=parms["fg"]; parms["fg"]=exp(fg0); if(fg0<=0) parms["fg"]=0;
     if(is.null(datFit)) datFit=datFitMake2(c(parms,parmsG)); #print(parms);
     Preds=stepIterv(datFit,c(parms,parmsG),PE=0);
     pASS=NA; LL=sum((Preds-datFit[,"K"])^2);
     # if(fitgive==1){ datASS=rbind(datFit,datFit); datASS[,3]=rep(c(0.01,1.5),each=nrow(datFit));
     #      pASS=mean(apply(matrix(stepIterv(datASS,c(parms,parmsG),NT=50),ncol=2),1,diff)>0.05); }
     resm=cbind(clspred=Preds,datFit[,c("K","KI")]); sresm=resm>0.5;
     dfa=head(serialAgg(cbind(datFit,pred=resm[,"clspred"]),4,FUN=mean),-1)
     sr2m=cor(sresm)^2; r2shy=cor(resm[,"clspred"],resm[,"KI"])^2;
     pEQo=mean(sresm[,"clspred"]==sresm[,"K"]); pEQi=mean(sresm[,"clspred"]==sresm[,"KI"]); iEQo=mean(sresm[,"KI"]==sresm[,"K"]);
     stats=c(pASS=pASS,r2r=round(cor(dfa)["pred","K"]^2,2),sr2shy=sr2m["KI","clspred"],sr2=sr2m["K","clspred"],icsr2=sr2m["K","KI"],pEQi=pEQi,pEQo=pEQo,iEQo=iEQo,LL=LL); 
     if(fitgive==105){ return(PEsetto0=Preds); }
     if(fitgive>0){ print(stats); return(list(stats,resm[,"clspred"])); }
     parms["fg"]=fg0; assign("REPORT",rbind(REPORT,c(parms,stats)),envir=.GlobalEnv); LL;
}
fitrun=function(pfit,parmsG,Ntries=300,LB=lowerBex,UB=upperBex,parmInit=parmStr,GA=1,LLfun=LoglikFunLS,ptryLoc=0.15,plotReg=TRUE,datFit=NULL){
     if(is.null(datFit)) datFit=datFitMake2(parmsG)
     g=fit=nloptr(parmsFit=parmnames[pfit], opts=list("algorithm"=c("NLOPT_GN_DIRECT","NLOPT_GN_ISRES","NLOPT_GN_BFGS")[GA],"xtol_rel"=1e-5,"maxeval"=ceiling((1-ptryLoc)*Ntries),"print_level"=2), x0=parmInit[pfit],eval_f=LLfun,lb=LB[pfit],ub=UB[pfit],fitgive=FALSE,parmsG=parmsG,parmsUse=parmInit,datFit=datFit);
     if(ptryLoc>0) fit=nloptr(parmsFit=parmnames[pfit], opts=list("algorithm"="NLOPT_LN_COBYLA","xtol_rel"=1e-5,"maxeval"=ceiling(ptryLoc*Ntries),"print_level"=2), x0=g$solution,eval_f=LLfun,lb=LB[pfit],ub=UB[pfit],fitgive=FALSE,parmsG=parmsG,parmsUse=parmInit,datFit=datFit)
     parmsBest=parmInit; parmsBest[pfit]=fit$solution; parmsBest["fg"]=exp(parmsBest["fg"]);
     preds=LLfun(parmsTry=fit$solution,parmsFit=parmnames[pfit],parmsG=parmsG,fitgive=(1+plotReg),parmsUse=parmInit,datFit=datFit)
     
     dfa=head(serialAgg(cbind(datFit,pred=preds[[2]]),which(colnames(datFit)=="Y"),FUN=mean),-1); 
     CR=cor(dfa)["K",ncol(dfa)]; r2r=sign(CR)*CR^2; rSSE=sum((dfa[,"K"]-dfa[,"pred"])^2);
     if(plotReg) matplot(dfa[,c("K","pred")],type="l",main=r2r)
     list(pfit=pfit,short=round(c(preds[[1]],rSSE=rSSE,r2r=r2r,rmin=min(dfa[,"pred"])),2),parms=c(parmsBest,parmsG),Set=c(GA,Ntries),B=rbind(LB,UB),dfa=dfa);
}





#Differences in parameter names between paper and code:
#r=r_0, RecEat=delta_R, Eff=beta, fg=f, and EXP=theta.
#reps and sigma deprecated

parmsG=c(NSA=2,reps=10,EXP=1.88,a=0.32,b=2.2); datFit=datFitMake2(parmsG);
#varType says what PE affects: 1=r; 2=delta; 3=fg; 4=tau; 5=mu; NOT FITTED.
parms=parmStr=c(r=10,delta=8.6,fg=0,  Eff=3,  sigma=0.345,RecEat=1.6);
upperBex=c(15,25,3,  30,  0.5,3); lowerBex=c(2,8,-0.01,  0,  0.001,1); parmnames=names(lowerBex)=names(upperBex)=names(parms);
REPORT=matrix(nrow=0,ncol=length(parms)+9) #Track results: ncol = #parms + # stats in LoglikFun


ODEfit=fitrun(pfit=c(1:4,6),Ntries=2500,LLfun=LoglikFunLS,ptryLoc=0.2,plotReg=FALSE,LB=lowerBex,UB=upperBex,parmInit=parmStr,parmsG=parmsG); 
lowerBex[c("r","delta")]=c(8,1); upperBex["delta"]=15; 
NULLfit=fitrun(pfit=c(1,2,4,6),Ntries=1800,LLfun=LoglikFunLS,ptryLoc=0.15,plotReg=FALSE,LB=lowerBex,UB=upperBex,parmInit=parmStr,parmsG=parmsG);
lowerBex[c("r","delta")]=c(2,8); upperBex["delta"]=25;

#HERE'S THE KEY METRIC: HOW MUCH BETTER DO WE PREDICT 2-YR-AHEAD REGIONAL CHANGES VS A RANDOM WALK?
rp=ODEfit$dfa; (cor(apply(rp[,c("K","pred")],2,diff))[1,2]^2) - (cor(apply(rp[,c("K","KI")],2,diff))[1,2]^2)
rp=NULLfit$dfa; (cor(apply(rp[,c("K","pred")],2,diff))[1,2]^2) - (cor(apply(rp[,c("K","KI")],2,diff))[1,2]^2)


fiti=ODEfit
fitinull=NULLfit










#### ANALYSES ####
###### Role of TPs ####
fiteval2=function(fiti,plotgive=1){
     pfiti=fiti$pfit; pf=parmnames[pfiti]; parmsi=pg=fiti$parms;
     pt=parmsi[pf]; pt["fg"]=log(pt["fg"]); pg["NSA"]=1; if(is.na(pt["fg"])) pt["fg"]=0;
     datFiti=cbind(datFitMake2(c(NSA=1,pg),temp=1:5), pred=LoglikFunLS(parmsTry=pt,parmsFit=names(pt),parmsG=pg,parmsUse=parmStr,fitgive=99)[[2]])
     if(plotgive==0) return(datFiti)
     sfa=aggregate(datFiti[,c("K","pred")]~datFiti[,"Y"],FUN=mean,na.rm=TRUE)
     matplot(sfa[,1],1-sfa[,-1],type="l",lty=1,lwd=3,col=c(1,4),main="Barren state frequency"); datFiti;
}
stateThresh=0.4; giveSmoo=function(x,info,ic,kmx=5,Q=0.05){ 
     temp=info[,"T"]; tr=quantile(temp,c(Q,1-Q),na.rm=TRUE);
     us=ic<stateThresh & !is.na(temp); R=x[us]>stateThresh; md=glm(R~temp[us],binomial); #print(rsq(md)); 
     cbind(temp[us],1/md$fitted.values)[temp[us]>tr[1] & temp[us]<tr[2],];
}


layout(matrix(c(1,1,2,3),nrow=2,byrow=TRUE))
par(oma = c(1, 1, 1, 1),   # outer margin of whole figure
    mar = c(3, 3, 2, 1))   # inner margin around each plot
mdbBG=fiteval2(fiti); nllBG=fiteval2(fitinull,plotgive=0)

obs=mdbBG[,"K"]; ic=mdbBG[,"KI"]; pred=mdbBG[,"pred"]; pnull=nllBG[,"pred"]; info=mdbBG[,c("trans","Y","T")]; #pred=mdbo[,"pred"]
info[,"T"]=datFitMake2(c(a=0.32,EXP=1.88,b=2.2,NSA=1),temp=4)[,"T"]; 
info=cbind(info,rec=gamdat$recL1[match(codify(info[,1:2]),codify(gamdat[,c("trans","yr")]))])
info=cbind(info,npx=gamdat$nsatpix[match(codify(info[,1:2]),codify(gamdat[,c("trans","yr")]))])
plot(giveSmoo(obs,info,ic),ylim=c(1.5,10),las=1); points(giveSmoo(pred,info,ic),col=4); points(giveSmoo(pnull,info,ic),col=3); #points(giveSmoo(info[,"rec"]),col=4)
sfa=aggregate(cbind(obs,pred,pnull)>stateThresh~info[,"Y"],FUN=mean,na.rm=TRUE); 
boxplot(1-sfa[,2:4],horizontal=FALSE,col=c(8,4,3),las=1); #cor(sfa)["pred","obs"]^2
par(mfrow=c(1,1));








####### Figure 2,4B ####
odePot2=function(n=NA,parmsi=parmsPot,Var,type=1){ #Var=c(lusm,E)
     r=ltrim(parmsi["r"]*(1+parmsi["Eff"]*Var[2])*(Var[2]>0),potminr); f=parmsi["fg"]; m=0;
     d=parmsi["delta"]*Var[1]; dr=parmsi["RecEat"]*d;
     if(type==2) n=seq(0,1,len=1e3); if(sum(is.na(c(n,parmsi,Var)))>0) return(NA);
     # antiderivative rn(1-n)/(1+dr/(1+fn^2)) - nd/(1+fn^2) - mn with respect to n
     npot=(-3*f*(dr*r*log(dr+1+f*n^2) + d*log(1+f*n^2)) + n*r*f*(6*dr+n*f*(3-2*n) - 3*m/r) - 6*r*dr*sqrt(f*(1+dr))*atan(n*sqrt(f/(1+dr)))) / (6*f^2)
     
     if(type==2){ rpot=npot-min(npot); tip=which.min(rpot); return(sum(c(0,rpot[-(1:tip)]))/sum(rpot)); }
     ddv=n[diff(sign(diff(-npot)))==-2]; if(length(ddv)==0) if(min(diff(-npot))>0) ddv=1 else ddv=0; if(type==3) return(round(ddv,2)); npot
}
potplot2=function(Us=2.5,Ts=seq(15,20,len=5),parmsPot,axscl=c(1,1,0.15,1)){
     pn=100; pots=-apply(cbind(Us,e.t3(Ts)),1,function(x) odePot2(n=seq(0-axscl[3],1+axscl[3],len=pn),parms=parmsPot,Var=x))
     if(axscl[4]==1) potsn=apply(pots%*%diag(1/pots[1,]),2,rev)
     if(axscl[4]==2) potsn=apply(apply(pots,2,nrm,TRUE),2,rev)
     madj=axscl[2]*0.4*min(apply(potsn,2,function(x) diff(range(x)))); 
     radj=range(potsn)+0.5*madj*c(-1,1)/axscl[2];
     tn=ncol(potsn); typ=1+(length(Ts)==1); slot_h=diff(radj)+madj/2;
     matplot(potsn,col=0,ylim=c(0,slot_h*tn),xlim=c(pn,0),main=paste0("At ",list(round(exp(Us),1)-1,Ts)[[typ]],c(" urchins per m2"," degrees C")[typ]))
     shades=adjustcolor(colorRampPalette(c("dodgerblue","red"))(tn), alpha.f=0.5) #gray.colors(tn,start=0.7,end=0.95,gamma=2.2,rev=FALSE);
     for(i in 1:tn) polygon(c(1,1:pn,pn), (i-1)*slot_h-radj[1] + c(radj[1],potsn[,i],radj[1]),col=shades[i])
}


#Figure 2b: Plot potentials for a given urchin or temp level
parmsPot=fiti$parms[c("r","delta","fg","Eff","RecEat")]; potminr=0.01; #parmsPot[1];
datFiti=datFitMake2(fiti$parms,temp=1:5)
etg3=gam(E~T,data=data.frame(datFiti)); rsq(etg3); #plot(etg3)
e.t3=function(T) pmax(as.numeric(predict.gam(etg3,newdata=data.frame(T=T),type="response")),0.12)
# potplot2(Us=1.8,Ts=c(15,18,19,20.5),parmsPot,axscl=c(1,1,0.1,2)) #shows Heat stress collapse
potplot2(Us=1.5,Ts=c(10,18,19.5,20.5),parmsPot,axscl=c(1,1,0.1,2)) #shows TP creation



#Figure 4b: Plot attraction basin volume across urchins and temps
nt=200; Useq=seq(0.8,2.8,len=nt); Tseq=seq(15.5,20.5,len=nt); conds=as.matrix(expand.grid(U=Useq,T=e.t3(Tseq)));
pn=100; potsM=-apply(conds,1,function(x) odePot2(n=seq(0,1,len=pn),parmsi=parmsPot,Var=x))
resM=apply(potsM,2,function(x){ sg=sign(diff(x)); if(sd(sg)==0) return(sg[1]==-1); xn=1-nrm(x); 1-sum(xn[1:which(sg!=1)[1]])/sum(xn); })
resM=matrix(resM,nt,nt); SpatiotempPlot(resM,pal=3,Xax=Tseq,Yax=(Useq),Conts=c(0.05,0.95),cont=2,cexCont=0.5,las=1)
# points(datFiti[,c("T","U")],cex=0.2,col=6,pch=16); round(exp(c(1,1.5,2,2.5)),1)
points(cbind(0.5+c(17,19.5,19.5),c(1.65,1.65,log(1+(exp(1.65)-1)/2))),lwd=3,col=6,cex=1.5)
# for(i in unique(datFiti[,"trans"])){ di=datFiti[datFiti[,"trans"]==i & datFiti[,"Y"]<2027,]; xr=range(di[,"T"],na.rm=TRUE); y=di[1,"U"]; segments(lwd=1/2,xr[1],y,xr[2],col=6); }



#Figure 2c
ulv=exp(Useq)-1; filts=function(x,s=3) as.numeric(filter(x,rep(1/s,s)))
tHSC=filts(apply(resM,1,function(x) Tseq[which(x<0.001)[1]]))
tRLC1=filts(apply(resM,1,function(x) Tseq[which(x<0.33)[1]]))
tRLC2=filts(apply(resM,1,function(x) Tseq[which(x<0.48)[1]]))
# plot(ulv,tRLC2,type="l",las=1,lwd=3,ylim=c(15,21.5)); lines(ulv,tHSC,col=2,lwd=3)
plot(ulv,tRLC1,type="l",las=1,lwd=4,ylim=c(15.75,20.26),col=4,xlim=c(0,15)); lines(ulv,tHSC,col=2,lwd=4)
rug(jitter(utrim(exp(unique(datFiti[,"U"]))-1,14.75),200)); rug(datFiti[,"T"],side=2)


# Figure 2d
gdFR=function(gamdatFuti){
     SET=c(CritNO3 = 0.4, CritTemp = 200, Ufct = 1, reps = 30, Uyr = 2040)
     rngi=range(mdli$fitted.values); Epred=predict.gam(mdli,newdata=gamdatFuti,type="response");
     gamdatFuti$E=(pmax(Epred-rngi[1],0)/(rngi[2]-rngi[1]))^parmsi["EXP"]
     nmx=gamdatFuti$no3w; lo=nmx<gamdatFuti$no3s; nmx[lo]=gamdatFuti$no3s[lo]; gamdatFuti$E[nmx<SET["CritNO3"] | gamdatFuti$temp>SET["CritTemp"]]=0;
     gamdatFuti$U=log(1+((exp(gamdatFuti$U)-1)*pmin(pmax(1-(1-SET["Ufct"])*(gamdatFuti$yr-2019)/(SET["Uyr"]-2019),SET["Ufct"]),1)))
     
     pn=100; potsDat=-apply(gamdatFuti[,c("U","E")],1,function(x) odePot2(n=seq(0,1,len=pn),parmsi=parmsPot,Var=x))
     gamdatFuti$res=apply(potsDat,2,function(x){ sg=sign(diff(x)); if(sd(sg)==0) return(sg[1]==-1); xn=1-nrm(x); 1-sum(xn[1:which(sg!=1)[1]])/sum(xn); })
     t(matricize(gamdatFuti[,c("trans","yr","res")]))
}
parmsi=fiti$parms; pfit=fiti$pfit; mdli=datFitMake2(fiti$parms,gamGive=TRUE); datFiti=datFitMake2(c(NSA=1,fiti$parms),temp=1:5);
pf=parmnames[pfit]; pt=parmsi[pf]; pt["fg"]=log(pt["fg"]); pg=parmsi; pg[c("NSA","reps")]=c(1,2); datFiti=cbind(datFiti,pred=LoglikFunLS(parmsTry=pt,parmsG=pg,fitgive=99,datFit=datFiti)[[2]]); 
pn=100; potsDat=-apply(datFiti[,c("U","E")],1,function(x) odePot2(n=seq(0,1,len=pn),parmsi=parmsPot,Var=x))
datFiti=cbind(datFiti,res=apply(potsDat,2,function(x){ sg=sign(diff(x)); if(sd(sg)==0) return(sg[1]==-1); xn=1-nrm(x); 1-sum(xn[1:which(sg!=1)[1]])/sum(xn); }))
gdrdat=t(matricize(datFiti[,c("trans","Y","res")])); 


first_n_consecutive_after=function(x,n,y,m){
     make_runs=function(v,len){ r=rle(v); ends=cumsum(r$lengths);
     list(starts=ends - r$lengths + 1L, ends=ends, qual=which(r$values & r$lengths >= len)); }
     rx=make_runs(x,n); ry=make_runs(y,m);
     if(!length(rx$qual)) return(x*NA)
     
     x_search=rx$qual; y_ends=c(0L, if(length(ry$qual)) ry$ends[ry$qual]);
     out=unlist(lapply(y_ends, function(y_end) {
          hit=x_search[rx$starts[x_search] > y_end][1]
          x_search <<- x_search[x_search > hit]; rx$starts[hit]; }))
     head(c(out,NA*x),length(x))
}
gd=gamdatFutAll[[2]]; gdr=rbind(gdrdat,gdFR(gd[[1]])); for(i in 2:length(gd)) gdr=abind(gdr,rbind(gdrdat,gdFR(gd[[i]])),along=3);
lg=4; TH=1/3; recp=c(1-TH,lg); NB=500; ys=as.numeric(rownames(gdr)); # gdr=gdr2 #gdr=gdr3 
rlcy=t(apply(gdr,2,function(x) c(na.omit(as.numeric(apply(x,2,function(xi) first_n_consecutive_after(xi<TH,lg,xi>recp[1],recp[2])))),rep(NA,NB))[1:NB]))
hscy=t(apply(gdr,2,function(x) c(na.omit(as.numeric(apply(x,2,function(xi) first_n_consecutive_after(xi==0,lg,xi>recp[1],recp[2])))),rep(NA,NB))[1:NB]))
susi=serialAgg(datFiti,"trans",c("K","U"),FUN=mean); susi=susi[order(susi[,1]),];
sus=susi[,3]>0.85; hscy[!sus,]=NA; rlcy[!sus,]=NA
rlcy[rlcy==hscy]=NA; startdie=which(hscy[,1]==1); rlcy[startdie,]=NA; hscy[startdie,]=NA; mean(!sus)+length(startdie)/nrow(rlcy);
# par(mfrow=1:2); htp=function(tp) hist(ys[tp],breaks=50,main=mean(is.na(tp)),xlim=range(ys)); htp(rlcy); htp(hscy);
par(mfrow=c(1,1)); boxplot(cbind(ys[asvt(rlcy)],ys[asvt(hscy)]),col=c(4,2),ylim=c(1988,2100),main= 1 - (mean(is.na(rlcy[,1])) - mean(!sus) - length(startdie)/nrow(rlcy)));
# matplot(unique(gd[[1]][,"U"]),cbind(rowMeans(rlcy,na.rm=TRUE),rowMeans(hscy,na.rm=TRUE)),ylim=c(0,110)); points(cbind(unique(gd[[1]][,"U"]),10),col=4)







##### RLC vs HSC ####
prm=fiti$parms; prm["NSA"]=1; datFit=datFitMake2(prm,temp=1:5); VAR=datFit[,c("U","E")];
ns=seq(0,1,len=1e3); TPs=1-apply(VAR,1,function(x) odePot2(ns,parms=prm,Var=x,type=3))

TPtrans=aggregate(TPs~datFit[,"trans"],FUN=mean)[,2]; TPf=TPs; TPf[TPtrans<0.05 | TPtrans>0.95]=NA; #Omit sites where model never predicts state shifts #hist(TPf,breaks=1e2)
TPA=aggregate(TPs~datFit[,"Y"],FUN=mean)

#How many reefs experienced creation of ASS
ASScre=aggregate(cbind(datFit[,"U"],TPs==0,(TPs>0 & TPs<1))~datFit[,"trans"],FUN=max); sum(ASScre[,3]==1 & ASScre[,4]==1); sum(ASScre[,3]==1 & ASScre[,4]==1)/length(unique(datFit[,"trans"]));


#Next, ID which collapses fall in (a) ASS region and (b) urchin-barren-only region
TPsCol=cbind(TPs,gamdat[match(codify(datFit[,c("Y","trans")]),codify(gamdat[,c("yr","trans")])),c("nPixCol","nsatpix","trans","yr")])
TPsCol$col=(TPsCol$nPixCol/TPsCol$nsatpix)>=0.45; nc=sum(TPsCol$col);
TPsCol$HSC=TPsCol$col==1 & TPsCol$TPs==1; TPsCol$RLC=TPsCol$col==1 & TPsCol$TPs>0 & TPsCol$TPs<1;
nc; tail(colSums(TPsCol,na.rm=TRUE),2)/nc; #colMeans(aggregate(cbind(HSC,RLC)~trans,FUN=sum,TPsCol)[,-1]>0)
sum(TPsCol$col==1 & TPs==0)/nc; sum(TPsCol$col==1 & TPs==0 & TPsCol$yr>2014)/nc; #collapses not explained by RLC or HSC







###### Figure 3 ####
#Plotting temp effects on growth rate and resilience
par(mfrow=1:2); kmx=4; XLM=c(15.5,20.5); LW=5; CX=0.5; PH=1; S=0.7; ptCol=rgb(0.2,0.6,.9,alpha=0.25);
datFit=cbind(datFit,E2=prm["r"]*(1+prm["Eff"]*datFit[,"E"]));
plot(datFit[,c("T","E2")],cex=CX,col=ptCol,pch=PH,xlim=XLM,ylim=c(0,150),las=1,xlab="5yr temp",ylab="r(Et)"); 
lines(gampred(datFit[datFit[,"Y"]<2017,c("T","E2")],1,2,k=kmx,plotgive=0)$mu,lwd=LW)
potsE=-apply(VAR,1,function(x) odePot2(ns,parms=prm,Var=x))
resE=apply(potsE,2,function(x){ sg=sign(diff(x)); if(sd(sg)==0) return(sg[1]==-1); xn=1-nrm(x); 1-sum(xn[1:which(sg!=1)[1]])/sum(xn); })
plot(datFit[,"T"],ltrim(utrim(jitter(resE,amount=0.02),1),0),cex=CX,col=ptCol,pch=PH,xlim=XLM,las=1,xlab="5yr temp",ylab="Forest attraction basin volume"); 
lines(gampred(cbind(datFit[,"T"],resE),1,2,k=kmx,plotgive=0)$mu,lwd=LW)
par(mfrow=c(1,1))


mda=serialAgg(mdbBG,"Y","pred",FUN=mean);
gamdati=gamdat; gamdati$luc=nrm2(gamdat$lu,Lvl=gamdat$trans,center=TRUE); gamdati$luc[gamdati$trans>20]=NA;
gda=serialAgg(gamdati,"yr",c("res","luc"),FUN=mean); gda[,2]=1-gda[,2]; gda=cbind(gda,NA);
gda[gda[,1]>=1988,4]=1-mda[,2]; colnames(gda)=c("y","o","b","p");

plot(gda[,c("y","b")],col=2,lwd=0,type="l",ylim=c(-1,0.7),xlim=c(1984,2020),las=1)
bp=!is.na(gda[,"b"]); polygon(c(gda[bp,"y"],rev(gda[bp,"y"])),c(gda[bp,"b"],gda[bp,"b"]-500),col=rgb(0.87,0.32,0.42,alpha=0.25),border=NA)
plot2(gda[,c("y","o")],cex=0,col=1); points(gda[,c("y","o")],type="o",pch=16,col=8);
lines(gda[,"y"],.06+gda[,"p"],lwd=3,col=4)













##### Grazing index ####
par(mfrow=1:2,mgp=c(2.25,0.75,0)); X=jitter(gamdat$res,4); CL=adjustcolor(4,alpha.f=0.2);
Qylm=0.1; V=exp(gamdat$lu)-1; gampred(cbind(X,V),1,2,k=3,spanSignif=NA,seCol=1,seFill=0.4,Lcol=1,col=CL,pch=16,las=1,ylim=quantile(V,c(Qylm,1-Qylm),na.rm=TRUE),ylab="Detected urchin density",xlab="Forest cover")
legend("topright",expression("Mean trend" %+-%"95%CI"),col=1,lwd=4,cex=0.9,box.col=0,seg.len=0.8)

Qylm=0.25; V=nrm2(gamdat$lu,Lvl=gamdat$trans,center=TRUE); 
gampred(cbind(X,V),1,2,k=3,spanSignif=NA,seCol=2,seFill=0.4,Lcol=2,col=CL,pch=16,las=1,ylim=quantile(V,c(Qylm,1-Qylm),na.rm=TRUE),ylab="Reef-level anomaly in \n log deteced urchin density",xlab="Forest cover")
legend("topright",expression("Mean trend" %+-%"95%CI"),col=2,lwd=4,cex=0.9,box.col=0,seg.len=0.8)

xmd=gampred(cbind(X,V),1,2,plotgive=FALSE)
plot(xmd$mu[,1],nrm(xmd$mu[,2]),type="l",lwd=LW,col=2,las=1,xlab="Forest cover",ylab="Grazing activity");
rng=range(xmd$mu[,2]); polygon(xmd$polyse[,1],(xmd$polyse[,2]-rng[1])/diff(rng),border=NA,col=adjustcolor(2,alpha.f=0.4))
xs=(0:100)/100; lines(xs,1/(1+parmsi["fg"]*xs^2),col="purple",lwd=4)
par(mfrow=c(1,1)); legend("topright",c(expression("Mean trend" %+-%"95%CI"),"Assumed Type IV","functional response"),col=c("#DF536B","purple","white"),lwd=4,cex=0.9,box.col=0,seg.len=0.8)




##### Dens vs Beh ####
#####Parsing role of behavior vs recruitment pulses
#Role of behavior plots
Evars=c("temps","temp","tempL1","tempL2","tempL3","MEI","NPGO","PDO","n2","no3","no3L1")
gamdatFMD=na.omit(gamdat[gamdat$nsatpix>4,c("lusm","lu","kp","trans","yr","res","restm1",Evars)])
gamdatFMD$uresid2=exp(gamdatFMD$lu)-exp(gamdatFMD$lusm)
# gamdatFMD$uresid2=nrm2(exp(gamdatFMD$lu),gamdatFMD$trans) #z-scores don't work
gamdatFMD$lusmn=exp(gamdatFMD$lusm)-1
gamdatFMD$targetFactor=1/(1+parmsPot["fg"]*gamdatFMD$kp^2)
# gamdatFMD[,c("res","restm1",Evars)]=data.frame(apply(gamdatFMD[,c("res","restm1",Evars)],2,nrm2,Lvl=gamdatFMD$trans,center=FALSE))


gamdatFMD$pu=log(1+alldat0[match(codify(gamdatFMD[,c("trans","yr")]),codify(alldat0[,c("Site","Yr")])),"PrpU"])
gamdatFMD$ru=log(1+alldat0[match(codify(gamdatFMD[,c("trans","yr")]),codify(alldat0[,c("Site","Yr")])),"RedU"])


szsCT=15 #Maximum test diameter to consider urchin as a "recruit"
usf=usf[usf$Species%in%c(11005,11006) & !usf$SiteNumber%in%c(17:20,2001:3003,37),]; usfe=usf[rep(1:nrow(usf),usf$NoOfInd),];
xi=serialAgg(usfe,c("Species","SiteNumber","SurveyYear"),"NoOfInd",CatsCode=FALSE)
xi=cbind(xi,  szs=serialAgg(usfe,c("Species","SiteNumber","SurveyYear"),"Size_mm",CatsCode=FALSE,FUN=function(x){ if(length(x)<30) return(NA); mean(x<szsCT); })[,4])
qtrim=function(x,Q=0.99){ Qs=quantile(x,c(1-Q,Q),na.rm=TRUE); utrim(ltrim(x,Qs[1]),Qs[2]); } #Prevent outliers from wrecking recruitment-density relation
gamdatFMD$szsR=qtrim(xi[xi[,"Species"]==11005,][match(codify(gamdatFMD[,c("trans","yr")]),codify(xi[xi[,"Species"]==11005,c("SiteNumber","SurveyYear")])),"szs"],.95)
gamdatFMD$szsP=qtrim(xi[xi[,"Species"]==11006,][match(codify(gamdatFMD[,c("trans","yr")]),codify(xi[xi[,"Species"]==11006,c("SiteNumber","SurveyYear")])),"szs"],.95)


x=ARMs[ARMs[,"Species"]%in%c(11005,11006) & !ARMs[,"SiteNumber"]%in%c(17:20,2001:3003,37),]; narms=serialAgg(x,c("Species","SiteNumber","SurveyYear"),"ArmNo",FUN=function(xi) length(unique(xi)),CatsCode=FALSE)[,4];
serialAgg(x,c("SiteNumber"),"ArmNo",FUN=function(xi) length(unique(xi)),CatsCode=FALSE);
xa=serialAgg(cbind(x,sn=x[,"Size_mm"]*x[,"NoOfInd"]*(x[,"Size_mm"]<szsCT),rec=x[,"NoOfInd"]*(x[,"Size_mm"]<szsCT)),c("Species","SiteNumber","SurveyYear"),c("NoOfInd","sn","rec"),FUN=sum,CatsCode=FALSE)
xa=cbind(xa[,1:4],szs=xa[,"sn"]/xa[,"rec"],rec=xa[,"rec"]/xa[,"NoOfInd"],nrec=log(1+xa[,"rec"]/narms)); serialAgg(xa,"Species",c("szs","NoOfInd","rec","nrec"),FUN=mean);
getv=function(spp=11005,V="szs") xa[match(codify(cbind(spp,gamdatFMD[,c("trans","yr")])), codify(xa[,c("Species","SiteNumber","SurveyYear")])),V]
gamdatFMD$nrecParm=getv(11006,"nrec"); gamdatFMD$nrecRarm=getv(11005,"nrec"); gamdatFMD$nrecAarm=gamdatFMD$nrecParm+gamdatFMD$nrecRarm;


SpDensCVs=apply(serialAgg(gamdatFMD,"trans",c("ru","pu"),FUN=sd),2,mean,na.rm=TRUE)
Vp="nrecParm"; Vr="nrecRarm"; 
ulp=makelags(gamdatFMD,Vp,E=3,pop="trans",tau=1); sel=!is.na(ulp[,1]);
y=domin(gamdatFMD$pu[sel]~gamdatFMD$res[sel]+ulp[sel,1], lm, list(summary,"r.squared"),complete=F,conditional=F); c(y$Standardized[1],y$Fit_Statistic_Overall)
gampred(cbind(ulp[,1],gamdatFMD$pu)[sel,],1,2,k=3,spanSignif=1.2,las=1,col="purple",xlab="Recruit density last year",ylab="log Detected >25mm urchin density now"); 
gampred(gamdatFMD[sel,],"res","pu",k=3,spanSignif=1.2,las=1,col="purple",xlab="Proportion reef forested",ylab="log Detected >25mm urchin density now")


ulr=makelags(gamdatFMD,Vr,E=3,pop="trans",tau=1); sel=!is.na(ulr[,1]);
y=domin(gamdatFMD$ru[sel]~gamdatFMD$res[sel]+ulr[sel,1], lm, list(summary,"r.squared"),complete=F,conditional=F); c(y$Standardized[1],y$Fit_Statistic_Overall)
gampred(cbind(ulr[,1],gamdatFMD$ru)[sel,],1,2,k=3,spanSignif=1.2,las=1,col=2,xlab="Recruit density last year",ylab="log Detected >25mm urchin density now"); 
gampred(gamdatFMD[sel,],"res","ru",k=3,spanSignif=1.2,las=1,col=2,xlab="Proportion reef forested",ylab="log Detected >25mm urchin density now")


#Relative dominance of kelp cover in predicting sizes of both urchins together
ul=makelags(gamdatFMD,"nrecAarm",E=3,pop="trans",tau=1); sel=!is.na(ul[,1]);
y=domin(gamdatFMD$lu[sel]~gamdatFMD$res[sel]+ul[sel,1], lm, list(summary,"r.squared"),complete=F,conditional=F); c(y$Standardized[1],y$Fit_Statistic_Overall)









#### FUT CLIMATES ####
##### PROJECT MODELS UNDER FUTURE CLIMATES
# sd(serialAgg(cbind(datFiti,err=abs(devsn)),"Y","err",FUN=mean)[,2])
# sd(serialAgg(cbind(datFiti,err=abs(devsn)),"trans","err",FUN=mean)[,2])
#Errors vary more by reef than by year... Hence, we line up process error by reef
simIter=function(datFiti,gamdatFutiAll,parmsi,simPE,seed,gdF){
     gamdatFuti=gdF(seed); gamdatFuti$devsn=gamdatFuti$devPE=0;
     devsn=datFiti[,"K"]-datFiti[,"pred"]; devsn=devsn-mean(devsn)-0.015; devsnr=cbind(pr=round(datFiti[,"pred"],1),devsn);
     YPs=2020:2097; YRs=rep(1987:2019,length(YPs))[1:length(YPs)]; 
     set.seed(seed); if(simPE==2) YRs=sample(YRs);
     for(yr in YPs[-1]){
          Vari=as.matrix(gamdatFuti[gamdatFuti$yr==yr,c("U","E","PE","trans")])
          predi=stepIterv(Var=Vari[,1:2],parms=parmsi,PE=0,ICi=gamdatFuti[gamdatFuti$yr==yr-1,"K"],NT=1,NE=1)
          if(simPE<3){ devsi=devsn[match(paste0(Vari[,"trans"],"_",YRs[YPs==yr]),codify(datFiti[,c("trans","Y")]))]; devsi[is.na(devsi)]=0; }
          if(simPE==3){ set.seed(seed+yr); devsi=sapply(round(predi,1),function(x) sample(devsnr[devsnr[,1]%in%(x+c(-1:1)*(x>0 & x<0.8)/10),2],1)); }
          predid=pmin(pmax(predi+devsi*(simPE>1),0),0.99); gamdatFuti$K[gamdatFuti$yr==yr]=predid; #Trims ~14% of values
          gamdatFuti$devsn[gamdatFuti$yr==yr]=devsi; gamdatFuti$devPE[gamdatFuti$yr==yr]=predid-predi;
     }; gamdatFuti$K;
     # sa=serialAgg(as.matrix(gamdatFuti),"yr",c("K","devsn","devPE"),FUN=mean)[,-1]; matplot(sa,type="l"); colMeans(sa)
}
simFun=function(SET=c(CritNO3=0.4,sstCrit=200,Ufct=1,Uyr=2020,reps=1),gamdatFutiAll,fiti,mdli,datFiti,plotgive=1,Rfct=1,simPE=0,XL=c(1980,2100)){
     parmsi=fiti$parms; parmsi["NSA"]=1; if(!is.null(dim(gamdatFutiAll))) gamdatFutiAll=list(gamdatFutiAll);
     gdF=function(rep=1){
          gamdatFuti=gamdatFutiAll[[utrim(rep,length(gamdatFutiAll))]];
          rngi=range(mdli$fitted.values); Epred=predict.gam(mdli,newdata=gamdatFuti,type="response");
          gamdatFuti$E=Rfct*(pmax(Epred-rngi[1],0)/(rngi[2]-rngi[1]))^parmsi["EXP"]
          # gamdatFuti$E[gamdatFuti$no3s<SET["CritNO3"] | gamdatFuti$temp>SET["CritTemp"]]=-1;
          nmx=gamdatFuti$no3w; lo=nmx<gamdatFuti$no3s; nmx[lo]=gamdatFuti$no3s[lo]; gamdatFuti$E[nmx<SET["CritNO3"] | gamdatFuti$temp>SET["CritTemp"]]=-1;
          gamdatFuti$U=log(1+((exp(gamdatFuti$U)-1)*pmin(pmax(1-(1-SET["Ufct"])*(gamdatFuti$yr-2030)/(SET["Uyr"]-2030),SET["Ufct"]),1)))
          gamdatFuti
     }
     # plot(serialAgg(as.matrix(gdF(12)),"yr","E",FUN=function(x) mean(ltrim(x,0))),ylim=c(0,0.6))
     Ksim=NA; for(i in 1:SET["reps"]) Ksim=cbind(Ksim,simIter(datFiti,gamdatFutiAll,parmsi,simPE=simPE,i,gdF))
     gamdatFuti=0*gdF(); for(i in 1:SET["reps"]) gamdatFuti=gamdatFuti+gdF(i)/SET["reps"]; gamdatFuti$K=apply(Ksim,1,median,na.rm=TRUE);
     
     gamdatAll=rbind(as.matrix(gamdatFuti)[,c("yr","K","U","E","temp", "no3","no3s","n5")],cbind(datFiti[,c("Y","K","U","E","T")],NA,NA,NA))
     dfa=serialAgg(gamdatAll,"yr",FUN=mean); dfa[dfa[,"yr"]%in%(1989:2019),"K"]=fiti$dfa[,"pred"]; if(plotgive==0) return(dfa);
}
allScenRun=function(SET=c(CritNO3=1,CritTemp=200,Ufct=1/2,reps=50,Uyr=2050),PE=0,NC=6,FIT=fiti,MDL=mdli,DF=datFiti){
     cnd=cbind(gdf=c(1,1,2,2,3,3),Ufct=c(1,SET["Ufct"])); if(NC>1) registerDoParallel(cores=NC);
     out=foreach(i=1:nrow(cnd),.combine=cbind,.packages=c("abind","mgcv","MASS","deSolve","GPEDM"),.export=ls()) %dopar% 
          simFun(c(cnd[i,"Ufct"],SET),gamdatFutAll[[cnd[i,"gdf"]]],FIT,MDL,DF,0,simPE=PE)[,c("yr"[i==1],"K")]
     if(NC>1) stopImplicitCluster(); out;
}


#####REGIONAL GLMS
coastSim2=function(sstfut,R=6,reps=1,CritNO3=0.4,CritTemp=200,subsn=300){
     subs=which(edatAnnUP[,"Reg","info"]==R); set.seed(1); if(!is.na(subsn)) subs=sample(subs,subsn);
     train=data.frame(apply(edatAnnUP[subs,,c("UPrecover","UPdie","temps","temp","tempL1","tempL2","tempL3","tempL4","no3s","no3w","UPbarren")],3,asvt))
     Vdie=c("tempL4","tempL1","temp")[R-3]; Vrec=c("temp","tempL2")[1+(R==5)];
     mdie1=glm(paste0("UPdie~",Vdie),binomial,train); mrec1=glm(paste0("UPrecover~",Vrec),binomial,train);
     Ys=2021:2097; NS=length(subs); out=array(NA,c(NS*(1+length(Ys)),reps));
     for(i in 1:reps){
          dfut=sstfut[[i]]/1000; dfut=dfut[dfut$trans%in%subs,]; dfut$K0=dfut$K;
          for(yr in Ys){
               tpf=function(mdl) ltrim(utrim(predict(mdl,type="response",newdata=dfut[dfut$yr==yr,]),1),0); set.seed(1+yr);
               KI=dfut[dfut$yr==yr-1,"K0"]; K0=KI + rbinom(NS,1-KI,tpf(mrec1)) - rbinom(NS,KI,tpf(mdie1));
               dfut[dfut$yr==yr,"K0"]=out[dfut$yr==yr,i]=K0;
          }; };
     outa=cbind(dfut$yr,apply(out,1,median))
     od=cbind(yr=rep(unique(floor(eyrs)),length(subs)),K=1-asvt(edatAnnUP[subs,,"UPbarren"]))
     x=serialAgg(rbind(od[,c("yr","K")],as.matrix(outa)),1,FUN=mean)
}
allScenRunGLM=function(CritNO3=0.4,CritTemp=200,reps=30,Rs=5:6,NC=length(Rs)*3){
     cnd=as.matrix(expand.grid(R=Rs,gdf=1:3)); if(NC>1) registerDoParallel(cores=NC);
     out=foreach(i=1:nrow(cnd),.combine=cbind,.packages=c("abind","mgcv","MASS","deSolve","GPEDM"),.export=ls()) %dopar% 
          coastSim2(subsn=1e3,gamdatFutCWAll[[cnd[i,"gdf"]]],cnd[i,"R"],reps,CritNO3,CritTemp)[,c(1[i==1],2)]
     if(NC>1) stopImplicitCluster(); out;
}




parmsi=fiti$parms; pfit=fiti$pfit; mdli=datFitMake2(fiti$parms,gamGive=TRUE); datFiti=datFitMake2(c(NSA=1,fiti$parms),temp=1:5);
pf=parmnames[pfit]; pt=parmsi[pf]; pt["fg"]=log(pt["fg"]); pg=parmsi; pg[c("NSA","reps")]=c(1,2); datFiti=cbind(datFiti,pred=LoglikFunLS(parmsTry=pt,parmsG=pg,fitgive=99,datFit=datFiti)[[2]]); 
allScenGreatNOPE2_2025_4_2=allScenRun(c(CritNO3=0.4,CritTemp=200,Ufct=1/2,reps=30,Uyr=2050),PE=2) ##THIS ONE!!


##### Figure 4A,C ####
parmsin=fitinull$parms; pfitn=fitinull$pfit; mdlin=datFitMake2(parmsin,gamGive=TRUE); datFitin=datFitMake2(c(NSA=1,parmsin),temp=1:5);
pfn=c(parmnames[pfitn],"fg"); ptn=parmsi[pfn]; ptn["fg"]=0; pgn=parmsi; pgn[c("NSA","reps")]=c(1,2); 
datFitin=cbind(datFitin,pred=LoglikFunLS(parmsTry=ptn,parmsFit=pfn,parmsG=pgn,fitgive=99)[[2]]); 
allScenGreatNOTP=allScenRun(c(CritNO3=0.4,CritTemp=200,Ufct=1/2,reps=30,Uyr=2040),PE=0,FIT=fitinull,MDL=mdlin,DF=datFitin)
nllS=gam(allScenGreatNOTP[,4]~s(allScenGreatNOTP[,1],k=22))$fitted.values


allScen=allScenS=allScenGreatNOPE2_2025_4_2; KK=8; CS=2; for(i in 2:7) allScenS[,i]=gam(allScen[,i]~s(I(1:nrow(allScen)),k=KK))$fitted.values
trns=function(x,RF=mean(allScenS[1:30,CS*2])) 100*(x-RF)/RF; tail(allScenS[,(CS*2)+(0:1)],1)
plot(allScenS[,1],trns(allScenS[,CS*2]),type="l",lwd=3,col=1,lty=1,las=1,ylim=100*c(-0.75,0.3),yaxt="n"); axis(2,at=seq(25,-75,by=-25),las=2); 
lines(allScenS[,1],trns(allScenS[,CS*2+1]),type="l",lwd=3,col=4,lty=1); lines(allScenGreatNOTP[,1],trns(nllS),col=3,lwd=3); #legend(cex=0.85,"bottomleft",c("RCP 4.5 Best-fit prediction",paste0("Best-fit prediction, ",sprintf('\u2193'),"urchins 50%"),"Base case, no tipping points"),lwd=3,lty=1,col=c(1,4,3),box.col=rgb(0,0,0,alpha=0),bg=rgb(0,0,0,alpha=0)); 
lines(allScenS[1:30,1],trns(allScenS[,CS*2+1])[1:30],type="l",lwd=6,col=6,lty=1);
abline(0,0,col=8,lty=2)





#How often is critical SST threshold exceeded in the past? 2% of sites & yrs
mean(gamdat$no3s<0.4 | gamdat$temp>200)


allScenGreatGLM_CW=allScenRunGLM(CritNO3=0.4,CritTemp=200,Rs=4:6,NC=8)


dfun=function(x,nyi=20,nyf=10){ str=colMeans(head(x,nyi)); fin=colMeans(tail(x,nyf)); (str-fin)/str; }
dglm=dfun(allScenGreatGLM_CW[,-1][,c(FALSE,TRUE,TRUE)])
dode=dfun(allScenGreatNOPE2_2025_4_2[,-1]); tp=-100*c(dglm,dode); 
CL=c("darkorange",palette()[c(6)],"gray25",palette()[c(4)]);
barplot(tp[c(2,1,7,8,NA,4,3,9,10,NA,6,5,11,12,NA)],ylim=c(-120,10),las=1,col=c(CL,NA),yaxt="n"); axis(2,at=c(seq(0,-100,by=-25)),las=2); 
legend(bty="n","bottomleft",fill=CL,cex=0.8,c("Logistic, Baja CA","Logistic, Southern CA","Dynamical, LTM reefs","Dynamical, LTM reefs, 50% urchins")); 






#### VERIFICATION  ####
NSA=2; datFitUse=datFit[!datFit[,"trans"]%in%c(2,6),1:6] #Omitting 2 sites not sampled in some years
datFitMake3=function(parmsi,newdat){
     mdli=datFitMake2(parmsi,gamGive=TRUE); rngi=range(mdli$fitted.values); 
     Epred=predict.gam(mdli,newdata=newdat,type="response");
     newdat$E=(pmax(Epred-rngi[1],0)/(rngi[2]-rngi[1]))^parmsi["EXP"]
     newdat=cbind(r=round(parmsi["r"]*(1+newdat[,"E"]*parmsi["Eff"]),1),newdat); newdat;
}
BDtests3=function(parmsi){
     parmsii=c(r=6,parms["delta"],fg=as.numeric(exp(parms["fg"])),  Eff=0,  parmsi["sigma"],  parms["RecEat"],NSA=NSA,reps=1)
     datSim0=cbind(datFitUse[,c("KI","U")], PE=rlnormT(nrow(datFitUse),1,parmsii["sigma"])-1)
     dget=function(Yi) gamdat[match(codify(cbind(Yi,datSim0[,"U"])),codify(gamdat[,c("yr","lusm")])),  c("no3s","no3sL1","tempL1","tempL2","tempL3","n2","NPGO")]
     datTp1=datFitMake3(parmsi,dget(datFitUse[,"Y"])); datTp2=datFitMake3(parmsi,dget(datFitUse[,"Y"]+1));
     datSim1=cbind(rE=datTp1[,"r"],rEFut=datTp2[,"r"],Y=datFitUse[,"Y"],datSim0); colnames(datTp2)=paste0("F.",colnames(datTp2));
     
     ntp1sim=apply(datSim1,1,function(x){ parmsi["r"]=x["rE"]; lsoda(x["KI"],0:1,odeFull1,parms=parmsi,Vari=c(x["U"],0),PE=x["PE"])[2,-1]; })
     ntp2sim=NA; if(NSA==2) ntp2sim=apply(cbind(datSim1,ntp1sim),1,function(x){ parmsi["r"]=x["rEFut"]; lsoda(x["ntp1sim"],0:1,odeFull1,parms=parmsi,Vari=c(x["U"],0),PE=x["PE"])[2,-1]; });
     cbind(datSim1[,-2], ntp1j=ntp1sim, ntp2j=ntp2sim, rEFut=datSim1[,2], datTp1[,-1],  datTp2[,-1])
}
simdatFitMake2=function(fdat,parmsG,atheta="after2"){
     #!!!Remember time structure:  rE takes KI to ntp1,  rEFut takes ntp1 to ntp2 = K. ==>SO: rE=Elag, rEFut=E
     fdat=as.matrix(fdat); nt=fdat[,"KI"]; ntp1=fdat[,"ntp1j"]; C=parmsG["a"]; rpii = (C+(ntp1+C)/(nt+C)) / (parmsG["b"]-nt);
     gamRes=gam(rpii~s(no3s,k=6)+s(no3sL1,k=6)+s(tempL1,k=6)+s(tempL2,k=10)+s(tempL3,k=8)+s(n2,k=8)+s(NPGO,k=8),data=data.frame(fdat)); rsq(gamRes);
     envtTp1=fdat[,tail(1:ncol(fdat),8)]; colnames(envtTp1)=substr(colnames(envtTp1),3,nchar(colnames(envtTp1)));
     Etp1_0=predict.gam(gamRes,newdata=data.frame(envtTp1),type="response");
     if(atheta=="after"){ Et=nrm(gamRes$fitted.values^parmsG["EXP"]); EFut=nrm(Etp1_0^parmsG["EXP"]); }
     if(atheta=="after2"){ Et=nrm(gamRes$fitted.values)^parmsG["EXP"]; EFut=nrm(Etp1_0)^parmsG["EXP"]; }
     
     NS=nrow(fdat)/length(unique(fdat[,"Y"])); fdat=cbind(fdat,trans=1:NS,resEpred=round(Et,6)); rrr=range(fdat[,"rE"]);
     cr=round(cor(cbind(fdat[,c("rE","resEpred","rEFut")],EFut))[rbind(1:2,3:4)]^2,2); plot(fdat[,c("rE","resEpred")],main=cr[1]); abline(-rrr[1]/diff(rrr),1/diff(rrr),col=2); 
     fdatUE=unique(fdat[,c("rE","resEpred")]); Kp=fdat[,"ntp2j"]; KI=nt;
     datFit=cbind(U=fdat[,"U"], K=Kp,KI=KI,Y=fdat[,"Y"],E=EFut,ELag=fdat[,"resEpred"],fdat[,c("trans","rE","rEFut")]);
     as.matrix(datFit[complete.cases(datFit),]);
}


#Parameter recovery
pfiti=fiti$pfit; pf=parmnames[pfiti]; parmsi=fiti$parm; pt=parmsi[pf]; pt["fg"]=log(pt["fg"]); 
PS=parmStr; PS[pfiti]=colMeans(fiti$B)[pfiti]; parmsi["sigma"]=0; REPORT=matrix(nrow=0,ncol=length(parms)+9);
simDatnR3i=BDtests3(parmsi); DATFITi=simdatFitMake2(simDatnR3i,parmsG,atheta="after2");
fitiSim=fitrun(pfit=pfiti,Ntries=1400,ptryLoc=0,plotReg=FALSE,LB=fiti$B[1,],UB=fiti$B[2,],parmInit=PS,parmsG=parmsG,datFit=DATFITi)



rbind(fitiSim$parms,fiti$parms)
(fitiSim$parms-fiti$parms)/fiti$parms






