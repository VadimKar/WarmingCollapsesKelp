


#Script 4: build future climate projections


library(abind); library(ncf); library(ncdf4); library(mgcv);  library(doParallel); library(foreach);

setwd("/Users/aqeco/Dropbox/Projects/KelpAllee/CodebaseFinal/Data")
nc_e=nc_open(cpath('NEPackelpCanopyEnv_2023.nc'))
edatAnnUP=readRDS("ProcessedData/edatAnnUP.rds")
gamdat=readRDS("ProcessedData/gamdat.rds")


########## GET DATA
#get historic data to calibrate projections
latek=ncvar_get(nc_e,"lat"); lonek=ncvar_get(nc_e,"lon"); siteuse=latek>33.49 & latek<34.51 & lonek>-120.41 & lonek<(-118.71);
yus=5:148; edatF=abind(ncvar_get(nc_e,"temperature"),ncvar_get(nc_e,"nitrate"),along=3)[siteuse,yus,]; eyrs=ncvar_get(nc_e,"year")[yus]
no3.sstfun=function(sst,no3){ mdl=gam(no3~s(sst)); print(rsq(mdl)); function(sst) as.numeric(predict.gam(mdl,newdata=data.frame(sst=sst))); }
sbs=sample(1:nrow(edatF),2e3); no3.sst=no3.sstfun(asvt(edatF[sbs,,1]),asvt(edatF[sbs,,2]))
osst=colMeans(edatF[,,1],na.rm=TRUE); ono3=colMeans(edatF[,,2],na.rm=TRUE);


     

#Get long-term SST projections
wd0=getwd(); setwd(cpath("C:\\Users\\vkara\\Dropbox\\Projects\\KelpAllee\\Data\\climate\\PSL NOAA Climate projections"));
Toread=Sys.glob("*.nc"); spi=regexpr("_",Toread); infos=data.frame(Mdl=substr(Toread,5,spi-1),RCP=as.numeric(substr(Toread,spi+5,spi+6))/10);
LMEwant=3; getfun=function(x) t(apply(ncvar_get(nc_open(x),"SST")[LMEwant,1,1:12,], 2, binmean, 3));
sstA=getfun(Toread[1]); for(i in Toread[-1]) sstA=abind(sstA,getfun(i),along=3); setwd(wd0); 

#aggregate across 24/43 models for which we have projections at RCP 2.6, 4.5, 7, AND 8.5
sus=infos[,"RCP"]%in%c(2.6,4.5,8.5); mdsn=summary(as.factor(infos[sus,"Mdl"])); sum(mdsn==3); mdsUse=sus & infos[,1]%in%names(mdsn[mdsn==3]);
aggfun=function(x) serialAgg(cbind(infos[mdsUse,"RCP"],x),1,FUN=median)[,-1]
sstA0=aperm(sstA[,,mdsUse],3:1); sstA2=aggfun(sstA0[,,1]); for(i in 2:dim(sstA0)[3]) sstA2=abind(sstA2,aggfun(sstA0[,,i]),along=3); 
sstA2=aperm(sstA2,3:1)[tail(1:dim(sstA)[3],121),,] #Earlier we put models in dim 1, years in dim 3. Now reverse back. Dimension 2 is winter, spring, summer, fall temps
# matplot(sstA2[,3,],type="l",lty=1,lwd=3)







###### MAKE FUTURE CLIMATES
normper=function(X,Y,per=4){ x=colMeans(X); y=colMeans(Y);
     di=cbind(x,y); for(ml in 1:per) di=cbind(di,c(y[-(1:ml)],y[1:ml]));
     CR=cor(di)[1,-1]; max(CR); trm=which.max(CR)-1; YN=Y[,(trm+1):(ncol(Y)-(per-trm))];
     seasonDevs=binapply(colMeans(X),ncol(X)/per,clustered=FALSE) - binapply(colMeans(YN),ncol(YN)/per,clustered=FALSE)
     YN+matrix(seasonDevs,nrow(YN),ncol(YN),byrow=TRUE)
}
# 1 - get surrogate TS and put them together for years 2020-2100. surrog2 adapted from wsyn & preserves local diffs in var
surrog2=function(dati,nsurrogs,syncpres,phasekeepPer=NA,pkrng=0){
     wasvect=!is.matrix(dati); if(wasvect) dati=matrix(dati,1,length(dati));
     nloc=nrow(dati); nt=ncol(dati); freq=(0:(nt-1))/nt; fkeep=1:nt;
     if(!is.na(phasekeepPer)) fkeep=freq==(1/phasekeepPer) | freq==(1-(1/phasekeepPer))
     fftdat=t(apply(dati,1,fft)); fftmod=Mod(fftdat); fftarg=Arg(fftdat); 
     res=list(); for(n in 1:nsurrogs){
          h=matrix(rnorm(nt*nloc),nloc,nt)
          if(syncpres) h=h[rep(1,nloc),]; #h[,fkeep]=0; 
          randomizedphases=(fftarg + Arg(t(apply(h,1,fft))))%%(2*pi)
          randomizedphases[,fkeep]=0
          fftsurrog=matrix(complex(modulus=fftmod,argument=randomizedphases),nloc,nt)
          inverse=Re(t(apply(fftsurrog,1,fft,inverse=TRUE))/nt)
          inverse=normper(dati,inverse,phasekeepPer)
          res[[n]]=inverse
     }
     # plot(colMeans(dati[,1:36]),type="l"); lines(colMeans(res[[8]][,1:36]),col=2)
     if(wasvect) for(n in 1:length(res)) res[[n]]=as.vector(res[[n]]); res;
}
matricize=function(di){
     tl=range(di[,2],na.rm=TRUE); ts=sort(unique(di[,2])); if(length(unique(diff(ts)))>1){ print("goofy times!!"); ts=tl[1]+(0:diff(tl)); };
     pops=sort(unique(di[,1])); m=matrix(NA,length(pops),length(ts)); for(i in 1:nrow(di)) m[pops==di[i,1], ts==di[i,2]]=di[i,3]; m[1:nrow(m),]; 
     rownames(m)=pops; colnames(m)=ts; m;
}
rollmean=function(x,ny,binuse=c(F,F,T,F),binl=length(binuse)) filter(x,rep(1/(ny*binl),ny*binl))[binuse]
centf=function(x) t(apply(x,1,nrm2,center=TRUE));
surrogfun=function(climScen=3,yw=2020:2100,seed=1,yrsBase=yw[1]-(10:1)){
     dat=as.matrix(gamdat[gamdat$yr>1985,]); nyn=4*length(yw); #note: omitting yrs 1 and 2 cause lack temp for first 5 quarters
     d0=rbind(dat[,c("trans","yr","tempw")],dat[,c("trans","yr","temps")],dat[,c("trans","yr","temp")],dat[,c("trans","yr","tempf")])
     d0[,"yr"]=d0[,"yr"]+rep((0:3)/4,each=nrow(d0)/4); d0=d0[order(d0[,"yr"]),]; d1=matricize(d0);
     #This next step (1) subtracts local means, (2) plugs in region-wide averages for NAs, and (3) subtracts means again
     locdevs=centf(apply(centf(d1),2,function(x){ x[is.na(x)]=mean(x,na.rm=TRUE); x; }))
     set.seed(seed); eiter=do.call("cbind",surrog2(locdevs,1+ceiling(nyn/ncol(d1)),TRUE,4))[,1:nyn]
     eiter=eiter + matrix(rowMeans(d1,na.rm=TRUE),nrow(eiter),ncol(eiter))
     if(climScen==0) return(eiter);
     
     # ###Some diagnostics on aligning seasons and synchrony
     # plot(binapply(colMeans(eiter[,128+(1:128)]),16,clustered=FALSE),type="o",col=2,ylab="Deviation from global mean, deg C", xlab="quarter (2y total)",main="Data (black) vs surrog() (red)"); points(binapply(colMeans(d1[,1:120]),15,clustered=FALSE),type="o");
     # plot(d1[,c(FALSE,FALSE,TRUE,FALSE)],d1[,c(FALSE,TRUE,FALSE,FALSE)]); points(eiter[,c(FALSE,FALSE,TRUE,FALSE)],eiter[,c(FALSE,TRUE,FALSE,FALSE)],col=2)
     # surrog=round(cor(t(eiter[,c(F,T,F,F)]))[1,],2); obs=round(cor(t(d1[,c(F,T,F,F)]))[1,],2); plot(obs,surrog,ylim=c(0,1),xlim=c(0,1),main="comparing synchrony \n in annualized data"); abline(0,1)
     
     #add smoothed relative annual changes in SST from reference years
     trndGet=function(x) smooth.spline(x,nknots=4)$y - mean(x[(2101-(121:1))%in%yrsBase])
     trndS=tail(as.vector(t(apply(sstA2[,,climScen],2,trndGet))),nyn)
     efin=eiter+matrix(trndS,nrow(eiter),nyn,byrow=TRUE); colnames(efin)=rep(yw,each=4); efin;
}
gamdatFutfun2=function(climScen=3,seed=1){
     sstf=surrogfun(climScen,seed=seed); no3f=matrix(no3.sst(as.vector(sstf)),nrow(sstf),ncol(sstf));
     sstfs=sstf[,c(FALSE,FALSE,TRUE,FALSE)]; nt=ncol(sstfs); ns=nrow(sstfs);
     no3s=no3f[,c(FALSE,TRUE,FALSE,FALSE)]; no3=no3f[,c(FALSE,FALSE,TRUE,FALSE)]; no3w=no3f[,c(TRUE,FALSE,FALSE,FALSE)];
     no3L1=cbind(NA,no3[,-nt]); no3L2=cbind(NA,no3L1[,-nt]);
     sstL1=cbind(NA,sstfs[,-nt]); sstL2=cbind(NA,sstL1[,-nt]); sstL3=cbind(NA,sstL2[,-nt]); sstL4=cbind(NA,sstL3[,-nt]);
     
     YPs=2020:2100; gamdatFut=expand.grid(yr=YPs,trans=sort(unique(gamdat$trans)));
     gamdatFut=cbind(gamdatFut,K=NA,NPGO=NA,temp=asvt(sstfs),no3=asvt(no3),no3L1=asvt(no3L1),no3L2=asvt(no3L2),no3s=asvt(no3s),no3w=asvt(no3w),no3sL1=asvt(cbind(NA,no3s[,-nt])),
                     tempL1=asvt(sstL1),tempL2=asvt(sstL2),tempL3=asvt(sstL3),tempL4=asvt(sstL4),temps=asvt(sstf[,c(FALSE,TRUE,FALSE,FALSE)]),
                     n2=rep(rollmean(colMeans(no3f),2),ns),n5=rep(rollmean(colMeans(no3f),5),ns),
                     U=gamdat$lusm[match(gamdatFut$trans,gamdat$trans)],PE=0)
     gamdatFut$n2[gamdatFut$yr==2020]=gamdat$n2[gamdat$yr==2019][1]; gamdatFut$n5[gamdatFut$yr<2022]=gamdat$n5[gamdat$yr==2019][1];
     for(i in unique(gamdatFut$trans)){
          gsel=function(y) gamdat$trans==i & gamdat$yr==y; gfsel=function(y) gamdatFut$trans==i & gamdatFut$yr==y;
          gamdatFut$K[gfsel(2020)]=gamdat$res[gsel(2018)]; gamdatFut$no3sL1[gfsel(2020)]=gamdat$no3s[gsel(2019)];
          gamdatFut$no3L1[gfsel(2020)]=gamdatFut$no3L2[gfsel(2021)]=gamdat$no3[gsel(2019)]; gamdatFut$no3L2[gfsel(2020)]=gamdat$no3[gsel(2018)];
          gamdatFut$tempL1[gfsel(2020)]=gamdatFut$tempL2[gfsel(2021)]=gamdatFut$tempL3[gfsel(2022)]=gamdatFut$tempL4[gfsel(2023)]=gamdat$temp[gsel(2019)]
          gamdatFut$tempL2[gfsel(2020)]=gamdatFut$tempL3[gfsel(2021)]=gamdatFut$tempL4[gfsel(2022)]=gamdat$temp[gsel(2018)]
          gamdatFut$tempL3[gfsel(2020)]=gamdatFut$tempL4[gfsel(2021)]=gamdat$temp[gsel(2017)]
          gamdatFut$tempL4[gfsel(2020)]=gamdat$temp[gsel(2016)]
     }; gamdatFut[gamdatFut$yr<=2097,];
}
gamdatFutfunCW=function(climScen=3,seed=1){ 
     yw=2016:2100; nt=length(yw); yrsBase=yw[1]-(10:1);
     dat=edatAnnUP[,-(1:2),c("tempw","temps","temp","tempf")]; nyn=4*nt; ns=nrow(dat);
     #note: omitting yrs 1 and 2 cause lack temp for first 5 quarters
     nyd=ncol(edatAnnUP)-2; dy=rep(1985+(1:nyd),4)+rep((0:3)/4,each=nyd);
     dim(dat)=c(ns,nyd*4); dat=dat[,order(dy)]; locdevs=centf(dat);
     set.seed(seed); eiter=do.call("cbind",surrog2(locdevs,1+ceiling(nyn/ncol(dat)),TRUE,4))[,1:nyn]
     eiter=eiter + matrix(rowMeans(dat,na.rm=TRUE),ns,nyn)
     # ###Some diagnostics on aligning seasons and synchrony
     # plot(binapply(colMeans(eiter[,(16*8*1.25)+(1:128)]),16,clustered=FALSE),type="o",col=2,ylab="Deviation from global mean, deg C", xlab="quarter (2y total)",main="Data (black) vs surrog() (red)"); points(binapply(colMeans(dat[,1:120]),15,clustered=FALSE),type="o");
     # plot(dat[,c(FALSE,FALSE,TRUE,FALSE)],dat[,c(FALSE,TRUE,FALSE,FALSE)]); points(eiter[,c(FALSE,FALSE,TRUE,FALSE)],eiter[,c(FALSE,TRUE,FALSE,FALSE)],col=2)
     # ss=sample(1:ns,300); surrog=round(cor(t(eiter[ss,c(F,T,F,F)]))[1,],2); obs=round(cor(t(dat[ss,c(F,T,F,F)]))[1,],2); 
     # plot(obs,surrog,ylim=0:1,xlim=0:1,main="comparing synchrony of annual temp \n in data (X) vs surrogates (Y)"); abline(0,1)
     
     #add smoothed relative annual changes in SST from reference years
     trndGet=function(x) smooth.spline(x,nknots=4)$y - mean(x[(2101-(121:1))%in%yrsBase])
     trndS=0; if(climScen>0) trndS=tail(as.vector(t(apply(sstA2[,,climScen],2,trndGet))),nyn);
     sstf=eiter + matrix(trndS,ns,nyn,byrow=TRUE)
     
     no3fv=round(as.vector(sstf),3); no3fu=unique(no3fv); no3f=matrix(no3.sst(no3fu)[match(no3fv,no3fu)],ns,nt*4);
     no3f2=abind(lapply(split.data.frame(t(no3f),1:4),t),along=3); no3M=no3f2[,,2]; no3M[no3M<no3f2[,,1]]=no3f2[,,1][no3M<no3f2[,,1]]; #This effectively takes the max of winter and spring nitrate
     sstfs=sstf[,c(FALSE,FALSE,TRUE,FALSE)]; sstfspr=sstf[,c(FALSE,TRUE,FALSE,FALSE)]; sstL1=cbind(NA,sstfs[,-nt]); sstL2=cbind(NA,sstL1[,-nt]); sstL3=cbind(NA,sstL2[,-nt]); sstL4=cbind(NA,sstL3[,-nt]);
     no3sL1=cbind(NA,no3f2[,-nt,2]); no3=no3f2[,,3]; no3L1=cbind(NA,no3[,-nt]); no3L2=cbind(NA,no3L1[,-nt]);
     
     gdf=data.frame(yr=yw,trans=rep(1:nrow(dat),each=nt),K=NA,temp=asvt(sstfs),temps=asvt(sstfspr),tempL1=asvt(sstL1),tempL2=asvt(sstL2),tempL3=asvt(sstL3),tempL4=asvt(sstL4),
                    no3s=asvt(no3f2[,,2]),no3sL1=asvt(no3sL1),no3w=asvt(no3f2[,,1]),no3mx=asvt(no3M),no3=asvt(no3),no3L1=asvt(no3L1),no3L2=asvt(no3L2))
     gdf$K[gdf$yr==2020]=1-edatAnnUP[,35,"UPbarren"]; final1=gdf[gdf$yr%in%2020:2097,];
     final2=data.frame(matrix(as.integer(1000*as.matrix(final1)),nrow(final1),ncol(final1))); names(final2)=names(gdf); final2;
}
gamdatFutMake=function(climScen,reps=1,FUN=gamdatFutfun2,NC=8){ if(NC>1) registerDoParallel(NC)
     out=foreach(i=1:reps,.combine=c,.packages=c("abind","MASS"),.export=ls()) %dopar% list(FUN(climScen,i))
     if(NC>1) stopImplicitCluster(); out; }





NR=30; gamdatFut1=gamdatFutMake(1,NR); gamdatFut2=gamdatFutMake(2,NR); gamdatFut3=gamdatFutMake(3,NR);
saveRDS(list(gamdatFut1,gamdatFut2,gamdatFut3),"gamdatFutAll.rds")

NR=30; gamdatFutCW0=gamdatFutMake(0,NR,gamdatFutfunCW); gamdatFutCW1=gamdatFutMake(1,NR,gamdatFutfunCW); gamdatFutCW2=gamdatFutMake(2,NR,gamdatFutfunCW); gamdatFutCW3=gamdatFutMake(3,NR,gamdatFutfunCW);
saveRDS(list(gamdatFutCW1,gamdatFutCW2,gamdatFutCW3,gamdatFutCW0),"gamdatFutCWAll.rds")















