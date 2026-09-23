


#Script 3: Assemble region-wide remote sensing data and run analyses in Figure 1


library(ncdf4); library(abind); library(raster)


setwd("/Users/aqeco/Dropbox/Projects/KelpAllee/CodebaseFinal/Data")
nc_e=nc_open(cpath('NEPackelpCanopyEnv_2023.nc'))
nc_data=nc_open(cpath('kelpCanopyFromLandsat_2021.nc'))
gamdat=readRDS("ProcessedData/gamdat.rds")

#In gamdat file
#nsatpix - # pixels at transect (constant over time)
#nPixForL1 - # pixels present last year
#nPixBarL1 - # pixels absent last year
#nPixColL1 = # pixels present last year and absent this year
#nPixRecL1 = # pixels absent last year and present this year
#restm1 = proportion of pixels forested on the site (resilience) last year (t-1)
#no3 = min quarterly nitrate between t-1 and t. Will roll up into temp later.
#no3L1 = min quarterly nitrate between t-2 and t-1. Will roll up into temp later.
#temp = max quarterly temp between t-1 and t
#tempLX = max quarterly temp between t-1-X and t-X
#lusm = log urchin density at the site, averaged over time
#hs = wave stress between last year and now





scaleup=function(lonlat,SZ=4,distThresh=11e-4,plotgive=FALSE){ 
     cutsFine=c(60,47,42,38,36.2,35.5,33.6,34.4,32.4,31,28.17)
     RegH=rowSums(outer(lonlat[,2],cutsFine,"<"))
     dst=cbind(lonlat[,1:2],NA,NA,NA,RegH,NA,NA,NA)
     for(R in unique(RegH)){
          gps=lonlat[RegH==R,]; ds=c(0.0003160301,0.0002679593);
          lonu=seq(min(gps[,1]),max(gps[,1]),by=ds[1]); dlon=ds[1]/2;
          latu=seq(min(gps[,2]),max(gps[,2]),by=ds[2]); dlat=ds[2]/2;
          mlon=findInterval(gps[,1],c(lonu[1]-dlon,lonu+dlon))-1
          mlat=findInterval(gps[,2],c(latu[1]-dlat,latu+dlat))-1
          mlonR=match(mlon,sort(unique(mlon))); mlatR=match(mlat,sort(unique(mlat)))
          d=abs(gps[,1]-lonu[pmax(mlon,1)])+abs(gps[,2]-latu[pmax(mlat,1)])
          dst[RegH==R,3:5]=cbind(d,mlonR,mlatR)
     }
     print(mean(dst[dst[,3]<10e-4,3]*sqrt(2)/2)/0.0002703456)
     if(plotgive) hist(dst[dst[,3]<10e-4,3]*100*1000*sqrt(2)/2,xlab="meters from nearest pixel ID") 
     for(i in unique(dst[,6])){
          dsti=dst[dst[,6]==i,]; loncelr=range(dsti[,4]); latcelr=range(dsti[,5]);
          upcel=cbind(findInterval(dsti[,4],seq(loncelr[1],loncelr[2],by=SZ)),findInterval(dsti[,5],seq(latcelr[1],latcelr[2],by=SZ)))
          dst[dst[,6]==i,7:8]=upcel+i/10
     }
     miss=dst[,3]>distThresh; dst[miss,7:8]=NA; dst;
}
rowmeans=function(x,grp=rep(1,nrow(x))){ 
     n=rowsum(1+0*x,grp,reorder=FALSE,na.rm=TRUE); 
     avg=rowsum(x,grp,reorder=FALSE,na.rm=TRUE)/n; avg[n==0]=NA; 
     if(nrow(avg)==1) avg=as.vector(avg); avg; 
}
#Spatial analyses
clumpan=function(valmat,prop=FALSE){
     if(min(is.na(valmat))==1) return(NA); if(max(valmat,na.rm=T)==0) return(valmat*0)
     cr=range(which(colMeans(is.na(valmat))<1)); rr=range(which(rowMeans(is.na(valmat))<1));
     if(cr[1]!=cr[2] & rr[1]!=rr[2]) valmat=valmat[rr[1]:rr[2],cr[1]:cr[2]]
     cr=clump(raster(valmat)); cf=na.omit(freq(cr));
     cf[match(as.vector(as.matrix(cr)),cf[,1]),2] / c(1,sum(!is.na(valmat)))[1+prop]
}
clumpEval=function(vals,cordsYX,binlen=10,prop=FALSE){
     inreg=cordsYX[,1]>0 & cordsYX[,2]>0; vals=vals[inreg]; valsOut=0*vals-1; cordsYX=cordsYX[inreg,]; 
     m=ref=matrix(NA,max(cordsYX[,1]),max(cordsYX[,2])); m[cordsYX]=vals; ref[cordsYX]=1:length(vals);
     binr=ceiling((1:nrow(m))/(binlen/0.03)); lb=length(binr)+(0:-1); binr[lb[1]]=binr[lb[2]];	
     for(bi in unique(binr)){
          mi=m[binr==bi,]; mir=na.omit(as.vector(mi)); 
          refi=na.omit(as.vector(ref[binr==bi,]))
          valsOut[refi][mir==1]=na.omit(clumpan(mi,prop))
          valsOut[refi][mir==0]=na.omit(clumpan(!mi,prop))
     }
     out=NA*inreg; out[inreg]=valsOut; out;
}
BarSpatCalc2=function(datBarTemp,dst,minlen=0,maxlen=100,windsize=10,prop=FALSE){
     datBarSpat=NA*datBarTemp; for(Yr in 1:ncol(datBarTemp)) for(R in unique(dst[,6])){ sts=(dst[,6]==R & log(dst[,3])<(-6)); print(c(Yr,R));
     datBarSpat[sts,Yr]=clumpEval((datBarTemp[sts,Yr]>=minlen & datBarTemp[sts,Yr]<=minlen),dst[sts,5:4],windsize,prop); }; datBarSpat;
}
durtrim=function(x,minl=2,zrokeep=FALSE){ x[x<minl & (!zrokeep | x!=0)]=NA; x; }
spatrel=function(Y,YP=!is.na(rowSums(Y)),lat=dat[,2],span="cv",bass=0){ 
     smoo=function(yi) supsmu(order(lat[YP],decreasing=TRUE),yi,span=span,bass=bass)$y; 
     out=NA*Y[,1]; out[YP]=smoo(Y[YP,1])/smoo(Y[YP,2]); out; 
}
spatplot=function(Y,ynm=NULL,X=NULL,Y2=NULL,Y3=NULL,Y4=NULL,cls=c(1,3,2),lws=rep(5,4),labRel=0.9,alpha=0.02,us=AlleeUse,lat=dat0[us,2],regi=Reg0[us],bass=0,...){
     X0=order(lat,decreasing=TRUE); if(is.null(X)) X=X0; xnm="Patch #,   ordered North to South"; 
     smool=function(YI,X,cl=1,lw=5){ sm=supsmu(X,YI,bass=bass); lines(sm$x,sm$y,col=cl,lwd=lw); }
     plot(X,Y,xlab=xnm,ylab=ynm,col=rgb(0,0,1,alpha=alpha),pch=16,...); abline(v=cumsum(summary(as.factor(regi)))[c(1,3,4,5)],lwd=5,col=8,lty=1); 
     text(c("WA","Central CA","Southern CA","Baja CA"),x=c(1e3,3e4,6.5e4,9.5e4),y=max(Y,na.rm=TRUE)*rep(labRel,4))
     YP=list(Y,Y2,Y3,Y4); for(i in which(unlist(lapply(YP,length))>1)) smool(YP[[i]],X,cls[i],lws[i])
}





#Full region
pTH=0.2; cuts=c(60,47,42,38,34.65,32.4);
tsel=1:160; edatget=function(vnm,mod=1) cbind(ncvar_get(nc_e,"lon"), ncvar_get(nc_e,"lat"), mod*ncvar_get(nc_e,vnm)[,tsel])
edat0=abind(edatget("area",1/900),edatget("temperature"),edatget("nitrate"),along=3)
eyrs=(ncvar_get(nc_e,"year")+(ncvar_get(nc_e,"quarter")-1)/4)[tsel]; colnames(edat0)=c("lon","lat",eyrs); remove(nc_e); dimnames(edat0)[[3]]=c("area","temp","no3");
presFreq=rowMeans(edat0[,3:150,1]>0,na.rm=TRUE); edatAlleeUse=rowMeans(is.na(edat0[,3:150,1]))<0.3378 & presFreq>pTH;
eReg0=rowSums(outer(edat0[,"lat",1],cuts,"<")); eReg=eReg0[edatAlleeUse];
eupscaleIDs=scaleup(edat0[edatAlleeUse,1:2,1],SZ=4) #Match each 30m pixel to an id of 4*30=120m pixels


edatAnnNms=c("info","area","forT","barT","barren","barL1","recover","die","temp","temps","no3","no3s","no3f","no3w","tempf","tempw")
edatAnn=array(NA, c(sum(edatAlleeUse), (ncol(edat0)-2)/4, length(edatAnnNms))); dimnames(edatAnn)[[3]]=edatAnnNms;
colnm=c("Reg","presFreq","lon","lat","PxDstGPS","TiffLonPx","TiffLatPx","OutsideTiffGrpID"); coldiff=ncol(edatAnn)-length(colnm);
edatAnn[,,"info"]=cbind(eReg,presFreq[edatAlleeUse],eupscaleIDs[,-(4:6)],matrix(NA,nrow(edatAnn),coldiff)); 
#First line below calculates edatAnnSall, used in model fitting. Second line calculates edatAnnThall, used in Fig 1.
# edatAnn[,,"area"]=t(apply(edat0[edatAlleeUse,-(1:2),1][,c(FALSE,TRUE,TRUE,FALSE)],1,binmean,2)); bar=edatAnn[,,"area"]==0;
edatAnn[,,"area"]=t(apply(edat0[edatAlleeUse,-(1:2),1][,c(F,T,T,T)],1,binmean,3)); bar=edatAnn[,,"area"]<0.04;
edatAnn[,,"forT"]=t(apply(!bar,1,function(x){ r=rle(x); rep(r$lengths*r$values,r$lengths); }))
edatAnn[,,"barT"]=t(apply(bar,1,function(x){ r=rle(x); rep(r$lengths*r$values,r$lengths); }))
edatAnn[,,"barren"]=Bs=edatAnn[,,"barT"]>1; edatAnn[,,"barL1"]=cbind(NA,edatAnn[,-ncol(edatAnn),"barren"]);
erecover=t(apply(Bs==1,1,diff)==-1); erecover[Bs[,-ncol(edatAnn)]==0]=NA; edatAnn[,,"recover"]=cbind(erecover,NA);
edie=t(apply(Bs==0,1,diff)==-1); edie[Bs[,-ncol(edatAnn)]==1]=NA; edatAnn[,,"die"]=cbind(edie,NA);
edatAnn[,,c("tempw","temps","temp","tempf")]=abind(lapply(split.data.frame(t(edat0[edatAlleeUse,-(1:2),"temp"]),1:4),t),along=3)
edatAnn[,,c("no3w","no3s","no3","no3f")]=abind(lapply(split.data.frame(t(edat0[edatAlleeUse,-(1:2),"no3"]),1:4),t),along=3)
#Add lags of spring and summer temp and nitrate    #edatAnn[1:6,3:8,c("barren","die","barL1")]
matExpi=function(x,ncols=ncol(edatAnn)) abind(array(NA,c(nrow(x),ncols-ncol(x),dim(x)[3])),x,along=2)
for(L in 1:4) edatAnn=abind(edatAnn,matExpi(edatAnn[,head(1:ncol(edatAnn),-L),c("temp","temps","no3","no3s")]),along=3)
dimnames(edatAnn)[[3]][tail(1:dim(edatAnn)[3],16)]=paste0(c("temp","temps","no3","no3s"),"L",rep(1:4,each=4))
dimnames(edatAnn)[[2]]=c(colnm,rep("NA",coldiff))



dat0=cbind(ncvar_get(nc_data,"lon"), ncvar_get(nc_data,"lat"), cover=ncvar_get(nc_data,"area")/900); 
yrs=ncvar_get(nc_data,"year")+(ncvar_get(nc_data,"quarter")-1)/4; colnames(dat0)=c("lon","lat",yrs); remove(nc_data);
presFreq=rowMeans(dat0[,-(1:2)]>0,na.rm=TRUE); datAlleeUse=rowMeans(is.na(dat0[,-(1:2)]))<0.3378 & presFreq>pTH;
Reg0=rowSums(outer(dat0[,"lat"],cuts,"<")); Reg=Reg0[datAlleeUse];
upscaleIDs=scaleup(dat0[datAlleeUse,1:2],SZ=8) #Match each 30m pixel to an id of 4*30=120m pixels

datAnnNms=c("info","area","forT","barT","forE","barE","barren")
datAnn=array(NA, c(sum(datAlleeUse), (ncol(dat0)-2)/4, length(datAnnNms))); dimnames(datAnn)[[3]]=datAnnNms;
colnm=c("Reg","presFreq","lon","lat","PxDstGPS","TiffLonPx","TiffLatPx","OutsideTiffGrpID"); coldiff=ncol(datAnn)-length(colnm);
datAnn[,,"info"]=cbind(Reg,presFreq[datAlleeUse],upscaleIDs[,-(4:6)],matrix(NA,nrow(datAnn),coldiff)); 
#First line below calculates edatAnnSall, used in model fitting. Second line calculates edatAnnThall
# datAnn[,,"area"]=t(apply(dat0[datAlleeUse,-(1:2),1][,c(FALSE,TRUE,TRUE,FALSE)],1,binmean,2)); bar=datAnn[,,"area"]==0;
datAnn[,,"area"]=t(apply(dat0[datAlleeUse,-(1:2)][,c(F,T,T,T)],1,binmean,3)); bar=datAnn[,,"area"]<0.04;
datAnn[,,"forT"]=t(apply(!bar,1,function(x){ r=rle(x); rep(r$lengths*r$values,r$lengths); }))
datAnn[,,"barT"]=t(apply(bar,1,function(x){ r=rle(x); rep(r$lengths*r$values,r$lengths); }))
datAnn[,,"barren"]=datAnn[,,"barT"]>1; dimnames(datAnn)[[2]]=c(colnm,rep("NA",coldiff));
datBarSpat2_5km_propCalc=BarSpatCalc2(datAnn[,,"barT"],upscaleIDs,2,windsize=5,prop=TRUE)
datBarSpat=datForSpat=datBarSpat2_5km_propCalc; datForSpat[datAnn[,,"barT"]>0]=NA; datBarSpat[datAnn[,,"barT"]<2]=NA;
datAnn[,,"forE"]=datForSpat; datAnn[,,"barE"]=datBarSpat;




fill1NA=function(x){
     r=rle(is.na(x)); nal=rep(r$lengths*r$values,r$lengths);
     for(i in which(nal==1)) x[i]=sample(x[i+c(-1,1)],1)
     x
}
upscale=function(edatAnn,areaUPcalc=FALSE,SZcrit=4,stateThresh=0.4){
     regPres=!is.na(edatAnn[,"TiffLonPx","info"]); grpPres=!is.na(edatAnn[,"OutsideTiffGrpID","info"]);
     idUp=codify(edatAnn[,c("TiffLonPx","TiffLatPx"),"info"]); idUp[!regPres & grpPres]=edatAnn[!regPres & grpPres,"OutsideTiffGrpID","info"];
     us=regPres | grpPres; edatAnnCls=edatAnn[us,,]; idUp=idUp[us];
     nsamp=edatAnnUP=rowsum(matrix(1,nrow(edatAnnCls),ncol(edatAnnCls)),idUp,reorder=FALSE)
     for(i in 1:dim(edatAnnCls)[3]) edatAnnUP=abind(edatAnnUP,rowmeans(edatAnnCls[,,i],idUp),along=3)
     usUP=nsamp[,1]>SZcrit; edatAnnUP=edatAnnUP[usUP,,-1]; nsamp=nsamp[usUP,1]; remove(edatAnnCls)
     dimnames(edatAnnUP)[[3]]=dimnames(edatAnn)[[3]];
     edatAnnUP[,1,"info"]=round(edatAnnUP[,1,"info"]); Bs=(1-edatAnnUP[,,"barren"])<stateThresh;
     Bs=t(apply(Bs,1,fill1NA))
     
     euprecover=t(apply(Bs==1,1,diff)==-1); euprecover[Bs[,-ncol(edatAnnUP)]==0]=NA;
     eupdie=t(apply(Bs==0,1,diff)==-1); eupdie[Bs[,-ncol(edatAnnUP)]==1]=NA;
     edatAnnUP=abind(edatAnnUP,UPbarren=Bs,UPrecover=cbind(euprecover,NA),UPdie=cbind(eupdie,NA),
                     UPforT=t(apply(!Bs,1,function(x){ r=rle(x); rep(r$lengths*r$values,r$lengths); })),
                     UPbarT=t(apply(Bs,1,function(x){ r=rle(x); rep(r$lengths*r$values,r$lengths); })),along=3)
     if(areaUPcalc){
          datBarSpat2_5km_propCalcUP=BarSpatCalc2(edatAnnUP[,,"UPbarT"],edatAnnUP[,c(3:7,1),"info"],2,windsize=10,prop=TRUE)
          datBarSpat=datForSpat=datBarSpat2_5km_propCalcUP; datForSpat[edatAnnUP[,,"UPbarT"]>0]=NA; datBarSpat[edatAnnUP[,,"UPforT"]>0]=NA;
          edatAnnUP=abind(edatAnnUP,UPforE=datForSpat,UPbarE=datBarSpat,along=3)
     }; dimnames(edatAnnUP)[[2]]=dimnames(edatAnn)[[2]]; edatAnnUP;
}









edatAnnUP=upscale(edatAnn)
datUP=upscale(datAnn,TRUE)
saveRDS(edatAnnUP,"ProcessedData/edatAnnUP.rds")





#Fig. 1A
library(ggOceanMaps); library(sf); library(rnaturalearth)
lims=c(xmin=-131,xmax=-110.4,ymin=22,ymax=51); bbox=st_bbox(lims,crs=4326) |> st_as_sfc();
states=ne_states(country=c("united states of america", "mexico", "canada"), returnclass="sf") |> st_crop(bbox)
coastline=ne_coastline(scale="medium", returnclass="sf") |> st_crop(bbox)
gps=data.frame(unique(round(datAnn[,c("lon","lat","Reg"),"info"],1)),cl=NA); gps=gps[!(gps[,2]<46 & gps[,2]>44) & !(gps[,2]<38 & gps[,2]>37.2) & !(gps[,2]<42 & gps[,2]>41),];
gps[gps[,3]%in%(1:2),4]=cols[1]; gps[gps[,3]%in%(3:4),4]=cols[2]; gps[gps[,3]==5,4]=cols[3]; gps[gps[,3]==6,4]=cols[4]; 
basemap(limits=lims, bathy.style="rcb", grid.col=NA, land.col="ivory2") +
     geom_sf(data=states,fill=NA,color="gray70",linewidth=0.2,inherit.aes=FALSE) +
     geom_point(data=gps, aes(x=lon,y=lat,color=cl), cex=2.5) + scale_color_identity() +
     geom_sf(data=coastline,fill=NA,color="black",linewidth=0.35,inherit.aes=FALSE) + 
     coord_sf(xlim=lims[1:2],ylim=lims[3:4],expand=FALSE) + theme_void() + theme(legend.position="none")




staTP=data.frame(RG=rep(datUP[,"Reg","info"],each=ncol(datUP)), pB=asvt(datUP[,,"barren"]), dF=asvt(durtrim(datUP[,,"UPforT"],1)),
                 dB=asvt(durtrim(datUP[,,"UPbarT"],1)), aB=asvt(datUP[,,"UPbarE"]));
#lump small sub-regions with central CA for clarity
staTP[,1][staTP[,1]==2]=1; staTP[,1][staTP[,1]==3]=4;



#Compare regional stats
samd=serialAgg(cbind(S=as.numeric(staTP[,1]),staTP),1,FUN=function(x) round(median(x,na.rm=TRUE),3)); samd; samd[2,]/samd[1,]; summary(as.factor(staTP[,1]))/37
samd=serialAgg(cbind(S=as.numeric(staTP[,1]>4),staTP),1,FUN=function(x) round(median(x,na.rm=TRUE),3)); samd; samd[2,]/samd[1,]; summary(as.factor(staTP[,1]))/37


library(ggplot2)
cols=c("springgreen3","#11C8F0",palette()[c(6)],"darkorange"); cl=staTP[, 1]; cl[cl==1]=cols[1]; cl[cl==4]=cols[2]; cl[cl==5]=cols[3]; cl[cl==6]=cols[4]; 
df <- data.frame(xi = 1-staTP[,"pB"], grp = factor(staTP[, 1]), col = cl)
ggplot(df, aes(xi, color = col, fill = col)) + geom_density(bw = 0.055) + scale_color_identity() + scale_fill_identity() +
     facet_wrap(~grp, ncol = 1) + theme_void() + theme(panel.grid = element_blank(), legend.position = "none",strip.text = element_blank())












#Fig. 1B
mplot=function(m,cols=rep("springgreen3",nrow(m))){
     nr=nrow(m); nc=ncol(m)
     plot(0, type="n", xlim=c(0.5,nc+0.5), ylim=c(0.5,nr+0.5),xlab="Column", ylab="Row")
     abline(h=seq(0.5,nr+0.5), col="gray80", lwd=0.5)
     idx <- which(m == 1, arr.ind=TRUE)
     rect(idx[,2]-.5, nr-idx[,1]+.5, idx[,2]+.5, nr-idx[,1]+1.5, col=cols[idx[,1]], border=NA)
}
mplotadd=function(m) {
     nr <- nrow(m); nc <- ncol(m); idx <- which(m == 1, arr.ind=TRUE)
     rect(idx[,2]-.5, nr-idx[,1]+.5, idx[,2]+.5, nr-idx[,1]+1.5, col="gray", border=NA)
}
tsiUP=(1-edatAnnUP[4151,,"barren"])>0.4; edatAnnUP[4151,c("lat","lon"),"info"]

regs=datUP[,"Reg","info"]; regs[regs%in%(1:2)]=1; regs[regs%in%(3:4)]=4;
set.seed(1); n=80; nregs=round(n*summary(as.factor(regs))/nrow(datUP)); 
sels=unlist(sapply(nregs,function(ni) sample(which(regs==c(1,4:6)[nregs==ni]),ni)))
selr=regs[sels]; sels=sels[order(selr)]; selc=selr[order(selr)]; selc[selc==1]=cols[1]; selc[selc==4]=cols[2]; selc[selc==5]=cols[3]; selc[selc==6]=cols[4]; 
m=round(0.38*n); tp=1-datUP[sels,,"UPbarren"]; tp[m,]=tsiUP[1:37]; mplot(tp,cols=selc); abline(h=n-m+c(-0.75,0.75),lwd=3,col=2); mplotadd(is.na(tp))



# par(ask=TRUE); for(i in sample(which(edatAnnUP[,"Reg","info"]>4),1e2)) 
i=4151; plot(unique(floor(eyrs)),1-edatAnnUP[i,,"barren"],type="o",pch=16,main=i,las=1,lwd=2); par(ask=FALSE); #4140; 4151









#Fig. 1C
datget=function(reg,nsamp=4e3){
     if(reg==0){ dati=gamdat[gamdat$nsatpix>4 & !is.na(gamdat$tempL4),]; dati$recover=dati$rec; return(dati); }
     if(reg[1]>0){ pixids=sample(which(edatAnn[,"Reg","info"]%in%reg),nsamp); dati=edatAnn[pixids,,]; }
     dati=abind(dati,matrix(list(NA,transids)[[1+(reg==0)]],dim(dati)[1],dim(dati)[2]),along=3)
     if(filt | reg[1]>0){ us=rowMeans(dati[,1:37,"area"]>0,na.rm=TRUE)>pTH; dati=dati[us,,]; pixids=pixids[us]; }
     NY=ncol(dati); nsamp=nrow(dati); dim(dati)=c(nsamp*NY,dim(dati)[3]);
     dati=data.frame(dati); names(dati)=c(dimnames(edatAnn)[[3]],"trans");
     dati$yr=rep(unique(floor(eyrs)),each=nsamp); dati$pixid=rep(pixids,NY);
     dati$T5=rowMeans(dati[,c("temp","tempL1","tempL2","tempL3","tempL4")])
     dati[!is.na(dati$T5),]
}
datgetUP=function(reg,stateThresh=0.4){
     if(reg==0){
          datiUP=gamdat[gamdat$nsatpix>4 & !is.na(gamdat$tempL4),]; icB=datiUP$restm1<stateThresh; B=datiUP$res<stateThresh;
          datiUP$recover=icB & !B; datiUP$recover[!icB]=NA;
          datiUP$die=!icB & B; datiUP$die[icB]=NA; return(datiUP)
     }
     dati=edatAnnUP[edatAnnUP[,"Reg","info"]==reg,,-1]; pixids=1:nrow(dati); 
     NY=ncol(dati); nsamp=nrow(dati); dim(dati)=c(nsamp*NY,dim(dati)[3]);
     dati=data.frame(dati); names(dati)=dimnames(edatAnnUP)[[3]][-1];
     dati$yr=rep(unique(floor(eyrs)),each=nsamp); dati$pixid=rep(pixids,NY);
     dati$T5=rowMeans(dati[,c("temp","tempL1","tempL2","tempL3","tempL4")])
     dati$recover=dati$UPrecover; dati$die=dati$UPdie;
     dati#[!is.na(dati$T5),]
}
tplotfun=function(reg,Vdie="tempL1",Vrec="tempL2",dur=TRUE,cent=TRUE,Qcut=0.99,Cbar=2,Cfor=1,LW=2,plot=0,UP=TRUE,ylbL="Mean state duration",warmProj=NA,...){
     if(UP) dati=datgetUP(reg) else dati=datget(reg);
     glmDieA=glm(formula(paste0("die~",Vdie)),binomial,dati); glmRecA=glm(formula(paste0("recover~",Vrec)),binomial,dati)
     if(plot==0 & is.na(warmProj)) return(c(rsq(glmDieA),rsq(glmRecA)))
     fn=(1:2)[1+(LW>0)]; ptt=c("o","l")[1+(LW>0)]
     xdie=na.omit(dati[,c("die",Vdie)])[,Vdie]; diep=aggregate(glmDieA$fitted.values,by=list(round(xdie*fn)/fn),mean);
     xrec=na.omit(dati[,c("recover",Vrec)])[,Vrec]; recp=aggregate(glmRecA$fitted.values,by=list(round(xrec*fn)/fn),mean);
     xdier=quantile(xdie,c(1-Qcut,Qcut),na.rm=TRUE); xrecr=quantile(xrec,c(1-Qcut,Qcut),na.rm=TRUE);
     diep[diep[,1]<xdier[1] | diep[,1]>xdier[2],]=NA; recp[recp[,1]<xrecr[1] | recp[,1]>xrecr[2],]=NA;
     if(cent){ diep[,1]=diep[,1]-mean(xdier); recp[,1]=recp[,1]-mean(xrecr); }
     if(!is.na(warmProj)) return(c(diep[c(which.min(abs(diep[,1])),which.min(abs(diep[,1]-warmProj))),2],
          recp[c(which.min(abs(recp[,1])),which.min(abs(recp[,1]-warmProj))),2]))
     if(dur){ diep[,2]=1/diep[,2]; recp[,2]=1/recp[,2]; }
     if(plot==1) plot(diep,lwd=LW,type="l",col=Cfor,las=1,ylab=ylbL,...)
     points(diep,col=Cfor,lwd=LW,type="l",...); points(recp,col=Cbar,lwd=LW,type="l",...);
     points(diep[c(FALSE,TRUE),],col=Cfor,lwd=LW,type="p",...); points(recp[c(FALSE,TRUE),],col=Cbar,lwd=LW,type="p",...);
}
par(mfrow=c(1,1),mar=c(4,4,4,5)) #c(bottom, left, top, right)
tplotfun(5,"tempL1","tempL2",plot=1,pch=16,xlim=c(-4,4),ylim=c(2,12),xlab="Temperature anomaly above regional mean",main="Temperature drives presence of \n forests (black) and barrens (red)")
tplotfun(6,"temp","temp",plot=2,pch=17)
tplotfun(0,"tempL2","tempL4",plot=2,pch=18) #for UP=FALSE, best lags are L2 and L2
# tplotfun(4,"tempL4","temp",plot=2,pch=0,Cbar=2,Cfor=1)
legend("topright",c("South CA region","LTM sites","Baja CA region"),pch=c(16,18,17),box.col=NA,cex=0.8)
# legend("top",c("South CA region","LTM sites","Baja CA region","Central CA region"),pch=c(16,18,17,0),box.col=0,cex=0.7)
# legend("topright",c("Central CA region","South CA region","Baja CA region"),pch=c(0,16,17),col=8,box.col=0,pt.cex=1.25,cex=0.8)
labs=c(2,4,6,8,10,12,14); axis(side=4,at=labs,labels=round(1/labs,2),las=2); mtext("Exit probability",side=4,padj=5);



tplotfun(5,"tempL1","tempL2",plot=1,pch=16,cent=FALSE,Qcut=0.995,xlim=c(13,24),ylim=c(2,14),xlab="Summer temperature",main="Temperature drives presence of \n forests (black) and barrens (red)")
tplotfun(6,"temp","temp",plot=2,pch=17,cent=FALSE,Qcut=0.995)
tplotfun(0,"tempL2","tempL4",plot=2,pch=18,cent=FALSE,Qcut=0.995)
tplotfun(4,"tempL4","temp",plot=2,pch=0,cent=FALSE,Qcut=0.995)
legend("topright",c("South CA region","LTM sites","Baja CA region","Central CA region"),pch=c(16,18,17,0),box.col=0,cex=0.7)
labs=c(2,4,6,8,10,12,14); axis(side=4,at=labs,labels=round(1/labs,2),las=2); mtext("Exit probability",side=4,padj=5);






#Figure A2
yru=unique(floor(yrs)); ts=cbind(reg=rep(asvt(datUP[,"Reg","info"]),each=ncol(datUP)), yr=rep(yru,nrow(datUP)), res=1-asvt(datUP[,,"UPbarren"]))
ts[ts[,1]==2,1]=1; ts[ts[,1]==3,1]=4; tsa=serialAgg(ts,c("reg","yr"),"res",FUN=mean,CatsCode=FALSE)
tsa=rbind(tsa,cbind(0,serialAgg((gamdat),"yr","res",FUN=mean)))
matplot(yru,1-t(matricize(tsa))[,2:3],type="l",lty=1,lwd=3,ylim=0:1,las=1,xlab="",ylab="Proportion of region barren")
legend("topleft",c("Washington and Oregon","Northern CA and Central CA"),col=1:2,lwd=6,seg.len=0.5,box.col=rgb(1,1,1,alpha=0))
matplot(yru,1-t(matricize(tsa))[,-(2:3)],col=c(3,4,6),type="l",lty=1,lwd=3,ylim=0:1,las=1,xlab="",ylab="Proportion of region barren")
legend("topleft",c("LTM Reefs","Southern CA","Mexico"),col=c(3:4,6),lwd=6,seg.len=0.5,box.col=rgb(1,1,1,alpha=0))


#Figure A3
library(ncf)
spatcorUP=function(var="UPbarren",N=1e4){
     sel=which(edatAnnUP[,"Reg","info"]>4); if(var!="dbar") br=edatAnnUP[sel,,var] else {
          br=edatAnnUP[sel,,"UPrecover"]; bna=is.na(br); br[bna]=-edatAnnUP[sel,,"UPdie"][bna]; }
     cr=cor2(t(br)); gps=edatAnnUP[sel,c("lon","lat"),"info"]; xdist=gcdist(gps[,1],gps[,2]);
     tp=cbind(D=as.numeric(xdist[lower.tri(cr)]),C=as.numeric(cr[lower.tri(cr)]))
     set.seed(1); tpp=rbind(tp[sample(1:nrow(tp),N),],tp[sample(which(tp[,1]<2),N),]);
     tpp[,1]=log(tpp[,1]); gam(C~s(D,k=5),data=data.frame(tpp));
}
mstate=spatcorUP(); mdstate=spatcorUP("dbar"); msst=spatcorUP("temp"); mres=spatcorUP("barren"); mno3s=spatcorUP("no3s");
Ds=seq(-1.7,7,len=100); mtp=function(md,Di=Ds) predict.gam(md,data.frame(D=Di),"response");

matplot(Ds,cbind(mtp(mstate),mtp(mdstate),mtp(mres),mtp(msst),mtp(mno3s)-0.01),lwd=3,lty=1,type="l",col=c(1,4,3,6,7),las=1,ylab="Spatial correlation",xlab="Distance, km",xaxt="n")
axis(1,at=c(-1.7,0,2,4,7),labels=c(round(exp(c(-1.7,0,2,4)),1),1100)); abline(0,0,lty=2,lwd=2,col=8);
legend(x=0.72,y=0.77,c("Summer SST","Spring nitrate","Reef kelp cover","Reef state","Reef state change"),col=c(6,7,3,1,4),lwd=6,seg.len=0.5,box.col=rgb(1,1,1,alpha=0),cex=0.8)
round(mtp(mstate,log(c(0.18,1,5,50,100,300))),2)






