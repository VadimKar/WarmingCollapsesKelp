


#Script 2: Assemble diver and satellite data at LTM reefs


library(abind); library(ncf); library(ncdf4); library(raster); library(mgcv); #library(home)


setwd("/Users/aqeco/Dropbox/Projects/KelpAllee/CodebaseFinal/Data")
nc_e=nc_open('CAkelpCanopyEnv_2021.nc')
pdo0=import("WC_climate_annual.csv"); beuti=import("BEUTI_monthly.csv"); cuti=import("CUTI_monthly.csv");

KFM5m=import("LTM data\\KFM_5mQuadrat_Summary_1996-2021.csv")[,c("SiteNumber","SurveyYear","Species","MeanDensity_sqm")]
KFM1m=import("LTM data\\KFM_1mQuadrat_Summary_1982-2021.csv")[,c("SiteNumber","SurveyYear","Species","MeanDensity_sqm","SiteCode")]
KFMgps=import("LTM data\\KFM_transect_endpoints.csv")[,c("Site_Code","Long_DD","Lat_DD")]
LTERkelp=import("LTM data\\Annual_Kelp_All_Years_20211020.csv")
LTERq=import("LTM data\\Annual_Quad_Swath_All_Years_20211020.csv"); 
LTERgps=import("LTM data\\Benthic_Transect_Depths_coors.csv")[,c("SITE","TRANSECT","DEPTH_MLLW_M","LATITUDE","LONGITUDE")]


###### Import environment and kelp canopy data
#For spatial analyses: assign each patch a row,column index based on where it falls in a region-specific matrix of patches
#The regions here, RegionsNew, are build to minimize the size of the matrices (and save RAM)
getcords=function(ids){ 
     for(RegHenry in 8:1){
          gps=coordinates(raster(paste0("geomat\\R",RegHenry,"_kelpCanopy_test_202104.tif"),6))
          lonu=unique(gps[,1]); dlon=mean(diff(lonu))/2; latu=unique(gps[,2]); dlat=mean(diff(latu))/2; remove(gps);
          mlon=findInterval(ids[,1],c(lonu[1]-dlon,lonu+dlon))-1
          mlat=length(latu)-(findInterval(ids[,2],rev(c(latu[1]-dlat,latu+dlat))))
          ids=cbind(ids,mlon,mlat,lonu[pmax(mlon,1)],latu[pmax(mlat,1)],RegHenry)
     }
     dst=t(apply(ids,1,function(x){ i=matrix(x[-(1:2)],nrow=5); d=abs(x[1]-i[3,])+abs(x[2]-i[4,]); c(x[1:2],min(d),i[c(1,2,5),which.min(d)]); }))
     list(ids,dst)
}
tsel=1:148; edatget=function(vnm,mod=1) cbind(ncvar_get(nc_e,"lon"), ncvar_get(nc_e,"lat"), mod*ncvar_get(nc_e,vnm)[,tsel])
edat0=abind(edatget("area",1/900),edatget("temperature"),edatget("nitrate"),NA*edatget("nitrate"),NA*edatget("nitrate"),NA*edatget("nitrate"),along=3)
eyrs=(ncvar_get(nc_e,"year")+(ncvar_get(nc_e,"quarter")-1)/4)[tsel]; colnames(edat0)=c("lon","lat",eyrs); remove(nc_e); dimnames(edat0)[[3]]=c("area","temp","no3","hs","parmu","parmx");
edatAlleeUse=rowMeans(is.na(edat0[,-(1:2),1]))<0.3378 & rowMeans(edat0[,-(1:2),1]>0,na.rm=TRUE)>0.3









###### Import diver data
CD=KFM1m$Species; Y=KFM1m$SurveyYear;
KFM5m$Species=1.1*(KFM5m$Species==2002.25)
KFM1m$Species=1.1*(CD==2002 & Y<1996) + 1.2*(CD==2002.5) + 6.1*(CD==11006) + 6.2*(CD==11005)
KFMgps=aggregate(as.matrix(KFMgps[,-1])~KFM1m$SiteNumber[match(KFMgps$Site_Code,KFM1m$SiteCode)],FUN=mean)
KFMall=round(serialAgg(rbind(KFM5m,KFM1m[,-5]),1:3,4,CatsCode=FALSE),2); KFMdat=unique(KFMall[,1:2]);
KFMdat=cbind(KFMdat,KFMgps[match(KFMdat[,1],KFMgps[,1]),2:3]);
for(FG in c(1.1,1.2,6.1,6.2)){
     dati=KFMall[KFMall[,"Species"]==FG,]; KFMdat=cbind(KFMdat,NA); if(nrow(dati)==0) next;
     for(i in 1:nrow(dati)) KFMdat[KFMdat[,1]==dati[i,1] & KFMdat[,2]==dati[i,2],ncol(KFMdat)]=dati[i,4]
}
colnames(KFMdat)=c("Site","Yr","Lon","Lat", "AdKelp","JvKelp","PrpU","RedU")


#LTER: combine site and transect
DT=LTERkelp; mkLTERst=function(DT) 100+as.numeric(SerialReplace(DT$SITE,c("ABUR","AHND","AQUE","BULL","CARP","GOLB","IVEE","MOHK","NAPL","SCDI","SCTW"),1:11))+DT$TRANSECT/10
LTERkads=serialAgg(cbind(mkLTERst(DT),DT$YEAR,FG=1.1,(DT$FRONDS>0)/80)[DT$COMMON_NAME=="Giant Kelp",], 1:3, 4, CatsCode=FALSE)

DT=LTERq; DT$COUNT=DT$COUNT/DT$AREA; # a few species sampled in 20m^2 swaths
LTERquad0=serialAgg(cbind(ST=mkLTERst(DT),DT), c("ST","YEAR","TAXON_GENUS","GROUP","MOBILITY"), "COUNT", FUN=function(x) mean(x,na.rm=TRUE), CatsCode=FALSE)
GN=LTERquad0[,"TAXON_GENUS"]; Groups=1.2*(GN=="Macrocystis") + 6.1*(GN=="Strongylocentrotus") + 6.2*(GN=="Mesocentrotus")
LTERquad=serialAgg(apply(cbind(LTERquad0[,1:2],FG=Groups,N=LTERquad0[,6])[Groups>0,],2,as.numeric), c("ST","YEAR","FG"), "N", CatsCode=FALSE)

LTERdat0=rbind(LTERkads,LTERquad); LTERdat0=LTERdat0[!LTERdat0[,"FG"]%in%c(0,4.01),]; 
LTERdat0=LTERdat0[order(LTERdat0[,3]),]; LTERdat=round(cbind(LTERkads[,1:2], matrix(LTERdat0[,-(1:3)],nrow=nrow(LTERkads))),2);

LTERdat=cbind(LTERdat[,1:2], LTERgps[match(LTERdat[,1],mkLTERst(LTERgps)),5:4], LTERdat[,-(1:2)])
colnames(LTERdat)=c("Site","Yr","Lon","Lat","AdKelp","JvKelp","PrpU","RedU")

alldat0=as.matrix(rbind(KFMdat,LTERdat)); adlm=unique(alldat0[,c(1,3:4)]);
saveRDS(alldat0,"ProcessedData/alldat0.rds")






###### Subselect environment and canopy data that are near LTM reefs
adlm=unique(alldat0[,c(1,3:4)]); 
latSBCI=c(34.5,33.4); edat0SBCI=edat0[,2,1]<latSBCI[1] & edat0[,2,1]>latSBCI[2] & edat0[,1,1]<(-119)
RScands=cbind(1:nrow(edat0),edat0[,1:2,1])[edat0SBCI,]; dim(RScands);
Rad=70; CID70=apply(as.matrix(adlm[,3:2]),1,function(x) RScands[closest(x,RScands[,3:2],Rad/1e3),1]); quantile(unlist(lapply(CID70,length)),c(0.01,0.5,0.99)); length(unlist(CID70))
CID1=CID70; 
Ysat=unique(floor(eyrs))
CID2=lapply(CID1, function(x){
     if(length(x)==1) x=c(x,x); edati=edat0[x,-(1:2),];
     tsia=apply(abind(edati[,c(FALSE,TRUE,FALSE,FALSE),],edati[,c(FALSE,FALSE,TRUE,FALSE),],along=4),1:3,mean,na.rm=TRUE)
     tsia[,,"temp"]=edati[,c(FALSE,FALSE,TRUE,FALSE),"temp"]; tsia[,,"no3"]=edati[,c(FALSE,FALSE,TRUE,FALSE),"no3"];
     tsia[,,"hs"]=apply(abind(edati[,c(TRUE,FALSE,FALSE,FALSE),"hs"],cbind(NA,edati[,c(FALSE,FALSE,FALSE,TRUE),"hs"])[,1:length(Ysat)],along=3),1:2,max,na.rm=TRUE)
     tsia[abs(tsia)==Inf]=NA; tsia=abind(tsia,"recr"=NA*tsia[,,"temp"],
                "no3s"=edati[,c(FALSE,TRUE,FALSE,FALSE),"no3"],"temps"=edati[,c(FALSE,TRUE,FALSE,FALSE),"temp"],
                "tempw"=edati[,c(TRUE,FALSE,FALSE,FALSE),"temp"],"tempf"=edati[,c(FALSE,FALSE,FALSE,TRUE),"temp"]); 
     
     
     Ls=t(apply(tsia[,,1],1,function(x){ r=rle(x==0); rep(r$lengths*r$values,r$lengths); }))
     up=t(apply(Ls>1,1,diff)==-1); up[Ls[,-ncol(Ls)]<2]=NA; dwn=t(apply(Ls>1,1,diff)==1); dwn[Ls[,-ncol(Ls)]>1]=NA;
     
     out=cbind(Ysat, colMeans(Ls>1), colMeans(tsia[,,1]>0), colMeans(tsia[,,1]), c(colMeans(up,na.rm=TRUE),NA), c(colMeans(dwn,na.rm=TRUE),NA), apply(tsia[,,c(2:4,6:11)],2:3,mean,na.rm=TRUE))
     out=rbind(cbind(1982:1983,NA*out[1:2,])[,-2], out); if(!2021%in%Ysat) out=rbind(out,c(2021,NA*out[1,-1])); out;
})
CID2n=CID2; for(i in 1:length(CID2)) CID2n[[i]]=cbind(adlm[i,1],CID2[[i]],NA,NA,NA,NA)[CID2[[i]][,1]%in%(1982:2021),1:16] 
alldatF=round(do.call("rbind",CID2n),2); refs=match(codify(alldatF[,1:2]),codify(alldat0[,1:2]));
alldatF=cbind(alldatF[,1:2],alldat0[refs,-(1:2)],alldatF[,-(1:2)]); colnames(alldatF)[tail(1:ncol(alldatF),14)]=c("pBarren","pPres","mCvr","recover","extinct","temp","no3","hs","parmx","recr","no3s","temps","tempw","tempf");





#Transform into 3d array
alldati=alldatF; pdo=pdo0[pdo0[,1]<2022,];
YA=sort(unique(alldati[,2])); SA=sort(unique(alldati[,1])); ada=array(NA,c(length(SA),length(YA),ncol(alldati)));
for(i in 1:nrow(alldati)) ada[SA==alldati[i,1],YA==alldati[i,2],]=alldati[i,]
ada=abind(ada,NA*ada[,,1:3],along=3); for(Y in which(pdo[,1]>1981)) ada[,pdo[Y,1]-1981,dim(ada)[3]+(-2:0)]=matrix(as.numeric(pdo[Y,-1]),dim(ada)[1],3,byrow=TRUE);
dimnames(ada)[[3]]=c(colnames(alldati),"npgo","mei","pdo")





#Adding more indeces:
upIndxTrm=function(cdat,yrs=1988:2021,latbin=4,FUN=mean) binapply(cdat[cdat$year%in%yrs,2+latbin],12,FUN=FUN)
upIndx=cbind(yr=1988:2021,beutiAv=upIndxTrm(beuti),beutiMx=upIndxTrm(beuti,FUN=max),cutiAv=upIndxTrm(beuti),cutiMx=upIndxTrm(beuti,FUN=max))

rollmeanFill=function(x,ny,binkeeps=c(F,F,T,F),ARends=TRUE){
     # if(mean(is.na(x))>0.75) return(NA*x[binkeeps])
     binl=length(binkeeps); ra=filter(x,rep(1/(ny*binl),ny*binl))[binkeeps]; 
     if(!ARends) return(ra); mid=round(length(ra)/2); miss=which(is.na(ra));	
     #These settings will need refining, if we keep arma approach:
     ARnsset=c(1,0,0); ARsset=list(order=c((ny<8),1,1),period=ny);
     ARF=predict(arima(ra[-miss],ARnsset,ARsset,method="ML"),n.ahead=sum(miss>=mid))$pred
     ARB=rev(predict(arima(rev(ra[-miss]),ARnsset,ARsset,method="ML"),n.ahead=sum(miss<mid))$pred)
     c(ARB,ra[-miss],ARF)
}
mno3=colMeans(edat0[edatAlleeUse,-(1:2),][,,"no3"]); pdo=pdo0[pdo0[,1]<2021,];
pdo=cbind(pdo, rbind(NA,NA,NA,NA,NA,cbind(n2=rollmeanFill(mno3,2),n5=rollmeanFill(mno3,5),n10=rollmeanFill(mno3,10))))
par(mfrow=1:2); matplot(1984:2020,apply(na.omit(pdo[,c(2,5)]),2,nrm),col=c(1,4),main="R2=0.54:   NPGO (black) vs. \n 2y running regional mean NO3 (blue)",type="l",lty=1,lwd=2,ylab="Relative values",xlab="")
matplot(1984:2020,apply(na.omit(cbind(-pdo$PDO,pdo[,c(5)])),2,nrm),col=c(1,4),main="R2=0.62:   -PDO (black) vs. \n 2y running regional mean NO3 (blue)",type="l",lty=1,lwd=2,ylab="Relative values",xlab="")
plot(na.omit(pdo[,c(2,5)]),ylab="2y running regional mean NO3")
plot(pdo$PDO,pdo[,5],ylab="2y running regional mean NO3",xlab="PDO")







#Add lagged variables
gvl=function(L=0,V="temp"){ outi=ada[,,V][,c(0:L,1:(ncol(ada)-L))]; outi[,0:L]=NA; as.vector(outi); }
gamdat=data.frame(mc=gvl(0,"mCvr"),kp=gvl(0,"pPres"),res=1-gvl(0,"pBarren"),die=gvl(0,"extinct"),rec=gvl(0,"recover"),
                  lu=log(1+gvl(0,"PrpU")+gvl(0,"RedU")),temp=gvl(),tempL1=gvl(1),tempL2=gvl(2),tempL3=gvl(3),tempL4=gvl(4),tempL5=gvl(5),
                  no3=gvl(0,"no3"),no3L1=gvl(1,"no3"),no3L2=gvl(2,"no3"),no3L3=gvl(3,"no3"),
                  no3s=gvl(0,"no3s"),no3sL1=gvl(1,"no3s"),no3sL2=gvl(2,"no3s"),no3sL3=gvl(3,"no3s"),
                  temps=gvl(0,"temps"),tempw=gvl(0,"tempw"),tempf=gvl(0,"tempf"),hs=gvl(0,"hs"),yr=gvl(0,2))
gamdat=cbind(gamdat, upIndx[match(gamdat$yr,upIndx[,"yr"]),-1], pdo[match(gamdat$yr,pdo$year),-1])
gamdat$restm1=1-gvl(1,"pBarren"); gamdat$dres=gamdat$restm1-gamdat$res; gamdat$logdn=log(1+as.vector(cbind(1-ada[,2:40,"pBarren"],NA))/pmax(gamdat$res,0.02));
gamdat$logmn=as.vector(t(apply(1-ada[,,"pBarren"],1,function(x) log(1e-3 + x/mean(x,na.rm=TRUE)))))
gamdat$trans=as.vector(t(apply(ada[,,1],1,function(x) rep(na.omit(x)[1],length(x))))); gamdat$site=floor(gamdat$trans);
gamdat$resn=nrm2(gamdat$res,gamdat$trans); 
gamdat$nsatpix=unlist(lapply(CID1,length))[match(gamdat$trans,adlm[,1])];
gamdat$recr=log(pmax(as.vector(cbind(NA,ada[,,22])[,1:40]),0.85)); #round(cor(gamdat),2);
gamdat$nPixFor=round(gamdat$nsatpix*gamdat$res); gamdat$nPixBar=round(gamdat$nsatpix*(1-gamdat$res)); 
gamdat$nPixRec=round(gamdat$nPixBar*gamdat$rec); gamdat$nPixRec[is.na(gamdat$rec) & !is.na(gamdat$res)]=0;
gamdat$nPixCol=round(gamdat$nPixFor*gamdat$die); gamdat$nPixCol[is.na(gamdat$die) & !is.na(gamdat$res)]=0;
gamdat$recL1=gvl(1,"recover"); gamdat$dieL1=gvl(1,"extinct"); 
gamdat$nPixForL1=round(gamdat$nsatpix*gamdat$restm1); gamdat$nPixBarL1=round(gamdat$nsatpix*(1-gamdat$restm1)); 
gamdat$nPixRecL1=round(gamdat$nPixBarL1*gamdat$recL1); gamdat$nPixRecL1[is.na(gamdat$recL1) & !is.na(gamdat$restm1)]=0;
gamdat$nPixColL1=round(gamdat$nPixForL1*gamdat$dieL1); gamdat$nPixColL1[is.na(gamdat$dieL1) & !is.na(gamdat$restm1)]=0;
gamdat=gamdat[!is.na(gamdat[,1]),];
gamdat$lusm=gamdat$lun=gamdat$lu; Si=gamdat$trans; for(s in unique(Si)){ sis=Si==s; lus=gamdat$lu[sis]; gamdat$lusm[sis]=mean(lus,na.rm=TRUE); gamdat$lun[sis]=(lus - mean(lus[gamdat$yr[sis]>2001],na.rm=TRUE))/sd(lus,na.rm=TRUE); }



#Trim off a few outlying observations and poorly observed transects (<5 30x30m pixels)
gamdat[,c("no3sL1","no3sL2","no3sL3")]=apply(gamdat[,c("no3sL1","no3sL2","no3sL3")],2,utrim,9.5)
gamdat[,c("no3","no3L1")]=apply(gamdat[,c("no3","no3L1")],2,utrim,0.85)
gamdat=gamdat[gamdat$nsatpix>4,]

saveRDS(gamdat,"ProcessedData/gamdat.rds")







