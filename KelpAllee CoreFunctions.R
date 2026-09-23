

#Script 1: enter core functions to use throughout analyses


#utility functions
Sys.glob=function(pathnm,...) base::Sys.glob(cpath(pathnm),...)
readRDS=function(pathnm) base::readRDS(cpath(pathnm))
saveRDS=function(x,pathnm) base::saveRDS(x,cpath(pathnm))
setwd=function(pathnm) base::setwd(cpath(pathnm))
import=function(pathnm) rio::import(cpath(pathnm))
rsq=function(GAM) round(cor(GAM$y,GAM$fitted.values)^2,3)
rsq2=function(GAM) round(cor(GAM$y>0.5,GAM$fitted.values>0.5)^2,3)
asvt=function(x) as.vector(t(x))
cpath=function(pathnm,machine=Sys.info()["sysname"]){ #not yet edited to go from linux to windows/Mac
     if(machine=="Linux") pathnm=SerialReplace(pathnm,c("\\\\","C:","Users/vkara", "Users/aqeco"),c("/","","home/vadim", "home/vadim"))
     if(machine=="Darwin") pathnm=SerialReplace(pathnm,c("\\\\","C:","vkara"),c("/","","aqeco"))
     if(machine=="Windows") pathnm=SerialReplace(pathnm,c("/","aqeco"),c("\\\\","vkara")); pathnm;
}


#transformation functions
#' Normalize x. Default is global; if levels Lvl supplied will be local.
#' @param x vector of values to normalize
#' @param Lvl optional vector denoting category of each x value for local lormalization
#' @param xuse optional T/F vector, where FALSEs signal to omit corresponding entries in x when calculating mean and sd for normalization
nrm2=function (x, Lvl = rep(1, length(x)), xuse = rep(TRUE, length(x)), zeroNullSD = FALSE, center = FALSE){
     normcent = function(x, xuse) { SD = sd(x[xuse], na.rm = TRUE)
     if ((zeroNullSD & SD == 0) | center) SD = 1
     (x - mean(x[xuse], na.rm = TRUE))/SD
     }
     xout = x; for (i in unique(Lvl)) xout[Lvl == i] = normcent(x[Lvl == i], xuse[Lvl == i]); xout
}
nrm=function(x,minScl=min(x,na.rm=TRUE),Lvl=rep(1,length(x))){
     normfun=function(x) (x - min(x,na.rm=TRUE)) / (max(x,na.rm=TRUE) - min(x,na.rm=TRUE))
     xout=x; for(i in unique(Lvl)) xout[Lvl==i]=normfun(x[Lvl==i]); xout;
}
matricize=function(di){
     tl=range(di[,2],na.rm=TRUE); ts=sort(unique(di[,2])); if(length(unique(diff(ts)))>1){ print("goofy times!!"); ts=tl[1]+(0:diff(tl)); };
     pops=sort(unique(di[,1])); m=matrix(NA,length(pops),length(ts)); for(i in 1:nrow(di)) m[pops==di[i,1], ts==di[i,2]]=di[i,3]; m[1:nrow(m),]; 
     rownames(m)=pops; colnames(m)=ts; m;
}
binmean=function(x,binsize,NRM=FALSE,FUN=function(x) mean(x,na.rm=TRUE)){ out=apply(matrix(x,nrow=binsize),2,FUN); if(NRM) out=out/max(out,na.rm=TRUE); out; }
#' Apply a function to clusters grouped in a simple pattern
#' @param x vector of values to bin
#' @param nreps number of individuals in each bin
#' @param clustered whether groups are clustered: rep(1:nreps,each=5)=TRUE, rep(1:nreps,5)=FALSE
#' @param FUN function to apply; defaults to mean
#' @param ... other inputs to FUN
binapply=function(x,nreps,clustered=TRUE,FUN=function(x) mean(x,na.rm=TRUE),...) apply(matrix(x,nrow=nreps,byrow=!clustered),2,FUN,...)
SerialReplace=function(x,Targets,Replacements){
     Ntargets=length(Targets); NReplace=length(Replacements); maxN=max(Ntargets,NReplace);
     if(Ntargets!=NReplace) print("SerialReplace warning: unequal number of targets and replacements; recycling shorter argument")
     Index=data.frame(Targets,Replacements); for(i in 1:maxN) x=gsub(as.character(Index[i,1]), as.character(Index[i,2]), x); x;
}

#Changing some defaults
mean=function(x,na.rm=TRUE,...) base::mean(x,na.rm=na.rm,...)
sum=function(x,na.rm=TRUE,...) base::sum(x,na.rm=na.rm,...)
min=function(x,na.rm=TRUE,...) base::min(x,na.rm=na.rm,...)
max=function(x,na.rm=TRUE,...) base::max(x,na.rm=na.rm,...)
median=function(x,na.rm=TRUE,...) stats::median(x,na.rm=na.rm,...)
utrim=function(x,lim=quantile(x,q,na.rm=TRUE),q=0.95) base::pmin(x,lim)
ltrim=function(x,lim=quantile(x,q,na.rm=TRUE),q=0.05) base::pmax(x,lim)
cor=function(x,...,use="pairwise.complete.obs") stats::cor(x,use=use,...)



#Analysis functions
codify=function(x,cols=1:ncol(x),sep="_") as.matrix(cbind(Index=apply(x[,cols],1,paste0,collapse=sep),x[,-cols]))
#' aggregate potentially multiple values in columns of x by categories given by other columns of x
#' @param CatsCode whether to also return the composite category index
serialAgg=function(x,AggCats,AggTarg=NULL,CatsCode=length(AggCats)>1,FUN=function(x) sum(x,na.rm=TRUE)){
     if(is.null(AggTarg)){
          if(is.numeric(AggCats)) AggTarg=(1:ncol(x))[-AggCats]
          if(is.character(AggCats)) AggTarg=colnames(x)[!colnames(x)%in%AggCats]
     }
     Numbers=prod(apply(t(t(x[,AggCats])),2,is.numeric))
     ncat=length(AggCats); if(ncat==1) Cats=as.character(x[,AggCats]) else Cats=codify(x[,AggCats]);
     agged=as.matrix(aggregate(x[,AggTarg], by=list(Cats), FUN=FUN))
     if(!CatsCode){
          if(ncat>1) agged=cbind(matrix(unlist(strsplit(agged[,1],"_")),ncol=ncat,byrow=TRUE), matrix(agged[,-1],ncol=ncol(agged)-1))
          if(Numbers) agged=t(apply(agged,1,as.numeric))
          colnames(agged)=colnames(cbind(x[,c(AggCats[1],AggCats)],x[,c(AggTarg,AggTarg[1])]))[c(1,3:(ncol(agged)+1))]
     }; agged;
}
#' Function to fit, evaluate, and plot 1d splines in mgcv
#' @param dat data matrix or frame
#' @param varSeq dependent variable in spline. Can be column index or name in dat. If GAM supplied, defaults to all.vars(formula(GAM))[2].
#' @param varResp response variable in spline. Can be column index or name in dat. If GAM supplied, defaults to all.vars(formula(GAM))[1].
#' @param GAM an optional pre-fitted gam to use. NOTE: if variable names in GAM don't match those in dat, need to sperately specify varSeq, varResp. 
#' @param k if fitting gam, smoothness parameter. Can specify a range of k values, in which one will be selected based on BIC. Defaults to k=5:30.
#' @param fam if fitting gam, family of error distribution.
#' @param wt weights to use if fitting gam. Can be vector or index of dat.
#' @param seqRng range of x-values to project; defaults to 0.01 and 0.99 quantiles.
#' @param plotgive whether to return a plot or only associated data.
#' @param seFact width of mean confidence intervals to plot. Defaults to 1.96SE.
#' @param spanSignif min length of sign-of-change-runs to include when analyzing changes. Measured as proportion of seqRng, defauts to 10%. If NA trend signifinance not evaluated.
#' @param seFill color level of SE polygon shade; defaults to gray if col in ... unspecified.
#' @... optional inputs for plot()
gampred=function(dat,varSeq=NULL,varResp=NULL,k=5:30,fam=gaussian,wt=NULL,seqRng=NULL,GAM=NULL, plotgive=TRUE,plotfun=plot,seFact=1.96,spanSignif=0.1,seFill=0.2,seCol=NULL,Lcol=NULL,predType="link",...){
     #assumes that either (a) GAM is built on variables named in dat or (b) gam is unbuilt
     if("matrix"%in%class(dat)) dat=data.frame(dat);
     # if(length(varSeq)>1 & length(varResp)>1){ dat=data.frame(x=varSeq,y=varResp); varSeq="x"; varResp="y"; }
     require(mgcv); xvr=varSeq; yvr=varResp;
     if(!is.null(GAM)){ k="k given"
     if(is.numeric(varSeq)) varResp=colnames(dat)[varSeq]; if(is.numeric(varResp)) varResp=colnames(dat)[varResp];
     if(is.null(varSeq)) varSeq=all.vars(formula(GAM))[2]; if(is.null(varResp)) varResp=all.vars(formula(GAM))[1];
     } else {
          if(is.numeric(varSeq)) colnames(dat)[varSeq]="X" else colnames(dat)[colnames(dat)==varSeq]="X"; varSeq="X";
          if(is.numeric(varResp)) colnames(dat)[varResp]="Y" else colnames(dat)[colnames(dat)==varResp]="Y"; varResp="Y";
          if(length(k)){ Bs=1[-1]; for(ki in k) Bs=c(Bs,tryCatch(BIC(gam(Y~s(X,k=ki),data=dat)),error=function(e) NA)); k=k[c(tail(1+which(diff(Bs)<(-2)),1),1)[1]]; }
          if(length(wt)==1) wt=dat[,wt]; GAM=gam(Y~s(X,k=k), data=dat,weights=wt,family=fam);
          # if(length(k)==1 & k==2) GAM=gam(Y~X, data=dat,weights=wt,family=fam);
     }
     
     if(is.null(seqRng)) seqRng=quantile(dat[,varSeq],c(0.01,0.99),na.rm=TRUE);
     N=1e2; datNew=data.frame(dat[rep(c(which(!is.na(dat))[1],1)[1],N),]);
     varNum=unlist(lapply(dat,class))=="numeric"; datNew[,varNum]=colMeans(t(t(dat[,varNum])),na.rm=TRUE); 
     xseq=seq(seqRng[1],seqRng[2],len=N); datNew[,varSeq]=xseq;
     
     PRED=predict.gam(GAM,newdata=datNew,se.fit=TRUE,type=predType); Bnd=cbind(c(xseq,rev(xseq)),c(PRED$fit+seFact*PRED$se.fit,rev(PRED$fit-seFact*PRED$se.fit)));
     rsq=tryCatch(cor(GAM$fitted.values,GAM$y)^2,error=function(e) NA); pval=summary(GAM)$s.table[,"p-value"][1]; dSignif=1+0*xseq;
     if(!is.na(spanSignif)){
          require(tsgam); #remotes::install_github("gavinsimpson/tsgam") #tsgam not available on CRAN
          Deriv=fderiv(GAM,newdata=datNew); muDeriv=Deriv$derivatives[[varSeq]]$deriv; se2Deriv=1.96*Deriv$derivatives[[varSeq]]$se.deriv;
          dSignif[(muDeriv+se2Deriv)<0]=4; dSignif[(muDeriv-se2Deriv)>0]=2;
          dSrle=rle(dSignif); dSignif[rep(dSrle$lengths,dSrle$lengths)<(spanSignif*N)]=1;
     }
     OUT=list(mu=cbind(xseq,PRED$fit),  polyse=Bnd,   pts=cbind(dat[,varSeq],dat[,varResp]),  dSignif=dSignif,  rsq=rsq,  pval=pval,  k=k); if(!plotgive) return(OUT);
     
     xpt=list(...); chkin=(!c("main","ylab","xlab","col")%in%names(xpt)); rgbCol=c(0,0,0); if(!chkin[4]) rgbCol=col2rgb(xpt$col)[,1]/255;
     if(chkin[2]) xpt$ylab=yvr; if(chkin[3]) xpt$xlab=xvr; if(chkin[1]) xpt$main=paste0(c("R2=",round(rsq,3),"__Pval=",round(pval,3),"__k=",k),collapse="");
     if(!is.null(seCol)) seCol=col2rgb(seCol)/255 else seCol=rgbCol;
     do.call(plotfun, c(list(OUT[["pts"]]), xpt)); polygon(OUT[["polyse"]], border=NA, col=do.call(rgb, c(as.list(seCol),alpha=seFill))); 
     if(!is.na(spanSignif)){ segs=c(0,cumsum(diff(dSignif)!=0)); for(si in unique(segs)) lines(OUT[["mu"]][segs==si,],lwd=4,col=OUT[["dSignif"]][segs==si][1]); do.call(points, c(list(OUT[["pts"]]), xpt));
     } else {
          if(is.null(Lcol)) Lcol=do.call(rgb, as.list(rgbCol))
          lines(OUT[["mu"]],lwd=4,col=Lcol)
     }
}
#' distance in km between two GPS points in decimal degrees
#' @param v0 vector of long1, lat1, long2, lat2
distCalc=function(v0){
     v=v0*pi/180; long1=v[1]; lat1=v[2]; long2=v[3]; lat2=v[4];
     a = (sin((lat2-lat1)/2))^2 + cos(lat1)*cos(lat2)*(sin((long2-long1)/2))^2
     c = 2*atan2(sqrt(a), sqrt(1-a)); return(d=6378.145*c);
}
#' Find which pair of coordinates from a candidate list is closest to coordinates of target, and distance in km.
#' @param target length-2 vector of DD GPS latitude and longitude
#' @param candidates 2-column matrix of DD GPS latitudes and longitudes
#' @param thresh if not NULL, return candidantes (indexed by row #) within this distance (in km) from target
closest=function(target,candidates,thresh=NULL){ 
     dists=apply(cbind(target[1],target[2],candidates),1,function(x) distCalc(as.numeric(x[c(2:1,4:3)]))); 
     if(is.null(thresh)) return(c(which.min(dists),min(dists)));
     which(dists<thresh)
}




#Plotting functions
plot2=function(...,axcol=TRUE){
     par(new=TRUE); tp=list(...); tp$xlab=NULL; Col=c(1,tp[["col"]])[1+axcol];
     tav=which(names(tp)%in%c("ylab","at","outer","col.ticks","labels","pos","ylab"))
     ta=c(side=4,col=Col,col.axis=Col,tp[tav]); tp[tav]=NULL; 
     do.call(plot,c(tp,yaxt="n",xaxt="n",xlab="",ylab="")); do.call(axis,ta);
}
SpatiotempPlot=function(spacetime,Xax=1:dim(spacetime)[2],Yax=1:dim(spacetime)[1], XaxN="X",YaxN="Y",figtitle="Title",Zlim=c(0,max(spacetime,na.rm=TRUE)),
                        cont=NULL,cexAx=1,contPlot=spacetime,cexCont=1.5*cexAx,lwCont=cexCont,Conts=NULL,contSpan=1,palette=1,las=1){
     require(fields); spacetime[is.na(spacetime)]=Zlim[1]-1;
     COL=rev(rainbow(1e3,start=0,end=0.7)); if(palette>1){ require(pals); COL=list(parula(1e3),head(tail(parula(1e3),-50),-50),terrain.colors(1e3))[[palette-1]]; }
     if(length(Zlim)==1){
          Zlim=quantile(spacetime,c(Zlim,1-Zlim),na.rm=TRUE) #can provide just a quantile percentage to set limits
          Rng=range(spacetime); if(Zlim[1]<Rng[1]) Zlim[1]=Rng[1]; if(Zlim[2]>Rng[2]) Zlim[2]=Rng[2];
     }
     spacetime[which(is.na(spacetime),arr.ind=TRUE)]=max(Zlim)+1 #NAs painted white
     image.plot(x=Xax,y=Yax,z=t(spacetime), zlim=Zlim, xlab=XaxN, ylab=YaxN,axis.args=list(las=2),
                cex.axis=cexAx, cex.lab=cexAx, legend.cex=cexAx, main=figtitle, col=COL,las=las);  box(); 
     
     #Add contour lines to plot
     if(!is.null(cont)){
          if(abs(log(max(contPlot,na.rm=TRUE),10))>2)  options(scipen=-10)
          if(contSpan>1){
               smoo1=t(apply(contPlot, 1, function(x) supsmu(1:ncol(contPlot),x,span=contSpan/ncol(contPlot))$y))
               smoo2=t(apply(smoo1, 1, function(x) supsmu(1:ncol(smoo1),x,span=contSpan/ncol(smoo1))$y))
               contPlot=smoo2
          }
          if(is.null(Conts)) contour(x=Xax,y=Yax,z=t(contPlot),add=TRUE,col=cont,lwd=lwCont,labcex=cexCont);
          if(!is.null(Conts)) contour(x=Xax,y=Yax,z=t(contPlot),levels=Conts,add=TRUE,col=cont,lwd=lwCont,labcex=cexCont);
          if(abs(log(max(contPlot,na.rm=TRUE),10))>2)  options(scipen=0)
     }
}





