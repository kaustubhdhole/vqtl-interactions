zin = read.csv("QTL_and_vQTL.csv")
bins = read.csv("Bins.csv",stringsAsFactors=FALSE)
Z=NULL
for(r in 1:nrow(zin))
{
if(zin$only.vQTL[r] || zin$onlyQTL[r] || zin$X.1[r])
	{ 
	chr=as.numeric(strsplit(as.character(zin$Chromosome[r]),"chr")[[1]][2])
	print(chr);
	  bp=zin$Base.pair[r]
	  print(bp)
	  rb=intersect(which(bins$chr==chr),which(bins$pos==bp))
	  print(rb)
	  B=bins$bin[rb]
	  M=bins$marker[rb]
	  print(zin[r,])
	  print(B)
	  print(M)
	  Z=rbind(Z,cbind(zin[r,],B,M))
	}
}
FM = read.csv("First_markers_of_bins.csv",stringsAsFactors=F);
Zinb=Z
ubins=unique(Zinb$B)
zvect=NULL
for(b in ubins)
{ rb=which(b==Zinb$B)
  bnvnames=NULL
  for(r in rb) bnvnames = c(bnvnames,as.character(Zinb$Environment[r]))
  if(length(bnvnames)>1)
  {
  bnvpairs=t(combn(bnvnames,2))
  marker = FM[FM$bin==b,2];
  rr=which(genodata[,marker]=="RR")
  bb=which(genodata[,marker]=="BB")
  vect=NULL
	for(r in 1:nrow(bnvpairs))
	 {  e1=bnvpairs[r,1];e2=bnvpairs[r,2]
		mk=strsplit(marker,"X")[[1]][2]
		
		x=phenodata[,e1];y=phenodata[,e2]
		corRM = cor(x[rr],y[rr],use="complete.obs"); covRM = cov(x[rr],y[rr],use="complete.obs")
		corBY = cor(x[bb],y[bb],use="complete.obs"); covBY = cov(x[bb],y[bb],use="complete.obs")
		
		m_rm <- mcreg(x[rr],y[rr],method.reg="Deming", mref.name=e1,mtest.name=e2, na.rm=TRUE);z1=getResiduals(m_rm)$optimized;y1=abs(z1)
		m_by <- mcreg(x[bb],y[bb],method.reg="Deming", mref.name=e1,mtest.name=e2, na.rm=TRUE);z2=getResiduals(m_by)$optimized;y2=abs(z2)
		m_all <- mcreg(x,y,method.reg="Deming", mref.name=e1,mtest.name=e2, na.rm=TRUE);z3=getResiduals(m_all)$optimized;y3=abs(z3)
		
		ttest=t.test(y1,y2)[1]$statistic;
		
		m1=mean(y1,na.rm=T);m2=mean(y2,na.rm=T);m=mean(y3,na.rm=T)
		v1=var(y1,na.rm=T);v2=var(y2,na.rm=T);v=var(y3,na.rm=T)
		vz1=var(z1,na.rm=T);vz2=var(z2,na.rm=T)
		mRRe1 = mean(x[rr],na.rm=T); mRRe2 = mean(y[rr],na.rm=T); vRRe1 = var(x[rr],na.rm=T); vRRe2 = var(y[rr],na.rm=T); 
		mBBe1 = mean(x[bb],na.rm=T); mBBe2 = mean(y[bb],na.rm=T); vBBe1 = var(x[bb],na.rm=T); vBBe2 = var(y[bb],na.rm=T); 
		zvect = rbind(zvect,c(Bin=b,Environment1 = e1,Environment2 = e2,ttest=ttest,meanRM=m1,meanBY=m2,varRM=v1,varBY=v2,meanRM_env1=mRRe1,meanBY_env1=mBBe1,varRM_env1=vRRe1,varBY_env1=vBBe1,meanRM_env2=mRRe2,meanBY_env2=mBBe2,varRM_env2=vRRe2,varBY_env2=vBBe2));
	 }
	}
}
write.csv(zvect,file.path(paste("DataFile2.csv",sep="")),row.names=F)
