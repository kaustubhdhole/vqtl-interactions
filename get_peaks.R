
flodcutoff = 3  # Set LODF cutoff. Peaks above this cut-off will be considered.
bflodcutoff = 3 # Set LODBF cutoff.Peaks above this cut-off will be considered.
chrbegin=c(1,210,1113,1281,2547,3110,3445,4539,4981,5527,6384,7221,8433,9274,9908,11000,11623)

DataFile1 = NULL;
for(env in envnames){
etable=read.csv(file.path("QTL",paste(env,".csv",sep="")))
lodctff=flodcutoff;
phenocol=as.numeric(phenodata[,env])
pcmarkers=NULL
for(i in 1:16)
{ ctable = etable[chrbegin[i]:chrbegin[i+1],]
  k = 51
  while(k<(nrow(ctable)-49))
  { j=which.max(ctable$LOD[k:(k+49)])
      pos=j+k-1
	  if(ctable$LOD[pos]>lodctff)
	  {
		if(sum(ctable$LOD[pos]<ctable$LOD[(pos+1):(pos+50)],na.rm=T)==0 & sum(ctable$LOD[pos]<ctable$LOD[(pos-50):(pos-1)],na.rm=T)==0)
		{ pcmarkers = rbind(pcmarkers,ctable[pos,])
		  marker = ctable[pos,1];
		  genocol=genodata[,match(marker,colnames(genodata))];
		  y1 = phenocol[which(genocol=="BB")]; y2 = phenocol[which(genocol=="RR")];
		  DataFile1 = rbind(DataFile1,c(env,ctable[pos,1],"-",mean(y2,na.rm=T),var(y2,na.rm=T),mean(y1,na.rm=T),var(y1,na.rm=T)));
		}
	  }
	  k = k + 50
  }
}
write.csv(pcmarkers,file.path("QTL",paste(env,".csv",sep="")))
print(env)

}

for(env in envnames){
etable=read.csv(file.path("vQTL",paste(env,".csv",sep="")))
lodctff=bflodcutoff;
pcmarkers=NULL
phenocol=as.numeric(phenodata[,env])
for(i in 1:16)
{ ctable = etable[chrbegin[i]:chrbegin[i+1],]
  k = 51
  while(k<(nrow(ctable)-49))
  { j=which.max(ctable$LODBF[k:(k+49)])
      pos=j+k-1
	  if(ctable$LODBF[pos]>lodctff)
	  {
		if(sum(ctable$LODBF[pos]<ctable$LODBF[(pos+1):(pos+50)],na.rm=T)==0 & sum(ctable$LODBF[pos]<ctable$LODBF[(pos-50):(pos-1)],na.rm=T)==0)
		{ pcmarkers = rbind(pcmarkers,ctable[pos,])
		  marker = ctable[pos,1];
		  genocol=genodata[,match(marker,colnames(genodata))];
		  y1 = phenocol[which(genocol=="BB")]; y2 = phenocol[which(genocol=="RR")];
		  DataFile1 = rbind(DataFile1,c(env,"-",ctable[pos,1],mean(y2,na.rm=T),var(y2,na.rm=T),mean(y1,na.rm=T),var(y1,na.rm=T)));
		}
	  }
	  k = k + 50
  }
}
write.csv(pcmarkers,file.path("vQTL",paste(env,".csv",sep="")))
print(env)
}
colnames(DataFile1)=c("Environment","QTL","vQTL","mRM","varRM","mBY","varBY")
write.csv(DataFile1,file.path("Output","DataFile1_2.csv"))