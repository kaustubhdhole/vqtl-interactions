getBin <- function(m)
{ binTable = read.csv("C:\\Users\\Delll\\Desktop\\Map_QTL\\Bins\\Bins.csv",stringsAsFactors=FALSE);
  bin = binTable$bin[which(binTable$marker==m)]
  return(bin)
}

getBinRepresentatives <- function(b)
{ brepTable = read.csv("C:\\Users\\Delll\\Documents\\Educational\\2014 Thesis\\Rqtl\\QTL mapping in an Environment Gradient\\bins\\first_markers_of_bins.csv",stringsAsFactors=FALSE);
 m = brepTable$first_marker[which(brepTable$bin==b)]
 return(m)
}

data3 = read.csv("C:\\Users\\Delll\\Desktop\\Map_QTL\\Interactions\\DataFile3_withPvalues.csv",stringsAsFactors=FALSE);
data4 = NULL
for(r in 1:nrow(data3)){
	env = data3$Environment[r]
	m1 = data3$M1[r]; b1 = getBin(paste("X",m1,sep="")); rm1 = getBinRepresentatives(b1);
	m2 = data3$M2[r]; b2 = getBin(paste("X",m2,sep="")); rm2 = getBinRepresentatives(b2);
	q2 = data3[r,5]
	vq2 = data3[r,6]
	data4 = rbind(data4,c(env,rm1,rm2,q2,vq2))
	#data4 = rbind(data4,c(env,b1,b2,q2,vq2))
}
write.csv(data4,"C:\\Users\\Delll\\Desktop\\Map_QTL\\Bins\\DataFile3BinsFirstMarkers.csv")
