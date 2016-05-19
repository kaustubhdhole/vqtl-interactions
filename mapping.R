# Set the following two variables
no_of_env = 34 # Specify the Number of environments in the phenogeno file
nperm = 10     # Specify the desired Number of Permutations

library("parallel"); 
library("doParallel"); 
library("foreach");
source("getphenogeno.R");
source('BF_test.R')

##############################################
#           QTL and vQTL Mapping             #
##############################################
#               Kaustubh Dhole               #
##############################################
if (!file.exists("QTL"))  dir.create("QTL")
if (!file.exists("vQTL")) dir.create("vQTL")

Data=getphenogeno(file="phenogeno.csv", no_of_env, normalizePheno=TRUE); #### with normalization
phenodata=Data[[1]];
genodata=Data[[2]];
markers = colnames(genodata);
envnames=colnames(phenodata)

for(env in envnames)
{ 
phenocol=as.numeric(phenodata[,env])

perm_f = perm_bf = rep(0,nperm)
for(perm in 1:length(perm_f))
{ 	clust=makeCluster(detectCores()-2);
	registerDoParallel(clust,cores=detectCores()-2); ##Using 2 cores lesser than system capacity
	lod_vector = foreach (marker = markers,.combine=rbind) %dopar%
		{	genocol=genodata[,match(marker,colnames(genodata))];
			y1 = phenocol[which(genocol=="BB")]
			y2 = phenocol[which(genocol=="RR")]
			nmin=min(length(y1),length(y2))
			f=F_test(y1,y2,TRUE);
			bf=BF_test(y1,y2,TRUE);
			c(f,bf)
		}
	stopCluster(clust);
	perm_f[perm] = max(lod_vector[,1]); # Genome-wide Maximum
	perm_bf[perm] = max(lod_vector[,2]);# Genome-wide Maximum
}

### Find the LOD scores and the p-values
	clust=makeCluster(detectCores()-2);
	registerDoParallel(clust,cores=detectCores()-2); ##Using 20 out of 24 cores
lod_vector = foreach (marker = markers,.combine=rbind) %dopar%
{	genocol=genodata[,match(marker,colnames(genodata))]; # only take the column corresponding to the marker
	y1 = phenocol[which(genocol=="BB")]
	y2 = phenocol[which(genocol=="RR")]
	nmin=min(length(y1),length(y2))
	f=F_test(y1,y2,FALSE)
	bf=BF_test(y1,y2,FALSE)
	lod=(nmin)*log(f*(1/(2*nmin-2))+1,10)
	lodb=(nmin)*log(bf*(1/(2*nmin-2))+1,10)
	pvalf = sum(f<perm_f)/nperm 
	pvalbf= sum(bf<perm_bf)/nperm
	#c(f_lod=lod,bf_lod=lodb,mBY=mean(y1,na.rm=T),mRM=mean(y2,na.rm=T),vBY=var(y1,na.rm=T),vRM=var(y2,na.rm=T)))
	c(lod,lodb,pvalf,pvalbf)
}
	stopCluster(clust);
	rownames(lod_vector) = markers
	colnames(lod_vector) = c("LOD","LODBF","P-value(F)","P-value (BF)")

write.csv(lod_vector[,c(1,3)], file.path("QTL",paste(env,".csv",sep="")))
write.csv(lod_vector[,c(2,4)], file.path("vQTL",paste(env,".csv",sep="")))

if (!file.exists(file.path("QTL","Plots")))  dir.create(file.path("QTL","Plots"))
if (!file.exists(file.path("vQTL","Plots")))  dir.create(file.path("vQTL","Plots"))
if(TRUE){
	pdf(file.path("QTL","Plots",paste(env,'.pdf',sep="")))
	plot(lod_vector[,1],ylab='LOD_F',cex=0.5,type="no");lines(lod_vector[,1],lwd=1,col='blue')
	box(lwd=1,col='lightblue')
	dev.off()
	pdf(file.path("vQTL","Plots",paste(env,'.pdf',sep="")))
	plot(lod_vector[,2],ylab='LOD_BF',cex=0.5,type="no");lines(lod_vector[,2],lwd=1,col='blue');
	box(lwd=1,col='lightblue')
	dev.off()
	}
print(paste("Environment",env,"Done"))
}
print("Thank You.")
