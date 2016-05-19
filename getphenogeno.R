getphenogeno <- function (file, no_of_env, normalizePheno=FALSE, computeParallel=FALSE)
{ print(paste('No. of phenotypes =',no_of_env))
  crossf=read.csv(file,stringsAsFactors=F);
  phenof=data.matrix(crossf[3:nrow(crossf),1:no_of_env])
  
if (normalizePheno==TRUE)
{ phenof = data.matrix(phenof); # converts a data frame into a matrix (numeric)
  nphenof=phenof-rep(colMeans(phenof,na.rm=T),each=nrow(phenof));
 colVar=apply(phenof,2,var,na.rm=T);
 colSD=sqrt(colVar);
 phenof=nphenof/rep(colSD,each=nrow(phenof));
}

genof=crossf[3:nrow(crossf),(no_of_env+1):ncol(crossf)];
markers=colnames(genof);
print(paste('No. of markers =',length(markers)))
return(list(phenof,genof));
}