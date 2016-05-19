BF_test <- function(y1,y2,maxNULL=FALSE)
{ # BF_test for two groups
  N1 = length(y1); N2 = length(y2);
  N = N1 + N2;
  y=c(y1,y2)
  if(maxNULL==TRUE){rr=sample(1:N,N); y1=y[rr[1:N1]];y2=y[rr[(N1+1):N]];}
  y1med = median(y1,na.rm=T); y2med = median(y2,na.rm=T); # to use Levene's test, replace median by mean
  #y1med = mean(y1,na.rm=T); y2med = mean(y2,na.rm=T);
  z1 = abs(y1-y1med); z2 = abs(y2-y2med);
  zmean = mean(c(z1,z2),na.rm=T);
  z1mean = mean(z1,na.rm=T); z2mean = mean(z2,na.rm=T);
  num = (N-2)*(N1*(z1mean - zmean)^2 + N2*(z2mean - zmean)^2)
  den = (sum((z1 - z1mean)^2,na.rm=T)+sum((z1 - z1mean)^2,na.rm=T))
  W = num/den
  return(W)
}

#### Kaustubh Dhole || Nov 2014 BITS Pilani

F_test<-function(y1,y2,maxNULL=FALSE)
{ # F_test for two groups
  y = c(y1,y2);
  N1 = length(y1); N2 = length(y2);
  N = N1 + N2;
  
  if(maxNULL){rr=sample(1:N,N); y1=y[rr[1:N1]];y2=y[rr[(N1+1):N]];}
  
  y1m = mean(y1,na.rm=T); y2m = mean(y2,na.rm=T); 
  ym = mean(y,na.rm=T)
  
  z1 = sum((y1 - y1m)^2,na.rm=T);
  z2 = sum((y2 - y2m)^2,na.rm=T);
  z = sum((y - ym)^2,na.rm=T);
  
  SS0=sum(z,na.rm=T);
  SS1=sum(z1+z2,na.rm=T);
  
  F_stat = (N*log((SS0/SS1),10))
  return(F_stat)

}