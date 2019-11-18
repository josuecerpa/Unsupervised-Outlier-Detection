#Libraries for parallel programming
library(parallel)
library(doParallel)
library(foreach)

#Cleaning NA's

sanitize_data = function(data){
  t(apply(data, 1, function(x) {
    x[which(is.na(x))] <- mean(x, na.rm=T)
    x
  }))
}

#Shape, magnitude and amplitude indexes for detecting the different
#types of outliers

is.biweight =  function(mtx,cores=detectCores(),partition=ncol(mtx)/cores){
  #Number of clusters in terms of the cores
  cluster = makeCluster(cores)
  #Necessary to register the parallel backendn
  registerDoParallel(cluster)
  #Number of splits for the parallel processing
  splits <- parallel:::splitList(1:ncol(mtx), max(1, as.integer(partition)))
  #Users' median
  MEDIAN=foreach(i=1:length(splits),.combine='c') %dopar% {
    apply(mtx[,splits[[i]]],2,median)
  }
  #Users' mad
  MAD=foreach(i=1:length(splits),.combine='c') %dopar% {
    apply(mtx[,splits[[i]]],2,mad)
  }
  #Computation biweight midcovariance#######################
  ui=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    numerator=t(mtx[,splits[[i]]])-MEDIAN[splits[[i]]]
    denominator=9*qnorm(0.75)*MAD[splits[[i]]]
    result=t(numerator/denominator)
  }
  ai=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    I1=ui[,splits[[i]]]
    I1[(-1<=I1) & (I1<=1)]=1
    I1[(I1< -1) | I1 > 1]=0
    I1
  }
  
  minusmedian=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    #In ordert to compute user-median(user) we need to make the transpose
    #of the user's matrix and the end we transpose again.
    t(t(mtx[,splits[[i]]])-MEDIAN[splits[[i]]])
  }
  
  biwmidcov=foreach(i=1:ncol(mtx),.combine='rbind') %:%
    foreach(j=1:length(splits),.combine='c') %dopar% {
      #Elements numerator############  
      ai.x=matrix(ai[,i],ncol = 1)
      minusmedian.x=matrix(minusmedian[,i],ncol=1)
      ui.x=(1-matrix(ui[,i],ncol=1)^2)
      bi.y=data.frame(ai[,splits[[j]]])
      minusmedian.y=data.frame(minusmedian[,splits[[j]]])
      vi.y=data.frame((1-ui[,splits[[j]]]^2))
      
      #Elements denominator###########
      ui5.x=(1-5*matrix(ui[,i],ncol=1)^2)
      vi5.y=data.frame((1-5*ui[,splits[[j]]]^2))
      
      #biweight midcovariance
      numerator.1=ai.x * minusmedian.x * ui.x^2 * bi.y * minusmedian.y *vi.y^2
      numerator.1.1=nrow(mtx)*apply(numerator.1,2,sum)
      
      denominator.1=sum(ai.x * ui.x * ui5.x)
      denominator.2=bi.y * vi.y * vi5.y 
      denominator.2.1=apply(denominator.2,2,sum)
      
      sbxy=numerator.1.1/(denominator.1 * denominator.2.1)
    }
  
  ai.midvariance=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    I1=abs(ui[,splits[[i]]])
    I1[I1<1]=1
    I1[!(I1=1)]=0
    I1
  }
  
  biwmidvar=foreach(i=1:length(splits),.combine='c') %dopar% {
    numerator1=ai.midvariance[,splits[[i]]]*minusmedian[,splits[[i]]]^2 *(1-ui[,splits[[i]]]^2)^4
    numerator1.1=sqrt(apply(numerator1,2,sum))
    
    denominator1=ai.midvariance[,splits[[i]]] * (1-ui[,splits[[i]]]^2) * (1-5*ui[,splits[[i]]]^2)
    denominator1.1=abs(apply(denominator1,2,sum))
    
    ((sqrt(nrow(mtx))*numerator1.1)/denominator1.1)^2
  }
  
  biwmidcor=foreach(i=1:ncol(mtx),.combine='rbind') %:%
    foreach(j=1:ncol(mtx),.combine='c') %dopar% {
      biwmidcov[i,j]/sqrt(biwmidvar[i]*biwmidvar[j])
    }
  is=foreach(i=1:length(splits),.combine='c') %dopar% {
    unname(abs(colMeans(biwmidcor[,splits[[i]]]-1)))
  }
  stopCluster(cluster)
  return(is)
}
iamp.biweight =  function(mtx,cores=detectCores(),partition=ncol(mtx)/cores){
  #Number of clusters in terms of the cores
  cluster = makeCluster(cores)
  #Necessary to register the parallel backendn
  registerDoParallel(cluster)
  #Number of splits for the parallel processing
  splits <- parallel:::splitList(1:ncol(mtx), max(1, as.integer(partition)))
  #Users' median
  MEDIAN=foreach(i=1:length(splits),.combine='c') %dopar% {
    apply(mtx[,splits[[i]]],2,median)
  }
  #Users' mad
  MAD=foreach(i=1:length(splits),.combine='c') %dopar% {
    apply(mtx[,splits[[i]]],2,mad)
  }
  #Computation biweight midcovariance#######################
  ui=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    numerator=t(mtx[,splits[[i]]])-MEDIAN[splits[[i]]]
    denominator=9*qnorm(0.75)*MAD[splits[[i]]]
    result=t(numerator/denominator)
  }
  ai=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    I1=ui[,splits[[i]]]
    I1[(-1<=I1) & (I1<=1)]=1
    I1[(I1< -1) | I1 > 1]=0
    I1
  }
  minusmedian=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    #In ordert to compute user-median(user) we need to make the transpose
    #of the user's matrix and the end we transpose again.
    t(t(mtx[,splits[[i]]])-MEDIAN[splits[[i]]])
  }
  biwmidcov=foreach(i=1:ncol(mtx),.combine='rbind') %:%
    foreach(j=1:length(splits),.combine='c') %dopar% {
      #Elements numerator############  
      ai.x=matrix(ai[,i],ncol = 1)
      minusmedian.x=matrix(minusmedian[,i],ncol=1)
      ui.x=(1-matrix(ui[,i],ncol=1)^2)
      bi.y=data.frame(ai[,splits[[j]]])
      minusmedian.y=data.frame(minusmedian[,splits[[j]]])
      vi.y=data.frame((1-ui[,splits[[j]]]^2))
      
      #Elements denominator###########
      ui5.x=(1-5*matrix(ui[,i],ncol=1)^2)
      vi5.y=data.frame((1-5*ui[,splits[[j]]]^2))
      
      #biweight midcovariance
      numerator.1=ai.x * minusmedian.x * ui.x^2 * bi.y * minusmedian.y *vi.y^2
      numerator.1.1=nrow(mtx)*apply(numerator.1,2,sum)
      
      denominator.1=sum(ai.x * ui.x * ui5.x)
      denominator.2=bi.y * vi.y * vi5.y 
      denominator.2.1=apply(denominator.2,2,sum)
      
      sbxy=numerator.1.1/(denominator.1 * denominator.2.1)
    }
  ai.midvariance=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    I1=abs(ui[,splits[[i]]])
    I1[I1<1]=1
    I1[!(I1=1)]=0
    I1
  }
  biwmidvar=foreach(i=1:length(splits),.combine='c') %dopar% {
    numerator1=ai.midvariance[,splits[[i]]]*minusmedian[,splits[[i]]]^2 *(1-ui[,splits[[i]]]^2)^4
    numerator1.1=sqrt(apply(numerator1,2,sum))
    
    denominator1=ai.midvariance[,splits[[i]]] * (1-ui[,splits[[i]]]^2) * (1-5*ui[,splits[[i]]]^2)
    denominator1.1=abs(apply(denominator1,2,sum))
    
    ((sqrt(nrow(mtx))*numerator1.1)/denominator1.1)^2
  }
  
  biwmidcor=foreach(i=1:ncol(mtx),.combine='rbind') %:%
    foreach(j=1:ncol(mtx),.combine='c') %dopar% {
      biwmidcov[i,j]/sqrt(biwmidvar[i]*biwmidvar[j])
    }
  is=foreach(i=1:length(splits),.combine='c') %dopar% {
    colMeans(biwmidcor[,splits[[i]]])
  }
  ia=foreach(i=1:length(splits),.combine = 'c')  %dopar% {
    abs(is[splits[[i]]] * (MAD[splits[[i]]]))
  }
  stopCluster(cluster)
  names(ia)=colnames(mtx)
  return(ia)
}
imag.biweight =  function(mtx,cores=detectCores(),partition=ncol(mtx)/cores){
  #Number of clusters in terms of the cores
  cluster = makeCluster(cores)
  #Necessary to register the parallel backendn
  registerDoParallel(cluster)
  #Number of splits for the parallel processing
  splits <- parallel:::splitList(1:ncol(mtx), max(1, as.integer(partition)))
  #Users' median
  MEDIAN=foreach(i=1:length(splits),.combine='c') %dopar% {
    apply(mtx[,splits[[i]]],2,median)
  }
  #Users' mad
  MAD=foreach(i=1:length(splits),.combine='c') %dopar% {
    apply(mtx[,splits[[i]]],2,mad)
  }
  #Computation biweight midcovariance#######################
  ui=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    numerator=t(mtx[,splits[[i]]])-MEDIAN[splits[[i]]]
    denominator=9*qnorm(0.75)*MAD[splits[[i]]]
    result=t(numerator/denominator)
  }
  ai=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    I1=ui[,splits[[i]]]
    I1[(-1<=I1) & (I1<=1)]=1
    I1[(I1< -1) | I1 > 1]=0
    I1
  }
  minusmedian=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    #In ordert to compute user-median(user) we need to make the transpose
    #of the user's matrix and the end we transpose again.
    t(t(mtx[,splits[[i]]])-MEDIAN[splits[[i]]])
  }
  biwmidcov=foreach(i=1:ncol(mtx),.combine='rbind') %:%
    foreach(j=1:length(splits),.combine='c') %dopar% {
      #Elements numerator############  
      ai.x=matrix(ai[,i],ncol = 1)
      minusmedian.x=matrix(minusmedian[,i],ncol=1)
      ui.x=(1-matrix(ui[,i],ncol=1)^2)
      bi.y=data.frame(ai[,splits[[j]]])
      minusmedian.y=data.frame(minusmedian[,splits[[j]]])
      vi.y=data.frame((1-ui[,splits[[j]]]^2))
      
      #Elements denominator###########
      ui5.x=(1-5*matrix(ui[,i],ncol=1)^2)
      vi5.y=data.frame((1-5*ui[,splits[[j]]]^2))
      
      #biweight midcovariance
      numerator.1=ai.x * minusmedian.x * ui.x^2 * bi.y * minusmedian.y *vi.y^2
      numerator.1.1=nrow(mtx)*apply(numerator.1,2,sum)
      
      denominator.1=sum(ai.x * ui.x * ui5.x)
      denominator.2=bi.y * vi.y * vi5.y 
      denominator.2.1=apply(denominator.2,2,sum)
      
      sbxy=numerator.1.1/(denominator.1 * denominator.2.1)
    }
  ai.midvariance=foreach(i=1:length(splits),.combine='cbind') %dopar% {
    I1=abs(ui[,splits[[i]]])
    I1[I1<1]=1
    I1[!(I1=1)]=0
    I1
  }
  biwmidvar=foreach(i=1:length(splits),.combine='c') %dopar% {
    numerator1=ai.midvariance[,splits[[i]]]*minusmedian[,splits[[i]]]^2 *(1-ui[,splits[[i]]]^2)^4
    numerator1.1=sqrt(apply(numerator1,2,sum))
    
    denominator1=ai.midvariance[,splits[[i]]] * (1-ui[,splits[[i]]]^2) * (1-5*ui[,splits[[i]]]^2)
    denominator1.1=abs(apply(denominator1,2,sum))
    
    ((sqrt(nrow(mtx))*numerator1.1)/denominator1.1)^2
  }
  
  biwmidcor=foreach(i=1:ncol(mtx),.combine='rbind') %:%
    foreach(j=1:ncol(mtx),.combine='c') %dopar% {
      biwmidcov[i,j]/sqrt(biwmidvar[i]*biwmidvar[j])
    }
  #computes the index
  imag=foreach(i=1:length(splits),.combine = 'c') %dopar%{
    #median's users
    median=MEDIAN[splits[[i]]]
    #transpose of biweight correlation
    biweight=t(biwmidcor[,splits[[i]]])
    #computation of the MAD
    MADs=MAD[splits[[i]]]
    #biweight(x,x_j)*MAD(X)
    biwepermads=data.frame(t(biweight*MADs))
    #biweight(x,x_j)*MAD(x)*median(x_j)
    biw.mads.medianj=t(biwepermads*matrix(MEDIAN,ncol=1))
    #median(x)-kendall(x,x_j)*MAD(x)*median(x_j)
    index=t(median-biw.mads.medianj)
    #absolute value of the index
    fresults=abs(colMeans(index))
  }
  stopCluster(cluster)
  names(imag)=colnames(mtx)
  return(imag)
}

#Final function

#Note:-In the function getOutliers, each column is a multivariate data point.
#     -Matrices with strings are not allowed
getOutliers=function(data){
  #Different indexes
  shape.index=is.biweight(mtx=data)
  amplitude.index=iamp.biweight(mtx=data)
  magnitude.index=imag.biweight(mtx=data)
  
  #Outliers positions
  out.shape=unname(which(shape.index>boxplot.stats(shape.index)$stats[5]))
  out.amp=unname(which(amplitude.index>boxplot.stats(amplitude.index)$stats[5]))
  out.mag=unname(which(magnitude.index>boxplot.stats(magnitude.index)$stats[5]))
  
  return(list(position_shape_outliers=out.shape,
              position_amplitude_outliers=out.amp,
              position_magnitude_outliers=out.mag))
}






