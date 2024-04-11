

 
SCMboot<-function(data, 
                  indices,
                  setting,
                  lambda){
   
  df_b=data[indices,]
    
  n0=sum(df_b$S==1 & df_b$A==0) # number of rct control subjects
  m=sum(df_b$S==0) # number of external control subjects
  
  
  # create data matrices: attributes by row and subject by column
  if(setting %in% c(1,2,3,4)){
    X10 = t(as.matrix(df_b%>%filter(S==1&A==0)%>%dplyr::select(X1:X5,Y1:Y4)))
    X00 = t(as.matrix(df_b%>%filter(S==0)%>%dplyr::select(X1:X5,Y1:Y4)))
    
  }else{
    X10 = t(as.matrix(df_b%>%filter(S==1&A==0)%>%dplyr::select(X1:X5,Y1:Y6)))
    X00 = t(as.matrix(df_b%>%filter(S==0)%>%dplyr::select(X1:X5,Y1:Y6)))
    
  }
  
  # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
  #n_cores <- detectCores()
  ans=lapply(1:n0, subject_SC, setting=setting, X10=X10,X00=X00,lambda=lambda)
   
  # estimated outcome matrix: n0*T (combining by rows)
  y.est.mat=do.call(rbind,lapply(ans, function(x) as.vector(x[[2]])))

  # Aggregate group level synthetic control estiamte
  if(setting %in% c(1,2,3,4)){
    Y.trt=filter(df_b, S==1&A==1)%>%
      summarise(Y3obs=mean(Y3),Y4obs=mean(Y4))%>%
      dplyr::select(Y3obs, Y4obs)%>%unlist()
    tau=Y.trt-colMeans(y.est.mat)
  }else{
    Y.trt=filter(df_b, S==1&A==1)%>%
      summarise(Y5obs=mean(Y5),Y6obs=mean(Y6))%>%
      dplyr::select(Y5obs, Y6obs)%>%unlist()
    tau=Y.trt-colMeans(y.est.mat)
  }
   
  return(tau)
  
}


 

