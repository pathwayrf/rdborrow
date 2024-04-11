

 
SCMboot<-function(data, 
                  indices,
                  setting,
                  lambda){
  #browser()
  df_b=data[indices,]
   
  n0=sum(df_b$S==1 & df_b$A==0) # number of rct control subjects
  m=sum(df_b$S==0) # number of external control subjects
  
  # standardize numeric vars
  columns_to_standardize=c("SMN2_Copy_Number","Age_Enrollment","Y0")
  df_b[columns_to_standardize] <- lapply(df_b[columns_to_standardize], scale)
  
  
  # create data matrices: attributes by row and subject by column
  if(setting==1){
    X10 = t(as.matrix(df_b%>%filter(S==1&A==0)%>%dplyr::select(SMA_Type, SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    X00 = t(as.matrix(df_b%>%filter(S==0)%>%dplyr::select(SMA_Type, SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    
  }else{
    X10 = t(as.matrix(df_b%>%filter(S==1&A==0)%>%dplyr::select(SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    X00 = t(as.matrix(df_b%>%filter(S==0)%>%dplyr::select(SMN2_Copy_Number,Scoliosis, Age_Enrollment,Y0,Y1,Y2,Y3,Y4)))
    
  }
  
  # # For each rct control subject, find the synthetic control weight and estimate for all time points, return is a list: each element is a list (weight vector and estimated outcome vector)
  ans=lapply(1:n0, subject_SC,
             setting=setting,
             X10=X10,
             X00=X00,
             lambda=lambda)


  # estimated outcome matrix: n0*T (combining by rows)
  y.est.mat=do.call(rbind,lapply(ans, function(x) as.vector(x[[2]])))

  # Aggregate group level synthetic control estiamte
  Y.trt=filter(df_b, S==1&A==1)%>%
    summarise(Y3obs=mean(Y3),Y4obs=mean(Y4))%>%
    dplyr::select(Y3obs, Y4obs)%>%unlist()
  Y.trt-colMeans(y.est.mat)

   
  
}


