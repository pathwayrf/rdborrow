### RBMI
rbmi.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/Update R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/Update R/RBMI.R",sep=""))
   
  
  
  
  # setting 1
  result=cal.simXY(path=mydir,setting=1,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=RBMI(df=mysimXY,setting=1)
  res1=cbind(1,"RBMI",1:2, do.call(rbind,ans))
  
  # ans=RBMI(df=mysimXY,setting=1, mis=TRUE)
  # res1_1=cbind(1.1,"RBMI",1:2, do.call(rbind,ans))
    
  # setting 2
  result=cal.simXY(path=mydir,setting=2,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=RBMI(df=mysimXY,setting=2)
  res2=cbind(2,"RBMI",1:2, do.call(rbind,ans))
  
  # ans=RBMI(df=mysimXY,setting=2, mis=TRUE)
  # res2_1=cbind(2.1,"RBMI",1:2, do.call(rbind,ans))
  # 
  # setting 3: study bias
  result=cal.simXY(path=mydir,setting=3,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=RBMI(df=mysimXY,setting=3)
  res3=cbind(3,"RBMI",1:2, do.call(rbind,ans))
  
  # ans=RBMI(df=mysimXY,setting=3, mis=TRUE)
  # res3_1=cbind(3.1,"RBMI",1:2, do.call(rbind,ans))
  # 
  
  # setting 4
  result=cal.simXY(path=mydir,setting=4,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=RBMI(df=mysimXY,setting=4)
  res4=cbind(4,"RBMI",1:2, do.call(rbind,ans))
  
  # ans=RBMI(df=mysimXY,setting=4, mis=TRUE)
  # res4_1=cbind(4.1,"RBMI",1:2, do.call(rbind,ans))
  # 
  # setting 5
  result=cal.simXY(path=mydir,setting=5,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=RBMI(df=mysimXY,setting=5)
  res5=cbind(5,"RBMI",1:2, do.call(rbind,ans))
  
  # ans=RBMI(df=mysimXY,setting=5, mis=TRUE)
  # res5_1=cbind(5.1,"RBMI",1:2, do.call(rbind,ans))
  # 
  
  rbind(res1, res2, res3, res4, res5)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp', "RhpcBLASctl","rbmi","RefBasedMI","mice"))
ids=batchMap(fun = rbmi.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)
 
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
getStatus(reg = reg)
getErrorMessages(reg = reg)

rbmi.results=reduceResults(rbind, reg=reg)

killJobs(reg=reg)
clearRegistry(reg = reg)

 

 

 




### DID OR
did.or.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/Update R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/Update R/DiD.R",sep=""))
  source(paste(mydir,"myproj2/Update R/DiDboot.R",sep=""))
  
  # setting 1
  result=cal.simXY(path=mydir,setting=1,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-OR",
          out.model="Y~X1+X2+X3+X4+X2*X4+X5+S*time+A*time")
  res1=cbind(1,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-OR",
          out.model="Y~X1+X2+X3+X4+X5+S*time+A*time")
  res1_1=cbind(1.1,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  # setting 2
  result=cal.simXY(path=mydir,setting=2,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time")
  res2=cbind(2,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X5+S*time+A*time")
  res2_1=cbind(2.1,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  # setting 2
  result=cal.simXY(path=mydir,setting=3,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time")
  res3=cbind(3,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X5+S*time+A*time")
  res3_1=cbind(3.1,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  
  
  # setting 4
  result=cal.simXY(path=mydir,setting=4,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time")
  res4=cbind(4,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X5+S*time+A*time")
  res4_1=cbind(4.1,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  # setting 5
  result=cal.simXY(path=mydir,setting=5,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time")
  res5=cbind(5,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-OR",
          out.model="Y~X2+X3+X4+X5+S*time+A*time")
  res5_1=cbind(5.1,"DID-EC-OR",1:2, do.call(rbind,ans))
  
  
  rbind(res1, res1_1, res2, res2_1, res3, res3_1, res4, res4_1,res5, res5_1)
  
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp',  "RhpcBLASctl","rbmi"))
ids=batchMap(fun = did.or.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)

start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
getStatus(reg = reg)
getErrorMessages(reg = reg)

did.or.results=reduceResults(rbind, reg=reg)

killJobs(reg=reg)
clearRegistry(reg = reg)

 
 

### DID IPW
did.ipw.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/Update R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/Update R/DiD.R",sep=""))
  source(paste(mydir,"myproj2/Update R/DiDboot.R",sep=""))
  
  # setting 1
  result=cal.simXY(path=mydir,setting=1,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-IPW",
          ps.model="S~X1+X2+X3+X4+X2*X4+X5",
          piA.model="A~X1+X2+X3+X4+X5")
  res1=cbind(1,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-IPW",
          ps.model="S~X1+X2+X3+X4+X5",
          piA.model="A~X1+X2+X3+X4+X5")
  res1_2=cbind(1.2,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  # setting 2
  result=cal.simXY(path=mydir,setting=2,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res2=cbind(2,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res2_2=cbind(2.2,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  # setting 3
  result=cal.simXY(path=mydir,setting=3,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res3=cbind(3,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res3_2=cbind(3.2,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  
  
  # setting 4
  result=cal.simXY(path=mydir,setting=4,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res4=cbind(4,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res4_2=cbind(4.2,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  # setting 5
  result=cal.simXY(path=mydir,setting=5,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res5=cbind(5,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-IPW",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res5_2=cbind(5.2,"DID-EC-IPW",1:2, do.call(rbind,ans))
  
  
  rbind(res1, res1_2, res2, res2_2, res3, res3_2, res4, res4_2, res5, res5_2)
  

}

reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp', "RhpcBLASctl","rbmi"))
ids=batchMap(fun = did.ipw.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)

start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
getStatus(reg = reg)
getErrorMessages(reg = reg)

did.ipw.results=reduceResults(rbind, reg=reg)

killJobs(reg=reg)
clearRegistry(reg = reg)




### DID AIPW
did.aipw.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/Update R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/Update R/DiD.R",sep=""))
  source(paste(mydir,"myproj2/Update R/DiDboot.R",sep=""))
  
  # setting 1
  result=cal.simXY(path=mydir,setting=1,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-AIPW",
          out.model="Y~X1+X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X1+X2+X3+X4+X2*X4+X5",
          piA.model="A~X1+X2+X3+X4+X5")
  res1=cbind(1,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-AIPW",
          out.model="Y~X1+X2+X3+X4+X5+S*time+A*time",
          ps.model="S~X1+X2+X3+X4+X2*X4+X5",
          piA.model="A~X1+X2+X3+X4+X5")
  res1_1=cbind(1.1,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=1,method="DID-EC-AIPW",
          out.model="Y~X1+X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X1+X2+X3+X4+X5",
          piA.model="A~X1+X2+X3+X4+X5")
  res1_2=cbind(1.2,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  
  # setting 2
  result=cal.simXY(path=mydir,setting=2,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res2=cbind(2,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res2_1=cbind(2.1,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=2,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res2_2=cbind(2.2,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  # setting 3
  result=cal.simXY(path=mydir,setting=3,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res3=cbind(3,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res3_1=cbind(3.1,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=3,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res3_2=cbind(3.2,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  # setting 4
  result=cal.simXY(path=mydir,setting=4,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res4=cbind(4,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res4_1=cbind(4.1,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=4,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res4_2=cbind(4.2,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  # setting 5
  result=cal.simXY(path=mydir,setting=5,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res5=cbind(5,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X2*X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res5_1=cbind(5.1,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  ans=DiD(df=mysimXY,setting=5,method="DID-EC-AIPW",
          out.model="Y~X2+X3+X4+X2*X4+X5+S*time+A*time",
          ps.model="S~X2+X3+X4+X5",
          piA.model="A~X2+X3+X4+X5")
  res5_2=cbind(5.2,"DID-EC-AIPW",1:2, do.call(rbind,ans))
  
  rbind(res1, res1_1, res1_2, res2, res2_1,res2_2, res3, res3_1,res3_2, res4, res4_1, res4_2, res5, res5_1, res5_2)
  
  
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp', "RhpcBLASctl","rbmi"))
ids=batchMap(fun = did.aipw.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)

start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
getStatus(reg = reg)
getErrorMessages(reg = reg)

did.aipw.results=reduceResults(rbind, reg=reg)

killJobs(reg=reg)
clearRegistry(reg = reg)

 
## SCM
SCM.sim=function(setting,B,mydir){
  
  # share package folder
  #.libPaths("/home/pangh3/R/focal/x86_64/4.2.1")
  
  source(paste(mydir,"myproj2/Update R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/Update R/SCM.R",sep=""))
  source(paste(mydir,"myproj2/Update R/SCMboot.R",sep=""))
  
   
  result=cal.simXY(path=mydir,setting=setting,size=220)
  list2env(setNames(result,c("mysimXY")), envir = .GlobalEnv)
  #table(mysimXY$S)
  ans=SCM(df=mysimXY,setting=setting)
  cbind(setting,"SCM",1:2, do.call(rbind,ans))
   
  
}

reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','parallel','nlme','multcomp','CVXR'))

# share head node liberary with workers, put template.tmpl in working directory
reg$cluster.functions=makeClusterFunctionsLSF(template = "template")

setting=c(1,2,3,4,5) 
B=rep(1,3000)  
parameters<- expand.grid(setting = setting,B=B)
ids=batchMap(fun = SCM.sim, B=parameters$B, setting=parameters$setting,mydir=mydir, reg = reg)
submitJobs(ids,reg = reg) #, resources =list(walltime = 10000, memory = 512)


start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
getStatus(reg = reg)
getErrorMessages(reg = reg)

scm.results=reduceResults(rbind, reg=reg)

killJobs(reg=reg)
clearRegistry(reg = reg)



# place-holder for final results, only run this once
# scm.results=NULL


# To simulate 3000: each for loop simulate 500 for 2 times= 1000. 
# To avoid stuck in one submission, repeat the following loop 3 successfuly times
# for(i in 1:2){
#   
#   reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
#                       make.default = FALSE,
#                       packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
#                                  'ggsci','gridExtra','knitr','scales',
#                                  'boot','parallel','nlme','multcomp','CVXR'),
#                       template = "myproj2/Update R/lsf.tmpl")
#   
#   
#   setting=c(1,2,3,4,5) 
#   B=rep(1,500)  
#   parameters<- expand.grid(setting = setting,B=B)
#   ids=batchMap(fun = SCM.sim, B=parameters$B, setting=parameters$setting,mydir=mydir, reg = reg)
#   submitJobs(ids,reg = reg, resources =list(walltime = 10000, memory = 512))
#   
#   
#   start_time <- Sys.time()
#   waitForJobs(reg = reg)
#   end_time <- Sys.time()
#   # when all is done, print out status
#   getStatus(reg = reg)
#   # if any error message?
#   getErrorMessages(reg = reg)
#   # time?
#   print(end_time - start_time)
#   
#   # get simulation reults
#   tmp=reduceResults(rbind, reg=reg)
#   
#   # add to final results
#   scm.results=rbind(scm.results, tmp)
#   
#   # kill and clear
#   killJobs(reg=reg) 
#   clearRegistry(reg = reg)
# }




sim.results=as.data.frame(rbind(rbmi.results,
                                did.or.results,
                                did.ipw.results,
                                did.aipw.results, 
                                scm.results))
names(sim.results)=c("Setting","Method","Time","Estimate","SE","Lower","Upper")

sim.results=sim.results%>%
  mutate_at(vars("Setting","Estimate","SE","Time","Lower","Upper"), as.numeric)

save(sim.results, file =paste(mydir,"myproj2/output/update.sim.results.RData",sep=""))


