
### RBMI
rbmi.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/R/cal.simXY.R", sep=""))
  source(paste(mydir,"myproj2/R/RBMI.R", sep=""))
  source(paste(mydir,"myproj2/R/DiD.R", sep=""))
  source(paste(mydir,"myproj2/R/DiDboot.R", sep=""))
   
  # param grid
  setting=c(1,2,3,4,5,6)
  
  dumm=function(setting){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=300)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    
    mi=RBMI(df=mysimXY,setting=setting)
    mi=do.call(rbind,mi)
    mi=cbind(setting,ATE[3:4],"RBMI",c(3,4),mi)
   
    return(mi)
  }
  results=lapply(c(1,2,3,4,5,6), dumm)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp', "RhpcBLASctl","rbmi"))
ids=batchMap(fun = rbmi.sim, B=rep(1,3),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)
getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)

getErrorMessages(reg = reg)

rbmi.results=reduceResults(rbind, reg=reg)

# ~1 hr

### DID OR
did.or.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/R/RBMI.R",sep=""))
  source(paste(mydir,"myproj2/R/DiD.R",sep=""))
  source(paste(mydir,"myproj2/R/DiDboot.R",sep=""))
  
  # param grid
  setting=1:6
  
  dumm=function(setting){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=300)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    
    did=DiD(df=mysimXY,setting=setting,
             method="DID-EC-OR")
    did=do.call(rbind,did)
    did=cbind(setting,ATE[3:4],"DID-EC-OR",c(3,4),did)
    
    return(did)
  }
  results=lapply(1:6, dumm)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp',  "RhpcBLASctl","rbmi"))
ids=batchMap(fun = did.or.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)
getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)

getErrorMessages(reg = reg)

did.or.results=reduceResults(rbind, reg=reg)

# >1 hr

### DID IPW
did.or.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/R/RBMI.R",sep=""))
  source(paste(mydir,"myproj2/R/DiD.R",sep=""))
  source(paste(mydir,"myproj2/R/DiDboot.R",sep=""))
  
  # param grid
  setting=1:6
  
  dumm=function(setting){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=300)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    
    did=DiD(df=mysimXY,setting=setting,
            method="DID-EC-IPW")
    did=do.call(rbind,did)
    did=cbind(setting,ATE[3:4],"DID-EC-IPW",c(3,4),did)
    
    return(did)
  }
  results=lapply(1:6, dumm)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp', "RhpcBLASctl","rbmi"))
ids=batchMap(fun = did.or.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)
getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)

getErrorMessages(reg = reg)

did.ipw.results=reduceResults(rbind, reg=reg)
 
# >1hr
 
### DID AIPW
did.or.sim=function(B,mydir){
  
  source(paste(mydir,"myproj2/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/R/RBMI.R",sep=""))
  source(paste(mydir,"myproj2/R/DiD.R",sep=""))
  source(paste(mydir,"myproj2/R/DiDboot.R",sep=""))
  
  # param grid
  setting=1:6
  
  dumm=function(setting){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=300)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    
    did=DiD(df=mysimXY,setting=setting,
            method="DID-EC-AIPW")
    did=do.call(rbind,did)
    did=cbind(setting,ATE[3:4],"DID-EC-AIPW",c(3,4),did)
    
    return(did)
  }
  results=lapply(1:6, dumm)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                    make.default = FALSE,
                    packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','multcomp', "RhpcBLASctl","rbmi"))
ids=batchMap(fun = did.or.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)
getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)

getErrorMessages(reg = reg)

did.aipw.results=reduceResults(rbind, reg=reg)


# >1hr


## scm
sim=function(setting,B,mydir){
  source(paste(mydir,"myproj2/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj2/R/SCM.R",sep=""))
  source(paste(mydir,"myproj2/R/SCMboot.R",sep=""))
    
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=300)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    scm=SCM(df=mysimXY,setting=setting)
    scm=do.call(rbind,scm)
    scm=cbind(setting,ATE[3:4],"SCM",c(3,4),scm)

}

reg <- makeRegistry(file.dir=NA,work.dir = getwd(),  
                     make.default = FALSE,
                     packages=c('nnet','dplyr', 'tidyr', 'ggplot2','reshape2',
                                'ggsci','gridExtra','knitr','scales',
                                'boot','parallel','nlme','multcomp','CVXR'))

 
setting=c(1,2,3,4,5,6) 
B=rep(1,3000)  
parameters<- expand.grid(setting = setting,B=B)


ids=batchMap(fun = sim, B=parameters$B, setting=parameters$setting,mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)
getStatus(reg = reg)


start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)

getErrorMessages(reg = reg)
scm.results=reduceResults(rbind, reg=reg)

sim.results=as.data.frame(rbind(rbmi.results,did.or.results,did.ipw.results,did.aipw.results, scm.results))
names(sim.results)=c("Setting","ATE","Method","Time","Estimate","SE","Lower","Upper")

sim.results=sim.results%>%
  mutate_at(vars("Setting","ATE","Estimate","SE","Time","Lower","Upper"), as.numeric)
 
save(sim.results, file =paste(mydir,"myproj2/output/sim.results.RData",sep=""))
   



