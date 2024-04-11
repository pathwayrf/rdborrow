################### MMRM
MMRM.sim=function(B,mydir){
  
  source(paste(mydir,"myproj1/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj1/R/MMRM.rct.R",sep=""))
  source(paste(mydir,"myproj1/R/EC_(A)IPW_OPT.R",sep=""))
  source(paste(mydir,"myproj1/R/EC_(A)IPW_OPT_boot.R",sep=""))

  # param grid
  size=c(220)
  setting=c(1,2,3,4,5)
  true.effect=seq(0,1.5,by=0.05)
  parameters<- expand.grid(size = size, setting = setting, true.effect=true.effect)

  dumm=function(size,setting,true.effect){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=size,
                     power.cal=T,
                     tau=true.effect)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    ans=MMRM.rct(df=mysimXY,setting=setting)
    return(c(setting,true.effect,ans[[2]]))
  }
  results=mapply(dumm, size=parameters$size, setting=parameters$setting, true.effect=parameters$true.effect, SIMPLIFY = FALSE)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),
                    make.default = FALSE,
                    packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','lme4','lmerTest','multcomp','cmdstanr','posterior','nnet'))




ids=batchMap(fun = MMRM.sim, B=rep(1,3000), mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)

getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)

getErrorMessages(reg = reg)

clearRegistry(reg = reg)


MMRM.results=reduceResults(rbind, reg=reg)

# 16 min


###################### EC-IPW-OPT

IPW.sim=function(B,mydir){
   
  source(paste(mydir,"myproj1/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj1/R/MMRM.rct.R",sep=""))
  source(paste(mydir,"myproj1/R/EC_(A)IPW_OPT.R",sep=""))
  source(paste(mydir,"myproj1/R/EC_(A)IPW_OPT_boot.R",sep=""))

  # param grid
  size=c(220)
  setting=c(1,2,3,4,5)
  true.effect=seq(0,1.5,by=0.05)
  parameters<- expand.grid(size = size, setting = setting, true.effect=true.effect)

  dumm=function(size,setting,true.effect){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=size,
                     power.cal=T,
                     tau=true.effect)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    ans=EC_IPW_OPT(df=mysimXY,wt=0,method="EC-IPW-OPT",setting=setting,Bootstrap=F,optimal.weight=T)
    return(c(setting,true.effect,ans[[1]][2],ans[[2]][2],ans[[1]][2]-qnorm(0.975)*ans[[2]][2],ans[[1]][2]+qnorm(0.975)*ans[[2]][2],ans[[3]]))
  }
  results=mapply(dumm, size=parameters$size, setting=parameters$setting, true.effect=parameters$true.effect, SIMPLIFY = FALSE)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),
                    make.default = FALSE,
                    packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','lme4','lmerTest','multcomp','cmdstanr','posterior','nnet'))

ids=batchMap(fun =IPW.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)

getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)
getErrorMessages(reg = reg)
clearRegistry(reg = reg)
IPW.results=reduceResults(rbind, reg=reg)

# 15 min



###################### EC-AIPW-OPT

AIPW.sim=function(B,mydir){
  
  source(paste(mydir,"myproj1/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj1/R/MMRM.rct.R",sep=""))
  source(paste(mydir,"myproj1/R/EC_(A)IPW_OPT.R",sep=""))
  source(paste(mydir,"myproj1/R/EC_(A)IPW_OPT_boot.R",sep=""))

  # param grid
  size=c(220)
  setting=c(1,2,3,4,5)
  true.effect=seq(0,1.5,by=0.05)
  parameters<- expand.grid(size = size, setting = setting, true.effect=true.effect)

  dumm=function(size,setting,true.effect){
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=size,
                     power.cal=T,
                     tau=true.effect)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    ans=EC_IPW_OPT(df=mysimXY,wt=0,method="EC-AIPW-OPT",setting=setting,Bootstrap=F,optimal.weight=T)
    return(c(setting,true.effect,ans[[1]][2],ans[[2]][2],ans[[1]][2]-qnorm(0.975)*ans[[2]][2],ans[[1]][2]+qnorm(0.975)*ans[[2]][2],ans[[3]]))
  }
  results=mapply(dumm, size=parameters$size, setting=parameters$setting, true.effect=parameters$true.effect, SIMPLIFY = FALSE)
  do.call(rbind,results)
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),
                    make.default = FALSE,
                    packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','lme4','lmerTest','multcomp','cmdstanr','posterior','nnet'))

ids=batchMap(fun =AIPW.sim, B=rep(1,3000),mydir=mydir, reg = reg)
submitJobs(ids,reg = reg)

getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)
getErrorMessages(reg = reg)
clearRegistry(reg = reg)
AIPW.results=reduceResults(rbind, reg=reg)
# 15 min



################################## Power prior
# PP.sim=function(B,mydir){
# 
#   source(paste(mydir,"myproj1/R/cal.simXY.R",sep=""))
#   source(paste(mydir,"myproj1/R/PP.R",sep=""))
#   source(paste(mydir,"myproj1/R/CPP.R",sep=""))
# 
# 
#   check_cmdstan_toolchain(fix = TRUE, quiet = FALSE)
#   check_cmdstan_toolchain()
#   cmdstanr::set_cmdstan_path(paste(mydir,"myproj1/cmdstan-2.32.1",sep=""))
# 
#   # param grid
#   size=c(220)
#   setting=c(1,2,3,4,5)
#   true.effect=seq(0,1.5,by=0.05)
#   parameters<- expand.grid(size = size, setting = setting, true.effect=true.effect)
# 
#   dumm=function(size,setting,true.effect){
#     result=cal.simXY(path=mydir,
#                      setting=setting,
#                      size=size,
#                      power.cal=T,
#                      tau=true.effect)
#     list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
#     ans=PP(path=paste(mydir,"myproj1/R/",sep=""),df=mysimXY,setting=setting)
#     return(c(setting,true.effect,ans[[2]],ans[[3]][1]))
#   }
#   results=mapply(dumm, size=parameters$size, setting=parameters$setting, true.effect=parameters$true.effect, SIMPLIFY = FALSE)
#   do.call(rbind,results)
# 
# }
# reg <- makeRegistry(file.dir=NA,work.dir = getwd(),
#                     make.default = FALSE,
#                     packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
#                                'ggsci','gridExtra','knitr','scales',
#                                'boot','nlme','lme4','multcomp','cmdstanr','posterior','utils','nnet'))
# 
# ids=batchMap(fun =PP.sim, B=rep(1,1000),mydir=mydir, reg = reg)



PP.sim=function(B,mydir,size,setting,true.effect){

  source(paste(mydir,"myproj1/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj1/R/PP.R",sep=""))
  source(paste(mydir,"myproj1/R/CPP.R",sep=""))


  check_cmdstan_toolchain(fix = TRUE, quiet = FALSE)
  check_cmdstan_toolchain()
  cmdstanr::set_cmdstan_path(paste(mydir,"myproj1/cmdstan-2.32.1",sep=""))

    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=size,
                     power.cal=T,
                     tau=true.effect)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    ans=PP(path=paste(mydir,"myproj1/R/",sep=""),df=mysimXY,setting=setting)

    return(c(setting,true.effect,ans[[2]],ans[[3]][1]))
}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),
                         make.default = FALSE,
                         packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
                                    'ggsci','gridExtra','knitr','scales',
                                    'boot','nlme','lme4','multcomp','cmdstanr','posterior','utils','nnet'))
# param grid
B=rep(1,1000)
size=c(220)
setting=c(1,2,3,4,5)
true.effect=seq(0,1.5,by=0.05)
parameters<- expand.grid(B=B,size = size, setting = setting, true.effect=true.effect)

ids=batchMap(fun = PP.sim,B=parameters$B, mydir=mydir, size=parameters$size, setting=parameters$setting, true.effect=parameters$true.effect, reg = reg)

submitJobs(ids,reg = reg)

getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)
getErrorMessages(reg = reg)
getJobTable(reg = reg)

PP.results=reduceResults(rbind, reg=reg)


################################## Commensurate prior
CPP.sim=function(B,mydir,size,setting,true.effect){
   
  source(paste(mydir,"myproj1/R/cal.simXY.R",sep=""))
  source(paste(mydir,"myproj1/R/PP.R",sep=""))
  source(paste(mydir,"myproj1/R/CPP.R",sep=""))
  
  
  check_cmdstan_toolchain(fix = TRUE, quiet = FALSE)
  check_cmdstan_toolchain()
  cmdstanr::set_cmdstan_path(paste(mydir,"myproj1/cmdstan-2.32.1",sep=""))
  
    result=cal.simXY(path=mydir,
                     setting=setting,
                     size=size,
                     power.cal=T,
                     tau=true.effect)
    list2env(setNames(result,c("mysimXY","ATE")), envir = .GlobalEnv)
    ans=CPP(path=paste(mydir,"myproj1/R/",sep=""),df=mysimXY,setting=setting)
    return(c(setting,true.effect,ans[[2]],ans[[3]][1]))


}


reg <- makeRegistry(file.dir=NA,work.dir = getwd(),
                    make.default = FALSE,
                    packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','lme4','multcomp','cmdstanr','posterior','utils','nnet'))


# param grid
B=rep(1,1000)
size=c(220)
setting=c(1,2,3,4,5)
true.effect=seq(0,1.5,by=0.05)
parameters<- expand.grid(B=B,size = size, setting = setting, true.effect=true.effect)

ids=batchMap(fun = CPP.sim,B=parameters$B,mydir=mydir, size=parameters$size, setting=parameters$setting, true.effect=parameters$true.effect, reg = reg)

submitJobs(ids,reg = reg)
getStatus(reg = reg)
start_time <- Sys.time()
waitForJobs(reg = reg)
end_time <- Sys.time()
print(end_time - start_time)
killJobs(reg=reg)
clearRegistry(reg = reg)
getErrorMessages(reg = reg)

CPP.results=reduceResults(rbind, reg=reg)



sim.results=as.data.frame(rbind(cbind("Within Trial MMRM",MMRM.results,NA),
                                cbind("EC-IPW-OPT",IPW.results), 
                                cbind("EC-AIPW-OPT",AIPW.results),
                                cbind("Power Prior",PP.results),
                                cbind("Commensurate Prior",CPP.results)))
names(sim.results)=c("Method","Setting","ATE","Estimate","SE","Lower","Upper", "Weight")

sim.results=sim.results%>%
  mutate_at(vars("ATE","Estimate","SE","Lower","Upper", "Weight"), as.numeric)


save(sim.results, file =paste(mydir,"myproj1/output/sim.results.RData",sep=""))



