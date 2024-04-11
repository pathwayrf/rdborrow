
# Source the toy.simXY() function 
# source(paste(mydir, "myproj1/R/toy.simXY.R",sep=""))
source("EC_primaryanalysis_original/scripts/preload packages.R")
source("EC_primaryanalysis_original/R/toy.simXY.R")
# Simulate toy data 
mysimXY=toy.simXY()
mydir = "EC_primaryanalysis_original"

## To produce Fig. 2. (a) Distribution shift of a measured confounder X between R and E: Individuals with lower levels of a confounder (X) are over-represented in
# external controls PE compared to the trial PR. (b) Distribution shift of unmeasured confounder U between R and E: Individuals with higher levels of
# a U are over-represented in E compared to R. The existence of U would violate the Conditional Ignorability Assumption 3. (c) The result of selection
# bias is that PR(Y(0)) ̸= PE ( (0)), and the within-trial evidence (difference between the yellow and blue dashed lines) is unbiased, however, the external
# control evidence is biased (as shown by the difference between the green and blue dashed lines) due to either X or U.

# X distributions
p1=mysimXY%>%mutate(g = case_when(S==1& A==1 ~ "trial treated",
                                  S==1& A==0 ~ "trial control",
                                  S==0 ~ "external control"))%>%
  ggplot(aes(x = X, y = g, group = g)) +
  geom_density_ridges(aes(fill = g,col=g),alpha=0.5)+
  scale_fill_manual(values = c("trial treated"="darkgoldenrod1","trial control"="deepskyblue3","external control"="darkgreen"))+
  scale_color_manual(values = c("trial treated"="darkgoldenrod1","trial control"="deepskyblue3","external control"="darkgreen"))+
  guides(group="none",shape="none",color="none",fill="none")+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(color='black'),
        legend.position = "bottom",
        text = element_text(size = 4),
        strip.background = element_rect(
          color="grey",size=4
        ),
        strip.text = element_text(colour = 'red',size=3),
        plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(3, "line"),
        plot.caption = element_text(hjust = 0.5))+
  xlab("X")+ylab("Density")+
  annotate("text", x = -5, y = 3.5, label = TeX("$X |S=1,A=1$"),size=2,col="darkgoldenrod1")+
  annotate("text", x = -5, y = 2.5, label = TeX("$X |S=1,A=0$"),size=2,col="deepskyblue3")+
  annotate("text", x = -5, y = 1.5, label = TeX("$X |S=0$"),size=2,col="darkgreen")+
  labs(title = "Distributions of a measured confounder X")

# U distributions
p2=mysimXY%>%mutate(g = case_when(S==1& A==1 ~ "trial treated",
                                  S==1& A==0 ~ "trial control",
                                  S==0 ~ "external control"))%>%
  ggplot(aes(x = U, y = g, group = g)) +
  geom_density_ridges(aes(fill = g,col=g),alpha=0.5)+
  scale_fill_manual(values = c("trial treated"="darkgoldenrod1","trial control"="deepskyblue3","external control"="darkgreen"))+
  scale_color_manual(values = c("trial treated"="darkgoldenrod1","trial control"="deepskyblue3","external control"="darkgreen"))+
  guides(group="none",shape="none",color="none",fill="none")+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(color='black'),
        legend.position = "bottom",
        text = element_text(size = 4),
        strip.background = element_rect(
          color="grey",size=4
        ),
        strip.text = element_text(colour = 'red',size=3),
        plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(3, "line"),
        plot.caption = element_text(hjust = 0.5))+
  xlab("U")+ylab("Density")+
  annotate("text", x = -5, y = 3.5, label = TeX("$U |S=1,A=1$"),size=2,col="darkgoldenrod1")+
  annotate("text", x = -5, y = 2.5, label = TeX("$U |S=1,A=0$"),size=2,col="deepskyblue3")+
  annotate("text", x = -5, y = 1.5, label = TeX("$U |S=0$"),size=2,col="darkgreen")+
  labs(title = "Distributions of an unmeasured confounder U")
 
# Y distributions

# calcualte averages by groups to use in plot
temp=mysimXY%>%mutate(g = case_when(S==1& A==1 ~ "trial treated",
                                    S==1& A==0 ~ "trial control",
                                    S==0 ~ "external control"))%>%
  group_by(g)%>%summarise(ave=mean(Y))%>%ungroup()

p3=mysimXY%>%mutate(g = case_when(S==1& A==1 ~ "trial treated",
                                  S==1& A==0 ~ "trial control",
                                  S==0 ~ "external control"))%>%
  ggplot(aes(x = Y, y = g, group = g)) +
  geom_density_ridges(aes(fill = g,col=g),alpha=0.5)+
  scale_fill_manual(values = c("trial treated"="darkgoldenrod1","trial control"="deepskyblue3","external control"="darkgreen"))+
  scale_color_manual(values = c("trial treated"="darkgoldenrod1","trial control"="deepskyblue3","external control"="darkgreen"))+
  guides(group="none",shape="none",color="none",fill="none")+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(color='black'),
        legend.position = "bottom",
        text = element_text(size = 4),
        strip.background = element_rect(
          color="grey",size=4
        ),
        strip.text = element_text(colour = 'red',size=1),
        plot.title = element_text(hjust = 0.5),
        legend.key.width = unit(3, "line"),
        plot.margin = unit(c(0,0,0,0), "cm"))+
  xlab("Y")+ylab("Density")+
  annotate("text", x = -10, y = 3.5, label = TeX("$Y(1) |S=1,A=1$"),size=2,col="darkgoldenrod1")+
  annotate("text", x = -10, y = 2.5, label = TeX("$Y(0) |S=1,A=0$"),size=2,col="deepskyblue3")+
  annotate("text", x = -10, y = 1.5, label = TeX("$Y(0) |S=0$"),size=2,col="darkgreen")+
  geom_vline(aes(xintercept =ave,col=g), data=temp, linetype="dashed")+
  labs(title = "Distributions of observed/potential outcomes")
 
# combine 3 subfigures
(p1 | p2)/p3 +
  plot_annotation(tag_levels = c('A'))+
  plot_layout(guides='collect') &
  theme(legend.position='bottom')

# save it
# ggsave(paste(mydir, "myproj1/output/Fig1.png",sep=""), width = 10, height = 10, units = "cm",dpi = 1000)



## To simulate across different range of weights from 0 to 1, borrowing strength


sim=function(B,mydir){
  # specify the path to the R functions
  
  source(paste(mydir,"/R/toy.simXY.R",sep=""))
  source(paste(mydir,"/R/EC_(A)IPW_OPT for toy example.R",sep=""))
  
  
    mysimXY=toy.simXY(size=300,tau=0,beta.X=0.2,beta.U=0)
    ans1=EC_IPW_OPT_toy(df=mysimXY,wt=0,method="within trial")
    tmp2=lapply(seq(0,1,by=0.05), EC_IPW_OPT_toy, df=mysimXY,method="EC-IPW-OPT", optimal.weight=FALSE, form_x="X")
    ans2=do.call(rbind, tmp2)
    tmp3=lapply(seq(0,1,by=0.05), EC_IPW_OPT_toy, df=mysimXY,method="EC-AIPW-OPT", optimal.weight=FALSE, form_x="X")
    ans3=do.call(rbind, tmp3)
    case1=rbind(cbind("(a) Good overlap in X, No U", "Within Trial IPTW",matrix(ans1,nrow=1)),
              cbind("(a) Good overlap in X, No U", "EC-IPW",ans2),
              cbind("(a) Good overlap in X, No U","EC-AIPW",ans3))
    
    mysimXY=toy.simXY(size=300,tau=0,beta.X=1,beta.U=0)
    ans1=EC_IPW_OPT_toy(df=mysimXY,wt=0,method="within trial") 
    tmp2=lapply(seq(0,1,by=0.05), EC_IPW_OPT_toy, df=mysimXY,method="EC-IPW-OPT", optimal.weight=FALSE, form_x="X")
    ans2=do.call(rbind, tmp2)
    tmp3=lapply(seq(0,1,by=0.05), EC_IPW_OPT_toy, df=mysimXY,method="EC-AIPW-OPT", optimal.weight=FALSE, form_x="X")
    ans3=do.call(rbind, tmp3)
    case2=rbind(cbind("(b) Non-overlap in X, No U", "Within Trial IPTW",matrix(ans1,nrow=1)),
                cbind("(b) Non-overlap in X, No U", "EC-IPW",ans2),
                cbind("(b) Non-overlap in X, No U","EC-AIPW",ans3))
    
    
    mysimXY=toy.simXY(size=300,tau=0,beta.X=0.2,beta.U=-0.5)
    ans1=EC_IPW_OPT_toy(df=mysimXY,wt=0,method="within trial") 
    tmp2=lapply(seq(0,1,by=0.05), EC_IPW_OPT_toy, df=mysimXY,method="EC-IPW-OPT", optimal.weight=FALSE, form_x="X")
    ans2=do.call(rbind, tmp2)
    tmp3=lapply(seq(0,1,by=0.05), EC_IPW_OPT_toy, df=mysimXY,method="EC-AIPW-OPT", optimal.weight=FALSE, form_x="X")
    ans3=do.call(rbind, tmp3)
    case3=rbind(cbind("(c) Good overlap in X, strong unmeasured confounding U","Within Trial IPTW", matrix(ans1,nrow=1)),
                cbind("(c) Good overlap in X, strong unmeasured confounding U","EC-IPW",ans2),
                cbind("(c) Good overlap in X, strong unmeasured confounding U","EC-AIPW",ans3))
     
    results=rbind(case1,case2,case3)
  
  return(results)

}


reg.toy <- makeRegistry(file.dir=NA,work.dir = getwd(),
                    make.default = FALSE,
                    packages=c('dplyr', 'tidyr', 'ggplot2','reshape2',
                               'ggsci','gridExtra','knitr','scales',
                               'boot','nlme','lme4','lmerTest','multcomp','cmdstanr','posterior','nnet'))




ids=batchMap(fun = sim, B=rep(1,10),mydir=mydir,reg = reg.toy)
submitJobs(ids,reg = reg.toy)

getStatus(reg = reg.toy)
start_time <- Sys.time()
waitForJobs(reg = reg.toy)
end_time <- Sys.time()
print(end_time - start_time)
# killJobs(reg=reg.toy)
# clearRegistry(reg = reg.toy)
# getErrorMessages(reg = reg.toy)
# clearRegistry(reg = reg.toy)

toy.results=reduceResults(rbind, reg=reg.toy)
toy.results=as.data.frame(toy.results)
names(toy.results)=c("balance","method","est","se","wt")
toy.results=toy.results%>%
  mutate_at(vars("est","se","wt"), as.numeric)
# save the simulation results
# save(simple.results, file =paste(mydir,"/myproj1/output/toy.results.RData",sep=""))

## To produce Fig. 3. Bias (A), SE (b) and MSE (C) of EC-IPSW and EC-AIPSW with different weights.

# load if already have the simulation resutls
# load(paste(mydir,"myproj1/output/toy.results.RData",sep=""))
 
# define color and line styles
cols=c("Within Trial IPTW"="red",
       "EC-IPW"="forestgreen",
       "EC-AIPW"="purple"
)
lines=c("Within Trial IPTW"="dotted",
        "EC-IPW"="solid",
        "EC-AIPW"="solid"
)

 
temp=toy.results%>%
  mutate(method=factor(method,levels=c("Within Trial IPTW",
                                       "EC-IPW",
                                       "EC-AIPW" )))%>%
  group_by(balance, method,wt)%>%
  summarise(Bias=abs(mean(est)),
            `Standard Error`=sd(est),
            MSE=Bias^2+`Standard Error`^2)%>%
  ungroup()%>%
  gather(group, measurement, Bias,`Standard Error`,MSE, factor_key=TRUE)

within=filter(temp,method=="Within Trial IPTW")
temp=filter(temp,method!="Within Trial IPTW")


temp%>%
  ggplot(aes())+
  geom_line(aes(x=wt,y=measurement,col=method,linetype=method),size=1)+
  geom_ribbon(aes(x=wt, ymin = 0, ymax = measurement,fill = method), alpha = 0.05) +
  geom_hline(data=within,aes(yintercept=measurement,linetype=method,col=method))+
  scale_color_manual(values =cols)+
  scale_fill_manual(values=cols)+
  scale_linetype_manual(values=lines)+
  guides(fill="none",group="none",shape="none",color=guide_legend(title=""),linetype=guide_legend(title=""))+
  theme_bw() +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        text = element_text(size = 8),
        strip.background = element_rect(
          color="grey",size=20
        ),
        strip.text = element_text(colour = 'red',size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 30, vjust = 0.4, hjust=0.6),
        panel.spacing=unit(0, "lines"))+
  xlab("Weight w")+ylab("")+ggtitle("")+
  ggh4x::facet_grid2(group~balance,scales = "free_y", independent = "y", labeller = label_wrap_gen(multi_line = TRUE))

# ggsave(paste(mydir, "myproj1/output/Fig2.png",sep=""), width = 15, height = 15, units = "cm",dpi = 1000)

   
