library(readr)
require(ggplot2)
require(ggpubr)


aegis<- read_csv("C:/Users/wilso/Box/Proposals/Completed/Aegis/Copy of Resuspension - fomite contamination from activities on floor_Aegis.csv")


#----------QMRA regarding concentration of virus on floor by infection risk by activity/height/seeding material

  #ratio of concentration on floor to adjust concentrations at each height
  
  aegis$ratio<-aegis$`PFU/100cm2`/aegis$Floor
  
  mean<-rep(NA,max(aegis$Scenario))
  sd<-rep(NA,max(aegis$Scenario))
  height<-rep(NA,max(aegis$Scenario))
  floortype<-rep(NA,max(aegis$Scenario))
  activity<-rep(NA,max(aegis$Scenario))
  virussuspend<-rep(NA,max(aegis$Scenario))

  for (i in 1:max(aegis$Scenario)){
    
    mean[i]<-mean(as.numeric(aegis$ratio[aegis$Scenario==i]))
    sd[i]<-sd(as.numeric(aegis$ratio[aegis$Scenario==i]))
              height[i]<-aegis$`Height (cm)`[aegis$Scenario==i][1]
    floortype[i]<-aegis$`Floor type`[aegis$Scenario==i][1]
    activity[i]<-aegis$Activity[aegis$Scenario==i][1]
    virussuspend[i]<-aegis$`Virus suspended in Dust`[aegis$Scenario==i][1]
  }

aegis2<-data.frame(mean,sd,height,floortype,activity,virussuspend)
aegis2$activity[aegis2$activity=="walk"]<-"Walking"
aegis2$activity[aegis2$activity=="vacuum"]<-"Vacuuming"
aegis2$floortype[aegis2$floortype=="carpet"]<-"Carpet"
aegis2$floortype[aegis2$floortype=="Floor"]<-"Hard Flooring"

scenarios<-data.frame(id=c(1:12),
                      activities=rep(rep(c("Vacuuming","Walking"),2),3),
                      floortype=rep(rep(c("Carpet","Carpet","Hard Flooring","Hard Flooring")),3),
                      heights=c(rep("< 30cm",4),rep("55-105 cm",4),rep(">122 cm",4)),
                      ratiodust=NA,ratiotripartite=NA)
require(truncdist)

for (i in 1:12){
  tempdata<-aegis2[aegis2$floortype==scenarios$floortype[i] &
                     aegis2$activity==scenarios$activities[i] & 
                     aegis2$height==scenarios$heights[i] & aegis2$virussuspend=="Dust",]
  
  scenarios$ratiodust[scenarios$floortype==scenarios$floortype[i] &
                        scenarios$activities==scenarios$activities[i]& scenarios$heights==scenarios$heights[i]]<-
    tempdata$mean
  
  tempdata2<-aegis2[aegis2$floortype==scenarios$floortype[i] &
                     aegis2$activity==scenarios$activities[i] & 
                     aegis2$height==scenarios$heights[i] & aegis2$virussuspend=="Tripartite",]
  
  scenarios$ratiotripartite[scenarios$floortype==scenarios$floortype[i] &
                              scenarios$activities==scenarios$activities[i]& scenarios$heights==scenarios$heights[i]]<-
    tempdata2$mean
}

#Table S1
View(aegis2)


#Figure S1
require(ggplot2)
windows()
aegis2$height<-factor(aegis2$height,levels=c("< 30cm", "55-105 cm", ">122 cm"))
ggplot(aegis2)+geom_bar(aes(x=height,y=mean,fill=virussuspend,group=interaction(virussuspend,height)),color="black",stat="identity",alpha=0.2,position="dodge")+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd,x=height,group=interaction(virussuspend,height)),position="dodge",color="black")+
  facet_grid(activity~floortype)+
  scale_fill_discrete(name="Virus Application")+
  scale_x_discrete(name="Height")+
  scale_y_continuous(name=expression("Fraction of Floor Concentration on Surfaces"))+
  theme(axis.text = element_text(size=13),axis.title=element_text(size=13),legend.text = element_text(size=13),
        legend.title = element_text(size=13),strip.text = element_text(size=13))



#-------------- QMRA portion ------------------------------------

# log difference between dust/non-dust

#example with norovirus (noro /cm^2)

concentration.noro<-10^(sample(c(4,5,6,7,8),10000,replace=TRUE))

#MS2 surface to fingerpad from Anderson & Boehm (sd value corrected at revision stage)
TE.SH<-rtrunc(10000,"norm",mean=0.34,sd=0.12,a=0,b=1)

#ASTM tripartite soil load (Abney paper)
TE.HF<-rtrunc(10000,"norm",mean=0.41,sd=0.1098,a=0,b=1)

#AuYeung front partial fingers/5 to max of full front palm w/ fingers
S.H<-runif(10000,0.008,0.25) 

#AuYeung front partial fingers/5
S.F<-runif(10000,0.008,0.012)

#Beamer office paper/2 so it's for a single hand
A.hand<-runif(10000,890/2,1070/2)

risk.reduction.mean<-rep(NA,length(1:12))
risk.reduction.sd<-rep(NA,length(1:12))

for (i in 1:12){
  
  #loading on hands
  hands.baseline<-concentration.noro*scenarios$ratiodust[i]*TE.SH*S.H
  hands.intervention<-(concentration.noro*scenarios$ratiotripartite[i])*TE.SH*S.H
  
  #transfer to face
  dose.baseline<-hands.baseline*TE.HF*S.F*A.hand
  dose.intervention<-hands.intervention*TE.HF*S.F*A.hand
  
  #risks (fractional Poisson from van abel)
  P<-0.72
  u.a<-1106
  
  risk.baseline<-P*(1-exp(-dose.baseline/u.a))
  risk.intervention<-P*(1-exp(-dose.intervention/u.a))
  
  risk.reduction.mean[i]<-mean(log10(risk.baseline/risk.intervention))
  risk.reduction.sd[i]<-sd(log10(risk.baseline/risk.intervention))
  
  #save outputs
  if (i==1){
    allframe<-data.frame(risks=c(risk.baseline,risk.intervention),
                         type=c(rep("Dust",length(risk.baseline)),rep("Tripartite",length(risk.intervention))),
                         concentration.noro,
                          TE.SH=TE.SH,
                          TE.HF=TE.HF,
                          A.hand=A.hand,
                          activities=scenarios$activities[i],
                          heights=scenarios$heights[i],
                          floortype=scenarios$floortype[i])
  }else{
    tempframe<-data.frame(
                          risks=c(risk.baseline,risk.intervention),
                          type=c(rep("Dust",length(risk.baseline)),rep("Tripartite",length(risk.intervention))),
                          concentration.noro,
                          TE.SH=TE.SH,
                          TE.HF=TE.HF,
                          A.hand=A.hand,
                          activities=scenarios$activities[i],
                          heights=scenarios$heights[i],
                          floortype=scenarios$floortype[i])
    allframe<-rbind(allframe,tempframe)
  }
  
}

risk.reduction.frame<-data.frame(risk.reduction.mean,risk.reduction.sd,scenarios$activities,scenarios$heights,scenarios$floortype)
View(risk.reduction.frame)

conclist<-c(10^4,10^5,10^6,10^7,10^8)
heights<-c("< 30cm", ">122 cm", "55-105 cm")

for (j in 1:3){
  mean.dust.vacuum.carpet<-rep(NA,length(conclist))
  mean.dust.vacuum.hardflooring<-rep(NA,length(conclist))
  mean.dust.walking.carpet<-rep(NA,length(conclist))
  mean.dust.walking.hardflooring<-rep(NA,length(conclist))
  mean.tripartite.vacuum.carpet<-rep(NA,length(conclist))
  mean.tripartite.vacuum.hardflooring<-rep(NA,length(conclist))
  mean.tripartite.walking.carpet<-rep(NA,length(conclist))
  mean.tripartite.walking.hardflooring<-rep(NA,length(conclist))
  
  sd.dust.vacuum.carpet<-rep(NA,length(conclist))
  sd.dust.vacuum.hardflooring<-rep(NA,length(conclist))
  sd.dust.walking.carpet<-rep(NA,length(conclist))
  sd.dust.walking.hardflooring<-rep(NA,length(conclist))
  sd.tripartite.vacuum.carpet<-rep(NA,length(conclist))
  sd.tripartite.vacuum.hardflooring<-rep(NA,length(conclist))
  sd.tripartite.walking.carpet<-rep(NA,length(conclist))
  sd.tripartite.walking.hardflooring<-rep(NA,length(conclist))
  
  
  
  
  for (i in 1:length(conclist)){
    
    #dust----------------------------------
    mean.dust.vacuum.carpet[i]<-mean(allframe$risks[allframe$type=="Dust" & 
                                                      allframe$activities=="Vacuuming" & 
                                                      allframe$floortype=="Carpet" &
                                                      concentration.noro==conclist[i]&
                                                      allframe$heights==heights[j]])
    
    sd.dust.vacuum.carpet[i]<-sd(allframe$risks[allframe$type=="Dust" & 
                                                  allframe$activities=="Vacuuming" & 
                                                  allframe$floortype=="Carpet" &
                                                  concentration.noro==conclist[i]&
                                   allframe$heights==heights[j]])
    
    mean.dust.vacuum.hardflooring[i]<-mean(allframe$risks[allframe$type=="Dust" & 
                                                            allframe$activities=="Vacuuming" & 
                                                            allframe$floortype=="Hard Flooring" &
                                                            concentration.noro==conclist[i]&
                                             allframe$heights==heights[j]])
    sd.dust.vacuum.hardflooring[i]<-sd(allframe$risks[allframe$type=="Dust" & 
                                                        allframe$activities=="Vacuuming" & 
                                                        allframe$floortype=="Hard Flooring" &
                                                        concentration.noro==conclist[i]&
                                         allframe$heights==heights[j]])
    
    mean.dust.walking.carpet[i]<-mean(allframe$risks[allframe$type=="Dust" & 
                                                       allframe$activities=="Walking" & 
                                                       allframe$floortype=="Carpet" &
                                                       concentration.noro==conclist[i]&
                                        allframe$heights==heights[j]])
    sd.dust.walking.carpet[i]<-sd(allframe$risks[allframe$type=="Dust" & 
                                                   allframe$activities=="Walking" & 
                                                   allframe$floortype=="Carpet" &
                                                   concentration.noro==conclist[i]&
                                    allframe$heights==heights[j]])
    
    mean.dust.walking.hardflooring[i]<-mean(allframe$risks[allframe$type=="Dust" & 
                                                             allframe$activities=="Walking" & 
                                                             allframe$floortype=="Hard Flooring" &
                                                             concentration.noro==conclist[i]&
                                              allframe$heights==heights[j]])
    sd.dust.walking.hardflooring[i]<-sd(allframe$risks[allframe$type=="Dust" & 
                                                         allframe$activities=="Walking" & 
                                                         allframe$floortype=="Hard Flooring" &
                                                         concentration.noro==conclist[i]&
                                          allframe$heights==heights[j]])
    #-----------tripartite
    mean.tripartite.vacuum.carpet[i]<-mean(allframe$risks[allframe$type=="Tripartite" & 
                                                            allframe$activities=="Vacuuming" & 
                                                            allframe$floortype=="Carpet" &
                                                            concentration.noro==conclist[i]&
                                             allframe$heights==heights[j]])
    sd.tripartite.vacuum.carpet[i]<-sd(allframe$risks[allframe$type=="Tripartite" & 
                                                        allframe$activities=="Vacuuming" & 
                                                        allframe$floortype=="Carpet" &
                                                        concentration.noro==conclist[i]&
                                         allframe$heights==heights[j]])
    
    mean.tripartite.vacuum.hardflooring[i]<-mean(allframe$risks[allframe$type=="Tripartite" & 
                                                                  allframe$activities=="Vacuuming" & 
                                                                  allframe$floortype=="Hard Flooring" &
                                                                  concentration.noro==conclist[i]&
                                                   allframe$heights==heights[j]])
    sd.tripartite.vacuum.hardflooring[i]<-sd(allframe$risks[allframe$type=="Tripartite" & 
                                                              allframe$activities=="Vacuuming" & 
                                                              allframe$floortype=="Hard Flooring" &
                                                              concentration.noro==conclist[i]&
                                               allframe$heights==heights[j]])
    
    mean.tripartite.walking.carpet[i]<-mean(allframe$risks[allframe$type=="Tripartite" & 
                                                             allframe$activities=="Walking" & 
                                                             allframe$floortype=="Carpet" &
                                                             concentration.noro==conclist[i]&
                                              allframe$heights==heights[j]])
    sd.tripartite.walking.carpet[i]<-sd(allframe$risks[allframe$type=="Tripartite" & 
                                                         allframe$activities=="Walking" & 
                                                         allframe$floortype=="Carpet" &
                                                         concentration.noro==conclist[i]&
                                          allframe$heights==heights[j]])
    
    mean.tripartite.walking.hardflooring[i]<-mean(allframe$risks[allframe$type=="Tripartite" & 
                                                                   allframe$activities=="Walking" & 
                                                                   allframe$floortype=="Hard Flooring" &
                                                                   concentration.noro==conclist[i]&
                                                    allframe$heights==heights[j]])
    sd.tripartite.walking.hardflooring[i]<-sd(allframe$risks[allframe$type=="Tripartite" & 
                                                               allframe$activities=="Walking" & 
                                                               allframe$floortype=="Hard Flooring" &
                                                               concentration.noro==conclist[i]&
                                                allframe$heights==heights[j]])
    
  }
  
  if (j==1){
    plot.frame.floor<-data.frame(mean=c(mean.dust.vacuum.carpet,
                                        mean.dust.vacuum.hardflooring,
                                        mean.dust.walking.carpet,
                                        mean.dust.walking.hardflooring,
                                        mean.tripartite.vacuum.carpet,
                                        mean.tripartite.vacuum.hardflooring,
                                        mean.tripartite.walking.carpet,
                                        mean.tripartite.walking.hardflooring),
                                 sd=c(sd.dust.vacuum.carpet,
                                      sd.dust.vacuum.hardflooring,
                                      sd.dust.walking.carpet,
                                      sd.dust.walking.hardflooring,
                                      sd.tripartite.vacuum.carpet,
                                      sd.tripartite.vacuum.hardflooring,
                                      sd.tripartite.walking.carpet,
                                      sd.tripartite.walking.hardflooring),
                                 activity=rep(c(rep("Vacuuming",length(conclist)*2),rep("Walking",length(conclist)*2)),2),
                                 surface=rep(c(rep("Carpet",length(conclist)),rep("Hard Flooring",length(conclist))),4),
                                 type=c(rep("Dust",length(conclist)*4),rep("Tripartite",length(conclist)*4)),
                                 conc=rep(conclist,8),
                                 height=heights[j])
  }else{
    temp<-data.frame(mean=c(mean.dust.vacuum.carpet,
                                        mean.dust.vacuum.hardflooring,
                                        mean.dust.walking.carpet,
                                        mean.dust.walking.hardflooring,
                                        mean.tripartite.vacuum.carpet,
                                        mean.tripartite.vacuum.hardflooring,
                                        mean.tripartite.walking.carpet,
                                        mean.tripartite.walking.hardflooring),
                                 sd=c(sd.dust.vacuum.carpet,
                                      sd.dust.vacuum.hardflooring,
                                      sd.dust.walking.carpet,
                                      sd.dust.walking.hardflooring,
                                      sd.tripartite.vacuum.carpet,
                                      sd.tripartite.vacuum.hardflooring,
                                      sd.tripartite.walking.carpet,
                                      sd.tripartite.walking.hardflooring),
                                 activity=rep(c(rep("Vacuuming",length(conclist)*2),rep("Walking",length(conclist)*2)),2),
                                 surface=rep(c(rep("Carpet",length(conclist)),rep("Hard Flooring",length(conclist))),4),
                                type=c(rep("Dust",length(conclist)*4),rep("Tripartite",length(conclist)*4)),
                                 conc=rep(conclist,8),
                                  height=heights[j])
    plot.frame.floor<-rbind(plot.frame.floor,temp)
  }

}


# Figure 8
windows()
ggplot(plot.frame.floor)+geom_ribbon(aes(x=conc,ymin=mean-sd,ymax=mean+sd,fill=height,group=interaction(type,height)),alpha=0.1)+
  geom_line(aes(x=conc,y=mean,color=height,linetype=type),size=1)+
  facet_grid(activity~surface)+
  scale_x_continuous(trans="log10",name=expression("Floor Concentration (viral particles/cm"^2*")"))+
  scale_y_continuous(trans="log10",name="Infection Risk (Mean +/- SD)")+
  scale_linetype_discrete(name="")+
  scale_fill_discrete(name="")+
  scale_color_discrete(name="")+
  theme(axis.text=element_text(size=13),axis.title = element_text(size=13),
        strip.text=element_text(size=13),legend.text = element_text(size=13))
