
#Code to process SODA output
#and create SODA model figures

#set working directory
setwd("G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA")

#load data
total<- as.data.frame(read.table("soda.total.txt", sep = "", header = F))
total<-total[,1:10]

total.mean<-apply(total,2,mean)
total.sort<-apply(total,2,sort)
total.lower<-data.frame(total.sort[13,])
total.upper<-data.frame(total.sort[488,])
soda.total<-data.frame(total.mean, total.lower,total.upper)
colnames(soda.total)[1]<-"mean"
colnames(soda.total)[2]<-"lower"
colnames(soda.total)[3]<-"upper"
soda.total$x<-1:10
soda.total

#2014n
n14<- as.data.frame(read.table("soda.2014n.output.txt", sep = "", header = F))
n14<-n14[,1:10]

n14.mean<-apply(n14,2,mean)
n14.sort<-apply(n14,2,sort)
n14.lower<-data.frame(n14.sort[13,])
n14.upper<-data.frame(n14.sort[488,])
soda.n14<-data.frame(n14.mean, n14.lower,n14.upper)
colnames(soda.n14)[1]<-"mean"
colnames(soda.n14)[2]<-"lower"
colnames(soda.n14)[3]<-"upper"
soda.n14$group<-"2014 Northern Region"

soda.n14

#2015n
n15<- as.data.frame(read.table("soda.2015n.output.txt", sep = "", header = F))
n15<-n15[,1:10]

n15.mean<-apply(n15,2,mean)
n15.sort<-apply(n15,2,sort)
n15.lower<-data.frame(n15.sort[13,])
n15.upper<-data.frame(n15.sort[488,])
soda.n15<-data.frame(n15.mean, n15.lower,n15.upper)
colnames(soda.n15)[1]<-"mean"
colnames(soda.n15)[2]<-"lower"
colnames(soda.n15)[3]<-"upper"
soda.n15$group<-"2015 Northern Region"

soda.n15

#2014s
s14<- as.data.frame(read.table("soda.2014s.output.txt", sep = "", header = F))
s14<-s14[,1:10]

s14.mean<-apply(s14,2,mean)
s14.sort<-apply(s14,2,sort)
s14.lower<-data.frame(s14.sort[13,])
s14.upper<-data.frame(s14.sort[488,])
soda.s14<-data.frame(s14.mean, s14.lower,s14.upper)
colnames(soda.s14)[1]<-"mean"
colnames(soda.s14)[2]<-"lower"
colnames(soda.s14)[3]<-"upper"
soda.s14$group<-"2014 Southern Region"

soda.s14

#2015s
s15<- as.data.frame(read.table("soda.2015s.output.txt", sep = "", header = F))
s15<-s15[,1:10]

s15.mean<-apply(s15,2,mean)
s15.sort<-apply(s15,2,sort)
s15.lower<-data.frame(s15.sort[13,])
s15.upper<-data.frame(s15.sort[488,])
soda.s15<-data.frame(s15.mean, s15.lower,s15.upper)
colnames(soda.s15)[1]<-"mean"
colnames(soda.s15)[2]<-"lower"
colnames(soda.s15)[3]<-"upper"
soda.s15$group<-"2015 Southern Region"
soda.s15

soda.groups<-rbind(soda.n14,soda.n15,soda.s14,soda.s15)
soda.groups$x<-rep(1:10)

#use above for figures
#one for total stopover duration estimates
#one for stopover duration separated by year and region

library(ggplot2)

pd<-position_dodge(0.4)

group <- ggplot() + 
  geom_errorbar(data=soda.groups, aes(x=x, ymin=lower, ymax=upper, group=group), position=pd,  width=.1) +
  geom_point(data=soda.groups, aes(x=x, y=mean,  fill=group), color="black", shape=21,position=pd,size=4.5)+
  scale_fill_grey(aesthetics = c("fill"), breaks=c("2014 Northern Region",
                                                   "2015 Northern Region",
                                                   "2014 Southern Region",
                                                   "2015 Southern Region"), #labels=c("2014 Northern Region",
                                                                                     #"2015 Northern Region",
                                                                                     #"2014 Southern Region",
                                                                                     #"2015 Southern Region"),
                  name=NULL,start=0,end=1)+
  theme_bw() +
  xlab("Week of Staging Season") +
  ylab("Staging Duration at CCNS (weeks)") +
  scale_x_continuous(breaks=seq(1,10,1),
                     labels=c("16-22 July", "23-29 July", "30 July-5 Aug", "6-12 Aug", "13-19 Aug",
                     "20-26 Aug", "27 Aug-2 Sep", "3-9 Sep", "10-16 Sep", "17-24 Sep"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_y_continuous(breaks=seq(2,12,2))+
  theme(axis.title.x = element_text(size = 19, vjust=-.5, face="bold")) + #adjust font for axis titles and labels
  theme(axis.title.y = element_text(size = 19, vjust=1.5, face="bold")) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(panel.grid.minor = element_blank()) + #this tells it not to plot gridlines
  theme(panel.grid.major = element_blank()) +
  theme(legend.justification=c(0.98,0.98),
        legend.position=c(0.98,0.3))+
  theme(legend.text=element_text(size=16))  
#command to print plot
group
#this code creates a higher resolution image and saves as png to filepath below
jpeg(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//Fig4.jpg", width = 480 * 13, height = 480 * 10,  pointsize = 12 * 1.5, res = 600)
group
dev.off()


total <- ggplot() + 
  geom_errorbar(data=soda.total, aes(x=x, ymin=lower, ymax=upper),  width=.1) +
  geom_point(data=soda.total, aes(x=x, y=mean), fill="gray", color="black", shape=21,position=pd,size=4.5)+
  theme_bw() +
  xlab(NULL) +
  ylab("Staging Duration at CCNS (weeks)") +
  scale_x_continuous(breaks=seq(1,10,1))+
  #scale_y_continuous(breaks=seq(2,8,2))+
  theme(axis.title.x = element_text(size = 19, vjust=-.5, face="bold")) + #adjust font for axis titles and labels
  theme(axis.title.y = element_text(size = 19, vjust=1.5, face="bold")) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(panel.grid.minor = element_blank()) + #this tells it not to plot gridlines
  theme(panel.grid.major = element_blank()) +
  theme(legend.justification=c(0.98,0.98),
        legend.position=c(0.98,0.2))+
  theme(legend.text=element_text(size=16))  
#command to print plot
total
#this code creates a higher resolution image and saves as png to filepath below
png(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//total.png", width = 480 * 16, height = 480 * 12,  pointsize = 12 * 1.5, res = 600)
total
dev.off()

#group the previous plots together in grid 
library(cowplot)
soda<-plot_grid(total, group,  
                 labels=c("A", "B"),nrow=2, ncol=1)
soda
png(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK//soda.png", width = 480 * 16, height = 480 * 20,  pointsize = 12 * 1.5, res = 600)
soda
dev.off()