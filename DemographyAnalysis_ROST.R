#################################################################################################################################
#                    ROST Demography Analysis w/ RMARK                            ###############################################
#                                Author: Kayla Davis                              ###############################################
#                                                                                 ###############################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# set working directory
setwd("G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK")

# load in capture history files
ch.all<- as.matrix(read.csv("HYResights_AllProofed.csv", sep = ",", header = T))
ch<- as.matrix(read.csv("ch_ROST.csv", sep = ",", header = T))
chat<- as.matrix(read.csv("ch_reversed.csv", sep = ",", header = T)) 

# load libraries
library(plyr)
library(RMark)
library(ggplot2)
library(cowplot)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# First investigate naive stopover duration, we don't use this because it doesn't
# incorporate detection probability and is therefore kind of useless, but it is
# interesting as a first look at the data.

# This is mean minimum stopover, which is the difference in days/weeks
# between first and last days resighted.

# calculate mean minimum stopover duration
# for individuals resighted >=2 times
# this is the interval (7days=1 interval) between first and last resighting
ch.all<-as.data.frame(ch.all)

#change to correct date format
ch.all$DATE<-as.Date(ch.all$DATE, format="%m/%d/%Y")

#create a dataframe with only ID and date
ind<-ch.all$ID
day<-ch.all$DATE

stopover<-data.frame(ind, day)

# Find first and last occurrence by date
g1 <- aggregate(stopover$day, list(stopover$ind), min)
colnames(g1)[1] = "ID"
colnames(g1)[2] = "first"

g2 <- aggregate(stopover$day, list(stopover$ind), max)
colnames(g2)[1] = "ID"
colnames(g2)[2] = "last"

# Merge data frames by ID
a1 <- merge(g1, g2, by=c("ID"))

#calculate date difference in weeks between first and last resighting
a1$date.diff<-difftime(a1$last,a1$first, units=c("weeks"))

#merge date dataframe with original data frame with info for further modeling
names(group.df)[1]<-paste("ID")
stopover<-merge(a1,group.df,by="ID")

#eliminate individuals resighted only once
stopover<-subset(stopover, date.diff>0)

#find mean min stop duration in weeks
mean(stopover$date.diff)
sd(stopover$date.diff)
sqrt(var(stopover$date.diff)/length(stopover$date.diff)) #se

#look at mean min stop duration for all groups
ddply(stopover, c("region","year"), summarise,
      N = length(date.diff),
      mean = mean(date.diff),
      sd   = sd(date.diff),
      se   = sd/sqrt(N))

#test for differences in mean min stopover duration between years and regions     
t.test(stopover$date.diff~stopover$year)

t.test(stopover$date.diff~stopover$region)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# Now we will start formatting the capture history data for our RMark analysis 

# Convert capture history file into a data frame 
group.df<-as.data.frame(ch)

# function to concatenate capture history into strings
conc<-function(x)
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {ch<-(x[i,]>0)*1
  out[i]<-paste(ch[1],ch[2],ch[3],ch[4],ch[5],ch[6],ch[7],ch[8],ch[9],ch[10],sep="")}

  return (out)}

# new data frame with covariates and group effects
capt.hist<-data.frame(ch=conc(ch[,2:11]), year= group.df$year, region=group.df$region)

# Set up Pradel analysis
# Note the difference in time interval for the last capture occasion
# All capture windows were 7 days long (in both years) except for the last window
# which was 8 days. We have corrected for that by including 1.14 as the last interval.
capt.hist.process=process.data(capt.hist, model="Pradrec", groups=c("year","region"),
                               time.intervals = c(1,1,1,1,1,1,1,1,1.14))

#Set up the design data
pradel.ddl<-make.design.data(capt.hist.process, parameters = list(p=list(time.bins=unique(c(0,2,3,4,5,6,7,8,10.14)))))
names(pradel.ddl)
summary(pradel.ddl)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# Create the models to test 

# For the first step, we use full model on phi and f and find best structure for p
pradel.models1<- function()
{
  Phi.1=list(formula=~time*year*region)
  f.1=list(formula=~time*year*region)
  
  p.1=list(formula=~1)
  p.2=list(formula=~year)
  p.3=list(formula=~time)
  p.4=list(formula=~region)
  p.5=list(formula=~time*year*region)
  p.6=list(formula=~time*region)
  p.7=list(formula=~year*region)
  p.8=list(formula=~time*region+year)
  p.9=list(formula=~time*year+region)
  p.10=list(formula=~year+region+time)
  p.11=list(formula=~year*time)
  p.12=list(formula=~Time)
  p.13=list(formula=~I(Time^2))
  

# Compile a list of models
  pradel.cml<-create.model.list("Pradrec")
  pradel.cml
  results=mark.wrapper(pradel.cml, data=capt.hist.process,ddl=pradel.ddl, adjust=T,
                       invisible=F)
  return(results)
}

pradel.results1<-pradel.models1()
save.image()


# Create an AIC table to compare models
pradel.results1$model.table$Likelihood = exp(-0.5*pradel.results1$model.table$DeltaAICc)
pradel.results1$model.table

####################################################################################
# Sidebar to compute c-hat in RELEASE

# Do GOF testing on forward and reversed capture histories to adjust chat appropriately
# bring up program release to calculate c hat
# test 2+test3/df
release.gof(capt.hist.process,invisible = TRUE, title = "Release-gof", view = TRUE)

# adjust c-hat to account for extra binomial variation detected by release
# release results
# Goodness of Fit Results (TEST 2 + TEST 3) by Group

#Group  Chi-square   df   P-level
#-----  ----------  ----  -------
#1      52.4368    22    0.0003
#2      70.3737    26    0.0000
#3      38.2917    28    0.0930
#4      44.0631    26    0.0149
#Total     205.1652   102    0.0000
#c hat = 2.01

# make reversed capture history file a df to run through RELEASE
chat.df<-as.data.frame(chat)

# Do same procedure as above to create the correct capture history format
conc<-function(x)
{
  k<-ncol(x)
  n<-nrow(x)
  out<-array(dim=n)
  for (i in 1:n)
  {chat<-(x[i,]>0)*1
  out[i]<-paste(chat[1],chat[2],chat[3],chat[4],chat[5],chat[6],chat[7],chat[8],chat[9],chat[10],sep="")}
  
  return (out)}

# new data frame with covs and group effects
capt.hist.chat<-data.frame(chat=conc(chat[,2:11]),region=chat.df$region,year=chat.df$year,combo= chat.df$group)

capt.hist.chat.process=process.data(capt.hist.chat, model="CJS", groups=c("region","year","combo"),
                                    time.intervals = c(1,1,1,1,1,1,1,1,1.14))


# Run release on reversed capture history
release.gof(capt.hist.chat.process,invisible = TRUE, title = "Release-gof", view = TRUE)

# End sidebar
####################################################################################

# Adjust the c-hat using the highest c-hat from the forward and reversed RELEASE test
pradel.results1.adj<-adjust.chat(chat=2.03, pradel.results1)


# Compute matrices of model weights, number of parameters and Delta AICc values
p.model.weight<-pradel.results1.adj$model.table$weight
                          
p.model.npar<-pradel.results1.adj$model.table$npar
                         
p.model.DeltaAICc<-pradel.results1.adj$model.table$DeltaQAICc
                              
p.model.names<-pradel.results1.adj$model.table$model

p.model.AICc<-pradel.results1.adj$model.table$QAICc

p.model.deviance<-pradel.results1.adj$model.table$QDeviance
                            
# Make the above matrix a data frame for easy viewing
p.results<-data.frame(p.model.names,p.model.AICc, p.model.DeltaAICc,p.model.weight,p.model.npar,p.model.deviance)
p.results$Likelihood=exp(-0.5*p.results$p.model.DeltaAICc)

# And name the columns
colnames(p.results)[1] = "Model"
colnames(p.results)[2] = "QAICc"
colnames(p.results)[3] = "Delta QAICc"
colnames(p.results)[4]="weight"
colnames(p.results)[5]="K"
colnames(p.results)[6]="Deviance"
colnames(p.results)[7]="Likelihood"

# Tidy up with rounding
p.results[,2:7]<-round(p.results[,2:7],2)
p.results<-p.results[,c(1,2,3,4,7,5,6)]

# Voila, see the results
p.results

# And save to csv
write.csv(p.results,file="pradel.results.p.csv") # produces a csv file of above results

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#For the second step, use best structure from the above analysis for p and fit phi
pradel.models2<- function()
{
  
  f.1=list(formula=~time*year*region)
  
  Phi.1=list(formula=~1)
  Phi.2=list(formula=~year*time)
  Phi.3=list(formula=~time)
  Phi.4=list(formula=~region)
  Phi.5=list(formula=~time*year*region)
  Phi.6=list(formula=~time*region)
  Phi.7=list(formula=~year)
  Phi.8=list(formula=~time*region+year)
  Phi.9=list(formula=~time*year+region)
  Phi.10=list(formula=~Time)
  Phi.11=list(formula=~I(Time^2))
  Phi.12=list(formula=~time+year+region)
  Phi.13=list(formula=~year*region)
  
  
  p.5=list(formula=~time*year+region)
  
  
# Compile a list of models
  pradel.cml<-create.model.list("Pradrec")
  pradel.cml
  results=mark.wrapper(pradel.cml, data=capt.hist.process,ddl=pradel.ddl, adjust=T,
                       invisible=F)
  return(results)
}
pradel.results2<-pradel.models2()
save.image()

# Create the AIC table
pradel.results2$model.table$Likelihood = exp(-0.5*pradel.results2$model.table$DeltaAICc)
pradel.results2$model.table

# Adjust c-hat
pradel.results2.adj<-adjust.chat(chat=2.03, pradel.results2)

# Compute matrices of model weights, number of parameters and Delta AICc values
phi.model.weight<-pradel.results2.adj$model.table$weight

phi.model.npar<-pradel.results2.adj$model.table$npar

phi.model.DeltaAICc<-pradel.results2.adj$model.table$DeltaQAICc

phi.model.names<-pradel.results2.adj$model.table$model

phi.model.AICc<-pradel.results2.adj$model.table$QAICc

phi.model.deviance<-pradel.results2.adj$model.table$QDeviance

# Convert to a data frame
phi.results<-data.frame(phi.model.names,phi.model.AICc, phi.model.DeltaAICc,phi.model.weight,phi.model.npar,phi.model.deviance)
phi.results$Likelihood=exp(-0.5*phi.results$phi.model.DeltaAICc)

# Name the columns
colnames(phi.results)[1] = "Model"
colnames(phi.results)[2] = "QAICc"
colnames(phi.results)[3] = "Delta QAICc"
colnames(phi.results)[4]="weight"
colnames(phi.results)[5]="K"
colnames(phi.results)[6]="Deviance"
colnames(phi.results)[7]="Likelihood"

# Tidy up
phi.results[,2:7]<-round(phi.results[,2:7],2)
phi.results<-phi.results[,c(1,2,3,4,7,5,6)]

# And presto! Results!
phi.results

# Write results to csv
write.csv(phi.results,file="pradel.results.phi.csv")

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#For the third step, use best structure for p and phi, fit f
pradel.models3<- function()
{
  
  f.1=list(formula=~1)
  f.2=list(formula=~year*time)
  f.3=list(formula=~time)
  f.4=list(formula=~region)
  f.5=list(formula=~time*year*region)
  f.6=list(formula=~time*region)
  f.7=list(formula=~year)
  f.8=list(formula=~time*year+region)
  f.9=list(formula=~time+year+region)
  f.10=list(formula=~year*region)
  f.11=list(formula=~Time)
  f.12=list(formula=~I(Time^2))
  f.13=list(formula=~time*region+year)
  
  
  Phi.10=list(formula=~Time)
  Phi.11=list(formula=~I(Time^2))
  
  p.5=list(formula=~time*year+region)
  
  
# Compile a list of models
  pradel.cml<-create.model.list("Pradrec")
  pradel.cml
  results=mark.wrapper(pradel.cml, data=capt.hist.process,ddl=pradel.ddl, adjust=T,
                       invisible=F)
  return(results)
}
pradel.results3<-pradel.models3()
save.image()

# Create the AIC table
pradel.results3$model.table$Likelihood = exp(-0.5*pradel.results3$model.table$DeltaAICc)
pradel.results3$model.table

# Adjust with c-hat
pradel.results3.adj<-adjust.chat(chat=2.03, pradel.results3)

# Compute matrices of model weights, number of parameters and Delta AICc values
f.model.weight<-pradel.results3.adj$model.table$weight

f.model.npar<-pradel.results3.adj$model.table$npar

f.model.DeltaAICc<-pradel.results3.adj$model.table$DeltaQAICc

f.model.names<-pradel.results3.adj$model.table$model

f.model.AICc<-pradel.results3.adj$model.table$QAICc

f.model.deviance<-pradel.results3.adj$model.table$QDeviance

# Convert to data frame
f.results<-data.frame(f.model.names,f.model.AICc, f.model.DeltaAICc,f.model.weight,f.model.npar,f.model.deviance)
f.results$Likelihood=exp(-0.5*f.results$f.model.DeltaAICc)

# Name the columns
colnames(f.results)[1] = "Model"
colnames(f.results)[2] = "QAICc"
colnames(f.results)[3] = "Delta QAICc"
colnames(f.results)[4]="weight"
colnames(f.results)[5]="K"
colnames(f.results)[6]="Deviance"
colnames(f.results)[7]="Likelihood"

# Round it nicely
f.results[,2:7]<-round(f.results[,2:7],2)
f.results<-f.results[,c(1,2,3,4,7,5,6)]

# Take a look
f.results

# Write to csv
write.csv(f.results,file="pradel.results.f.csv")#produces an excel file of above results

# Now look at top model
pradel.results3$Phi.11.p.5.f.9

# See the results on the real scale
pradel.results3$Phi.11.p.5.f.9$results$real


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# this is to compute results for best model from above analyses
# instead of doing the entire analysis each time

pradel.models<- function()
{
  Phi.11=list(formula=~I(Time^2))
  
  p.5=list(formula=~time*year+region)
  
  f.9=list(formula=~time+year+region)
  
  #Compile a list of models
  pradel.cml<-create.model.list("Pradrec")
  pradel.cml
  results=mark.wrapper(pradel.cml, data=capt.hist.process,ddl=pradel.ddl, adjust=T,
                       invisible=F)
  return(results)
}
pradel.results<-pradel.models()


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
# Make some pretty figures

# Make a position nudge value so our plots don't have ugly overlapping points
position<-position_nudge(0.5,0)

# Pull our real results for phi into a data frame
# This will be easier to work with than the model name$$$$$
results.phi<-pradel.results$Phi.11.p.5.f.9$results$real[1:9,1:5]
results.phi$x<-(c("1","2","3","4","5","6","7","8","9"))

# This section creates the points with 95% CIs 
phigraph <- ggplot(data=results.phi, aes(x=x, y=estimate)) + 
  geom_errorbar(data=results.phi, aes(x=x, ymin=lcl, ymax=ucl), width=.1, color="black", position=position) +
  geom_point(data=results.phi, aes(x=x, y=estimate), size=4.5, color="black",fill="gray", shape=21, position=position) +
  theme_bw() +
  xlab(NULL) +
  ylab("Residency Rate") +
  theme(axis.title.x = element_text(size = 19, vjust=-.5, face="bold")) + #adjust font for axis titles and labels
  theme(axis.title.y = element_text(size = 19, vjust=1.5, face="bold")) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +

  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10"), labels=c("16-22 July", "23-29 July", "30 July-5 Aug", "6-12 Aug", "13-19 Aug",
                                                                               "20-26 Aug", "27 Aug-2 Sep", "3-9 Sep", "10-16 Sep", "17-24 Sep"))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  theme(panel.grid.minor = element_blank()) + #this tells it not to plot gridlines
  theme(panel.grid.major = element_blank()) +
  theme(legend.position="none")

# command to print plot
phigraph

#this code creates a higher resolution image and saves as png to filepath below
tiff(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK//phi.jpg", width = 480 * 16, height = 480 * 12,  pointsize = 12 * 1.5, res = 1200)
phigraph
dev.off()

#I used this bit of code to quickly figure out the appropriate date ranges for the x-axis
#first.day<-as.Date(c("07/16"),"%m/%d")
#last.day<-as.Date(c("09/23"),"%m/%d")
#dates<-seq(min(first.day),max(last.day),by=7)
#dates<-format(dates,format="%m-%d")
#dates

##########################################################################################################################

# Pull out the real results for recruitment and make dataframe for f
results.f<-pradel.results$Phi.11.p.5.f.9$results$real[42:77,1:5]
results.f$x<-rep(c("1","2","3","4","5","6","7","8","9"))
results.f$group<-as.factor(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
                   3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4))

# Make an appropriate value for dodging
pd<-position_dodge(0.4)

# next plot for recruitment
fgraph <- ggplot(data=results.f, aes(x=as.numeric(x)+0.5, y=estimate,  fill=group)) + 
  geom_errorbar(data=results.f, aes(x=as.numeric(x)+0.5, ymin=lcl, ymax=ucl, group=group), position=pd, width=.1) +
  geom_point(data=results.f, aes(x=as.numeric(x)+0.5, y=estimate,  fill=group), color="black", shape=21,position=pd, size=4.5)+
  scale_fill_grey(aesthetics = c("fill"), breaks=c("1", "2", "3", "4"), labels=c( "2014 Northern Region","2015 Northern Region",
                                                                                "2014 Southern Region","2015 Southern Region"), 
                                                                                  name=NULL,start=0,end=1)+
  theme_bw() +
  xlab(NULL) +
  ylab("Recruitment Rate") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10"), labels=c("16-22 July", "23-29 July", "30 July-5 Aug", "6-12 Aug", "13-19 Aug",
                                                                                "20-26 Aug", "27 Aug-2 Sep", "3-9 Sep", "10-16 Sep", "17-24 Sep"))+
  theme(axis.title.x = element_text(size = 19, vjust=-.5, face="bold")) + #adjust font for axis titles and labels
  theme(axis.title.y = element_text(size = 19, vjust=1.5, face="bold")) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(panel.grid.minor = element_blank()) + #this tells it not to plot gridlines
  theme(panel.grid.major = element_blank()) +
  theme(legend.justification=c(0.98,0.98),
        legend.position=c(0.98,0.98))+
  theme(legend.text=element_text(size=16))  

#command to print plot
fgraph

#this code creates a higher resolution image and saves as png to filepath below
png(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK//f.png", width = 480 * 16, height = 480 * 12,  pointsize = 12 * 1.5, res = 600)
fgraph
dev.off()

##################################################################################################################################

# Pull the real results and make dataframe for lambda
results.lambda<-pradel.results$Phi.11.p.5.f.9$results$derived$`Lambda Population Change`
results.lambda$x<-rep(c("1","2","3","4","5","6","7","8","9"))
results.lambda$group<-as.factor(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
                             3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4))
# next plot for lambda
lambdagraph <- ggplot() + 
  geom_abline(slope=0, intercept = 0, linetype = 2)+
  geom_errorbar(data=results.lambda, aes(x=as.numeric(x)+0.5, ymin=log(lcl), ymax=log(ucl), group=group), position=pd,  width=.1) +
  geom_point(data=results.lambda, aes(x=as.numeric(x)+0.5, y=log(estimate),  fill=group), color="black", shape=21,position=pd,size=4.5)+
  scale_fill_grey(aesthetics = c("fill"), breaks=c("1", "2", "3", "4"), labels=c( "2014 Northern Region","2015 Northern Region",
                                                                                  "2014 Southern Region","2015 Southern Region"), 
                                                                                    name=NULL,start=0,end=1)+
  theme_bw() +
  xlab("Week of Staging Season") +
  ylab("Log of Lambda (r)") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))+
  scale_x_discrete(limits=c("1","2","3","4","5","6","7","8","9","10"), labels=c("16-22 July", "23-29 July", "30 July-5 Aug", "6-12 Aug", "13-19 Aug",
                                                                                "20-26 Aug", "27 Aug-2 Sep", "3-9 Sep", "10-16 Sep", "17-24 Sep"))+
  theme(axis.title.x = element_text(size = 19, vjust=-.5, face="bold")) + #adjust font for axis titles and labels
  theme(axis.title.y = element_text(size = 19, vjust=1.5, face="bold")) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(panel.grid.minor = element_blank()) + #this tells it not to plot gridlines
  theme(panel.grid.major = element_blank()) +
  theme(legend.justification=c(0.98,0.98),
        legend.position=c(0.98,0.98))+
  theme(legend.text=element_text(size=16))

#command to print plot
lambdagraph

#this code creates a higher resolution image and saves as png to filepath below
png(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK//lambda.png", width = 480 * 16, height = 480 * 12,  pointsize = 12 * 1.5, res = 600)
lambdagraph
dev.off()

##################################################################################################################################

# group the previous plots together in grid 
parms<-plot_grid(phigraph,fgraph, lambdagraph,  
                 labels=c("A", "B", "C"),nrow=3, ncol=1)
parms
jpeg(filename = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK//Figure3.jpg", width = 480 * 10, height = 480 * 18,  pointsize = 12 * 1.5, res = 600)
parms
dev.off()

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#make the input files for SODA analysis

#file needs to have this format:
# no. birds caught(1630)  no. cap. occ.(10)
# $
# $
# capture histories w/ occ. separated by at least one blank

# For all of the files below, we exported to a csv to insert the leading rows
# Then we saved the files as txt files to run in SODA gui

# We ended up only using the total soda file and the 2014N, 2014S, 2015N, and 2015S files for the paper

# Create the total SODA file
soda<-matrix(c(group.df$wk1,group.df$wk2,group.df$wk3,group.df$wk4,
               group.df$wk5,group.df$wk6,group.df$wk7,group.df$wk8,
               group.df$wk9,group.df$wk10),nrow=1630,ncol = 10)

#replace 1s and 2s from R language with 0s and 1s
soda[soda<2]<-0
soda[soda>1]<-1
soda

# Export to csv to add the dollar sign to the leading rows as shown above (couldn't get this to work in R, so had to do manually for all files below
# which is a major bummer)
write.csv(soda,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//RMARK//soda_master.csv")

# capture histories w/ occ. separated by at least one blank for 2014 data
soda.2014<-subset(group.df,year=="2014")
soda.2014<-matrix(c(soda.2014$wk1,soda.2014$wk2,soda.2014$wk3,soda.2014$wk4,
                    soda.2014$wk5,soda.2014$wk6,soda.2014$wk7,soda.2014$wk8,
                    soda.2014$wk9,soda.2014$wk10),nrow=754,ncol = 10)

#replace 1s and 2s from R language with 0s and 1s
soda.2014[soda.2014<2]<-0
soda.2014[soda.2014>1]<-1

# Export to csv for manual manipulation
write.table(soda.2014,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.2014.txt")

# Same procedure as above for 2015 data
soda.2015<-subset(group.df,year=="2015")
soda.2015<-matrix(c(soda.2015$wk1,soda.2015$wk2,soda.2015$wk3,soda.2015$wk4,
                    soda.2015$wk5,soda.2015$wk6,soda.2015$wk7,soda.2015$wk8,
                    soda.2015$wk9,soda.2015$wk10),nrow=876,ncol = 10)
soda.2015[soda.2015<2]<-0
soda.2015[soda.2015>1]<-1
write.table(soda.2015,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.2015.txt")

# Do same for birds from northern region colonies 
soda.n<-subset(group.df, region=="north")
soda.n<-matrix(c(soda.n$wk1,soda.n$wk2,soda.n$wk3,soda.n$wk4,
                    soda.n$wk5,soda.n$wk6,soda.n$wk7,soda.n$wk8,
                    soda.n$wk9,soda.n$wk10),nrow=515,ncol = 10)
soda.n[soda.n<2]<-0
soda.n[soda.n>1]<-1
write.table(soda.n,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.n.txt")

# Again, same as above for birds from southern region colonies
soda.s<-subset(group.df, region=="south")
soda.s<-matrix(c(soda.s$wk1,soda.s$wk2,soda.s$wk3,soda.s$wk4,
                    soda.s$wk5,soda.s$wk6,soda.s$wk7,soda.s$wk8,
                    soda.s$wk9,soda.s$wk10),nrow=1115,ncol = 10)
soda.s[soda.s<2]<-0
soda.s[soda.s>1]<-1
write.table(soda.s,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.s.txt")

# Now do northern region in 2014
soda.2014n<-subset(group.df, year=="2014"& region=="north")
soda.2014n<-matrix(c(soda.2014n$wk1,soda.2014n$wk2,soda.2014n$wk3,soda.2014n$wk4,
                    soda.2014n$wk5,soda.2014n$wk6,soda.2014n$wk7,soda.2014n$wk8,
                    soda.2014n$wk9,soda.2014n$wk10),nrow=239,ncol = 10)
soda.2014n[soda.2014n<2]<-0
soda.2014n[soda.2014n>1]<-1
write.table(soda.2014n,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.2014n.txt")

# Northern region in 2015
soda.2015n<-subset(group.df, year=="2015"& region=="north")
soda.2015n<-matrix(c(soda.2015n$wk1,soda.2015n$wk2,soda.2015n$wk3,soda.2015n$wk4,
                     soda.2015n$wk5,soda.2015n$wk6,soda.2015n$wk7,soda.2015n$wk8,
                     soda.2015n$wk9,soda.2015n$wk10),nrow=276,ncol = 10)
soda.2015n[soda.2015n<2]<-0
soda.2015n[soda.2015n>1]<-1
write.table(soda.2015n,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.2015n.txt")

# Southern region in 2014
soda.2014s<-subset(group.df, year=="2014"& region=="south")
soda.2014s<-matrix(c(soda.2014s$wk1,soda.2014s$wk2,soda.2014s$wk3,soda.2014s$wk4,
                     soda.2014s$wk5,soda.2014s$wk6,soda.2014s$wk7,soda.2014s$wk8,
                     soda.2014s$wk9,soda.2014s$wk10),nrow=515,ncol = 10)
soda.2014s[soda.2014s<2]<-0
soda.2014s[soda.2014s>1]<-1
write.table(soda.2014s,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.2014s.txt")

#Southern region in 2015
soda.2015s<-subset(group.df, year=="2015"& region=="south")
soda.2015s<-matrix(c(soda.2015s$wk1,soda.2015s$wk2,soda.2015s$wk3,soda.2015s$wk4,
                     soda.2015s$wk5,soda.2015s$wk6,soda.2015s$wk7,soda.2015s$wk8,
                     soda.2015s$wk9,soda.2015s$wk10),nrow=600,ncol = 10)
soda.2015s[soda.2015s<2]<-0
soda.2015s[soda.2015s>1]<-1
write.table(soda.2015s,file = "G://My Drive//R Working Directory//ROST//Resights//DemographyAnalyses//ROST_Demography//SODA//soda.2015s.txt")

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

#did not plot p for the paper
#but here it is for kicks

#dataframe for p
results.p<-pradel.results3$Phi.11.p.5.f.9$results$real[10:41,1:5]

results.p$x<-rep(c("1","2","3","4","5","6","7","8"))
results.p$group<-as.factor(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,
                             3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4))

#next plot for detection
#next plot for recruitment
pd<-position_dodge(0.25)

pgraph <- ggplot() + 
  geom_errorbar(data=results.p, aes(x=x, ymin=lcl, ymax=ucl, group=group), position=pd,  width=.1) +
  geom_point(data=results.p, aes(x=x, y=estimate,  fill=group), color="black", shape=group,position=pd,size=4.5)+
  scale_fill_grey(aesthetics = c("fill"), breaks=c("1", "2", "3", "4"), labels=c( "2014 Northern Region","2015 Northern Region",
                                                                                  "2014 Southern Region","2015 Southern Region"),
                  name=NULL,start=0.1,end=0.9)+ 
  
  theme_bw() +
  xlab("Week of Staging Season") +
  ylab("Detection Probability") +
  #scale_x_continuous(breaks=seq(0,8,1)) +
  theme(axis.title.x = element_text(size = 19, vjust=-.5, face="bold")) + #adjust font for axis titles and labels
  theme(axis.title.y = element_text(size = 19, vjust=1.5, face="bold")) +
  theme(axis.text.x = element_text(size=16)) +
  theme(axis.text.y = element_text(size=16)) +
  theme(panel.grid.minor = element_blank()) + #this tells it not to plot gridlines
  theme(panel.grid.major = element_blank()) +
  theme(legend.justification=c(1,1),
        legend.position=c(1,1))  
#command to print plot
pgraph
#this code creates a higher resolution image and saves as png to filepath below
png(filename = "Z:\\Individual Files\\Kayla\\R Working Directory\\ROST\\Resights\\Ch3\\FigureData\\detection.png", width = 480 * 16, height = 480 * 12,  pointsize = 12 * 1.5, res = 600)
pgraph
dev.off()