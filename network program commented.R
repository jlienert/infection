

#Make the chemotherapy dataset from the large .rda file
#the first two lines are actually done in linux.  I am not sure of the corresponding commands in .
#subset to chemotherapy ward
cd /data/lienertjp/Data

grep CHEMO pnpas.rpt >chemo.csv
#remove extra white spaces
sed 's/ * / /g' <chemo.csv >chemo2.csv

module load R
R

#Begining of R code

#install.packages('statnet')
install.packages("statnet", repos="http://cran.rstudio.org")
#install.packages('Matrix')
install.packages("Matrix", repos="http://cran.rstudio.org")
library(statnet)
library(Matrix)
options(width=as.integer(Sys.getenv("COLUMNS")))

setwd('C:\\Users\\lienertjp\\Desktop')
#load into R
chemo=read.table("chemo2.csv",sep=' ')
names(chemo)=c("idsetuid","source","spellid","AdmissionDate","DischargeDate","WardStartDate","WardEndDate","WardName","bed","Bay","Building","Facility","AdminCategoryCode","AdmissionMethodCode","AdmissionSourceCode","DischargeDestinationCode","DischargeMethodCode","MainSpecialtyCode")
chemo2=with(chemo,data.frame(idsetuid=idsetuid,source=source,spellid=paste(spellid,AdmissionDate,sep=" "),AdmissionDate=paste(DischargeDate,WardStartDate,sep=" "),DischargeDate=paste(WardEndDate,WardName,sep=" "),WardStartDate=paste(bed,Bay,sep=" "),WardEndDate=paste(Building,Facility,sep=" "),WardName=AdminCategoryCode,bed=AdmissionMethodCode,Bay=AdmissionSourceCode,Building=DischargeDestinationCode,Facility=DischargeMethodCode,AdminCategoryCode=MainSpecialtyCode,AdmissionMethodCode=chemo[,19],AdmissionSourceCode=chemo[,20],DischargeDestinationCode=chemo[,21],DischargeMethodCode=chemo[,22],MainSpecialtyCode=chemo[,23]))

#add in demographic variables and sort by date
demo=read.table('pndemo.csv',sep=',')
chemodemo=demo[demo$idsetUID%in%chemo$idsetuid,]
chemodemo$idsetuid=chemodemo$idsetUID
names(chemodemo)=c("patientid","birthdate","deathdate","sex",'one',"idsetuid")
chemodemo=chemodemo[,c(-5)]

chemo2=merge(chemo2,chemodemo)
chemo=chemo2[order(chemo2$WardStartDate),]
chemo$death=0
chemo$deathdate=as.character(chemo$deathdate)
for(i in 1:nrow(chemo)){
  if(is.na(chemo$deathdate[i])){chemo$deathdate[i]="2015-06-30 00:00:000"}
  else chemo$death[i]=1
}

#Makes variables for the start and end times, in hours from Jan 1 2000 for each visit to the ward
chemo$WardStart=as.numeric(difftime(as.POSIXct(chemo$WardStartDate),as.POSIXct('2000-01-01 00:00:00.000'),units="hours"))
chemo$WardEnd=as.numeric(difftime(as.POSIXct(chemo$WardEndDate),as.POSIXct('2000-01-01 00:00:00.000'),units="hours"))

#Makes a dataset with only the first entry for each person
chemoFirst=chemo[!duplicated(chemo$idsetuid),]

#METHODS TO MAKE THE NETWORK

#Creates a list where each entry is a vector of the hour numbers for each patient in the ward
hourList=list()
for (i in 1:length(chemoFirst$idsetuid)) hourList[[i]]=c(0)
names(hourList)=chemoFirst$idsetuid

for (i in 1:nrow(chemo)){
  times=floor(chemo$WardStart[i]):floor(chemo$WardEnd[i])
  hourList[[as.character(chemo$idsetuid[i])]]=c(hourList[[as.character(chemo$idsetuid[i])]],times)
}

#removes the initial 0 in each entry in hourList, which was used as a placeholder
for (i in 1:length(hourList)) hourList[[i]]=hourList[[i]][2:length(hourList[[i]])]

#Make a matrix of the same thing
#Less efficient, but seful later when doing row and column sums.
p=length(unique(chemo$idsetuid))
ps=as.character(unique(chemo$idsetuid))
t=87600
ts=paste('T',1:t,sep='')
sMat=matrix(0,nrow=p,ncol=t,dimnames=list(ps,ts))
for(i in 1:length(hourList)){
  for (j in hourList[[i]]){
    sMat[chemoFirst$idsetuid[i],j]=1 
  }
}

#Sparse Matrix
#More efficient than regular matrix, given the large number of 0s
#sMat2=Matrix(sMat,sparse=T)

#Creates a matrix where each i,j is the hours of overlap between patient i and patient j
#diagonal is the total number of hours patient i was in the ward
ps=as.character(unique(chemo$idsetuid))
hMat=matrix(0,nrow=nrow(chemoFirst),ncol=nrow(chemoFirst),dimnames=list(ps,ps))
for(i in 1:nrow(chemoFirst)){
  for(j in i:nrow(chemoFirst)){
    num=length(intersect(hourList[[i]],hourList[[j]]))
    #The matrix is reciprocal
    hMat[i,j]=num
    hMat[j,i]=num
    #printout to make sure it is working, and how fast it is going
    print(paste(i,j,sep='-'))
    flush.console()
  }
}

#Sparse matrix of hours in the ward
hMat2=Matrix(hMat,sparse=T)

#Makes a matrix of Jaccard's Indices
jMat=matrix(0,nrow=nrow(hMat),ncol=ncol(hMat),dimnames=list(ps,ps))

riskset=list()
for(i in 1:length(hourList)) {
  riskset[[i]]=c(min(hourList[[i]]),max(hourList[[i]]))
}
names(riskset)=names(hourList)
for(i in 1:nrow(hMat)){
  for(j in (i+1):ncol(hMat)){
    jMat[i,j]=length(intersect(hourList[[ps[i]]],hourList[[ps[j]]]))/max(1,length(which((union(hourList[[ps[i]]],hourList[[ps[j]]])>=max(riskset[[ps[i]]][1],riskset[[ps[j]]][1]))&(union(hourList[[ps[i]]],hourList[[ps[j]]])<=min(riskset[[ps[i]]][2],riskset[[ps[j]]][2])))))
    jMat[j,i]=jMat[i,j]
  }
  print(i)
}
diag(jMat)=1
jacNet=network(jMat,matrix.type="adjacency",directed=F)
set.edge.value(jacNet,"weight",jMat)

#Function to make cutoffs for every time, to reduce computation.  
#This si used to determine who can potentially overlap with whom
#Essentially, for each possible hour, it lists patients that had been to the chemotherapy ward anytime within +/- 32 hours
#This is used as the risk set of overlap for patients in the simulation function below
breaker=function(){
  breaks=list()
  for(i in 1:(ncol(sMat))){
    if (length(which(rowSums(sMat[,max(0,(i-32)):min(ncol(sMat),(i+32))])>0))==0) breaks[[i]]=0
    else breaks[[i]]=which(rowSums(sMat[,max(0,(i-32)):min(ncol(sMat),(i+32))])>0)
    print(i)
    flush.console()
  }
  breaks
}

breaks=breaker()

popHours=colSums(sMat)
#Determine the distribution of overlap times for every person

#Do a simulation with 10000 trials, which takes about 2 days.  100 takes about 20 minutes, 1000 takes about 2 hours
#enc=encounters(10)

#Create vector of cutoff numbers for each person
#vec99=sapply(enc,pct,perct=0.99)

#Create the networks and add basic attributes
#el99=overlapEdgelist(hMat2,vec99)
#net99=as.network(el99,mode="edgelist",directed=F)
#delete.vertices(net99,4692)
#set.edge.value(net99,"hours",hMat2)
#set.vertex.attribute(net99,"start",chemoFirst$WardStart)
#network.vertex.names(net99)=as.character(chemoFirst$idsetuid)

#METHODS FOR DETERMINING OUTCOMES

#In Seconds, 6 months, and 5 years
tdif=15552000
tsurv=31556900*5
events=list()

#Looking at this is pretty painful.  
#Basically, it takes the dataframe of observations for each person, and determines good and bad outcomes.
#The first row is always a "start" ovservation, followed by an outcome
#If a person has more than one outcome, they have multiple sets of "start" observations followed by outcomes.
#An outcome occurs when there is at least 6 months between consecutive chemotherapy visits, and the otcome is determined by 5-year survival
#Most of the code is to account for all possibilities of the above

for(i in 1:nrow(chemoFirst)){
  patient=chemo[chemo$idsetuid==chemoFirst$idsetuid[i],]
  pevent=data.frame(id=patient$idsetuid[1],date=as.POSIXct(patient$WardStartDate[1]),type="start")
  if (nrow(patient)>2){
    for (j in 2:(nrow(patient)-1)){
      #At least 6 months between chemotherapy visits
      if(as.numeric(as.POSIXct(patient$WardStartDate[(j+1)])>(tdif+as.POSIXct(patient$WardEndDate[j])))) {
        if(patient$deathdate[1]=="NULL") pevent=rbind(pevent,data.frame(id=patient$idsetuid[j],date=as.POSIXct(patient$WardEndDate[j]),type="goodend"))
        else if(as.numeric(difftime(as.POSIXct(patient$deathdate[1]),as.POSIXct(patient$WardEndDate[j]),units="days"))>(365*5))  pevent=rbind(pevent,data.frame(id=patient$idsetuid[j],date=as.POSIXct(patient$WardEndDate[j]),type="goodend"))
        else pevent=rbind(pevent,data.frame(id=patient$idsetuid[j],date=as.POSIXct(patient$WardEndDate[j]),type="badend"))
      }
      if(as.numeric(as.POSIXct(patient$WardStartDate[j])>(tdif+as.POSIXct(patient$WardEndDate[(j-1)])))) pevent=rbind(pevent,data.frame(id=patient$idsetuid[j],date=as.POSIXct(patient$WardEndDate[j]),type="start"))
    }
  }
  if (nrow(patient)==2){
    if(as.numeric(as.POSIXct(patient$WardStartDate[(2)])>(tdif+as.numeric(as.POSIXct(patient$WardEndDate[1]))))){ 
      if(patient$deathdate[1]=="NULL") pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=c(as.POSIXct(patient$WardEndDate[1])),type=c("goodend")))
      else if(difftime(as.POSIXct(patient$deathdate[1]),as.POSIXct(patient$WardEndDate[1]),units="days")>(365*5)) pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=c(as.POSIXct(patient$WardEndDate[1])),type=c("goodend")))
      else if(difftime(as.POSIXct(patient$deathdate[1]),as.POSIXct(patient$WardEndDate[1]),units="days")<(365*5)) pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=c(as.POSIXct(patient$WardEndDate[1])),type=c("badend")))
    }
    else if (patient$deathdate[1]=="NULL") pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=as.POSIXct(patient$WardEndDate[2]),type=c("goodend")))
    else if (difftime(as.POSIXct(patient$deathdate[1]),as.POSIXct(patient$WardEndDate[2]))>tsurv) pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=as.POSIXct(patient$WardEndDate[2]),type=c("goodend")))
    else pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=as.POSIXct(patient$WardEndDate[2]),type=c("badend")))
  }
  if (patient$deathdate[1]!="NULL") {
    if ((as.POSIXct(patient$WardEndDate[nrow(patient)])+tsurv)<as.POSIXct(patient$deathdate[1])) pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=(as.POSIXct(patient$WardEndDate[nrow(patient)])),type="goodend")) 
    else pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=(as.POSIXct(patient$WardEndDate[nrow(patient)])),type="badend"))
  }
  else pevent=rbind(pevent,data.frame(id=patient$idsetuid[1],date=(as.POSIXct(patient$WardEndDate[nrow(patient)])),type="goodend"))
  events[[i]]=pevent
  print(i)
  flush.console()
}

#Puts data in data frame suitable for GEE regression
#Essentially removes the "start" observations, and includes some variables like the total time in the chemotherapy ward
make.logistic=function(x,df=chemo){
  data=df[df$idsetuid==x$id[1],]
  out=data.frame(id=x$id[1],stime=as.POSIXct(x$date[1]),etime=as.POSIXct(x$date[2]),time=(as.numeric(difftime(as.POSIXct(x$date[2]),as.POSIXct(x$date[1]),units="days"))/365),sexM=ifelse(data$sex[1]=="M",1,0),age=(as.numeric(difftime(as.POSIXct(x$date[1]),as.POSIXct(data$birthdate[1])),units="days")/365),outcome.Bad=ifelse(x$type[2]=="badend",1,0))
  if (nrow(x)>2){
    for(i in 3:nrow(x)){
      if (x$type[i]=="start") out=rbind(out,data.frame(id=x$id[i],stime=as.POSIXct(x$date[i]),etime=as.POSIXct(x$date[i+1]),time=(as.numeric(difftime(as.POSIXct(x$date[i+1]),as.POSIXct(x$date[i]),units="days"))/365),sexM=ifelse(data$sex[i]=="M",1,0),age=(as.numeric(difftime(as.POSIXct(x$date[i]),as.POSIXct(data$birthdate[i])),units="days")/365),outcome.Bad=ifelse(x$type[i+1]=="badend",1,0)))
    }
  }
  out
}

#Make one dataset of all outcome observations
dfs=lapply(events,make.logistic)
events.df=do.call(rbind,dfs)
events.df$id=as.character(events.df$id)

for (i in 1:nrow(events.df)){
  events.df$totHours[i]=diag(hMat)[events.df$id[i]]
  events.df$olHours[i]=rowSums(hMat)[events.df$id[i]]-events.df$totHours[i]
}

#Adds in dripping time and number of trips as variables
driptrip=function(x){
  for (i in 1:nrow(x)){
    person=chemo[chemo$idsetuid==x$id[i],]
    person=person[as.POSIXct(person$WardStartDate)>=as.POSIXct(x$stime[i]),]
    person=person[as.POSIXct(person$WardEndDate)<=as.POSIXct(x$etime[i]),]
    x$tripcount[i]=nrow(person)
    x$avgdrip[i]=sum(person$WardEnd-person$WardStart)/nrow(person)          
  }
  x
} 

driptrip(events.df)->events.df
#run through here prior to network simulation

