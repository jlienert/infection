# LINUX: cd /data/lienertjp/Data

options(width=as.integer(Sys.getenv("COLUMNS")))

library(dplyr)
library(magrittr)
library(lubridate)

#Code to fix pnmicro dataset
test=scan("pnmicro.txt",skip=2,sep=",",what="character")

tmat2=matrix(0,ncol=10,nrow=length(test))
colnames(tmat2)=c("idsetuid","accessionnumber","collectiondatetime","batchtestcode","bugcode","bugname","drugcode","drugname","susceptability","mic")

for(i in 1:length(test)){
  t=cbind(substr(test[i],1,50),substr(test[i],52,101),substr(test[i],103,125),substr(test[i],127,176),substr(test[i],178,185),substr(test[i],187,442),substr(test[i],444,451),substr(test[i],453,472),substr(test[i],474,487),substr(test[i],489,518))
  tmat2[i,]=t
  print(i)
}

write.table(tmat2,file="pnmic.csv",sep=",")

#Code to fix pnpas dataset
con=file("pnpas.txt")
open(con)
test=readLines(con)

tmat2=matrix(0,ncol=18,nrow=length(test))
colnames(tmat2)=c("idsetuid","source","spellid","AdmissionDate","DischargeDate","WardStartDate","WardEndDate","WardName","Bed","Bay","Building","Facility","AdminCategoryCode","AdmissionMethodCode","AdmissionSourceCode","DischargeDestinationCode","DischargeMethodCode","MainSpecialtyCode")

for(i in 1:length(test)){
  t=cbind(substr(test[i],1,50),substr(test[i],52,69),substr(test[i],71,326),substr(test[i],328,350),substr(test[i],352,374),substr(test[i],376,398),substr(test[i],400,422),substr(test[i],424,463),substr(test[i],465,720),substr(test[i],722,761),substr(test[i],763,802),substr(test[i],804,843),substr(test[i],845,944),substr(test[i],946,1045),substr(test[i],1047,1146),substr(test[i],1148,1247),substr(test[i],1249,1348),substr(test[i],1350,1449))
  tmat2[i,]=t
  print(i)
}

write.table(tmat2,file="pnpas.csv",sep=",",row.names=F)

#Code to fix pndemographics dataset
con=file("pndemographics.txt")
open(con)
test=readLines(con)

tmat2=matrix(0,ncol=5,nrow=length(test))
colnames(tmat2)=c("patientid","birthdate","deathdate","sex","idsetuid")

for(i in 1:length(test)){
  t=cbind(substr(test[i],1,11),substr(test[i],13,35),substr(test[i],37,59),substr(test[i],61,110),substr(test[i],112,143))
  tmat2[i,]=t
  print(i)
}

write.table(tmat2,file="pndemo.csv",sep=",",row.names=F)

#Code to fix wards dataset
test=scan("wards.csv",skip=1,sep=",",what="character")

tmat2=matrix(0,ncol=8,nrow=length(test))
colnames(tmat2)=c("WardName","ward","group","speciality","directorate","hospital","wardcat","nonoxf")

write.table(tmat2,file="pnmic.csv",sep=",")

#Code to fix pnlims dataset
test=scan("pnlims.rpt",skip=2,sep=",",what="character")

tmat2=matrix(0,ncol=8,nrow=length(test))
colnames(tmat2)=c("idsetuid","collectiondatetime","testname","minrange","maxrangeunits","value","numbervalue")

for(i in 1:length(test)){
  t=cbind(substr(test[i],1,50),substr(test[i],52,73),substr(test[i],75,124),substr(test[i],126,147),substr(test[i],149,170),substr(test[i],172,221),substr(test[i],223,272),substr(test[i],274,295))
  tmat2[i,]=t
  print(i)
}

write.table(tmat2,file="pnlims.csv",sep=",")



##########LINUX #################
sed -i -e 's/  //g' pnmic.csv
sed -i -e 's/" /"/g' pnmic.csv
sed -i -e 's/ "/"/g' pnmic.csv
sed -i -e 's/  //g' pnpas.csv
sed -i -e 's/" /"/g' pnpas.csv
sed -i -e 's/ "/"/g' pnpas.csv
sed -i -e 's/  //g' pndemo.csv
sed -i -e 's/" /"/g' pndemo.csv
sed -i -e 's/ "/"/g' pndemo.csv
sed -i -e 's/  //g' pnlims.csv
sed -i -e 's/" /"/g' pnlims.csv
sed -i -e 's/ "/"/g' pnlims.csv
##########END LINUX##############

#CODE FOR FINALIZING DATASETS, including variable recoding

pnmic=read.table("pnmic.csv",sep=",",stringsAsFactors=F,head=T)
pnpas=read.table("pnpas.csv",sep=",",stringsAsFactors=F,head=T)
pndemo=read.table("pndemo.csv",sep=",",stringsAsFactors=F,head=T)
icd=read.table("icd.csv",sep=",",stringsAsFactors=F,head=T)
pnlims=read.table("pnlims.csv",sep=",",stringsAsFactors=F,head=T)
names(pnlims)=c("idsetuid","collectiondatetime","testname","minrange","maxrangeunits","value","numbervalue")

dpas=merge(pnpas,pndemo)
dpas=rename(dpas,idsetuid.pas=idsetuid,birthdate=brithdate)
dpas=dpas[grepl("[a-zA-Z0-9]",dpas$idsetuid),]
dpas=dpas[-which(dpas$WardStartDate=="NULL"),]
dpas=dpas[-which(dpas$WardEndDate=="NULL"),]
dpas=dpas[-which(dpas$birthdate=="NULL"),]
dpas=dpas[-which(dpas$DischargeDate=="NULL"),]
dpas$WardStartDate %<>% as.POSIXct
dpas$WardEndDate %<>% as.POSIXct
dpas$AdmissionDate %<>% as.POSIXct
dpas$DischargeDate %<>% as.POSIXct
dpas$birthdate %<>% as.POSIXct
dpas<-transform(dpas,interv=as.interval(WardStartDate,WardEndDate))
dpas<-transform(dpas,interv2=as.interval(WardStartDate-hours(1),WardEndDate+hours(1)))

dmic=merge(pndemo,pnmic)
dmic=rename(dmic,idsetuid.mic=idsetuid,birthdate=brithdate)
dmic=dmic[-which(dmic$birthdate=="NULL"),]
dmic=dmic[-which(dmic$bugname=="NULL"),]
#dmic=dmic[-which(dmic$collectiondatetime=="NULL"),]
dmic$birthdate %<>% as.POSIXct
dmic$collectiondatetime %<>% as.POSIXct
dmic=dmic[dmic$patientid %in% dpas$patientid,]
dmic$spellid="0"

dlims=merge(pnlims,pndemo)
dlims=rename(dlims,idsetuid.lim=idsetuid,birthdate=brithdate)
dlims=dlims[grepl("[a-zA-Z0-9]",dlims$idsetuid.lim),]
dlims=dlims[-which(dlims$birthdate=="NULL"),]
dlims$birthdate %<>% as.POSIXct
dlims$collectiondatetime %<>% as.POSIXct
dlims=dlims[dlims$patientid %in% dpas$patientid,]
dlims$spellid="0"

dmic2=select(dmic,patientid,accessionnumber,collectiondatetime,batchtestcode,bugcode,bugname,drugcode,drugname,susceptability,mic,spellid)
dlims2=select(dlims,collectiondatetime,testname,minrange,maxrangeunits,value,numbervalue,patientid)

write.table(dpas,file="dpas.csv",sep=",")
write.table(dmic2,file="dmic.csv",sep=",")
write.table(dlims2,file="dlims.csv",sep=",")

######################CURRENTLY UNUSED
#Split pas and mic by unique patientd
people=split(dpas,as.factor(dpas$patientid))
for(i in 1:length(people)) names(people)[i]=people[[i]]$patientid[1]
tests=split(dmic,as.factor(dmic$patientid))
for(i in 1:length(tests)) names(tests)[i]=tests[[i]]$patientid[1]

#For each mic test, assign spellids to tests
for(i in 1:length(tests)) {
  temp.pas=people[names(people)==names(tests)[i]][[1]]
  for (j in 1:nrow(tests[[i]])){
    tests[[i]]$spellid[j]=temp.pas$spellid[tests[[i]]$collectiondatetime[j] %within% temp.pas$interv][1]
  }
  print(i)
}

dmic.sid=rbindlist(tests)

write.table(dmic.sid,file="dmspell.csv",sep=",",row.names=F)

#Extend window by an hour on each side
for(i in 1:length(tests)) {
  temp.pas=people[names(people)==names(tests)[i]][[1]]
  for (j in 1:nrow(tests[[i]])){
    tests[[i]]$spellid[j]=temp.pas$spellid[tests[[i]]$collectiondatetime[j] %within% temp.pas$interv2][1]
  }
  print(i)
}

dmic.sid2=rbindlist(tests)

write.table(dmic.sid2,file="dmspell2.csv",sep=",",row.names=F)
