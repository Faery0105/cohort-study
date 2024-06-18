library(WeightIt)
library(Survival)
library(cmprsk)
library(InformativeCensoring)

cancer_dx<-dxall[grepl("^1[4-9][0-9]|^20[0-9]",dxall$dx.code),]
check_cohort<-merge(demographic,cancer_dx,by="reference.key",all = F)
colnames(check_cohort)
halfyear<-check_cohort[which(check_cohort$date.dx<=check_cohort$date.vte&check_cohort$date.dx>check_cohort$date.vte-180),]
beforehalf<-check_cohort[which(check_cohort$date.dx<=check_cohort$date.vte-180),]

me1<-halfyear[grepl("^19[6-9]",halfyear$dx.code),]
rk1<-unique(me1$reference.key)

halfyear$check<-paste(halfyear$reference.key,halfyear$dx.code)
beforehalf$check<-paste(beforehalf$reference.key,beforehalf$dx.code)
check<-setdiff(halfyear$check,beforehalf$check)
rk2<-as.numeric(unique(word(check,1)))

cancer_drug<-rxall[grepl("^8",rxall$drug.code),]
drug_cohort<-merge(demographic,cancer_drug,by="reference.key",all = F)
treatment1<-drug_cohort[which(drug_cohort$date.rx.start<=drug_cohort$date.vte&drug_cohort$date.rx.start>=drug_cohort$date.vte-180),]
treatment2<-drug_cohort[which(drug_cohort$date.rx.start<drug_cohort$date.vte-180&drug_cohort$date.rx.end>=drug_cohort$date.vte-180),]
rk3_1<-unique(treatment1$reference.key)
rk3_2<-unique(treatment2$reference.key)
rk3<-union(rk3_1,rk3_2)

dx_tr<-dxall[grepl("^V58.1|^V66.2|^V67.2|^V58.0|^V66.1|^V67.1",dxall$dx.code),]
dx_tr_cohort<-merge(demographic,dx_tr,by="reference.key",all = F)
dx<-dx_tr_cohort[which(dx_tr_cohort$date.dx<=dx_tr_cohort$date.vte&
                         dx_tr_cohort$date.dx>dx_tr_cohort$date.vte-180),]
rk4<-unique(dx$reference.key)

px_tr_1<-pxall[grepl("^99.25|^99.28|^92.2|^92.3",pxall$px.code),]
px_tr_2<-pxall[grepl("chemotherapy|BRM|radiotherapy",pxall$detail,ignore.case = T),]
px_tr<-rbind(px_tr_1,px_tr_2)
px_tr_cohort<-merge(demographic,px_tr,by="reference.key",all = F)
px<-px_tr_cohort[which(px_tr_cohort$date.px<=px_tr_cohort$date.vte&
                         px_tr_cohort$date.px>px_tr_cohort$date.vte-180),]
rk5<-unique(px$reference.key)

ip_ca<-ipall[grepl("^1[4-9][0-9]|^20[0-9]",ipall$dxpx.code),]
ip_cat<-merge(demographic,ip_ca,by="reference.key",all = F)
ip1<-ip_cat[which(ip_cat$date.ip.admission<=ip_cat$date.vte&ip_cat$date.ip.admission>=ip_cat$date.vte-180),]
ip2<-ip_cat[which(ip_cat$date.ip.admission<ip_cat$date.vte-180&ip_cat$date.ip.discharge>=ip_cat$date.vte-180),]
rk6_1<-unique(ip1$reference.key)
rk6_2<-unique(ip2$reference.key)
rk6<-union(rk6_1,rk6_2)

rk_active<-union(union(union(rk1,rk2),union(rk3,rk4)),union(rk5,rk6))
cat_cohort<-demographic[which(demographic$reference.key %in% rk_active),]

cat_drug<-rxall[grepl("^2.8",rxall$drug.code),]
anticoagulant<-cat_drug[which(cat_drug$reference.key %in% rk_active),]
length(unique(anticoagulant$reference.key))

lmwh<-anticoagulant[grepl("ENOXAPARIN|NADROPARIN|TINZAPARIN",anticoagulant$drug.name),c(1:3)]
lmwh<-arrange(lmwh,lmwh$reference.key,lmwh$date.rx.start,lmwh$date.rx.end)
length(unique(lmwh$reference.key))
noac<-anticoagulant[grepl("APIXABAN|DABIGATRAN|EDOXABAN|RIVAROXABAN",anticoagulant$drug.name),c(1:3)]
noac<-arrange(noac,noac$reference.key,noac$date.rx.start,noac$date.rx.end)
length(unique(noac$reference.key))
rk_list<-unique(lmwh$reference.key)

lmwh_new<-data.frame(matrix(NA,ncol = ncol(lmwh)))
colnames(lmwh_new)<-colnames(lmwh)

for (i in c(1:length(rk_list))) {
  rk<-rk_list[i]
  lmwh_i<-lmwh[which(lmwh$reference.key==rk),]
  b<-NA
  if(nrow(lmwh_i)>1)for (a in c(nrow(lmwh_i):2)) {
    if(lmwh_i[a-1,"date.rx.end"]+30>=lmwh_i[a,"date.rx.start"]){
      lmwh_i[a-1,"date.rx.end"]<-lmwh_i[a,"date.rx.end"]
      b<-union(b,a)
    }
  }
  if(nrow(lmwh_i)>1&length(b)>1){
    b<-b[-1]
    lmwh_new_i<-lmwh_i[-b,]
  }else{lmwh_new_i<-lmwh_i}
  lmwh_new<-rbind(lmwh_new,lmwh_new_i)
}

lmwh_new<-lmwh_new[-1,]
lmwh_new$date.rx.start<-as.Date(lmwh_new$date.rx.start,origin = "1970-01-01")
lmwh_new$date.rx.end<-as.Date(lmwh_new$date.rx.end,origin = "1970-01-01")
colnames(lmwh_new)[c(2,3)]<-c("date.lmwh.s","date.lmwh.e")
setwd("C:/Users/LabPC12CSMPR/Desktop")
save(lmwh_new,file = "lmwh_new.RData")

noac_new<-data.frame(matrix(NA,ncol = ncol(noac)))
colnames(noac_new)<-colnames(noac)

for (i in c(1:length(rk_list))) {
  rk<-rk_list[i]
  noac_i<-noac[which(noac$reference.key==rk),]
  b<-NA
  if(nrow(noac_i)>1)for (a in c(nrow(noac_i):2)) {
    if(noac_i[a-1,"date.rx.end"]+30>=noac_i[a,"date.rx.start"]){
      noac_i[a-1,"date.rx.end"]<-noac_i[a,"date.rx.end"]
      b<-union(b,a)
    }
  }
  if(nrow(noac_i)>1&length(b)>1){
    b<-b[-1]
    noac_new_i<-noac_i[-b,]
  }else{noac_new_i<-noac_i}
  noac_new<-rbind(noac_new,noac_new_i)
}

noac_new<-noac_new[-1,]
noac_new$date.rx.start<-as.Date(noac_new$date.rx.start,origin = "1970-01-01")
noac_new$date.rx.end<-as.Date(noac_new$date.rx.end,origin = "1970-01-01")
colnames(noac_new)[c(2,3)]<-c("date.noac.s","date.noac.e")
setwd("C:/Users/LabPC12CSMPR/Desktop")
save(noac_new,file = "noac_new.RData")

check_lmwh<-merge(cat_cohort[,c("reference.key","date.vte")],lmwh_new,by="reference.key",all.x = T)
lmwh_cohort<-arrange(check_lmwh,reference.key,date.lmwh.s)
lmwh_cohort<-lmwh_cohort[which(lmwh_cohort$date.lmwh.s>=lmwh_cohort$date.vte),]
lmwh_cohort<-lmwh_cohort[!duplicated(lmwh_cohort$reference.key),]
cat_cohort<-merge(cat_cohort,lmwh_cohort[,c("reference.key","date.lmwh.s","date.lmwh.e")],by="reference.key",all.x = T)
colnames(cat_cohort)[ncol(cat_cohort)-1]<-"date.lmwh"

noac_new<-arrange(noac_new,noac_new$reference.key,noac_new$date.noac.s)
noac.index<-noac_new[!duplicated(noac_new$reference.key),]
rk.noac<-unique(noac.index$reference.key)
colnames(noac.index)[2]<-"date.noac"
cat_cohort<-merge(cat_cohort,noac.index,by="reference.key",all.x = T)
cat_cohort[which(cat_cohort$date.noac<=cat_cohort$date.lmwh),]$date.noac.e<-NA
cat_cohort[which(cat_cohort$date.noac<=cat_cohort$date.lmwh),]$date.noac<-NA
nrow(cat_cohort[which(is.na(cat_cohort$date.noac)==F),])

setwd("C:/Users/LabPC12CSMPR/Desktop")
save(cat_cohort,file = "cat_cohort.RData")

study_cohort<-cat_cohort

study_cohort$group<-0
study_cohort[which(is.na(study_cohort$date.noac)==F),]$group<-1
study_cohort$group<-as.factor(study_cohort$group)

study_cohort$age<-as.numeric(floor((as.Date(study_cohort$date.vte)-as.Date(study_cohort$date.birth))/365.25))
ex1<-study_cohort[which(is.na(study_cohort$sex)==T|is.na(study_cohort$age)==T),]$reference.key

ex2<-study_cohort[which(study_cohort$age<18),]$reference.key

ex3<-study_cohort[which(is.na(study_cohort$date.lmwh)==T),]$reference.key

check_lmwh<-merge(study_cohort[,c("reference.key","date.vte")],lmwh,by="reference.key",all.x = T)
ex4_1<-check_lmwh[which(check_lmwh$date.rx.start<check_lmwh$date.vte&check_lmwh$date.rx.start>=check_lmwh$date.vte-180),]
ex4_2<-check_lmwh[which(check_lmwh$date.rx.start<check_lmwh$date.vte-180&check_lmwh$date.rx.end>=check_lmwh$date.vte-180),]
ex4_cohort<-rbind(ex4_1,ex4_2)
ex4<-unique(ex4_cohort$reference.key)

check_anti<-merge(study_cohort[,c("reference.key","date.lmwh")],cat_drug,by="reference.key",all.x = T)
ex5_1<-check_anti[which(check_anti$date.rx.start<check_anti$date.lmwh&check_anti$date.rx.start>=check_anti$date.lmwh-180),]
ex5_2<-check_anti[which(check_anti$date.rx.start<check_anti$date.lmwh-180&check_anti$date.rx.end>=check_anti$date.lmwh-180),]
ex5_cohort<-rbind(ex5_1,ex5_2)
ex5<-unique(ex5_cohort$reference.key)

ex6<-study_cohort[which(study_cohort$date.lmwh>study_cohort$date.vte+90),]$reference.key

ex<-union(union(union(ex1,ex2),union(ex3,ex4)),union(ex5,ex6))
study_cohort<-study_cohort[-which(study_cohort$reference.key %in% ex),]

nrow(study_cohort[which(study_cohort$group==1),])
nrow(study_cohort[which(study_cohort$group==0),])
save(study_cohort,file = "study_cohort(after 1st exclusion).RData")

study_cohort<-merge(study_cohort[,c(1:15)],old[,c("reference.key","switch.days")],by="reference.key",all.x = T)
study_cohort[which(is.na(study_cohort$switch.days)==T),]$switch.days<-as.numeric(as.Date(study_cohort[which(is.na(study_cohort$switch.days)==T),]$date.noac))-as.numeric(as.Date(study_cohort[which(is.na(study_cohort$switch.days)==T),]$date.lmwh))
study_cohort[which(study_cohort$switch.days<=0),]$switch.days<-NA
mean(study_cohort[which(study_cohort$group==1),]$switch.days,na.rm = T)
sd(study_cohort[which(study_cohort$group==1),]$switch.days,na.rm = T)
median(study_cohort[which(study_cohort$group==1),]$switch.days,na.rm = T)

set.seed(10)
study_cohort[which(is.na(study_cohort$switch.days)==T),]$switch.days<-sample(study_cohort[which(study_cohort$group==1),]$switch.days,size =nrow(study_cohort[which(is.na(study_cohort$switch.days)==T),]),replace = T )
mean(study_cohort[which(study_cohort$group==0),]$switch.days,na.rm = T)
sd(study_cohort[which(study_cohort$group==0),]$switch.days,na.rm = T)
median(study_cohort[which(study_cohort$group==0),]$switch.days,na.rm = T)
study_cohort$index.date<-study_cohort$date.noac
study_cohort[which(study_cohort$group==0),]$index.date<-study_cohort[which(study_cohort$group==0),]$date.lmwh+study_cohort[which(study_cohort$group==0),]$switch.days
study_cohort$index.date<-as.Date(study_cohort$index.date,origin = "1970-01-01")
study_cohort$index.year<-year(study_cohort$index.date)
mean(study_cohort$switch.days)
sd(study_cohort$switch.days)
study_cohort$dur.rx<-NA
study_cohort[which(study_cohort$group==0),]$dur.rx<-as.numeric(study_cohort[which(study_cohort$group==0),]$date.lmwh.e-study_cohort[which(study_cohort$group==0),]$date.lmwh)
study_cohort[which(study_cohort$group==1),]$dur.rx<-as.numeric(study_cohort[which(study_cohort$group==1),]$date.noac.e-study_cohort[which(study_cohort$group==1),]$date.lmwh)
study_cohort$dur.rx.1<-NA
study_cohort[which(study_cohort$group==0),]$dur.rx.1<-as.numeric(study_cohort[which(study_cohort$group==0),]$index.date-study_cohort[which(study_cohort$group==0),]$date.lmwh)
study_cohort[which(study_cohort$group==1),]$dur.rx.1<-as.numeric(study_cohort[which(study_cohort$group==1),]$index.date-study_cohort[which(study_cohort$group==1),]$date.lmwh)
study_cohort$dur.rx.2<-NA
study_cohort[which(study_cohort$group==0),]$dur.rx.2<-as.numeric(study_cohort[which(study_cohort$group==0),]$date.lmwh.e-study_cohort[which(study_cohort$group==0),]$index.date)
study_cohort[which(study_cohort$group==1),]$dur.rx.2<-as.numeric(study_cohort[which(study_cohort$group==1),]$date.noac.e-study_cohort[which(study_cohort$group==1),]$index.date)

ex1<-study_cohort[which(study_cohort$index.date>as.Date("2022-12-31")|study_cohort$index.date<study_cohort$date.vte),]$reference.key

dx_his<-dxall[grepl("^286|^287.1|^287.[3-5]|^44[45]",
                    dxall$dx.code),]
check_dx<-merge(study_cohort[,c("reference.key","index.date")],dx_his,by="reference.key",all = F)
ex2_cohort<-check_dx[which(check_dx$date.dx<check_dx$index.date),]
ex2<-unique(ex2_cohort$reference.key)

ivcf<-pxall[grepl("^38.7|^39.74",pxall$px.code),]
ivcf<-merge(study_cohort[,c("reference.key","index.date")],ivcf,by="reference.key",all = F)
ivcf<-ivcf[which(ivcf$index.date<ivcf$date.px),]
ex3<-unique(ivcf$reference.key)

pregnant<-dxall[grepl("^6[3-7][0-9]",dxall$dx.code),]
pregnant_study<-merge(study_cohort[,c("reference.key","index.date")],pregnant,by="reference.key",all=F)
pregnant_study<-pregnant_study[which(pregnant_study$date.dx>=pregnant_study$index.date),]
ex4<-pregnant_study$reference.key

ex5<-study_cohort[which(study_cohort$date.death<study_cohort$index.date+7),]$reference.key

ex6<-study_cohort[which(study_cohort$date.lmwh.e<study_cohort$index.date&study_cohort$group==0),]$reference.key

ex<-union(union(union(ex1,ex2),union(ex3,ex4)),union(ex5,ex6))
study_cohort<-study_cohort[-which(study_cohort$reference.key %in% ex),]

nrow(study_cohort[which(study_cohort$group==1),])
nrow(study_cohort[which(study_cohort$group==0),])
summary(study_cohort[which(study_cohort$group==0),]$dur.rx.1)
summary(study_cohort[which(study_cohort$group==1),]$dur.rx.1)
summary(study_cohort[which(study_cohort$group==0),]$dur.rx.2)
summary(study_cohort[which(study_cohort$group==1),]$dur.rx.2)
summary(study_cohort[which(study_cohort$group==0),]$dur.rx)
summary(study_cohort[which(study_cohort$group==1),]$dur.rx)

study_cohort$start.tr<-as.numeric(study_cohort$date.lmwh-study_cohort$date.vte)
summary(study_cohort[which(study_cohort$group==0),]$start.tr)
summary(study_cohort[which(study_cohort$group==1),]$start.tr)

cohort<-merge(study_cohort[,c("reference.key","index.date")],ipall,by="reference.key",all.x = T)
cohort_after<-cohort[which(cohort$date.ip.admission>=cohort$index.date+14),]
cohort_re_vte<-cohort_after[grepl("^415.1|^453.[1245789]",cohort_after$dxpx.code),]
cohort_re_vte<-cohort_re_vte[order(cohort_re_vte$reference.key,cohort_re_vte$date.ip.admission),]
cohort_ip1<-cohort_re_vte[grepl("^D1",cohort_re_vte$dxpx.rank),]
cohort_ip1<-arrange(cohort_ip1,reference.key,date.ip.admission)
cohort_ip1_vte<-cohort_ip1[!duplicated(cohort_ip1$reference.key),]
study_cohort<-merge(study_cohort,cohort_ip1_vte[,c("reference.key","date.ip.admission")],by="reference.key",all.x = T)
colnames(study_cohort)[ncol(study_cohort)]<-"outcome.vte.date"

cohort<-merge(study_cohort[,c("reference.key","index.date")],ipall,by="reference.key",all.x = T)
cohort_after<-cohort[which(cohort$date.ip.admission>=cohort$index.date+14),]
cohort_re_dvt<-cohort_after[grepl("^453.[1245789]",cohort_after$dxpx.code),]
cohort_re_dvt<-cohort_re_dvt[order(cohort_re_dvt$reference.key,cohort_re_dvt$date.ip.admission),]
cohort_ip1<-cohort_re_dvt[grepl("^D1",cohort_re_dvt$dxpx.rank),]
cohort_ip1<-arrange(cohort_ip1,reference.key,date.ip.admission)
cohort_ip1_dvt<-cohort_ip1[!duplicated(cohort_ip1$reference.key),]
study_cohort<-merge(study_cohort,cohort_ip1_dvt[,c("reference.key","date.ip.admission")],by="reference.key",all.x = T)
colnames(study_cohort)[ncol(study_cohort)]<-"outcome.dvt.date"

cohort<-merge(study_cohort[,c("reference.key","index.date")],ipall,by="reference.key",all.x = T)
cohort_after<-cohort[which(cohort$date.ip.admission>=cohort$index.date+14),]
cohort_re_pe<-cohort_after[grepl("^415.1",cohort_after$dxpx.code),]
cohort_re_pe<-cohort_re_pe[order(cohort_re_pe$reference.key,cohort_re_pe$date.ip.admission),]
cohort_ip1<-cohort_re_pe[grepl("^D1",cohort_re_pe$dxpx.rank),]
cohort_ip1<-arrange(cohort_ip1,reference.key,date.ip.admission)
cohort_ip1_pe<-cohort_ip1[!duplicated(cohort_ip1$reference.key),]
study_cohort<-merge(study_cohort,cohort_ip1_pe[,c("reference.key","date.ip.admission")],by="reference.key",all.x = T)
colnames(study_cohort)[ncol(study_cohort)]<-"outcome.pe.date"

mb<-dxall[grepl("^43[0-2]|^53[1-4].[0246]|^535.[0-9]1|^537.83|^562.0[23]|^562.1[23]|^568.81|^569.3|^569.85|
                      ^578|^423.0|^459.0|^593.81|^599.7|^719.1|^784.7|^784.8|^786.3",dxall$dx.code),]
check_mb<-merge(study_cohort[,c("reference.key","index.date")],mb,by="reference.key",all = F)
outcome_mb<-check_mb[which(check_mb$date.dx>check_mb$index.date),]
mb2<-pxall[grepl("99.0",pxall$px.code),]
check_mb2<-merge(study_cohort[,c("reference.key","index.date")],mb2,by="reference.key",all = F)
outcome_mb2<-check_mb2[which(check_mb2$date.px>check_mb2$index.date),]
length(unique(outcome_mb2$reference.key))
rk2<-unique(outcome_mb2$reference.key)
rk1<-unique(outcome_mb$reference.key)
check<-intersect(rk1,rk2)

outcome_mb<-arrange(outcome_mb,reference.key,date.dx)
outcome_mb<-outcome_mb[!duplicated(outcome_mb$reference.key),]
colnames(outcome_mb)[3]<-"outcome.mb.date"
study_cohort<-merge(study_cohort,outcome_mb[,c(1,3)],by="reference.key",all.x = T)

ich<-dxall[grepl("^43[0-2]",dxall$dx.code),]
check_ich<-merge(study_cohort[,c("reference.key","index.date")],ich,by="reference.key",all = F)
outcome_ich<-check_ich[which(check_ich$date.dx>check_ich$index.date),]
outcome_ich<-arrange(outcome_ich,reference.key,date.dx)
outcome_ich<-outcome_ich[!duplicated(outcome_ich$reference.key),]
colnames(outcome_ich)[3]<-"outcome.ich.date"
study_cohort<-merge(study_cohort,outcome_ich[,c(1,3)],by="reference.key",all.x = T)

gib<-dxall[grepl("^53[1-4].[0246]|^535.[0-9]1|^537.83|^562.0[23]|^562.1[23]|^568.81|^569.3|^569.85|
                      ^578",dxall$dx.code),]
check_gib<-merge(study_cohort[,c("reference.key","index.date")],gib,by="reference.key",all = F)
outcome_gib<-check_gib[which(check_gib$date.dx>check_gib$index.date),]
outcome_gib<-arrange(outcome_gib,reference.key,date.dx)
outcome_gib<-outcome_gib[!duplicated(outcome_gib$reference.key),]
colnames(outcome_gib)[3]<-"outcome.gib.date"
study_cohort<-merge(study_cohort,outcome_gib[,c(1,3)],by="reference.key",all.x = T)

omb<-dxall[grepl("^423.0|^459.0|^593.81|^599.7|^719.1|^784.7|^784.8|^786.3",dxall$dx.code),]
check_omb<-merge(study_cohort[,c("reference.key","index.date")],omb,by="reference.key",all = F)
outcome_omb<-check_omb[which(check_omb$date.dx>check_omb$index.date),]
outcome_omb<-arrange(outcome_omb,reference.key,date.dx)
outcome_omb<-outcome_omb[!duplicated(outcome_omb$reference.key),]
colnames(outcome_omb)[3]<-"outcome.omb.date"
study_cohort<-merge(study_cohort,outcome_omb[,c(1,3)],by="reference.key",all.x = T)

study_cohort$outcome.death.date<-study_cohort$date.death

df_cox<-study_cohort
df_cox<- df_cox %>% 
  mutate(censor.date = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up = censor.date - index.date)
mean(as.numeric(df_cox$follow.up))
sd(as.numeric(df_cox$follow.up))
summary(as.numeric(df_cox[which(df_cox$group==0),]$follow.up))
summary(as.numeric(df_cox[which(df_cox$group==1),]$follow.up))

mean(study_cohort$age)
sd(study_cohort$age)

setwd("C:/Users/LabPC12CSMPR/Desktop")
save(study_cohort,file = "study_cohort.RData")

cohort<-study_cohort
all<-colnames(cohort)
a<-c(colnames(select(cohort, starts_with("outcome"))))
k<-sapply(a,grep,all)

cohort[which(cohort$outcome.vte.date>cohort$index.date+365),][,"outcome.vte.date"]<-NA
cohort[which(cohort$outcome.dvt.date>cohort$index.date+365),][,"outcome.dvt.date"]<-NA
cohort[which(cohort$outcome.pe.date>cohort$index.date+365),][,"outcome.pe.date"]<-NA
cohort[which(cohort$outcome.mb.date>cohort$index.date+365),][,"outcome.mb.date"]<-NA
cohort[which(cohort$outcome.ich.date>cohort$index.date+365),][,"outcome.ich.date"]<-NA
cohort[which(cohort$outcome.gib.date>cohort$index.date+365),][,"outcome.gib.date"]<-NA
cohort[which(cohort$outcome.omb.date>cohort$index.date+365),][,"outcome.omb.date"]<-NA
cohort[which(cohort$outcome.death.date>cohort$index.date+365),][,"outcome.death.date"]<-NA

save(cohort,file = "study_cohort(outcome events in 1 year).RData")

for (i in k) {
  cohort[,i]<-as.character(cohort[,i])
  if(nrow(cohort[which(is.na(cohort[,i])==F),])>0){cohort[which(is.na(cohort[,i])==F),][,i]<-"1"}
  if(nrow(cohort[which(is.na(cohort[,i])==T),])>0){cohort[which(is.na(cohort[,i])==T),][,i]<-"0"}
  cohort[,i]<-as.numeric(cohort[,i])
  cohort[,i]<-as.factor(cohort[,i])
}
vars<-c(colnames(select(cohort, starts_with("outcome"))))
table<-CreateTableOne(vars = vars, strata = c("group"), data = cohort, test = F)
table<-data.frame(print(table, smd = T))
write.xlsx(table, file = "study_cohort(outcome events in 1 year).xlsx",rowNames=T)

bmi_base<-merge(study_cohort[,c("reference.key","index.date")],BMI,by="reference.key",all = F)
bmi_base$bmi.dur<-as.numeric(as.numeric(bmi_base$index.date)-as.numeric(bmi_base$date.bmi))
bmi_base$bmi.dur.abs<-abs(as.numeric(bmi_base$index.date)-as.numeric(bmi_base$date.bmi))
colnames(bmi_base)
bmi_base<-arrange(bmi_base,reference.key,bmi.dur.abs)
bmi_base<-bmi_base[which(bmi_base$bmi.dur.abs<=365),]
bmi_base<-bmi_base[!duplicated(bmi_base$reference.key),]
study_cohort<-merge(study_cohort,bmi_base[,c(1,4)],by="reference.key",all.x = T)

aim.d<-merge(study_cohort,dxall,by="reference.key",all.x = T)
colnames(aim.d)
colnames(aim.d)[c((ncol(aim.d)-1):ncol(aim.d))]<-c("d.date","d.code")
aim.d<-aim.d[which(aim.d$d.date<aim.d$index.date),]
aim.d<-arrange(aim.d,reference.key,d.date)

colnames(study_cohort)
base<-data.frame(matrix(0,ncol = nrow(code.list),nrow = nrow(study_cohort)))
baseline_cohort<-cbind(study_cohort,base)
colnames(baseline_cohort)<-c(colnames(study_cohort),code.list$bs.dx.names)

for (i in c(1:8)) {
  baseline_cohort[grepl(code.list$diagnosis.codes[i],baseline_cohort$cancer.code),][,(ncol(study_cohort)+i)]<-1
}
for (i in c(9:(length(code.list$diagnosis.codes)))) {
  data<-aim.d[grepl(code.list$diagnosis.codes[i],aim.d$d.code),]
  rk<-unique(data$reference.key)
  baseline_cohort[which(baseline_cohort$reference.key %in% rk),][,(ncol(study_cohort)+i)]<-1
}

cancer_drug<-rxall[grepl("^8",rxall$drug.code),]
drug_cohort<-merge(study_cohort,cancer_drug,by="reference.key",all = F)
treatment1<-drug_cohort[which(drug_cohort$date.rx.start<=drug_cohort$index.date&drug_cohort$date.rx.start>=drug_cohort$index.date-90),]
treatment2<-drug_cohort[which(drug_cohort$date.rx.start<drug_cohort$index.date-90&drug_cohort$date.rx.end>=drug_cohort$index.date-90),]
rk_1<-unique(treatment1$reference.key)
rk_2<-unique(treatment2$reference.key)
dx_drug<-dxall[grepl("V58.1|V66.2|V67.2",dxall$dx.code),]
dx_drug_cohort<-merge(study_cohort,dx_drug,by="reference.key",all = F)
treatment3<-dx_drug_cohort[which(dx_drug_cohort$date.dx<=dx_drug_cohort$index.date&dx_drug_cohort$date.dx>dx_drug_cohort$index.date-90),]
rk_3<-unique(treatment3$reference.key)
px_drug_1<-pxall[grepl("^99.25|^99.28",pxall$px.code),]
px_drug_2<-pxall[grepl("chemotherapy|BRM",pxall$detail,ignore.case = T),]
px_drug<-rbind(px_drug_1,px_drug_2)
px_drug_cohort<-merge(study_cohort,px_drug,by="reference.key",all = F)
treatment4<-px_drug_cohort[which(px_drug_cohort$date.px<=px_drug_cohort$index.date&px_drug_cohort$date.px>px_drug_cohort$index.date-90),]
rk_4<-unique(treatment4$reference.key)
rk_drug<-union(union(rk_1,rk_2),union(rk_3,rk_4))
baseline_cohort$bs.tr.Drug.therapy<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_drug),]$bs.tr.Drug.therapy<-1

dx_radiotherapy<-dxall[grepl("V58.0|V66.1|V67.1",dxall$dx.code),]
dx_radiotherapy_cohort<-merge(study_cohort,dx_radiotherapy,by="reference.key",all = F)
treatment<-dx_radiotherapy_cohort[which(dx_radiotherapy_cohort$date.dx<=dx_radiotherapy_cohort$index.date&dx_radiotherapy_cohort$date.dx>dx_radiotherapy_cohort$index.date-90),]
rk_1<-unique(treatment$reference.key)
px_radiotherapy_1<-pxall[grepl("^92.[23]",pxall$px.code),]
px_radiotherapy_2<-pxall[grepl("radiotherapy",pxall$detail,ignore.case = T),]
px_radiotherapy<-rbind(px_radiotherapy_1,px_radiotherapy_2)
px_radiotherapy_cohort<-merge(study_cohort,px_radiotherapy,by="reference.key",all = F)
treatment<-px_radiotherapy_cohort[which(px_radiotherapy_cohort$date.px<=px_radiotherapy_cohort$index.date&px_radiotherapy_cohort$date.px>px_radiotherapy_cohort$index.date-90),]
rk_2<-unique(treatment$reference.key)
rk_radiotherapy<-union(rk_1,rk_2)
baseline_cohort$bs.tr.Radiotherapy<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_radiotherapy),]$bs.tr.Radiotherapy<-1

cancer_drug<-rxall[grepl("^8",rxall$drug.code),]
drug_cohort<-merge(study_cohort,cancer_drug,by="reference.key",all = F)
treatment1<-drug_cohort[which(drug_cohort$date.rx.start<=drug_cohort$index.date&drug_cohort$date.rx.start>=drug_cohort$index.date-90),]
treatment2<-drug_cohort[which(drug_cohort$date.rx.start<drug_cohort$index.date-90&drug_cohort$date.rx.end>=drug_cohort$index.date-90),]
rk_1<-unique(treatment1$reference.key)
rk_2<-unique(treatment2$reference.key)
dx_drug<-dxall[grepl("V58.1|V66.2|V67.2",dxall$dx.code),]
dx_drug_cohort<-merge(study_cohort,dx_drug,by="reference.key",all = F)
treatment3<-dx_drug_cohort[which(dx_drug_cohort$date.dx<=dx_drug_cohort$index.date&dx_drug_cohort$date.dx>dx_drug_cohort$index.date-90),]
rk_3<-unique(treatment3$reference.key)
px_drug_1<-pxall[grepl("^99.25|^99.28",pxall$px.code),]
px_drug_2<-pxall[grepl("chemotherapy|BRM",pxall$detail,ignore.case = T),]
px_drug<-rbind(px_drug_1,px_drug_2)
px_drug_cohort<-merge(study_cohort,px_drug,by="reference.key",all = F)
treatment4<-px_drug_cohort[which(px_drug_cohort$date.px<=px_drug_cohort$index.date&px_drug_cohort$date.px>px_drug_cohort$index.date-90),]
rk_4<-unique(treatment4$reference.key)
rk_drug<-union(union(rk_1,rk_2),union(rk_3,rk_4))
baseline_cohort$bs.tr.Drug.therapy<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_drug),]$bs.tr.Drug.therapy<-1

dx_radiotherapy<-dxall[grepl("V58.0|V66.1|V67.1",dxall$dx.code),]
dx_radiotherapy_cohort<-merge(study_cohort,dx_radiotherapy,by="reference.key",all = F)
treatment<-dx_radiotherapy_cohort[which(dx_radiotherapy_cohort$date.dx<=dx_radiotherapy_cohort$index.date&dx_radiotherapy_cohort$date.dx>dx_radiotherapy_cohort$index.date-90),]
rk_1<-unique(treatment$reference.key)
px_radiotherapy_1<-pxall[grepl("^92.[23]",pxall$px.code),]
px_radiotherapy_2<-pxall[grepl("radiotherapy",pxall$detail,ignore.case = T),]
px_radiotherapy<-rbind(px_radiotherapy_1,px_radiotherapy_2)
px_radiotherapy_cohort<-merge(study_cohort,px_radiotherapy,by="reference.key",all = F)
treatment<-px_radiotherapy_cohort[which(px_radiotherapy_cohort$date.px<=px_radiotherapy_cohort$index.date&px_radiotherapy_cohort$date.px>px_radiotherapy_cohort$index.date-90),]
rk_2<-unique(treatment$reference.key)
rk_radiotherapy<-union(rk_1,rk_2)
baseline_cohort$bs.tr.Radiotherapy<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_radiotherapy),]$bs.tr.Radiotherapy<-1

aspirin<-rxall[grepl("aspirin",rxall$drug.name,ignore.case=T),]
unique(aspirin$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],aspirin,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.aspirin<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.aspirin<-1}

Antiplatelet<-rxall[grepl("^2.9",rxall$drug.code),]
Antiplatelet<-Antiplatelet[!grepl("aspirin",Antiplatelet$drug.name,ignore.case=T),]
unique(Antiplatelet$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],Antiplatelet,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.Antiplatelet<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.Antiplatelet<-1}

NSAIDs<-rxall[grepl("^10.1.1",rxall$drug.code),]
NSAIDs<-NSAIDs[!grepl("aspirin",NSAIDs$drug.name,ignore.case=T),]
unique(NSAIDs$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],NSAIDs,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.NSAIDs<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.NSAIDs<-1}

Erythropoietin<-rxall[grepl("^9.1.3",rxall$drug.code),]
unique(Erythropoietin$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],Erythropoietin,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.Erythropoietin<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.Erythropoietin<-1}

EGFR<-rxall[grepl("gefitinib|erlotinib|afatinib|brigatinib|icotinib|cetuximab|Osimertinib",rxall$drug.name,ignore.case=T),]
unique(EGFR$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],EGFR,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.EGFR.inhibitor<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.EGFR.inhibitor<-1}

VEGF<-rxall[grepl("pazopanib|sunitinib|bevacizumab|sorafenib|regorafenib|cabozantinib|
                   lenvatinib|ponatinib|aflibercept|axitinib|tivozanib|ramucirumab|vandetanib",
                  rxall$drug.name,ignore.case=T),]
unique(VEGF$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],VEGF,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.VEGF.inhibitor<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.VEGF.inhibitor<-1}

CYP3A4<-rxall[grepl("Amiodarone|Dronederone|Quinidine|Verapamil|ritonavir|darunavir|fosamprenavir|indinavir|lopinavir|
                     nelfinavir|saquinavir|Itraconazole|Voriconazole|Posaconazole|Ketoconazole|Fluconazole|
                     Clarithromycin|Erythromycin|Rifampicin|Ciclosporin|Tacrolimus|Phenytoin|Carbamazepine|
                     Phenobarbitone",rxall$drug.name,ignore.case=T),]
unique(CYP3A4$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],CYP3A4,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.CYP3A4.Pgp<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.CYP3A4.Pgp<-1}

SSRIs<-rxall[grepl("BUPROPION|CITALOPRAM|DESVENLAFAXINE|DULOXETINE|ESCITALOPRAM|FLUOXETINE|FLUVOXAMINE|
                    MILNACIPRAN|PAROXETINE|SERTRALINE|VENLAFAXINE",rxall$drug.name,ignore.case=T),]
unique(SSRIs$drug.name)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],SSRIs,by="reference.key",all=F)
as_baseline1<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-90),]
as_baseline2<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-90&as_baseline$date.rx.end>=as_baseline$index.date-90),]
as_baseline<-rbind(as_baseline1,as_baseline2)
rk<-unique(as_baseline$reference.key)
baseline_cohort$bs.rx.SSRIs<-0
if(nrow(baseline_cohort[which(baseline_cohort$reference.key %in% rk),])>0){baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.rx.SSRIs<-1}

platelet<-laball[grepl("platelet",laball$lab.name,ignore.case=T),]
platelet_base<-merge(baseline_cohort[,c("reference.key","index.date")],platelet,by="reference.key",all = F)
length(unique(platelet_base$reference.key))
platelet_base$dur<-as.numeric(platelet_base$index.date-platelet_base$date.lab)
check1<-platelet_base[which(platelet_base$dur>=0&platelet_base$dur<=90),]
length(unique(check1$reference.key))
check2<-check1[which(check1$lab.result>=350),]
length(unique(check2$reference.key))
check3<-check2[duplicated(check2$reference.key),]
length(unique(check3$reference.key))
rk_platelet<-unique(check3$reference.key)
baseline_cohort$lab.platelet<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_platelet),]$lab.platelet<-1

haemoglobin<-laball[grepl("haemoglobin",laball$lab.name,ignore.case=T),]
haemoglobin<-haemoglobin[grepl("blood",haemoglobin$lab.name,ignore.case=T),]
haemoglobin_base<-merge(baseline_cohort[,c("reference.key","index.date")],haemoglobin,by="reference.key",all = F)
length(unique(haemoglobin_base$reference.key))
haemoglobin_base$dur<-as.numeric(haemoglobin_base$index.date-haemoglobin_base$date.lab)
check1<-haemoglobin_base[which(haemoglobin_base$dur>=0&haemoglobin_base$dur<=90),]
length(unique(check1$reference.key))
check2<-check1[which(check1$lab.result<10),]
length(unique(check2$reference.key))
check3<-check2[duplicated(check2$reference.key),]
length(unique(check3$reference.key))
rk_haemoglobin<-unique(check3$reference.key)
baseline_cohort$lab.haemoglobin<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_haemoglobin),]$lab.haemoglobin<-1

wbc<-laball[which(laball$lab.name=="WBC"),]
wbc_base<-merge(baseline_cohort[,c("reference.key","index.date")],wbc,by="reference.key",all = F)
length(unique(wbc_base$reference.key))
wbc_base$dur<-as.numeric(wbc_base$index.date-wbc_base$date.lab)
check1<-wbc_base[which(wbc_base$dur>=0&wbc_base$dur<=90),]
length(unique(check1$reference.key))
check2<-check1[which(check1$lab.result>11),]
length(unique(check2$reference.key))
check3<-check2[duplicated(check2$reference.key),]
length(unique(check3$reference.key))
rk_wbc<-unique(check3$reference.key)
baseline_cohort$lab.wbc<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_wbc),]$lab.wbc<-1

bmi_dx<-dxall[grepl("V85",dxall$dx.code),]
cancer2<-baseline_cohort[grepl("^151|^185|^19[12]",baseline_cohort$cancer.code),]
cancer1<-baseline_cohort[grepl("^162|^20[0-8]|^179|^18[01234689]|^172",baseline_cohort$cancer.code),]
baseline_cohort$bs.Khorana<-NA
for (i in c(1:nrow(baseline_cohort))) {
  rk<-baseline_cohort$reference.key[i]
  if(nrow(cancer2[which(cancer2$reference.key==rk),])>0){a<-2}else{a<-0}
  if(nrow(cancer1[which(cancer1$reference.key==rk),])>0){a<-1}
  if(baseline_cohort[i,"lab.platelet"]==1){b<-1}else{b<-0}
  if(baseline_cohort[i,"lab.haemoglobin"]==1){c<-1}else{c<-0}
  if(baseline_cohort[i,"bs.rx.Erythropoietin"]==1){c<-1}
  if(baseline_cohort[i,"lab.wbc"]==1){d<-1}else{d<-0}
  if(is.na(baseline_cohort[i,"bmi"])==F&baseline_cohort[i,"bmi"]>=35){e<-1}else{e<-0}
  score<-a+b+c+d+e
  if(score<2){baseline_cohort[i,"bs.Khorana"]<-1}
  if(score>=2){baseline_cohort[i,"bs.Khorana"]<-2}
}

aim.d<-merge(study_cohort,dxall,by="reference.key",all.x = T)
colnames(aim.d)
colnames(aim.d)[c((ncol(aim.d)-1):ncol(aim.d))]<-c("d.date","d.code")
aim.d<-aim.d[which(aim.d$d.date<=aim.d$index.date),]
aim.d<-arrange(aim.d,reference.key,d.date)
baseline_cohort$bs.CCI<-NA

for (z in c(1:nrow(baseline_cohort))) {
  rk<-baseline_cohort$reference.key[z]
  data<-baseline_cohort[which(baseline_cohort$reference.key==rk),]
  dx<-aim.d[which(aim.d$reference.key==rk),]
  if(data$age>=50&data$age<60){a<-1}else{a<-0}
  if(data$age>=60&data$age<70){a<-2}
  if(data$age>=70&data$age<80){a<-3}
  if(data$age>=80){a<-4}
  if(nrow(dx[grepl("^41[02]",dx$d.code),])>0){b<-1}else{b<-0}
  if(nrow(dx[grepl("^398.91|^402.[019]1|^404.[019][13]|^428",dx$d.code),])>0){c<-1}else{c<-0}
  if(nrow(dx[grepl("^441|^443.9|^785.4|^V43.4",dx$d.code),])>0){d<-1}else{d<-0}
  if(nrow(dx[grepl("^43[0-8]|^435",dx$d.code),])>0){e<-1}else{e<-0}
  if(nrow(dx[grepl("^290|^294.[0-2]|^331",dx$d.code),])>0){f<-1}else{f<-0}
  if(nrow(dx[grepl("^49[0-6]|^50[0-5]|^506.4",dx$d.code),])>0){g<-1}else{g<-0}
  if(nrow(dx[grepl("^71[04]|^725",dx$d.code),])>0){h<-1}else{h<-0}
  if(nrow(dx[grepl("^53[1-4]",dx$d.code),])>0){i<-1}else{i<-0}
  if(nrow(dx[grepl("^571.[2456]",dx$d.code),])>0){j<-1}else{j<-0}
  if(nrow(dx[grepl("^456.[0-2]|^572.[2-8]",dx$d.code),])>0){j<-3}
  if(nrow(dx[grepl("^250.[01237]",dx$d.code),])>0){k<-1}else{k<-0}
  if(nrow(dx[grepl("^250.[4-6]",dx$d.code),])>0){k<-2}
  if(nrow(dx[grepl("^342|^344",dx$d.code),])>0){l<-2}else{l<-0}
  if(nrow(dx[grepl("^582|^583.[0-7]|^58[568]",dx$d.code),])>0){m<-2}else{m<-0}
  if(nrow(dx[grepl("^1[4-8][0-9]|^19[0-4]",dx$d.code),])>0){n<-2}else{n<-0}
  if(nrow(dx[grepl("^19[6-9]",dx$d.code),])>0){n<-6}
  if(nrow(dx[grepl("^20[3-8]",dx$d.code),])>0){o<-2}else{o<-0}
  if(nrow(dx[grepl("^20[0-2]",dx$d.code),])>0){p<-2}else{p<-0}
  if(nrow(dx[grepl("^042|^079.53",dx$d.code),])>0){q<-6}else{q<-0}
  score<-a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q
  if(score<6){baseline_cohort[z,"bs.CCI"]<-1}
  if(score>=6){baseline_cohort[z,"bs.CCI"]<-2}
}

catheter<-pxall[grepl("central venous catheter",pxall$detail,ignore.case = T),]
catheter_base<-merge(baseline_cohort[,c("reference.key","index.date")],catheter,by="reference.key",all=F)
catheter_base<-catheter_base[which(catheter_base$date.px<=catheter_base$index.date),]
rk<-unique(catheter_base$reference.key)
baseline_cohort$bs.px.catheter<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.px.catheter<-1

surgery<-pxall[grepl("^0[1-9]|^[1-7][0-9]|^8[0-6]",pxall$px.code),]
px_surgery<-surgery[!grepl("scopy|care|biopsy|diagnostic",surgery$detail,ignore.case = T),]
px_surgery_cohort<-merge(study_cohort,px_surgery,by="reference.key",all = F)
treatment<-px_surgery_cohort[which(px_surgery_cohort$date.px<=px_surgery_cohort$index.date&px_surgery_cohort$date.px>px_surgery_cohort$index.date-90),]
rk<-unique(treatment$reference.key)
baseline_cohort$bs.px.Surgery<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.px.Surgery<-1

transfusion<-pxall[grepl("^99.0",pxall$px.code,ignore.case = T),]
transfusion_base<-merge(baseline_cohort[,c("reference.key","index.date")],transfusion,by="reference.key",all=F)
transfusion_base<-transfusion_base[which(transfusion_base$date.px<=transfusion_base$index.date&transfusion_base$date.px>transfusion_base$index.date-90),]
rk<-unique(transfusion_base$reference.key)
baseline_cohort$bs.px.transfusion<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$bs.px.transfusion<-1

ip_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],ipall,by="reference.key",all=F)
ip_baseline1<-ip_baseline[which(ip_baseline$date.ip.admission<=ip_baseline$index.date&ip_baseline$date.ip.admission>=ip_baseline$index.date-90),]
ip_baseline2<-ip_baseline[which(ip_baseline$date.ip.admission<ip_baseline$index.date-90&ip_baseline$date.ip.discharge>=ip_baseline$index.date-90),]
ip_baseline<-rbind(ip_baseline1,ip_baseline2)
ip_baseline<-distinct(ip_baseline,reference.key,date.ip.admission,date.ip.discharge)
baseline_cohort$bs.ip<-0
for (i in c(1:length(ip_baseline$reference.key))) {
  rk<-ip_baseline$reference.key[i]
  baseline_cohort[which(baseline_cohort$reference.key==rk),]$bs.ip<-nrow(ip_baseline[which(ip_baseline$reference.key==rk),])
}

ae_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],aeall,by="reference.key",all=F)
ae_baseline<-ae_baseline[which(ae_baseline$date.ae<=ae_baseline$index.date&ae_baseline$date.ae>=ae_baseline$index.date-90),]
ae_baseline<-distinct(ae_baseline,reference.key,date.ae)
baseline_cohort$bs.ae<-0
for (i in c(1:length(ae_baseline$reference.key))) {
  rk<-ae_baseline$reference.key[i]
  baseline_cohort[which(baseline_cohort$reference.key==rk),]$bs.ae<-nrow(ae_baseline[which(ae_baseline$reference.key==rk),])
}

renal<-aim.d[grepl("^40[34]|^58[0-6]|^590.0|^753.[1-3]|^V42.0|^45.1|^V56",aim.d$d.code),]
rk1<-unique(renal$reference.key)

liver<-aim.d[grepl("^571.[2456]|^456.[0-2]|^572.[2-8]",aim.d$d.code),]
rk2<-unique(liver$reference.key)

bleeding<-aim.d[grepl("^43[0-2]|^53[1-4].[0246]|^535.[0-9]1|^537.83|^562.0[23]|^562.1[23]|^568.81|^569.3|^569.85|
                        ^578|^423.0|^459.0|^593.81|^599.7|^719.1|^784.7|^784.8|^786.3",aim.d$d.code),]
b3<-bleeding[which(bleeding$d.date>bleeding$index.date-90),]
rk3<-unique(b3$reference.key)

surgery<-pxall[grepl("^[0-7][0-9]|^8[0-6]",pxall$px.code),]
surgery<-surgery[!grepl("scopy|care|biopsy|diagnostic",surgery$detail,ignore.case=T),]
check_surgery<-merge(baseline_cohort,surgery,by="reference.key",all.x = T)
b4<-check_surgery[which(check_surgery$date.px>check_surgery$index.date-14&check_surgery$date.px<check_surgery$index.date),]
rk4<-unique(b4$reference.key)

aspirin<-rxall[grepl("aspirin",rxall$drug.name,ignore.case=T),]
antiplatelet<-rxall[grepl("^2.9",rxall$drug.code),]
NSAIDs<-rxall[grepl("^10.1.1",rxall$drug.code),]
EGFR<-rxall[grepl("gefitinib|erlotinib|afatinib|brigatinib|icotinib|cetuximab|Osimertinib",rxall$drug.name,ignore.case=T),]
SSRIs<-rxall[grepl("BUPROPION|CITALOPRAM|DESVENLAFAXINE|DULOXETINE|ESCITALOPRAM|FLUOXETINE|FLUVOXAMINE|
                    MILNACIPRAN|PAROXETINE|SERTRALINE|VENLAFAXINE",rxall$drug.name,ignore.case=T),]
thrombolytics<-rxall[grepl("^2.10",rxall$drug.code),]
drugs<-rbind(aspirin,antiplatelet,NSAIDs,EGFR,SSRIs,thrombolytics)
as_baseline<-merge(baseline_cohort[,c("reference.key","index.date")],drugs,by="reference.key",all=F)
as_baseline3<-as_baseline[which(as_baseline$date.rx.start<=as_baseline$index.date&as_baseline$date.rx.start>=as_baseline$index.date-42),]
as_baseline4<-as_baseline[which(as_baseline$date.rx.start<as_baseline$index.date-42&as_baseline$date.rx.end>=as_baseline$index.date-42),]
as_baseline2<-rbind(as_baseline3,as_baseline4)
rk5<-unique(as_baseline2$reference.key)

metastatic<-aim.d[grepl("^19[6-9]",aim.d$d.code),]
rk6<-unique(metastatic$reference.key)

rk_bleed<-union(union(union(rk1,rk2),union(rk3,rk4)),union(rk5,rk6))
baseline_cohort$bleeding.risk<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk_bleed),]$bleeding.risk<-1
baseline_cohort$bleeding.risk<-as.factor(baseline_cohort$bleeding.risk)
nrow(baseline_cohort[which(baseline_cohort$bleeding.risk==1),])
nrow(baseline_cohort[which(baseline_cohort$bleeding.risk==0),])

colnames<-setdiff(colnames(select(baseline_cohort,starts_with("bs"))),c("bs.ip","bs.ae"))
for (i in colnames) {
  baseline_cohort[,i]<-as.factor(baseline_cohort[,i])
}

baseline_cohort[which(baseline_cohort$sex=="M"),]$sex<-0
baseline_cohort[which(baseline_cohort$sex=="F"),]$sex<-1
baseline_cohort$sex<-as.factor(baseline_cohort$sex)
str(baseline_cohort)

vars <- c("age", "sex",colnames(select(baseline_cohort, starts_with("bs"))))
table <- CreateTableOne(vars = vars, strata = c("group"), data = baseline_cohort, test = F)
table <- data.frame(print(table, smd = T))
setwd("C:/Users/LabPC12CSMPR/Desktop")
write.xlsx(table, file = "baseline.xlsx",rowNames=T)

save(baseline_cohort,file = "baseline.RData")

df_cox<-baseline_cohort

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1

df_cox<- df_cox %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.mb.date)%>%
  mutate(censor.date.2 = pmin(outcome.mb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1

df_cox<- df_cox %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.death.date)%>%
  mutate(censor.date.3 = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
df_cox[which(df_cox$outcome.date.3>df_cox$censor.date.3),]$outcome.date.3<-NA
df_cox[which(is.na(df_cox$outcome.date.3)==F),]$outcome.3<-1

formula<-str_c(c("group ~ age+sex",colnames(select(df_cox,starts_with("bs")))), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

vars <- c("age", "sex",colnames(select(baseline_cohort, starts_with("bs"))))
table <- svyCreateTableOne(vars=vars,strata = "group", data = design,test=F)
table <- data.frame(print(table, smd = T))
setwd("C:/Users/LabPC12CSMPR/Desktop")
write.xlsx(table, file = "baseline2.xlsx",rowNames=T)

fit <- survfit(Surv(follow.up.1,outcome.1) ~ group,data = df_cox,weights = W.glm$weights)
e<-ggsurvplot(fit, data = df_cox,conf.int = T,fun = "event",risk.table = T,xlim=c(0,180),ylim=c(0,0.2),break.x.by=30,pval=T,
              legend.labs=c("Non-switchers","Switchers"),xlab="Days",ylab="Cumulative incidence") 
f<-e+ggtitle("Hospitalization due to VTE")
jpeg("Hospitalization due to VTE.jpeg",width=4000,height = 2500,res = 300)
print(f)
dev.off()

fit <- survfit(Surv(follow.up.1,outcome.1) ~ group,data = df_cox,weights = W.glm$weights)
surv_pvalue(fit)$pval
pdf(file = "C:/Users/LabPC12CSMPR/Desktop/Figure 2A.pdf")
ggsurvplot(fit, data = df_cox,conf.int = T,fun = "event",risk.table = T,xlim=c(0,180),break.x.by=30,pval=paste("p","=","2.23e-6"),
           legend.labs=c("Non-switchers","Switchers"),xlab="Days",ylab="Cumulative incidence")+ggtitle("Hospitalization due to VTE")
dev.off()

fit <- survfit(Surv(follow.up.2,outcome.2) ~ group,data = df_cox,weights = W.glm$weights)
e<-ggsurvplot(fit, data = df_cox,conf.int = T,fun = "event",risk.table = T,xlim=c(0,180),ylim=c(0,0.2),break.x.by=30,pval=T,
              legend.labs=c("Non-switchers","Switchers"),xlab="Days",ylab="Cumulative incidence") 
f<-e+ggtitle("Major bleeding")
jpeg("Major Bleeding.jpeg",width=4000,height = 2500,res = 300)
print(f)
dev.off()

fit <- survfit(Surv(follow.up.2,outcome.2) ~ group,data = df_cox,weights = W.glm$weights)
surv_pvalue(fit)$pval
pdf(file = "C:/Users/LabPC12CSMPR/Desktop/Figure 2B.pdf") 
ggsurvplot(fit, data = df_cox,conf.int = T,fun = "event",risk.table = T,xlim=c(0,180),break.x.by=30,pval=paste("p","=","0.931"),
           legend.labs=c("Non-switchers","Switchers"),xlab="Days",ylab="Cumulative incidence")+ggtitle("Major bleeding")
dev.off()

fit <- survfit(Surv(follow.up.3,outcome.3) ~ group,data = df_cox,weights = W.glm$weights)
e<-ggsurvplot(fit, data = df_cox,conf.int = T,fun = "event",risk.table = T,xlim=c(0,180),break.x.by=30,pval=T,
              legend.labs=c("Non-switchers","Switchers"),xlab="Days",ylab="Cumulative incidence") 
f<-e+ggtitle("All-cause mortality")
jpeg("All-cause Mortality.jpeg",width=4000,height = 2500,res = 300)
print(f)
dev.off()
fit <- survfit(Surv(follow.up.3,outcome.3) ~ group,data = df_cox,weights = W.glm$weights)
surv_pvalue(fit)$pval
pdf(file = "C:/Users/LabPC12CSMPR/Desktop/Figure 2C.pdf") 
ggsurvplot(fit, data = df_cox,conf.int = T,fun = "event",risk.table = T,xlim=c(0,180),break.x.by=30,pval=paste("p","=","1.94e-28"),
           legend.labs=c("Non-switchers","Switchers"),xlab="Days",ylab="Cumulative incidence")+ggtitle("All-cause mortality")
dev.off()

df_cox<-baseline_cohort
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
df_cox$tr.end<-NA
df_cox$tr.end<-as.Date(df_cox$tr.end)
df_cox[which(df_cox$group==1),]$tr.end<-df_cox[which(df_cox$group==1),]$date.noac.e
df_cox[which(df_cox$group==0),]$tr.end<-df_cox[which(df_cox$group==0),]$date.lmwh.e

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1

df_cox<- df_cox %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.dvt.date)%>%
  mutate(censor.date.2 = pmin(outcome.dvt.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1

df_cox<- df_cox %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.pe.date)%>%
  mutate(censor.date.3 = pmin(outcome.pe.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
df_cox[which(df_cox$outcome.date.3>df_cox$censor.date.3),]$outcome.date.3<-NA
df_cox[which(is.na(df_cox$outcome.date.3)==F),]$outcome.3<-1

df_cox<- df_cox %>%
  mutate(outcome.4 = 0)%>%
  mutate(outcome.date.4 = outcome.mb.date)%>%
  mutate(censor.date.4 = pmin(outcome.mb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.4 = censor.date.4 - index.date)
df_cox[which(df_cox$outcome.date.4>df_cox$censor.date.4),]$outcome.date.4<-NA
df_cox[which(is.na(df_cox$outcome.date.4)==F),]$outcome.4<-1

df_cox<- df_cox %>%
  mutate(outcome.5 = 0)%>%
  mutate(outcome.date.5 = outcome.ich.date)%>%
  mutate(censor.date.5 = pmin(outcome.ich.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.5 = censor.date.5 - index.date)
df_cox[which(df_cox$outcome.date.5>df_cox$censor.date.5),]$outcome.date.5<-NA
df_cox[which(is.na(df_cox$outcome.date.5)==F),]$outcome.5<-1

df_cox<- df_cox %>%
  mutate(outcome.6 = 0)%>%
  mutate(outcome.date.6 = outcome.gib.date)%>%
  mutate(censor.date.6 = pmin(outcome.gib.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.6 = censor.date.6 - index.date)
df_cox[which(df_cox$outcome.date.6>df_cox$censor.date.6),]$outcome.date.6<-NA
df_cox[which(is.na(df_cox$outcome.date.6)==F),]$outcome.6<-1

df_cox<- df_cox %>%
  mutate(outcome.7 = 0)%>%
  mutate(outcome.date.7 = outcome.omb.date)%>%
  mutate(censor.date.7 = pmin(outcome.omb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.7 = censor.date.7 - index.date)
df_cox[which(df_cox$outcome.date.7>df_cox$censor.date.7),]$outcome.date.7<-NA
df_cox[which(is.na(df_cox$outcome.date.7)==F),]$outcome.7<-1

df_cox<- df_cox %>%
  mutate(outcome.8 = 0)%>%
  mutate(outcome.date.8 = outcome.death.date)%>%
  mutate(censor.date.8 = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.8 = censor.date.8 - index.date)
df_cox[which(df_cox$outcome.date.8>df_cox$censor.date.8),]$outcome.date.8<-NA
df_cox[which(is.na(df_cox$outcome.date.8)==F),]$outcome.8<-1

df_cox<- df_cox %>%
  mutate(outcome.9 = 0)%>%
  mutate(outcome.date.9 = outcome.vte.date)%>%
  mutate(censor.date.9 = pmin(outcome.vte.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.9 = censor.date.9 - index.date)
df_cox[which(df_cox$outcome.date.9>df_cox$censor.date.9),]$outcome.date.9<-NA
df_cox[which(is.na(df_cox$outcome.date.9)==F),]$outcome.9<-1

df_cox<- df_cox %>%
  mutate(outcome.10 = 0)%>%
  mutate(outcome.date.10 = outcome.dvt.date)%>%
  mutate(censor.date.10 = pmin(outcome.dvt.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.10 = censor.date.10 - index.date)
df_cox[which(df_cox$outcome.date.10>df_cox$censor.date.10),]$outcome.date.10<-NA
df_cox[which(is.na(df_cox$outcome.date.10)==F),]$outcome.10<-1

df_cox<- df_cox %>%
  mutate(outcome.11 = 0)%>%
  mutate(outcome.date.11 = outcome.pe.date)%>%
  mutate(censor.date.11 = pmin(outcome.pe.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.11 = censor.date.11 - index.date)
df_cox[which(df_cox$outcome.date.11>df_cox$censor.date.11),]$outcome.date.11<-NA
df_cox[which(is.na(df_cox$outcome.date.11)==F),]$outcome.11<-1

df_cox<- df_cox %>%
  mutate(outcome.12 = 0)%>%
  mutate(outcome.date.12 = outcome.mb.date)%>%
  mutate(censor.date.12 = pmin(outcome.mb.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.12 = censor.date.12 - index.date)
df_cox[which(df_cox$outcome.date.12>df_cox$censor.date.12),]$outcome.date.12<-NA
df_cox[which(is.na(df_cox$outcome.date.12)==F),]$outcome.12<-1

df_cox<- df_cox %>%
  mutate(outcome.13 = 0)%>%
  mutate(outcome.date.13 = outcome.ich.date)%>%
  mutate(censor.date.13 = pmin(outcome.ich.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.13 = censor.date.13 - index.date)
df_cox[which(df_cox$outcome.date.13>df_cox$censor.date.13),]$outcome.date.13<-NA
df_cox[which(is.na(df_cox$outcome.date.13)==F),]$outcome.13<-1

df_cox<- df_cox %>%
  mutate(outcome.14 = 0)%>%
  mutate(outcome.date.14 = outcome.gib.date)%>%
  mutate(censor.date.14 = pmin(outcome.gib.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.14 = censor.date.14 - index.date)
df_cox[which(df_cox$outcome.date.14>df_cox$censor.date.14),]$outcome.date.14<-NA
df_cox[which(is.na(df_cox$outcome.date.14)==F),]$outcome.14<-1

df_cox<- df_cox %>%
  mutate(outcome.15 = 0)%>%
  mutate(outcome.date.15 = outcome.omb.date)%>%
  mutate(censor.date.15 = pmin(outcome.omb.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.15 = censor.date.15 - index.date)
df_cox[which(df_cox$outcome.date.15>df_cox$censor.date.15),]$outcome.date.15<-NA
df_cox[which(is.na(df_cox$outcome.date.15)==F),]$outcome.15<-1

df_cox<- df_cox %>%
  mutate(outcome.16 = 0)%>%
  mutate(outcome.date.16 = outcome.death.date)%>%
  mutate(censor.date.16 = pmin(outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.16 = censor.date.16 - index.date)
df_cox[which(df_cox$outcome.date.16>df_cox$censor.date.16),]$outcome.date.16<-NA
df_cox[which(is.na(df_cox$outcome.date.16)==F),]$outcome.16<-1

result<-data.frame(matrix(NA,nrow = 16,ncol = 8))
rownames(result)<-c(paste0(colnames(df_cox)[22:29]," 6 months"),paste0(colnames(df_cox)[22:29]," 1 year"))

for (k in c(1:16)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,72+2+4*(k-1)])==F),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),72+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,72+2+4*(k-1)])==F),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),72+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,72+4+4*(k-1)], df_cox[,72+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.9, outcome.9) ~ group, design=design)
x<-summary(cox_result)
result[9,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[9,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.10, outcome.10) ~ group, design=design)
x<-summary(cox_result)
result[10,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[10,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.11, outcome.11) ~ group, design=design)
x<-summary(cox_result)
result[11,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[11,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.12, outcome.12) ~ group, design=design)
x<-summary(cox_result)
result[12,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[12,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.13, outcome.13) ~ group, design=design)
x<-summary(cox_result)
result[13,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[13,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.14, outcome.14) ~ group, design=design)
x<-summary(cox_result)
result[14,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[14,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.15, outcome.15) ~ group, design=design)
x<-summary(cox_result)
result[15,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[15,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.16, outcome.16) ~ group, design=design)
x<-summary(cox_result)
result[16,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[16,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

df_cox$pstat.1<-df_cox$outcome.1
df_cox$etime.1<-df_cox$follow.up.1
df_cox$event.1<-0
df_cox[which(df_cox$outcome.1==1),]$event.1<-1
df_cox[which(df_cox$outcome.1==0&df_cox$outcome.8==1),]$event.1<-2

df_cox$pstat.2<-df_cox$outcome.2
df_cox$etime.2<-df_cox$follow.up.2
df_cox$event.2<-0
df_cox[which(df_cox$outcome.2==1),]$event.2<-1
df_cox[which(df_cox$outcome.2==0&df_cox$outcome.8==1),]$event.2<-2

df_cox$pstat.3<-df_cox$outcome.3
df_cox$etime.3<-df_cox$follow.up.3
df_cox$event.3<-0
df_cox[which(df_cox$outcome.3==1),]$event.3<-1
df_cox[which(df_cox$outcome.3==0&df_cox$outcome.8==1),]$event.3<-2

df_cox$pstat.4<-df_cox$outcome.4
df_cox$etime.4<-df_cox$follow.up.4
df_cox$event.4<-0
df_cox[which(df_cox$outcome.4==1),]$event.4<-1
df_cox[which(df_cox$outcome.4==0&df_cox$outcome.8==1),]$event.4<-2

df_cox$pstat.5<-df_cox$outcome.5
df_cox$etime.5<-df_cox$follow.up.5
df_cox$event.5<-0
df_cox[which(df_cox$outcome.5==1),]$event.5<-1
df_cox[which(df_cox$outcome.5==0&df_cox$outcome.8==1),]$event.5<-2

df_cox$pstat.6<-df_cox$outcome.6
df_cox$etime.6<-df_cox$follow.up.6
df_cox$event.6<-0
df_cox[which(df_cox$outcome.6==1),]$event.6<-1
df_cox[which(df_cox$outcome.6==0&df_cox$outcome.8==1),]$event.6<-2

df_cox$pstat.7<-df_cox$outcome.7
df_cox$etime.7<-df_cox$follow.up.7
df_cox$event.7<-0
df_cox[which(df_cox$outcome.7==1),]$event.7<-1
df_cox[which(df_cox$outcome.7==0&df_cox$outcome.8==1),]$event.7<-2

df_cox$pstat.9<-df_cox$outcome.9
df_cox$etime.9<-df_cox$follow.up.9
df_cox$event.9<-0
df_cox[which(df_cox$outcome.9==1),]$event.9<-1
df_cox[which(df_cox$outcome.9==0&df_cox$outcome.8==1),]$event.9<-2

df_cox$pstat.10<-df_cox$outcome.10
df_cox$etime.10<-df_cox$follow.up.10
df_cox$event.10<-0
df_cox[which(df_cox$outcome.10==1),]$event.10<-1
df_cox[which(df_cox$outcome.10==0&df_cox$outcome.8==1),]$event.10<-2

df_cox$pstat.11<-df_cox$outcome.11
df_cox$etime.11<-df_cox$follow.up.11
df_cox$event.11<-0
df_cox[which(df_cox$outcome.11==1),]$event.11<-1
df_cox[which(df_cox$outcome.11==0&df_cox$outcome.8==1),]$event.11<-2

df_cox$pstat.12<-df_cox$outcome.12
df_cox$etime.12<-df_cox$follow.up.12
df_cox$event.12<-0
df_cox[which(df_cox$outcome.12==1),]$event.12<-1
df_cox[which(df_cox$outcome.12==0&df_cox$outcome.8==1),]$event.12<-2

df_cox$pstat.13<-df_cox$outcome.13
df_cox$etime.13<-df_cox$follow.up.13
df_cox$event.13<-0
df_cox[which(df_cox$outcome.13==1),]$event.13<-1
df_cox[which(df_cox$outcome.13==0&df_cox$outcome.8==1),]$event.13<-2

df_cox$pstat.14<-df_cox$outcome.14
df_cox$etime.14<-df_cox$follow.up.14
df_cox$event.14<-0
df_cox[which(df_cox$outcome.14==1),]$event.14<-1
df_cox[which(df_cox$outcome.14==0&df_cox$outcome.8==1),]$event.14<-2

df_cox$pstat.15<-df_cox$outcome.15
df_cox$etime.15<-df_cox$follow.up.15
df_cox$event.15<-0
df_cox[which(df_cox$outcome.15==1),]$event.15<-1
df_cox[which(df_cox$outcome.15==0&df_cox$outcome.8==1),]$event.15<-2

k<-df_cox[,colnames(select(df_cox,starts_with("bs"),"sex","age","group"))]

for (i in c(1:7)) {
  fit <- crr((df_cox[,136+2+3*(i-1)]), df_cox[,136+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
for (i in c(9:15)) {
  fit <- crr((df_cox[,136+2+3*(i-2)]), df_cox[,136+3+3*(i-2)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
result[c(8,16),c(7,8)]<-"-"
setwd("C:/Users/LabPC12CSMPR/Desktop")
write.xlsx(result,file="HR result-without treatment end nor index year.xlsx")

platin<-rxall[grepl("platin",rxall$drug.name,ignore.case=T),]
check<-merge(baseline_cohort,platin,by="reference.key",all = F)
platin_check1<-check[which(check$date.rx.start<=check$index.date&check$date.rx.start>=check$index.date-90),]
platin_check2<-check[which(check$date.rx.start<check$index.date-90&check$date.rx.end>=check$index.date-90),]
rk1<-unique(platin_check1$reference.key)
rk2<-unique(platin_check2$reference.key)
rk<-union(rk1,rk2)
baseline_cohort$platinum<-0
baseline_cohort[which(baseline_cohort$reference.key %in% rk),]$platinum<-1

df_cox<-baseline_cohort
df_cox$tr.end<-NA
df_cox$tr.end<-as.Date(df_cox$tr.end)
df_cox[which(df_cox$group==1),]$tr.end<-df_cox[which(df_cox$group==1),]$date.noac.e
df_cox[which(df_cox$group==0),]$tr.end<-df_cox[which(df_cox$group==0),]$date.lmwh.e

df_cox$group.switch<-0
df_cox[which(df_cox$group==1&df_cox$switch.days>7),]$group.switch<-2
df_cox[which(df_cox$group==1&df_cox$switch.days<=7),]$group.switch<-1

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1

df_cox<- df_cox %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.mb.date)%>%
  mutate(censor.date.2 = pmin(outcome.mb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1

df_cox<- df_cox %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.death.date)%>%
  mutate(censor.date.3 = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
df_cox[which(df_cox$outcome.date.3>df_cox$censor.date.3),]$outcome.date.3<-NA
df_cox[which(is.na(df_cox$outcome.date.3)==F),]$outcome.3<-1

sg1<-df_cox[which(df_cox$sex==0),]
sg2<-df_cox[which(df_cox$sex==1),]
sg3<-df_cox[which(df_cox$age<=65),]
sg4<-df_cox[which(df_cox$age>65),]
sg5<-df_cox[which(df_cox$bs.dx.Metastasis==0),]
sg6<-df_cox[which(df_cox$bs.dx.Metastasis==1),]
sg7<-df_cox[which(df_cox$bs.dx.Digestive.organs==0),]
sg8<-df_cox[which(df_cox$bs.dx.Digestive.organs==1),]
sg9<-df_cox[which(df_cox$platinum==0),]
sg10<-df_cox[which(df_cox$platinum==1),]
sg11<-df_cox[which(df_cox$bs.CCI==1),]
sg12<-df_cox[which(df_cox$bs.CCI==2),]
sg13<-df_cox[which(df_cox$bs.Khorana==1),]
sg14<-df_cox[which(df_cox$bs.Khorana==2),]
sg15<-df_cox[which(df_cox$bleeding.risk==0),]
sg16<-df_cox[which(df_cox$bleeding.risk==1),]
sg17<-df_cox[which(df_cox$group.switch!=2),]
sg18<-df_cox[which(df_cox$group.switch!=1),]

result<-data.frame(matrix(NA,nrow = 18,ncol = 4))
kk<-sg1

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[1,1]<-nrow(kk[which(kk$group==0),])
result[1,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[1,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg2

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[2,1]<-nrow(kk[which(kk$group==0),])
result[2,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[2,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[2,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg3

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[3,1]<-nrow(kk[which(kk$group==0),])
result[3,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[3,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg4

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[4,1]<-nrow(kk[which(kk$group==0),])
result[4,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[4,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[4,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg5

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[5,1]<-nrow(kk[which(kk$group==0),])
result[5,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[5,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[5,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg6
k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[6,1]<-nrow(kk[which(kk$group==0),])
result[6,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[6,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[6,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg7

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[7,1]<-nrow(kk[which(kk$group==0),])
result[7,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[7,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[7,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg8

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[8,1]<-nrow(kk[which(kk$group==0),])
result[8,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[8,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[8,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg9

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[9,1]<-nrow(kk[which(kk$group==0),])
result[9,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[9,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[9,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[9,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[9,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[9,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[9,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg10

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[10,1]<-nrow(kk[which(kk$group==0),])
result[10,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[10,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[10,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[10,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[10,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[10,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[10,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg11

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[11,1]<-nrow(kk[which(kk$group==0),])
result[11,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[11,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[11,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[11,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[11,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[11,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[11,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg12

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[12,1]<-nrow(kk[which(kk$group==0),])
result[12,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[12,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[12,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[12,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[12,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[12,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[12,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg13

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[13,1]<-nrow(kk[which(kk$group==0),])
result[13,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[13,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[13,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[13,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[13,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[13,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[13,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg14

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[14,1]<-nrow(kk[which(kk$group==0),])
result[14,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[14,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[14,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[14,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[14,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[14,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[14,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg15

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[15,1]<-nrow(kk[which(kk$group==0),])
result[15,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[15,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[15,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[15,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[15,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[15,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[15,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg16

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[16,1]<-nrow(kk[which(kk$group==0),])
result[16,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[16,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[16,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[16,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[16,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[16,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[16,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg17

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[17,1]<-nrow(kk[which(kk$group==0),])
result[17,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[17,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[17,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[17,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[17,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[17,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[17,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

kk<-sg18

k<-colnames(select(kk,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(kk[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(kk,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = kk, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = kk, weights = ~ W.glm$weights)
result[18,1]<-nrow(kk[which(kk$group==0),])
result[18,2]<-nrow(kk[which(kk$group==1),])
cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[18,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[18,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[18,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[18,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[18,7]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[18,8]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

setwd("C:/Users/LabPC12CSMPR/Desktop")
write.xlsx(result,file="HR result-subgroup-without treatment end nor index year.xlsx")

data_source<-baseline_cohort

data_source<- data_source %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
data_source[which(data_source$outcome.date.1>data_source$censor.date.1),]$outcome.date.1<-NA
data_source[which(is.na(data_source$outcome.date.1)==F),]$outcome.1<-1

data_source<- data_source %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.dvt.date)%>%
  mutate(censor.date.2 = pmin(outcome.dvt.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
data_source[which(data_source$outcome.date.2>data_source$censor.date.2),]$outcome.date.2<-NA
data_source[which(is.na(data_source$outcome.date.2)==F),]$outcome.2<-1

data_source<- data_source %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.pe.date)%>%
  mutate(censor.date.3 = pmin(outcome.pe.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
data_source[which(data_source$outcome.date.3>data_source$censor.date.3),]$outcome.date.3<-NA
data_source[which(is.na(data_source$outcome.date.3)==F),]$outcome.3<-1

data_source<- data_source %>%
  mutate(outcome.4 = 0)%>%
  mutate(outcome.date.4 = outcome.mb.date)%>%
  mutate(censor.date.4 = pmin(outcome.mb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.4 = censor.date.4 - index.date)
data_source[which(data_source$outcome.date.4>data_source$censor.date.4),]$outcome.date.4<-NA
data_source[which(is.na(data_source$outcome.date.4)==F),]$outcome.4<-1

data_source<- data_source %>%
  mutate(outcome.5 = 0)%>%
  mutate(outcome.date.5 = outcome.ich.date)%>%
  mutate(censor.date.5 = pmin(outcome.ich.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.5 = censor.date.5 - index.date)
data_source[which(data_source$outcome.date.5>data_source$censor.date.5),]$outcome.date.5<-NA
data_source[which(is.na(data_source$outcome.date.5)==F),]$outcome.5<-1

data_source<- data_source %>%
  mutate(outcome.6 = 0)%>%
  mutate(outcome.date.6 = outcome.gib.date)%>%
  mutate(censor.date.6 = pmin(outcome.gib.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.6 = censor.date.6 - index.date)
data_source[which(data_source$outcome.date.6>data_source$censor.date.6),]$outcome.date.6<-NA
data_source[which(is.na(data_source$outcome.date.6)==F),]$outcome.6<-1

data_source<- data_source %>%
  mutate(outcome.7 = 0)%>%
  mutate(outcome.date.7 = outcome.omb.date)%>%
  mutate(censor.date.7 = pmin(outcome.omb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.7 = censor.date.7 - index.date)
data_source[which(data_source$outcome.date.7>data_source$censor.date.7),]$outcome.date.7<-NA
data_source[which(is.na(data_source$outcome.date.7)==F),]$outcome.7<-1

data_source<- data_source %>%
  mutate(outcome.8 = 0)%>%
  mutate(outcome.date.8 = outcome.death.date)%>%
  mutate(censor.date.8 = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.8 = censor.date.8 - index.date)
data_source[which(data_source$outcome.date.8>data_source$censor.date.8),]$outcome.date.8<-NA
data_source[which(is.na(data_source$outcome.date.8)==F),]$outcome.8<-1

data_source$pstat.1<-data_source$outcome.1
data_source$etime.1<-data_source$follow.up.1
data_source$event.1<-0
data_source[which(data_source$outcome.1==1),]$event.1<-1
data_source[which(data_source$outcome.1==0&data_source$outcome.8==1),]$event.1<-2

data_source$pstat.2<-data_source$outcome.2
data_source$etime.2<-data_source$follow.up.2
data_source$event.2<-0
data_source[which(data_source$outcome.2==1),]$event.2<-1
data_source[which(data_source$outcome.2==0&data_source$outcome.8==1),]$event.2<-2

data_source$pstat.3<-data_source$outcome.3
data_source$etime.3<-data_source$follow.up.3
data_source$event.3<-0
data_source[which(data_source$outcome.3==1),]$event.3<-1
data_source[which(data_source$outcome.3==0&data_source$outcome.8==1),]$event.3<-2

data_source$pstat.4<-data_source$outcome.4
data_source$etime.4<-data_source$follow.up.4
data_source$event.4<-0
data_source[which(data_source$outcome.4==1),]$event.4<-1
data_source[which(data_source$outcome.4==0&data_source$outcome.8==1),]$event.4<-2

data_source$pstat.5<-data_source$outcome.5
data_source$etime.5<-data_source$follow.up.5
data_source$event.5<-0
data_source[which(data_source$outcome.5==1),]$event.5<-1
data_source[which(data_source$outcome.5==0&data_source$outcome.8==1),]$event.5<-2

data_source$pstat.6<-data_source$outcome.6
data_source$etime.6<-data_source$follow.up.6
data_source$event.6<-0
data_source[which(data_source$outcome.6==1),]$event.6<-1
data_source[which(data_source$outcome.6==0&data_source$outcome.8==1),]$event.6<-2

data_source$pstat.7<-data_source$outcome.7
data_source$etime.7<-data_source$follow.up.7
data_source$event.7<-0
data_source[which(data_source$outcome.7==1),]$event.7<-1
data_source[which(data_source$outcome.7==0&data_source$outcome.8==1),]$event.7<-2

s4<-dxall[grepl("^173|^20[0-8]",dxall$dx.code),]
s4<-merge(data_source,s4,by="reference.key",all.x = T)
s4<-s4[which(s4$date.dx<=s4$index.date+365),]
rk_s4<-unique(s4$reference.key)
data_source$s4<-0
data_source[which(data_source$reference.key %in% rk_s4),]$s4<-1
nrow(data_source[which(data_source$s4==1),])

df_cox<-data_source[which(data_source$s4==0),]
nrow(df_cox)
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

result<-data.frame(matrix(NA,nrow = 8,ncol = 8))
rownames(result)<-paste0(colnames(df_cox)[22:29]," 6 months")

for (k in c(1:8)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),71+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),71+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,71+4+4*(k-1)], df_cox[,71+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

k<-df_cox[,setdiff(colnames(select(data_source,starts_with("bs"),"sex","age","group")),rm)]

for (i in c(1:7)) {
  fit <- crr((df_cox[,103+2+3*(i-1)]), df_cox[,103+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
result[8,c(7,8)]<-"-"

setwd("C:/Users/LabPC12CSMPR/Desktop")
results4<-result
write.xlsx(results4,file="HR result-s4.xlsx")

data_source$s5<-0
data_source[which(as.integer(as.Date("2022-12-31")-data_source$index.date)>=180),]$s5<-1
nrow(data_source[which(data_source$s5==1),])

df_cox<-data_source[which(data_source$s5==1),]
nrow(df_cox)
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

result<-data.frame(matrix(NA,nrow = 8,ncol = 6))
rownames(result)<-paste0(colnames(df_cox)[22:29]," 6 months")

for (k in c(1:8)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),71+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),71+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,71+4+4*(k-1)], df_cox[,71+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

k<-df_cox[,setdiff(colnames(select(data_source,starts_with("bs"),"sex","age","group")),rm)]

for (i in c(1:7)) {
  fit <- crr((df_cox[,103+2+3*(i-1)]), df_cox[,103+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
result[8,c(7,8)]<-"-"

setwd("C:/Users/LabPC12CSMPR/Desktop")
results5<-result
write.xlsx(results5,file="HR result-s5.xlsx")

data_source$s6<-0
data_source[which(data_source$date.vte>=as.Date("2020-01-01")),]$s6<-1
nrow(data_source[which(data_source$s6==1),])

df_cox<-data_source[which(data_source$s6==0),]
nrow(df_cox)
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

result<-data.frame(matrix(NA,nrow = 8,ncol = 6))
rownames(result)<-paste0(colnames(df_cox)[22:29]," 6 months")

for (k in c(1:8)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),71+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),71+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,71+4+4*(k-1)], df_cox[,71+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

k<-df_cox[,setdiff(colnames(select(data_source,starts_with("bs"),"sex","age","group")),rm)]

for (i in c(1:7)) {
  fit <- crr((df_cox[,103+2+3*(i-1)]), df_cox[,103+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
result[8,c(7,8)]<-"-"

setwd("C:/Users/LabPC12CSMPR/Desktop")
results6<-result
write.xlsx(results6,file="HR result-s6.xlsx")

df_cox<-data_source
nrow(df_cox)
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}
result<-data.frame(matrix(NA,nrow = 8,ncol = 6))
rownames(result)<-paste0(colnames(df_cox)[22:29]," 6 months")

for (k in c(1:8)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),71+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,71+2+4*(k-1)])==F&df_cox[,71+2+4*(k-1)]<=df_cox[,71+3+4*(k-1)]),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),71+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,71+4+4*(k-1)], df_cox[,71+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex+index.year",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

k<-df_cox[,setdiff(colnames(select(data_source,starts_with("bs"),"sex","age","index.year","group")),rm)]

for (i in c(1:7)) {
  fit <- crr((df_cox[,103+2+3*(i-1)]), df_cox[,103+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
result[8,c(7,8)]<-"-"

setwd("C:/Users/LabPC12CSMPR/Desktop")
results7<-result
write.xlsx(results7,file="HR result-s7.xlsx")

df_cox<-baseline_cohort
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
df_cox$tr.end<-NA
df_cox$tr.end<-as.Date(df_cox$tr.end)
df_cox[which(df_cox$group==1),]$tr.end<-df_cox[which(df_cox$group==1),]$date.noac.e
df_cox[which(df_cox$group==0),]$tr.end<-df_cox[which(df_cox$group==0),]$date.lmwh.e

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1

df_cox<- df_cox %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.dvt.date)%>%
  mutate(censor.date.2 = pmin(outcome.dvt.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1

df_cox<- df_cox %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.pe.date)%>%
  mutate(censor.date.3 = pmin(outcome.pe.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
df_cox[which(df_cox$outcome.date.3>df_cox$censor.date.3),]$outcome.date.3<-NA
df_cox[which(is.na(df_cox$outcome.date.3)==F),]$outcome.3<-1

df_cox<- df_cox %>%
  mutate(outcome.4 = 0)%>%
  mutate(outcome.date.4 = outcome.mb.date)%>%
  mutate(censor.date.4 = pmin(outcome.mb.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.4 = censor.date.4 - index.date)
df_cox[which(df_cox$outcome.date.4>df_cox$censor.date.4),]$outcome.date.4<-NA
df_cox[which(is.na(df_cox$outcome.date.4)==F),]$outcome.4<-1

df_cox<- df_cox %>%
  mutate(outcome.5 = 0)%>%
  mutate(outcome.date.5 = outcome.ich.date)%>%
  mutate(censor.date.5 = pmin(outcome.ich.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.5 = censor.date.5 - index.date)
df_cox[which(df_cox$outcome.date.5>df_cox$censor.date.5),]$outcome.date.5<-NA
df_cox[which(is.na(df_cox$outcome.date.5)==F),]$outcome.5<-1

df_cox<- df_cox %>%
  mutate(outcome.6 = 0)%>%
  mutate(outcome.date.6 = outcome.gib.date)%>%
  mutate(censor.date.6 = pmin(outcome.gib.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.6 = censor.date.6 - index.date)
df_cox[which(df_cox$outcome.date.6>df_cox$censor.date.6),]$outcome.date.6<-NA
df_cox[which(is.na(df_cox$outcome.date.6)==F),]$outcome.6<-1

df_cox<- df_cox %>%
  mutate(outcome.7 = 0)%>%
  mutate(outcome.date.7 = outcome.omb.date)%>%
  mutate(censor.date.7 = pmin(outcome.omb.date, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.7 = censor.date.7 - index.date)
df_cox[which(df_cox$outcome.date.7>df_cox$censor.date.7),]$outcome.date.7<-NA
df_cox[which(is.na(df_cox$outcome.date.7)==F),]$outcome.7<-1

df_cox<- df_cox %>%
  mutate(outcome.8 = 0)%>%
  mutate(outcome.date.8 = outcome.death.date)%>%
  mutate(censor.date.8 = pmin(outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.8 = censor.date.8 - index.date)
df_cox[which(df_cox$outcome.date.8>df_cox$censor.date.8),]$outcome.date.8<-NA
df_cox[which(is.na(df_cox$outcome.date.8)==F),]$outcome.8<-1

result<-data.frame(matrix(NA,nrow = 8,ncol = 8))
rownames(result)<-c(paste0(colnames(df_cox)[22:29]," 6 months"),paste0(colnames(df_cox)[22:29]," 1 year"))

for (k in c(1:8)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,72+2+4*(k-1)])==F),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),72+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,72+2+4*(k-1)])==F),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),72+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,72+4+4*(k-1)], df_cox[,72+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

df_cox$pstat.1<-df_cox$outcome.1
df_cox$etime.1<-df_cox$follow.up.1
df_cox$event.1<-0
df_cox[which(df_cox$outcome.1==1),]$event.1<-1
df_cox[which(df_cox$outcome.1==0&df_cox$outcome.8==1),]$event.1<-2

df_cox$pstat.2<-df_cox$outcome.2
df_cox$etime.2<-df_cox$follow.up.2
df_cox$event.2<-0
df_cox[which(df_cox$outcome.2==1),]$event.2<-1
df_cox[which(df_cox$outcome.2==0&df_cox$outcome.8==1),]$event.2<-2

df_cox$pstat.3<-df_cox$outcome.3
df_cox$etime.3<-df_cox$follow.up.3
df_cox$event.3<-0
df_cox[which(df_cox$outcome.3==1),]$event.3<-1
df_cox[which(df_cox$outcome.3==0&df_cox$outcome.8==1),]$event.3<-2

df_cox$pstat.4<-df_cox$outcome.4
df_cox$etime.4<-df_cox$follow.up.4
df_cox$event.4<-0
df_cox[which(df_cox$outcome.4==1),]$event.4<-1
df_cox[which(df_cox$outcome.4==0&df_cox$outcome.8==1),]$event.4<-2

df_cox$pstat.5<-df_cox$outcome.5
df_cox$etime.5<-df_cox$follow.up.5
df_cox$event.5<-0
df_cox[which(df_cox$outcome.5==1),]$event.5<-1
df_cox[which(df_cox$outcome.5==0&df_cox$outcome.8==1),]$event.5<-2

df_cox$pstat.6<-df_cox$outcome.6
df_cox$etime.6<-df_cox$follow.up.6
df_cox$event.6<-0
df_cox[which(df_cox$outcome.6==1),]$event.6<-1
df_cox[which(df_cox$outcome.6==0&df_cox$outcome.8==1),]$event.6<-2

df_cox$pstat.7<-df_cox$outcome.7
df_cox$etime.7<-df_cox$follow.up.7
df_cox$event.7<-0
df_cox[which(df_cox$outcome.7==1),]$event.7<-1
df_cox[which(df_cox$outcome.7==0&df_cox$outcome.8==1),]$event.7<-2

k<-df_cox[,colnames(select(df_cox,starts_with("bs"),"sex","age","group"))]

for (i in c(1:7)) {
  fit <- crr((df_cox[,104+2+3*(i-1)]), df_cox[,104+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],e))
  print(x$coef[nrow(x$coef),5])
}
result[8,c(7,8)]<-"-"

setwd("C:/Users/LabPC12CSMPR/Desktop")
results8<-result
write.xlsx(results8,file="HR result-s8.xlsx")

df_cox<-baseline_cohort
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])

lmwh<-rxall[grepl("ENOXAPARIN|NADROPARIN|TINZAPARIN",rxall$drug.name),c(1:3)]
noac<-rxall[grepl("APIXABAN|DABIGATRAN|EDOXABAN|RIVAROXABAN",rxall$drug.name),c(1:3)]

check_group0<-merge(df_cox[which(df_cox$group==0),],noac,by="reference.key",all.x=T)
check_group0<-check_group0[which(check_group0$index.date<check_group0$date.rx.start),]
check_group0<-arrange(check_group0,check_group0$reference.key,check_group0$date.rx.start)
check_group0<-check_group0[!duplicated(check_group0$reference.key),]
check_group0<-check_group0[,c("reference.key","date.rx.start")]
colnames(check_group0)[2]<-"tr.change"

check_group1<-merge(df_cox[which(df_cox$group==1),],lmwh,by="reference.key",all.x=T)
check_group1<-check_group1[which(check_group1$index.date<check_group1$date.rx.start),]
check_group1<-arrange(check_group1,check_group1$reference.key,check_group1$date.rx.start)
check_group1<-check_group1[!duplicated(check_group1$reference.key),]
check_group1<-check_group1[,c("reference.key","date.rx.start")]
colnames(check_group1)[2]<-"tr.change2"

df_cox<-merge(df_cox,check_group0,by="reference.key",all.x = T)
df_cox<-merge(df_cox,check_group1,by="reference.key",all.x = T)
df_cox[which(is.na(df_cox$tr.change)==T),]$tr.change<-df_cox[which(is.na(df_cox$tr.change)==T),]$tr.change2
df_cox<-df_cox[,-ncol(df_cox)]

df_cox$tr.end<-NA
df_cox$tr.end<-as.Date(df_cox$tr.end)
df_cox[which(df_cox$group==1),]$tr.end<-df_cox[which(df_cox$group==1),]$date.noac.e
df_cox[which(df_cox$group==0),]$tr.end<-df_cox[which(df_cox$group==0),]$date.lmwh.e

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1

df_cox<- df_cox %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.dvt.date)%>%
  mutate(censor.date.2 = pmin(outcome.dvt.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1

df_cox<- df_cox %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.pe.date)%>%
  mutate(censor.date.3 = pmin(outcome.pe.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
df_cox[which(df_cox$outcome.date.3>df_cox$censor.date.3),]$outcome.date.3<-NA
df_cox[which(is.na(df_cox$outcome.date.3)==F),]$outcome.3<-1

df_cox<- df_cox %>%
  mutate(outcome.4 = 0)%>%
  mutate(outcome.date.4 = outcome.mb.date)%>%
  mutate(censor.date.4 = pmin(outcome.mb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.4 = censor.date.4 - index.date)
df_cox[which(df_cox$outcome.date.4>df_cox$censor.date.4),]$outcome.date.4<-NA
df_cox[which(is.na(df_cox$outcome.date.4)==F),]$outcome.4<-1

df_cox<- df_cox %>%
  mutate(outcome.5 = 0)%>%
  mutate(outcome.date.5 = outcome.ich.date)%>%
  mutate(censor.date.5 = pmin(outcome.ich.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.5 = censor.date.5 - index.date)
df_cox[which(df_cox$outcome.date.5>df_cox$censor.date.5),]$outcome.date.5<-NA
df_cox[which(is.na(df_cox$outcome.date.5)==F),]$outcome.5<-1

df_cox<- df_cox %>%
  mutate(outcome.6 = 0)%>%
  mutate(outcome.date.6 = outcome.gib.date)%>%
  mutate(censor.date.6 = pmin(outcome.gib.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.6 = censor.date.6 - index.date)
df_cox[which(df_cox$outcome.date.6>df_cox$censor.date.6),]$outcome.date.6<-NA
df_cox[which(is.na(df_cox$outcome.date.6)==F),]$outcome.6<-1

df_cox<- df_cox %>%
  mutate(outcome.7 = 0)%>%
  mutate(outcome.date.7 = outcome.omb.date)%>%
  mutate(censor.date.7 = pmin(outcome.omb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.7 = censor.date.7 - index.date)
df_cox[which(df_cox$outcome.date.7>df_cox$censor.date.7),]$outcome.date.7<-NA
df_cox[which(is.na(df_cox$outcome.date.7)==F),]$outcome.7<-1

df_cox<- df_cox %>%
  mutate(outcome.8 = 0)%>%
  mutate(outcome.date.8 = outcome.death.date)%>%
  mutate(censor.date.8 = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),tr.end,tr.change,na.rm = T)) %>% 
  mutate(follow.up.8 = censor.date.8 - index.date)
df_cox[which(df_cox$outcome.date.8>df_cox$censor.date.8),]$outcome.date.8<-NA
df_cox[which(is.na(df_cox$outcome.date.8)==F),]$outcome.8<-1

result<-data.frame(matrix(NA,nrow = 8,ncol = 8))
rownames(result)<-paste0(colnames(df_cox)[19:26]," 6 months")

colnames(df_cox)
for (k in c(1:8)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,70+2+4*(k-1)])==F),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),70+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,70+2+4*(k-1)])==F),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),70+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,70+4+4*(k-1)], df_cox[,70+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])

df_cox$pstat.1<-df_cox$outcome.1
df_cox$etime.1<-df_cox$follow.up.1
df_cox$event.1<-0
df_cox[which(df_cox$outcome.1==1),]$event.1<-1
df_cox[which(df_cox$outcome.1==0&df_cox$outcome.8==1),]$event.1<-2

df_cox$pstat.2<-df_cox$outcome.2
df_cox$etime.2<-df_cox$follow.up.2
df_cox$event.2<-0
df_cox[which(df_cox$outcome.2==1),]$event.2<-1
df_cox[which(df_cox$outcome.2==0&df_cox$outcome.8==1),]$event.2<-2

df_cox$pstat.3<-df_cox$outcome.3
df_cox$etime.3<-df_cox$follow.up.3
df_cox$event.3<-0
df_cox[which(df_cox$outcome.3==1),]$event.3<-1
df_cox[which(df_cox$outcome.3==0&df_cox$outcome.8==1),]$event.3<-2

df_cox$pstat.4<-df_cox$outcome.4
df_cox$etime.4<-df_cox$follow.up.4
df_cox$event.4<-0
df_cox[which(df_cox$outcome.4==1),]$event.4<-1
df_cox[which(df_cox$outcome.4==0&df_cox$outcome.8==1),]$event.4<-2

df_cox$pstat.5<-df_cox$outcome.5
df_cox$etime.5<-df_cox$follow.up.5
df_cox$event.5<-0
df_cox[which(df_cox$outcome.5==1),]$event.5<-1
df_cox[which(df_cox$outcome.5==0&df_cox$outcome.8==1),]$event.5<-2

df_cox$pstat.6<-df_cox$outcome.6
df_cox$etime.6<-df_cox$follow.up.6
df_cox$event.6<-0
df_cox[which(df_cox$outcome.6==1),]$event.6<-1
df_cox[which(df_cox$outcome.6==0&df_cox$outcome.8==1),]$event.6<-2

df_cox$pstat.7<-df_cox$outcome.7
df_cox$etime.7<-df_cox$follow.up.7
df_cox$event.7<-0
df_cox[which(df_cox$outcome.7==1),]$event.7<-1
df_cox[which(df_cox$outcome.7==0&df_cox$outcome.8==1),]$event.7<-2

k<-df_cox[,colnames(select(df_cox,starts_with("bs"),"sex","age","group"))]
colnames(df_cox)

for (i in c(1:7)) {
  fit <- crr((df_cox[,102+2+3*(i-1)]), df_cox[,102+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
result[8,c(7,8)]<-"-"

setwd("C:/Users/LabPC12CSMPR/Desktop")
write.xlsx(result,file="HR result-s8.xlsx",rowNames=T)

colnames(baseline_cohort)
df_cox<-baseline_cohort
df_cox$tr.end<-NA
df_cox$tr.end<-as.Date(df_cox$tr.end)
df_cox[which(df_cox$group==1),]$tr.end<-df_cox[which(df_cox$group==1),]$date.noac.e
df_cox[which(df_cox$group==0),]$tr.end<-df_cox[which(df_cox$group==0),]$date.lmwh.e
df_cox<-df_cox[which(df_cox$dur.rx>=30&df_cox$dur.rx<180),]
nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
summary(df_cox[which(df_cox$group==0),]$dur.rx)
summary(df_cox[which(df_cox$group==1),]$dur.rx)
k<-colnames(select(df_cox,starts_with("bs")))
rm<-NA
for(i in k){
  x<-unique(df_cox[,i])
  if(length(x)==1){
    print(i)
    rm<-c(rm,i)
  }
}

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.vte.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)
df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1

df_cox<- df_cox %>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.dvt.date)%>%
  mutate(censor.date.2 = pmin(outcome.dvt.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1

df_cox<- df_cox %>%
  mutate(outcome.3 = 0)%>%
  mutate(outcome.date.3 = outcome.pe.date)%>%
  mutate(censor.date.3 = pmin(outcome.pe.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.3 = censor.date.3 - index.date)
df_cox[which(df_cox$outcome.date.3>df_cox$censor.date.3),]$outcome.date.3<-NA
df_cox[which(is.na(df_cox$outcome.date.3)==F),]$outcome.3<-1

df_cox<- df_cox %>%
  mutate(outcome.4 = 0)%>%
  mutate(outcome.date.4 = outcome.mb.date)%>%
  mutate(censor.date.4 = pmin(outcome.mb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.4 = censor.date.4 - index.date)
df_cox[which(df_cox$outcome.date.4>df_cox$censor.date.4),]$outcome.date.4<-NA
df_cox[which(is.na(df_cox$outcome.date.4)==F),]$outcome.4<-1

df_cox<- df_cox %>%
  mutate(outcome.5 = 0)%>%
  mutate(outcome.date.5 = outcome.ich.date)%>%
  mutate(censor.date.5 = pmin(outcome.ich.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.5 = censor.date.5 - index.date)
df_cox[which(df_cox$outcome.date.5>df_cox$censor.date.5),]$outcome.date.5<-NA
df_cox[which(is.na(df_cox$outcome.date.5)==F),]$outcome.5<-1

df_cox<- df_cox %>%
  mutate(outcome.6 = 0)%>%
  mutate(outcome.date.6 = outcome.gib.date)%>%
  mutate(censor.date.6 = pmin(outcome.gib.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.6 = censor.date.6 - index.date)
df_cox[which(df_cox$outcome.date.6>df_cox$censor.date.6),]$outcome.date.6<-NA
df_cox[which(is.na(df_cox$outcome.date.6)==F),]$outcome.6<-1

df_cox<- df_cox %>%
  mutate(outcome.7 = 0)%>%
  mutate(outcome.date.7 = outcome.omb.date)%>%
  mutate(censor.date.7 = pmin(outcome.omb.date, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.7 = censor.date.7 - index.date)
df_cox[which(df_cox$outcome.date.7>df_cox$censor.date.7),]$outcome.date.7<-NA
df_cox[which(is.na(df_cox$outcome.date.7)==F),]$outcome.7<-1

df_cox<- df_cox %>%
  mutate(outcome.8 = 0)%>%
  mutate(outcome.date.8 = outcome.death.date)%>%
  mutate(censor.date.8 = pmin(outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.8 = censor.date.8 - index.date)
df_cox[which(df_cox$outcome.date.8>df_cox$censor.date.8),]$outcome.date.8<-NA
df_cox[which(is.na(df_cox$outcome.date.8)==F),]$outcome.8<-1

df_cox<- df_cox %>%
  mutate(outcome.9 = 0)%>%
  mutate(outcome.date.9 = outcome.vte.date)%>%
  mutate(censor.date.9 = pmin(outcome.vte.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.9 = censor.date.9 - index.date)
df_cox[which(df_cox$outcome.date.9>df_cox$censor.date.9),]$outcome.date.9<-NA
df_cox[which(is.na(df_cox$outcome.date.9)==F),]$outcome.9<-1

df_cox<- df_cox %>%
  mutate(outcome.10 = 0)%>%
  mutate(outcome.date.10 = outcome.dvt.date)%>%
  mutate(censor.date.10 = pmin(outcome.dvt.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.10 = censor.date.10 - index.date)
df_cox[which(df_cox$outcome.date.10>df_cox$censor.date.10),]$outcome.date.10<-NA
df_cox[which(is.na(df_cox$outcome.date.10)==F),]$outcome.10<-1

df_cox<- df_cox %>%
  mutate(outcome.11 = 0)%>%
  mutate(outcome.date.11 = outcome.pe.date)%>%
  mutate(censor.date.11 = pmin(outcome.pe.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.11 = censor.date.11 - index.date)
df_cox[which(df_cox$outcome.date.11>df_cox$censor.date.11),]$outcome.date.11<-NA
df_cox[which(is.na(df_cox$outcome.date.11)==F),]$outcome.11<-1

df_cox<- df_cox %>%
  mutate(outcome.12 = 0)%>%
  mutate(outcome.date.12 = outcome.mb.date)%>%
  mutate(censor.date.12 = pmin(outcome.mb.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.12 = censor.date.12 - index.date)
df_cox[which(df_cox$outcome.date.12>df_cox$censor.date.12),]$outcome.date.12<-NA
df_cox[which(is.na(df_cox$outcome.date.12)==F),]$outcome.12<-1

df_cox<- df_cox %>%
  mutate(outcome.13 = 0)%>%
  mutate(outcome.date.13 = outcome.ich.date)%>%
  mutate(censor.date.13 = pmin(outcome.ich.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.13 = censor.date.13 - index.date)
df_cox[which(df_cox$outcome.date.13>df_cox$censor.date.13),]$outcome.date.13<-NA
df_cox[which(is.na(df_cox$outcome.date.13)==F),]$outcome.13<-1

df_cox<- df_cox %>%
  mutate(outcome.14 = 0)%>%
  mutate(outcome.date.14 = outcome.gib.date)%>%
  mutate(censor.date.14 = pmin(outcome.gib.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.14 = censor.date.14 - index.date)
df_cox[which(df_cox$outcome.date.14>df_cox$censor.date.14),]$outcome.date.14<-NA
df_cox[which(is.na(df_cox$outcome.date.14)==F),]$outcome.14<-1

df_cox<- df_cox %>%
  mutate(outcome.15 = 0)%>%
  mutate(outcome.date.15 = outcome.omb.date)%>%
  mutate(censor.date.15 = pmin(outcome.omb.date, outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.15 = censor.date.15 - index.date)
df_cox[which(df_cox$outcome.date.15>df_cox$censor.date.15),]$outcome.date.15<-NA
df_cox[which(is.na(df_cox$outcome.date.15)==F),]$outcome.15<-1

df_cox<- df_cox %>%
  mutate(outcome.16 = 0)%>%
  mutate(outcome.date.16 = outcome.death.date)%>%
  mutate(censor.date.16 = pmin(outcome.death.date,index.date+365,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.16 = censor.date.16 - index.date)
df_cox[which(df_cox$outcome.date.16>df_cox$censor.date.16),]$outcome.date.16<-NA
df_cox[which(is.na(df_cox$outcome.date.16)==F),]$outcome.16<-1

result<-data.frame(matrix(NA,nrow = 16,ncol = 8))
rownames(result)<-c(paste0(colnames(df_cox)[22:29]," 6 months"),paste0(colnames(df_cox)[22:29]," 1 year"))

colnames(df_cox)
for (k in c(1:16)) {
  b<-nrow(df_cox[which(df_cox$group==0&is.na(df_cox[,72+2+4*(k-1)])==F),])
  c<-as.integer(sum(df_cox[which(df_cox$group==0),72+4+4*(k-1)])/365.25)
  d<-round(b/c,2)
  result[k,1]<-paste0(b,"/",d)
  e<-nrow(df_cox[which(df_cox$group==1&is.na(df_cox[,72+2+4*(k-1)])==F),])
  f<-as.integer(sum(df_cox[which(df_cox$group==1),72+4+4*(k-1)])/365.25)
  g<-round(e/f,2)
  result[k,2]<-paste0(e,"/",g)
  cox_result <- coxph(Surv(df_cox[,72+4+4*(k-1)], df_cox[,72+1+4*(k-1)]) ~ group, data = df_cox)
  x<-summary(cox_result)
  result[k,3]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
  result[k,4]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
  print(x$wald["pvalue"])
}

formula<-str_c(c("group ~ age+sex",setdiff(colnames(select(df_cox,starts_with("bs"))),rm)), collapse = "+")
W.glm <- weightit(formula(formula),data = df_cox, estimand="ATE", method="ps")
design <- svydesign(ids = ~ 1, data = df_cox, weights = ~ W.glm$weights)

cox_result <- svycoxph(Surv(follow.up.1, outcome.1) ~ group, design=design)
x<-summary(cox_result)
result[1,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[1,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.2, outcome.2) ~ group, design=design)
x<-summary(cox_result)
result[2,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[2,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.3, outcome.3) ~ group, design=design)
x<-summary(cox_result)
result[3,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[3,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.4, outcome.4) ~ group, design=design)
x<-summary(cox_result)
result[4,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[4,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.5, outcome.5) ~ group, design=design)
x<-summary(cox_result)
result[5,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[5,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.6, outcome.6) ~ group, design=design)
x<-summary(cox_result)
result[6,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[6,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.7, outcome.7) ~ group, design=design)
x<-summary(cox_result)
result[7,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[7,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.8, outcome.8) ~ group, design=design)
x<-summary(cox_result)
result[8,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[8,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
print(x$wald["pvalue"])
cox_result <- svycoxph(Surv(follow.up.9, outcome.9) ~ group, design=design)
x<-summary(cox_result)
result[9,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[9,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.10, outcome.10) ~ group, design=design)
x<-summary(cox_result)
result[10,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[10,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.11, outcome.11) ~ group, design=design)
x<-summary(cox_result)
result[11,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[11,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.12, outcome.12) ~ group, design=design)
x<-summary(cox_result)
result[12,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[12,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.13, outcome.13) ~ group, design=design)
x<-summary(cox_result)
result[13,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[13,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.14, outcome.14) ~ group, design=design)
x<-summary(cox_result)
result[14,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[14,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.15, outcome.15) ~ group, design=design)
x<-summary(cox_result)
result[15,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[15,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))
cox_result <- svycoxph(Surv(follow.up.16, outcome.16) ~ group, design=design)
x<-summary(cox_result)
result[16,5]<-paste0(round(x$coef[2],2)," (",round(x$conf.int[,"lower .95"],2),"-",round(x$conf.int[,"upper .95"],2),")")
result[16,6]<-ifelse(x$wald["pvalue"]<0.001,"<0.001",round(x$wald["pvalue"],3))

df_cox$pstat.1<-df_cox$outcome.1
df_cox$etime.1<-df_cox$follow.up.1
df_cox$event.1<-0
df_cox[which(df_cox$outcome.1==1),]$event.1<-1
df_cox[which(df_cox$outcome.1==0&df_cox$outcome.8==1),]$event.1<-2

df_cox$pstat.2<-df_cox$outcome.2
df_cox$etime.2<-df_cox$follow.up.2
df_cox$event.2<-0
df_cox[which(df_cox$outcome.2==1),]$event.2<-1
df_cox[which(df_cox$outcome.2==0&df_cox$outcome.8==1),]$event.2<-2

df_cox$pstat.3<-df_cox$outcome.3
df_cox$etime.3<-df_cox$follow.up.3
df_cox$event.3<-0
df_cox[which(df_cox$outcome.3==1),]$event.3<-1
df_cox[which(df_cox$outcome.3==0&df_cox$outcome.8==1),]$event.3<-2

df_cox$pstat.4<-df_cox$outcome.4
df_cox$etime.4<-df_cox$follow.up.4
df_cox$event.4<-0
df_cox[which(df_cox$outcome.4==1),]$event.4<-1
df_cox[which(df_cox$outcome.4==0&df_cox$outcome.8==1),]$event.4<-2

df_cox$pstat.5<-df_cox$outcome.5
df_cox$etime.5<-df_cox$follow.up.5
df_cox$event.5<-0
df_cox[which(df_cox$outcome.5==1),]$event.5<-1
df_cox[which(df_cox$outcome.5==0&df_cox$outcome.8==1),]$event.5<-2

df_cox$pstat.6<-df_cox$outcome.6
df_cox$etime.6<-df_cox$follow.up.6
df_cox$event.6<-0
df_cox[which(df_cox$outcome.6==1),]$event.6<-1
df_cox[which(df_cox$outcome.6==0&df_cox$outcome.8==1),]$event.6<-2

df_cox$pstat.7<-df_cox$outcome.7
df_cox$etime.7<-df_cox$follow.up.7
df_cox$event.7<-0
df_cox[which(df_cox$outcome.7==1),]$event.7<-1
df_cox[which(df_cox$outcome.7==0&df_cox$outcome.8==1),]$event.7<-2

df_cox$pstat.9<-df_cox$outcome.9
df_cox$etime.9<-df_cox$follow.up.9
df_cox$event.9<-0
df_cox[which(df_cox$outcome.9==1),]$event.9<-1
df_cox[which(df_cox$outcome.9==0&df_cox$outcome.8==1),]$event.9<-2

df_cox$pstat.10<-df_cox$outcome.10
df_cox$etime.10<-df_cox$follow.up.10
df_cox$event.10<-0
df_cox[which(df_cox$outcome.10==1),]$event.10<-1
df_cox[which(df_cox$outcome.10==0&df_cox$outcome.8==1),]$event.10<-2

df_cox$pstat.11<-df_cox$outcome.11
df_cox$etime.11<-df_cox$follow.up.11
df_cox$event.11<-0
df_cox[which(df_cox$outcome.11==1),]$event.11<-1
df_cox[which(df_cox$outcome.11==0&df_cox$outcome.8==1),]$event.11<-2

df_cox$pstat.12<-df_cox$outcome.12
df_cox$etime.12<-df_cox$follow.up.12
df_cox$event.12<-0
df_cox[which(df_cox$outcome.12==1),]$event.12<-1
df_cox[which(df_cox$outcome.12==0&df_cox$outcome.8==1),]$event.12<-2

df_cox$pstat.13<-df_cox$outcome.13
df_cox$etime.13<-df_cox$follow.up.13
df_cox$event.13<-0
df_cox[which(df_cox$outcome.13==1),]$event.13<-1
df_cox[which(df_cox$outcome.13==0&df_cox$outcome.8==1),]$event.13<-2

df_cox$pstat.14<-df_cox$outcome.14
df_cox$etime.14<-df_cox$follow.up.14
df_cox$event.14<-0
df_cox[which(df_cox$outcome.14==1),]$event.14<-1
df_cox[which(df_cox$outcome.14==0&df_cox$outcome.8==1),]$event.14<-2

df_cox$pstat.15<-df_cox$outcome.15
df_cox$etime.15<-df_cox$follow.up.15
df_cox$event.15<-0
df_cox[which(df_cox$outcome.15==1),]$event.15<-1
df_cox[which(df_cox$outcome.15==0&df_cox$outcome.8==1),]$event.15<-2

k<-df_cox[,colnames(select(df_cox,starts_with("bs"),"sex","age","group"))]
colnames(df_cox)

for (i in c(1:7)) {
  fit <- crr((df_cox[,136+2+3*(i-1)]), df_cox[,136+3+3*(i-1)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.001,"<0.001",round(x$coef[nrow(x$coef),5],3))
  print(x$coef[nrow(x$coef),5])
}
for (i in c(9:15)) {
  fit <- crr((df_cox[,136+2+3*(i-2)]), df_cox[,136+3+3*(i-2)], k,failcode=1,cencode=0)
  x<-summary(fit)
  result[i,7]<-paste0(round(x$conf.int[nrow(x$conf.int),1],2)," (",round(x$conf.int[nrow(x$conf.int),3],2),"-",round(x$conf.int[nrow(x$conf.int),4],2),")")
  result[i,8]<-ifelse(x$coef[nrow(x$coef),5]<0.05,"<0.05",round(x$coef[nrow(x$coef),5],2))
}
result[c(8,16),c(7,8)]<-"-"
setwd("C:/Users/LabPC12CSMPR/Desktop")
write.xlsx(result,file="HR result-s10.xlsx")

df_cox<-baseline_cohort

nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
df_cox$tr.end<-NA
df_cox$tr.end<-as.Date(df_cox$tr.end)
df_cox[which(df_cox$group==1),]$tr.end<-df_cox[which(df_cox$group==1),]$date.noac.e
df_cox[which(df_cox$group==0),]$tr.end<-df_cox[which(df_cox$group==0),]$date.lmwh.e

df_cox<- df_cox %>%
  mutate(outcome.1 = 0)%>%
  mutate(outcome.date.1 = outcome.vte.date)%>%
  mutate(censor.date.1 = pmin(outcome.date.1, outcome.death.date,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.1 = censor.date.1 - index.date)%>%
  mutate(outcome.2 = 0)%>%
  mutate(outcome.date.2 = outcome.vte.date)%>%
  mutate(censor.date.2 = pmin(outcome.date.2, outcome.death.date,tr.end,index.date+180,as.Date("2022-12-31"),na.rm = T)) %>% 
  mutate(follow.up.2 = censor.date.2 - index.date)

df_cox[which(df_cox$outcome.date.1>df_cox$censor.date.1),]$outcome.date.1<-NA
df_cox[which(is.na(df_cox$outcome.date.1)==F),]$outcome.1<-1
df_cox[which(df_cox$outcome.date.2>df_cox$censor.date.2),]$outcome.date.2<-NA
df_cox[which(is.na(df_cox$outcome.date.2)==F),]$outcome.2<-1
df_cox<-df_cox[which(df_cox$follow.up.2>0),]
df_cox$basegamma<-NA
df_cox[which(df_cox$follow.up.2<df_cox$follow.up.1),]$basegamma<-1

nrow(df_cox[which(df_cox$group==0),])
nrow(df_cox[which(df_cox$group==1),])
set.seed(1)
imputed.data.sets <- gammaImpute(formula=Surv(follow.up.2,outcome.2)~group+age+sex+bs.dx.Lip.oral.cavity.and.pharynx+bs.dx.Digestive.organs+
                                   bs.dx.Respiratory.system+bs.dx.Bone.skin.and.soft.tissue+bs.dx.Breast.and.genital.organs+bs.dx.Urinary.organs+
                                   bs.dx.Eye.brain.and.other.central.nervous.system.endocrine.glands+bs.dx.Lymphatic.and.hematopoietic.tissue+
                                   bs.dx.Metastasis+bs.dx.Obesity+bs.dx.Tobacco.use.disorder+bs.dx.Alcohol.use.disorder+bs.dx.Drug.use.disorder+
                                   bs.dx.Diabetes+bs.dx.Hypertension+bs.dx.Hyperlipidemia+bs.dx.Atrial.fibrillation+bs.dx.Congestive.heart.failure+
                                   bs.dx.Vascular.disease+bs.dx.Renal.disease+bs.tr.Drug.therapy+bs.tr.Radiotherapy+bs.rx.aspirin+bs.rx.Antiplatelet+
                                   bs.rx.NSAIDs+
                                   bs.rx.Erythropoietin+
                                   bs.rx.EGFR.inhibitor+
                                   bs.rx.VEGF.inhibitor+bs.rx.CYP3A4.Pgp+bs.rx.SSRIs+
                                   bs.Khorana+bs.CCI+bs.px.catheter+bs.px.Surgery+bs.px.transfusion+bs.ip+bs.ae,
                                 data = df_cox, m=20, gamma = "basegamma",
                                 gamma.factor = 0,DCO.time=180)
fits <- ImputeStat(imputed.data.sets)
y<-summary(fits)
paste0(round(exp(y[1,1]),2)," (",round(exp(y[1,6]),2),"-",round(exp(y[1,7]),2),")")
ifelse(y[1,5]<0.001,"<0.001",round(y[1,5],3))
print(y[1,5])