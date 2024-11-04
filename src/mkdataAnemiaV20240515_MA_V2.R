#------------------------------------------------------------------------------
# Data for HC Anemia Study (HCA)
# 
# THL data updated
# Created: 15-06-2023
# Updated: 23-02-2024, 15-05-2024, 13-08-2024 (manuscript V1)
# Authors: JKH, SE
#
#------------------------------------------------------------------------------


#OUTPUT FOLDER
out <- "W:/Veripalvelu/mikko/anemia_and_hc/out20240821/" # manuscript V1, fixed thrombosis
library(dplyr)


#------------------------
# 
# Load data
# 
#------------------------

# Load study population, only these individuals included in analyses
load(file="w:/HormonalC/HCdata2.RData")

# Load updated HILMO and Avo-HILMO data
load("W:/UpdateData/THL0224hilmo.RData")
load("W:/UpdateData/THL0224avohilmo.RData")

# Deliveries 1987-2022
load("W:/UpdateData/rawdeliv202311V2.RData")

# Abortions 1987-2022
load(file="w:/UpdateData/rawab202311V2.RData")

# Deaths
load(file="W:/UpdateData/rawdea2.RData")# until 31-12-2019

# Cancer
load("W:/Cancer/tmpcandt.RData")

load("W:/UpdateData/kelaerit2021.RData")

# Libraries

library(codeCollection)
library(Epi)
library(mapsFinland)

head(ICPCKoodit,2)
subset(ICPCKoodit,grepl("Anem|anem",fi))
View(subset(ICD10Koodit,grepl("Anem|anem",fi))[,1:3])
#------------------------



source('W:/multipleNCC/includefunctionsNCC.R', encoding = 'UTF-8')
source('W:/Anemia/mkHCdtV20240515.R', encoding = 'UTF-8')


#------------------------
# Make incidence data
# Here both HILMO and AvoHilmo
#------------------------

# lv.hilmo<-merge(x=HILMO.henk,
#                 y=lv.hilmo,
#                 by="shilmo_id",all.x=FALSE,all.y=TRUE)
# head(thl2023_5525_hilmo,2)


# cbind(substr(thl2023_5525_hilmo$TUPVA,1,10)[1:4],
# Epi::cal.yr(substr(thl2023_5525_hilmo$TUPVA,1,10)[1:4],format="%d.%m.%Y"))
# 
# cbind(substr(thl2023_5525_avohilmo$KAYNTI_ALKOI,1,10)[1:4],
#       Epi::cal.yr(substr(thl2023_5525_avohilmo$KAYNTI_ALKOI,1,10)[1:4],format="%d.%m.%Y"))

thl2023_5525_hilmo$TUPVA.pvm<-Epi::cal.yr(substr(thl2023_5525_hilmo$TUPVA,1,10),format="%d.%m.%Y")
thl2023_5525_avohilmo$KAYNTI_ALKOI.pvm<-Epi::cal.yr(substr(thl2023_5525_avohilmo$KAYNTI_ALKOI,1,10),format="%d.%m.%Y")

# SOF: 2019, EOF: 2020.999 (no drug purchases after this)
system.time(apu<-mk.HCdt(
  ICD10.dgn = "^D50",ICPC.code = "B80",
  popul=HC.data.2,
  HILMO.ICD=thl2023_5525_hilmo_icd10,
  HILMO.henk=thl2023_5525_hilmo,
  Avo.HILMO.ICD=thl2023_5525_ah_icd10,
  Avo.HILMO.ICPC=thl2023_5525_ah_icpc2,
  Avo.HILMO.henk=thl2023_5525_avohilmo,
  dea=raw.dea2,SOF=2019,EOF=2020.999)
)

save(apu,file=paste(out,"anemiaincidencedata.RData",sep = ""))
names(apu)

head(subset(apu$inci.dt,ep.type=="EOF"),4)
head(apu$inci.dt,3)
addmargins(table(apu$inci.dt$ep.type,apu$inci.dt$.i.prevalent))
# Flowchart for incidence
DiagrammeR::grViz(apu$flow)


# Nested case-control design with five controls
library(Epi)

apu.dt<-subset(apu$inci.dt)
with(apu.dt,tapply(first.event.pvm,ep.type,range))
with(apu.dt,tapply(SOF,ep.type,range))
with(apu.dt,tapply(EOF,ep.type,range))


tmp.ma<-match(apu.dt$shnro, HC.data.2$shnro)
apu.dt$kunta<-HC.data.2$kunta[tmp.ma]
apu.dt$ututkuR<-HC.data.2$ututkuR[tmp.ma]
apu.dt$soseR<-HC.data.2$soseR[tmp.ma]
apu.dt$syntyv<-HC.data.2$syntyv[tmp.ma]
apu.dt$kuntaC<-HC.data.2$kuntaC[tmp.ma]
apu.dt$sivsr<-HC.data.2$sivsr[tmp.ma]
apu.dt$ptoim1<-HC.data.2$ptoim1[tmp.ma]
apu.dt$seutukunta<-HC.data.2$seutukunta[tmp.ma]

# library(mapsFinland)
# cat(paste(names(HC.data.2),"=",names(HC.data.2),",",sep=""))
# head(subset(apu.dt,ep.type=="ICPC"),4)
# head(subset(apu.dt,ep.type%in%c("ICPC", "ICD10.hilmo", "ICD10.avohilmo")&(first.event.pvm>2019.99)))

# table(apu.dt$first.event.pvm==2019.999,apu.dt$ep.type)

set.seed(20231012+1)
apu.first.event<-with(apu.dt,ifelse(ep.type%in%c("ICPC", "ICD10.hilmo", "ICD10.avohilmo"),
                                    first.event.pvm+runif(nrow(apu.dt),0.0001,0.00011),
                                    first.event.pvm))
apu.dt$first.event.pvm.noise<-apu.first.event

set.seed(20231012)
system.time(
  apu.CC <- ccwc( entry=SOF, 
                  exit=first.event.pvm.noise, 
                  fail=ep.type%in%c("ICPC", "ICD10.hilmo", "ICD10.avohilmo"), 
                  origin=SOF,
                  controls=5, 
                  data=subset(apu.dt,!.i.prevalent),# Only non-prevalent
                  include=list(shnro=shnro,
                               pvm.ICPC=pvm.ICPC, 
                               ICPC2=ICPC2,
                               pvm.hilmo.ICD10=pvm.hilmo.ICD10,
                               hilmo.ICD10=hilmo.ICD10,
                               pvm.avohilmo.ICD10=pvm.avohilmo.ICD10,
                               avohilmo.ICD10=avohilmo.ICD10,
                               kunta=kunta, 
                               ptoim1=ptoim1, sivsr=sivsr,
                               soseR=soseR, ututkuR=ututkuR, 
                               kuntaC=kuntaC, 
                               seutukunta=seutukunta, 
                               syntyv=syntyv,
                               .i.prevalent=.i.prevalent
                  ), 
                  match=list(syntyv=syntyv,seutukunta=seutukunta))# matching variables
)



table(table(apu.CC$Set))


# Add NA as levels
apu.CC$soseR<-addNA(apu.CC$soseR)
apu.CC$ututkuR<-addNA(apu.CC$ututkuR)
levels(apu.CC$ututkuR)[7]<-"Unknown"
levels(Relevel(apu.CC$soseR,list(1,2,3,4,5,6,7,8:9)))
apu.CC$soseR<-Relevel(apu.CC$soseR,list(1,2,3,4,5,6,7,8:9))


# Assign class as factor
class(apu.CC$soseR)<-"factor"
class(apu.CC$ututkuR)<-"factor"

table(apu.CC$Fail)
length(unique(apu.CC$Set))
table(table(apu.CC$Set)) # Set size
table(table(apu.CC$shnro)) # Set size
range(apu.CC$Time)

save(apu,apu.dt,apu.CC,file=paste(out,"tempdata20240515.RData",sep=""))
#load(file="tempdata20240515.RData")

# Remove individuals pregnant in one year before EOF
# head(raw.deliv.2023.11.V2,2)
# head(raw.ab.2023.11.V2$TOIMENPIDE_PVM.pvm,2)
# summary(raw.ab.2023.11.V2[,c("RASK_KESTO_ARVIO" ,"RASK_KESTO_PV","RASK_KESTO_VK")])
raw.ab.2023.11.V2$preg.start<-with(raw.ab.2023.11.V2,TOIMENPIDE_PVM.pvm- (RASK_KESTO_PV+  RASK_KESTO_VK*7)/365.26)

tmp.poistetaan<-with(raw.ab.2023.11.V2,data.frame(shnro=shnro,alku=preg.start,loppu=TOIMENPIDE_PVM.pvm))
tmp.poistetaan<-with(raw.deliv.2023.11.V2,rbind(tmp.poistetaan,data.frame(shnro=shnro,alku=synnytyksen_pvm.pvm-9/12,loppu=synnytyksen_pvm.pvm)))
tmp.poistetaan<-subset(tmp.poistetaan,!is.na(alku)&!is.na(loppu)&(shnro%in%apu.CC$shnro))
tmp.poistetaan<-subset(tmp.poistetaan,loppu>2018)

tmp.poistetaan<-tmp.poistetaan[order(tmp.poistetaan$shnro,tmp.poistetaan$alku),]
apu.pois<-tapply(tmp.poistetaan$alku,tmp.poistetaan$shnro,function(x){seq(length(x))})
apu.pois.1<-do.call(what="c",apu.pois)

# First remove cases who were pregnant in one before EOF
apu.pois.lista<-vector(mode="list",length=length(unique(apu.pois.1)))
for(tmp.i in seq(apu.pois.lista)){
  apu.1<-merge(x=apu.CC,y=tmp.poistetaan[apu.pois.1==tmp.i,],by="shnro")
  apu.pois.lista[[tmp.i]]<-as.character(subset(apu.1,(Fail==1)&(alku<Time)&(loppu>Time))$shnro)
}
apu.pois.lista.1<-unique(unlist(apu.pois.lista))

length(apu.CC$Set[(apu.CC$shnro%in%apu.pois.lista.1)&(apu.CC$Fail==1)])# N=1113 pois
apu.CC.1<-subset(apu.CC,!(Set%in%apu.CC$Set[(apu.CC$shnro%in%apu.pois.lista.1)&(apu.CC$Fail==1)]))# N=1113 pois

# Then remove controls and cases who were pregnant in one year before EOF

apu.1<-merge(x=apu.CC,y=tmp.poistetaan[apu.pois.1==1,],by="shnro",all.x = TRUE)
apu.ulos<-subset(apu.1,ifelse(is.na(alku),TRUE,(alku<Time)&(loppu>Time)))
nrow(apu.CC)
nrow(apu.ulos)

for(tmp.i in seq(apu.pois.lista)[-1]){
  apu.1<-merge(x=apu.ulos[,!(names(apu.ulos)%in%c("alku","loppu"))],
               y=tmp.poistetaan[apu.pois.1==tmp.i,],
               by="shnro",all.x = TRUE,all.y = FALSE)
  apu.ulos<-subset(apu.1,ifelse(is.na(alku),TRUE,(alku<Time)&(loppu>Time)))
  print(nrow(apu.ulos))
}
 
apu<-tapply(apu.ulos$Fail,apu.ulos$Set,sum);table(apu)
apu1<-with(apu.ulos,duplicated(paste(Set,shnro,Fail)))
table(apu1)
apu.CC.1<-apu.ulos[!apu1,]
apu<-tapply(apu.CC.1$Fail,apu.CC.1$Set,sum);table(apu)
apu.CC.2<-subset(apu.CC.1,!(Set%in%names(apu)[apu==0]) )
apu<-tapply(apu.CC.2$Fail,apu.CC.2$Set,sum);table(apu)
table(apu.CC.2$Fail)
table(tapply(apu.CC.2$Fail,apu.CC.2$Set,length))

apu<-tapply(apu.CC.2$Fail==0,apu.CC.2$Set,sum);table(apu)
apu.CC.3<-subset(apu.CC.2,!(Set%in%names(apu)[apu==0]) )

table(tapply(apu.CC.3$Fail,apu.CC.3$Set,sum))
table(tapply(apu.CC.3$Fail==0,apu.CC.3$Set,sum))
table(tapply(apu.CC.3$Fail,apu.CC.3$Set,length))

 
# nrow(apu.CC);nrow(apu.CC.3)
# table(apu.CC$Fail);table(apu.CC.3$Fail)

anemia.flow<-mkflowchart(N=c(24726,17470), 
                      text.M="Pregnant", 
                      text.P=c("Anemia\\nCases (N=4121)\\nControls (N=20,605)",
                               "Anemia\\nCases (N=3316)\\nControls (N=14,154)"), 
                      type = 1)
DiagrammeR::grViz(anemia.flow)
save(anemia.flow,file=paste(out,"anemiaflow.RData",sep=""))

save(apu.CC.3,file=paste(out,"apuCC3.RData",sep=""))
#load(file="apuCC3.RData")

#-----------------------------------------------------------
# Adding covariates
#-----------------------------------------------------------

# Drug purchaces 2013-2020
load("W:/UpdateData/HCpurchV2021.RData") # 2018 and after, use this 

# ATC codes for HC
HC.atc.koodit<-paste("^",c("G02B","G03A","G03HB"),sep="",collapse="|")

tmp.purch<-subset(HC.purch.V2021,grepl(HC.atc.koodit,atc)&shnro%in%apu.CC.3$shnro)
tmp.purch$atc<-as.character(tmp.purch$atc)
dput(sort(unique((tmp.purch$atc))))

HC.atc.koodit.1<-c("G02BA03", "G02BB01", "G03AA", "G03AA07", "G03AA09", "G03AA10", 
  "G03AA11", "G03AA12", "G03AA13", "G03AA14", "G03AA16", "G03AB03", 
  "G03AB05", "G03AB06", "G03AB08", "G03AC01", "G03AC03", "G03AC06", 
  "G03AC08", "G03AC09", "G03AC10", "G03AD01", "G03AD02", "G03HB01")

 

apu.CC.3<-apu.CC.3[order(apu.CC.3$Set,apu.CC.3$Fail,apu.CC.3$shnro),]
apu.kerta<-tapply(apu.CC.3$Fail,apu.CC.3$shnro,function(x){seq(length(x))})
apu.kerta.1<-do.call(what="c",apu.kerta)

apu.2<-vector(mode="list",length = length(unique(apu.kerta.1)))
for(tmp.i in seq(apu.2)){apu.2[[tmp.i]]<-subset(apu.CC.3,apu.kerta.1==tmp.i)}
sapply(apu.2,nrow);table(apu.kerta.1)

apu.fun<-function(x,dt){
  lv.1<-subset(tmp.purch,(shnro%in%dt$shnro)&(atc==x))
  # print(nrow(lv.1))
  lv.1$Aika<-dt$Time[match(lv.1$shnro,dt$shnro)]
  lv.1<-subset(lv.1,(toim.pvm<Aika)&abs(toim.pvm-Aika)<=1) # ikkuna?
  if(nrow(lv.1)<1){return(rep(0,nrow(dt)))}
  lv.1<-table(lv.1$shnro)
  lv.1<-lv.1[match(dt$shnro,names(lv.1))]
  ifelse(is.na(lv.1),0,lv.1)
}

apu.3<-vector(mode="list",length = length(apu.2))
for(tmp.i in seq(apu.3)){apu.3[[tmp.i]]<-sapply(HC.atc.koodit.1,apu.fun,dt=apu.2[[tmp.i]])}
sapply(apu.2,dim)
sapply(apu.3,dim)

apply(apu.3[[1]],2,table)
tmp.atc.lkm<-do.call("rbind",apu.3)
apu<-(apply(tmp.atc.lkm,2,sum)>9)# ONly 10 or more
tmp.atc.lkm<-tmp.atc.lkm[,apu]
tmp.atc.i.<-(tmp.atc.lkm>0)*1
apply(tmp.atc.i.,2,sum)


# Special reimbursement in one year before EOF
apu.SK<-sort(unique(kela.erit.2021$SK1))

apu.fun1<-function(x,dt){
  lv.1<-subset(kela.erit.2021,(shnro%in%dt$shnro)&(SK1==x))
  lv.1$Aika<-dt$Time[match(lv.1$shnro,dt$shnro)]
  lv.1<-subset(lv.1,(alku.pvm<Aika)&(loppu.pvm>Aika))
  if(nrow(lv.1)<1){return(rep(0,nrow(dt)))}
  lv.1<-table(lv.1$shnro)
  lv.1<-lv.1[match(dt$shnro,names(lv.1))]
  ifelse(is.na(lv.1),0,lv.1)
}

apu.4<-vector(mode="list",length = length(apu.2))
for(tmp.i in seq(apu.4)){apu.4[[tmp.i]]<-sapply(apu.SK,apu.fun1,dt=apu.2[[tmp.i]])}
sapply(apu.2,dim);sapply(apu.4,dim)
tmp.erit<-do.call("rbind",apu.4)
tmp.erit<-(tmp.erit>0)*1
apply(tmp.erit,2,table)

save(tmp.erit,tmp.atc.lkm,tmp.atc.i.,file=paste(out,"erijatatc.RData",sep=''))

# Previous dgn, separately, (time checked by SE)
# table(grepl("D25|N80|N84|N85|N92|N93",thl2023_5525_hilmo_icd10$KOODI))
# table(grepl("G45",thl2023_5525_hilmo_icd10$KOODI))


# DIAGNOSIS
######
## Thrombosis, stroke, TIA from HILMO
######
tmp.thromb<-c("I80","I81","I82","I26","I21","I63", "G45")

table(grepl(paste(tmp.thromb,collapse = "|"),thl2023_5525_hilmo_icd10$KOODI))
# FALSE     TRUE 
# 13080243    42379  -> 42k koko datassa?

apu.hilmo<-subset(thl2023_5525_hilmo,shnro%in%apu.CC.3$shnro)
apu.hilmo<-merge(x=apu.hilmo,
                 y=subset(thl2023_5525_hilmo_icd10,grepl(paste(tmp.thromb,collapse = "|"),KOODI)),
                 by="HILMO_ID",all=FALSE)
table(apu.hilmo$KOODI)
# I2100 I2101 I2109 I2119 I2121 I2139 I2140 I2141 I2149 I2190 I2191   I26  I260  I269  I630  I631  I632  I634  I635  I636  I638  I639 
# 1   196     1     1     2    99     1    11     6     2     2     7   113   379     5     2     7    28    10    23    68   224 
# I80  I800  I801 I8020 I8029  I803  I808  I809   I82  I820  I821 I8220 I8221 I8280 I8288  I829   I84  I840  I841  I842  I843  I844 
# 5   103    61    29   121    42    53     9     1     1     2     1    29    28   109    34     3    11    53   186     9    30 
# I845  I846  I847  I848  I849 
# 52    39     5     8   101

# I2100 I2101 I2109 I2119 I2121 I2139 I2140 I2141 I2149 I2190 I2191   I26  I260  I269  I630  I631  I632  I634  I635  I636  I638  I639    
# 1   196     1     1     2    99     1    11     6     2     2     7   113   379     5     2     7    28    10    23    68   224        
#I80  I800  I801 I8020 I8029  I803  I808  I809   I81   I82   I820  I821 I8220 I8221 I8280 I8288  I829 
#5   103    61    29   121    42    53     9   137     1     1     2     1    29    28   109    34 
# TIA dg (G45) not available

sum(table(apu.hilmo$KOODI))
#1953 -> 2k kun rajataan apu.CC.3:n

apu.fun2<-function(x,dt){
#  lv.1<-subset(apu.hilmo,(shnro%in%dt$shnro)&(KOODI==x)) # This looks only exact matches
#> subset(apu.hilmo,(shnro%in%apu.2[[1]]$shnro)&(KOODI=="I80")) %>% nrow()
#[1] 5
  lv.1<-subset(apu.hilmo,(shnro%in%dt$shnro)&(grepl(pattern = x,x = KOODI)))
# While this finds all subcategories too
# > subset(apu.hilmo,(shnro%in%apu.2[[1]]$shnro)&(grepl(pattern = "I80",x = KOODI))) %>% nrow()
# [1] 423
  lv.1$Aika<-dt$Time[match(lv.1$shnro,dt$shnro)]
  lv.1<-subset(lv.1,(TUPVA.pvm<Aika)&abs(TUPVA.pvm-Aika)<5)# in last five years
  if(nrow(lv.1)<1){return(rep(0,nrow(dt)))}
  lv.1<-table(lv.1$shnro)
  lv.1<-lv.1[match(dt$shnro,names(lv.1))]
  ifelse(is.na(lv.1),0,(lv.1>0)*1)
}


#lv.1<-subset(apu.hilmo,(shnro%in%dt$shnro)&(KOODI==))




apu.5<-vector(mode="list",length = length(apu.2))
for(tmp.i in seq(apu.5)){apu.5[[tmp.i]]<-sapply(tmp.thromb,apu.fun2,dt=apu.2[[tmp.i]])}
sapply(apu.2,dim)
# [,1] [,2] [,3]
# [1,] 17183  255    3
# [2,]    24   24   24
sapply(apu.5,dim)
# [,1] [,2] [,3]
# [1,] 17183  255    3
# [2,]     6    6    6

tmp.thrombosis<-do.call("rbind",apu.5)
apply(tmp.thrombosis,2,sum)
#With (KOODI==x)
#I80 I84 I82 I26 I21 I63 
#4   2   0   3   0   0 
#With (grepl(pattern = x,x = KOODI)
# I80 I84 I82 I26 I21 I63 
# 84 136  23  28   3  21


######
## HMB
######
tmp.HMB<-c("D25","N80","N84","N85","N92","N93")
table(grepl(paste(tmp.HMB,collapse = "|"),thl2023_5525_hilmo_icd10$KOODI))
#FALSE     TRUE 
#12855438   267184 -> 267K koko datassa

apu.hilmo<-subset(thl2023_5525_hilmo,shnro%in%apu.CC.3$shnro)
apu.hilmo<-merge(x=apu.hilmo,
                # y=subset(thl2023_5525_hilmo_icd10,grepl(tmp.HMB,KOODI)),
                y=subset(thl2023_5525_hilmo_icd10,grepl(paste(tmp.HMB,collapse = "|"),KOODI)),
                 by="HILMO_ID",all=FALSE)
table(apu.hilmo$KOODI)
# D25  D250  D251  D252  D259   N80  N800  N801  N802  N803  N804  N805  N806 N8080 N8081 N8089  N809  N840  N841  N842  N843  N850 
# 20   362   426   151   921    19   230   935     9  1096   296   131    26    15    45    74   661   418    73     5     2   183 
# N851  N856  N857  N858  N859   N92  N920  N921  N922  N923  N924  N925  N926   N93  N930  N938  N939 
# 34     2    23     2    11    33  1470   963    49    12   154   296   562     9    87   101   143

#> apu.hilmo %>% filter(KOODI == 'N803' ) %>% inner_join(apu.2[[1]] %>% select(shnro), by=c("shnro" = "shnro")) %>% nrow()
#[1] 1082
# of the 1096 N803 persons 1082 are in apu.2[[1]]

# > apu.hilmo %>% filter(KOODI == 'N803' ) %>% inner_join(apu.2[[1]] %>% select(shnro), by=c("shnro" = "shnro")) %>% mutate(year = gsub(as.character(TUPVA.pvm),pattern = "\\..*",replacement = "") ) %>% count(year)
# year   n
# 1  2012  42
# 2  2013  54
# 3  2014  30
# 4  2015  68
# 5  2016  89
# 6  2017  79
# 7  2018  92
# 8  2019 120
# 9  2020 148
# 10 2021 169
# 11 2022 191


sum(table(apu.hilmo$KOODI))
#10049 -> 10k apu.CC.3:ssa

apu.5<-vector(mode="list",length = length(apu.2))
 for(tmp.i in seq(apu.5)){apu.5[[tmp.i]]<-sapply(tmp.HMB,apu.fun2,dt=apu.2[[tmp.i]])}
sapply(apu.2,dim)
# [,1] [,2] [,3]
# [1,] 17183  255    3
# [2,]    24   24   24
 
sapply(apu.5,dim)
# [,1] [,2] [,3]
# [1,] 17183  255    3
# [2,]     6    6    6

tmp.heavymb <- do.call("rbind",apu.5)
apply(tmp.heavymb,2,sum)
#With (KOODI==x)
# D25 N80 N84 N85 N92 N93 
# 2   1   0   0   5   1
#With (grepl(pattern = x,x = KOODI)
# D25 N80 N84 N85 N92 N93 
# 232 249 105  18 585  69 

###########################




######
## INFLAMMATION
######
tmp.infl<-c("E66")
table(grepl(paste(tmp.infl,collapse = "|"),thl2023_5525_hilmo_icd10$KOODI))
#FALSE     TRUE 
#12903209   219413 -> 219K koko datassa

apu.hilmo<-subset(thl2023_5525_hilmo,shnro%in%apu.CC.3$shnro)
apu.hilmo<-merge(x=apu.hilmo,
                 # y=subset(thl2023_5525_hilmo_icd10,grepl(tmp.infl,KOODI)),
                 y=subset(thl2023_5525_hilmo_icd10,grepl(paste(tmp.infl,collapse = "|"),KOODI)),
                 by="HILMO_ID",all=FALSE)
table(apu.hilmo$KOODI)
# E66  E660 E6600 E6601  E661  E662  E668  E669 
# 38    34   308  4556    10   221  3236   618 


#> apu.hilmo %>% filter(KOODI == 'E66' ) %>% inner_join(apu.2[[1]] %>% select(shnro), by=c("shnro" = "shnro")) %>% nrow()
#[1] 1082
# of the 1096 N803 persons 1082 are in apu.2[[1]]

# > apu.hilmo %>% filter(KOODI == 'N803' ) %>% inner_join(apu.2[[1]] %>% select(shnro), by=c("shnro" = "shnro")) %>% mutate(year = gsub(as.character(TUPVA.pvm),pattern = "\\..*",replacement = "") ) %>% count(year)
# year   n
# 1  2012  42
# 2  2013  54
# 3  2014  30
# 4  2015  68
# 5  2016  89
# 6  2017  79
# 7  2018  92
# 8  2019 120
# 9  2020 148
# 10 2021 169
# 11 2022 191


sum(table(apu.hilmo$KOODI))
#10049 -> 10k apu.CC.3:ssa

apu.5<-vector(mode="list",length = length(apu.2))
for(tmp.i in seq(apu.5)){apu.5[[tmp.i]]<-sapply(tmp.infl,apu.fun2,dt=apu.2[[tmp.i]])}
sapply(apu.2,dim)
# [,1] [,2] [,3]
# [1,] 17183  255    3
# [2,]    24   24   24

sapply(apu.5,dim)
# [,1] [,2] [,3]
# [1,] 17183  255    3
# [2,]     6    6    6

tmp.inflammation <- do.call("rbind",apu.5)
apply(tmp.inflammation,2,sum)
#With (KOODI==x)

#With (grepl(pattern = x,x = KOODI)


###########################





table(thl2023_5525_ah_icpc2$ICPC2)


# Combining data  into one data.frame

anemia.dt<-do.call(rbind,apu.2)
dim(anemia.dt);dim(apu.CC.3)
tmp.atc.lkm<-do.call("rbind",apu.3)
colnames(tmp.atc.lkm)<-paste0(colnames(tmp.atc.lkm),".lkm")
anemia.dt<-cbind(anemia.dt,tmp.atc.lkm)
apu<-(tmp.atc.lkm>1)*1
colnames(apu)<-gsub(x=colnames(apu),pattern = ".lkm",replacement = ".i.2",fixed=TRUE)
apply(apu,2,table)
anemia.dt<-cbind(anemia.dt,apu)


colnames(tmp.erit)<-paste0("SK",apu.SK)
apply(tmp.erit,2,table)
anemia.dt<-cbind(anemia.dt,tmp.erit)
anemia.dt<-cbind(anemia.dt,tmp.thrombosis)
anemia.dt<-cbind(anemia.dt,tmp.heavymb)
anemia.dt<-cbind(anemia.dt,tmp.inflammation)




# Combine substances

apu<-grepl(paste(tmp.atc.estradiaol,collapse="|"),colnames(tmp.atc.lkm))
dim(tmp.atc.lkm[,apu]);apu.estradiaol<-apply(tmp.atc.lkm[,apu],1,sum)

apu<-grepl(paste(tmp.atc.EE,collapse="|"),colnames(tmp.atc.lkm))
apu.EE<-apply(tmp.atc.lkm[,apu],1,sum)

apu<-grepl(paste(tmp.atc.progestine,collapse="|"),colnames(tmp.atc.lkm))
apu.progestine<-apply(tmp.atc.lkm[,apu],1,sum)

apu<-cbind(apu.estradiaol,apu.EE,apu.progestine)
apu1<-apply(apu,1,max)
apu2<-apply(apu,1,which.max)
table(apu1,apu2)
apu3<-ifelse(apu1<2,0,apu2)# At least two purchases, indicator
table(apu3)
apu4<-factor(apu3,levels=0:3,labels = c("no","estradiaol","EE","progestine"))



anemia.dt$HC.group<-apu4



anemia.dt$age<-(anemia.dt$Time-anemia.dt$syntyv)


apu<-ifelse(anemia.dt$SK103+anemia.dt$SK215>0,1,0)
anemia.dt$DM<-apu
save(anemia.dt,anemia.flow,file=paste(out,"anemiaFinal.RData",sep=''))

# Korjataan
# Näistä riittää yksi vuodessa
# c("G03AC08","G03AC03","G02BA03")

load(file=paste(out,"erijatatc.RData",sep=''))
apu<-(tmp.atc.lkm[,c("G03AC08","G03AC03","G02BA03")]>0)*1
# table(kaksi=anemia.dt$G03AC08.i.2,apu[,"G03AC08"])
# table(kaksi=anemia.dt$G03AC03.i.2,apu[,"G03AC03"])
# table(kaksi=anemia.dt$G02BA03.i.2,apu[,"G02BA03"])

anemia.dt$G03AC08.i.2<-apu[,"G03AC08"]
anemia.dt$G03AC03.i.2<-apu[,"G03AC03"]
anemia.dt$G02BA03.i.2<-apu[,"G02BA03"]


# Lisätään uudet mjat
tmp.atc.COC_EE <- c("G03AA07","G03AA09","G03AA10","G03AA12","G03AA16","G03AB03","G03HB01")
tmp.atc.COC_E2 <-c("G03AA14","G03AB08")
tmp.atc.ring <-c("G02BB01")#
tmp.atc.patch <- c("G03AA13")
tmp.atc.POC <-c("G03AC01","G03AC09","G03AC10")
tmp.atc.implant <- c("G03AC08","G03AC03")# yksi
tmp.atc.IUD <- c("G02BA03")# yksi

# Yksi luokittelu
apu1<-(apply(tmp.atc.lkm[,tmp.atc.COC_EE],1,sum)>1)*1
apu2<-(apply(tmp.atc.lkm[,tmp.atc.COC_E2],1,sum)>1)*1
apu3<-(apply(tmp.atc.lkm[,tmp.atc.POC],1,sum)>1)*1
apu4<-(apply(tmp.atc.lkm[,tmp.atc.implant],1,sum)>0)*1
apu5<-(tmp.atc.lkm[,tmp.atc.ring]>1)*1
apu6<-(tmp.atc.lkm[,tmp.atc.patch]>1)*1
apu7<-(tmp.atc.lkm[,tmp.atc.IUD]>0)*1

apu.1.7<-cbind(apu4,apu7,apu1,apu2,apu3,apu5,apu6)
apu.1.7.max<-apply(apu.1.7,1,which.max)
apu.HC.3<-ifelse(apply(apu.1.7,1,sum)==0,0,apu.1.7.max)
apu.HC.3<-factor(apu.HC.3,levels = 0:7,
                 labels = c("no","implant","IUD","COC_EE","COC_E2","POC","ring","patch"))

table(apu.HC.3,apu.1.7.max)

apu<-(apply(tmp.atc.lkm[,tmp.atc.COC_EE],1,sum)>1)*1
anemia.dt$COC_EE<-apu

apu<-(apply(tmp.atc.lkm[,tmp.atc.COC_E2],1,sum)>1)*1
anemia.dt$COC_E2<-apu

apu<-(apply(tmp.atc.lkm[,tmp.atc.POC],1,sum)>1)*1
anemia.dt$POC<-apu

apu<-(apply(tmp.atc.lkm[,tmp.atc.implant],1,sum)>0)*1
anemia.dt$implant<-apu
anemia.dt$HC.group.7<-apu.HC.3


# Add number of deliveries
# Deliveries 1987-2022
load("W:/UpdateData/rawdeliv202311V2.RData")
apu.deliv<-subset(raw.deliv.2023.11.V2,shnro%in%anemia.dt$shnro)
apu.deliv$Time<-anemia.dt$Time[match(apu.deliv$shnro,anemia.dt$shnro)]
apu.deliv<-subset(apu.deliv,Time>synnytyksen_pvm.pvm)
apu.deliv<-apu.deliv[order(apu.deliv$shnro,-apu.deliv$synnytyksen_pvm.pvm),]
apu.deliv<-apu.deliv[!duplicated(apu.deliv$shnro),]
apu<-apu.deliv$AIEMMATSYNNYTYKSET[match(anemia.dt$shnro,apu.deliv$shnro)]+1
table(apu,useNA = "always") 
apu<-ifelse(is.na(apu),0,apu) 

anemia.dt$N.deliveries<-apu
anemia.dt$N.deliveries.1<-factor(ifelse(anemia.dt$N.deliveries==0,1,2),levels = 1:2,labels = c("0",">0"))
table(anemia.dt$N.deliveries,anemia.dt$N.deliveries.1,useNA = "always")

# save(anemia.dt,anemia.flow,file=paste(out,"anemiaFinalV20240527.RData",sep=''))# Korjatut arvo


# mk.ncc(data=anemia.dt,formula=Fail~G02BB01.i.2+strata(Set)+SK103+SK112+soseR+sivsr+ututkuR,
#        tabvar="soseR",out.df=TRUE)


#############################################################
######### Relevel individual HC variables ###################
#############################################################

anemia.dt_relevel <- anemia.dt 

# COC_EE (G03AA07, G03AA09, G03AA10, G03AA12, G03AA16, G03AB03, G03HB01)
anemia.dt_relevel$G03AA07 <- factor(
  ifelse(anemia.dt_relevel$G03AA07.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA07.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AA09 <- factor(
  ifelse(anemia.dt_relevel$G03AA09.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA09.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AA10 <- factor(
  ifelse(anemia.dt_relevel$G03AA10.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA10.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AA12 <- factor(
  ifelse(anemia.dt_relevel$G03AA12.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA12.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AA16 <- factor(
  ifelse(anemia.dt_relevel$G03AA16.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA16.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AB03 <- factor(
  ifelse(anemia.dt_relevel$G03AB03.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AB03.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03HB01 <- factor(
  ifelse(anemia.dt_relevel$G03HB01.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03HB01.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

# COC_E2 (G03AA14, G03AB08)

anemia.dt_relevel$G03AA14 <- factor(
  ifelse(anemia.dt_relevel$G03AA14.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA14.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AB08 <- factor(
  ifelse(anemia.dt_relevel$G03AB08.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AB08.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

# POC (G03AC01, G03AC09, G03AC10)

anemia.dt_relevel$G03AC01 <- factor(
  ifelse(anemia.dt_relevel$G03AC01.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AC01.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AC09 <- factor(
  ifelse(anemia.dt_relevel$G03AC09.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AC09.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AC10 <- factor(
  ifelse(anemia.dt_relevel$G03AC10.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AC10.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)


# IUDs (G02BA03)
anemia.dt_relevel$G02BA03 <- factor(
  ifelse(anemia.dt_relevel$G02BA03.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G02BA03.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)


# ring (G02BB01)

anemia.dt_relevel$G02BB01 <- factor(
  ifelse(anemia.dt_relevel$G02BB01.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G02BB01.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

# patch (G03AA13)

anemia.dt_relevel$G03AA13 <- factor(
  ifelse(anemia.dt_relevel$G03AA13.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AA13.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)


# implant (G03AC08, G03AC03)

anemia.dt_relevel$G03AC08 <- factor(
  ifelse(anemia.dt_relevel$G03AC08.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AC08.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

anemia.dt_relevel$G03AC03 <- factor(
  ifelse(anemia.dt_relevel$G03AC03.i.2 == 1, "yes",
         ifelse(anemia.dt_relevel$G03AC03.i.2 == 0 & anemia.dt_relevel$HC.group.7 == "no", "no", "other" ))
)

save(anemia.dt_relevel,file=paste(out,"anemiaFinalV20240527_relevel.RData",sep=''))# Korjatut arvo




#####################################################################
################# ADD ATC CODES AS NEW LEVELS #######################
#####################################################################


load(file=paste(out,"erijatatc.RData",sep=''))
apu<-(tmp.atc.lkm[,c("G03AC08","G03AC03","G02BA03")]>0)*1
# table(kaksi=anemia.dt$G03AC08.i.2,apu[,"G03AC08"])
# table(kaksi=anemia.dt$G03AC03.i.2,apu[,"G03AC03"])
# table(kaksi=anemia.dt$G02BA03.i.2,apu[,"G02BA03"])

# tmp.atc.COC_EE <- c("G03AA07","G03AA09","G03AA10","G03AA12","G03AA16","G03AB03","G03HB01")
# tmp.atc.COC_E2 <-c("G03AA14","G03AB08")
# tmp.atc.ring <-c("G02BB01")#
# tmp.atc.patch <- c("G03AA13")
# tmp.atc.POC <-c("G03AC01","G03AC09","G03AC10")
# tmp.atc.implant <- c("G03AC08","G03AC03")# yksi
# tmp.atc.IUD <- c("G02BA03")# yksi

tmp.atc.G03AA07 <- c("G03AA07")
tmp.atc.G03AA09 <- c("G03AA09")
tmp.atc.G03AA10 <- c("G03AA10")
tmp.atc.G03AA12 <- c("G03AA12")
tmp.atc.G03AA16 <- c("G03AA16")
tmp.atc.G03AB03 <- c("G03AB03")
tmp.atc.G03HB01 <- c("G03HB01")
tmp.atc.G03AA14 <- c("G03AA14")
tmp.atc.G03AB08 <- c("G03AB08")
tmp.atc.G03AC01 <- c("G03AC01")
tmp.atc.G03AC09 <- c("G03AC09")
tmp.atc.G03AC10 <- c("G03AC10")
tmp.atc.G02BB01 <- c("G02BB01")
tmp.atc.G03AA13 <- c("G03AA13")
tmp.atc.G02BA03 <- c("G02BA03")
tmp.atc.G03AC08 <- c("G03AC08")
tmp.atc.G03AC03 <- c("G03AC03")



# Yksi luokittelu
apu1<-(tmp.atc.lkm[,tmp.atc.G03AA07]>1)*1
apu2<-(tmp.atc.lkm[,tmp.atc.G03AA09]>1)*1
apu3<-(tmp.atc.lkm[,tmp.atc.G03AA10]>1)*1
apu4<-(tmp.atc.lkm[,tmp.atc.G03AA12]>1)*1
apu5<-(tmp.atc.lkm[,tmp.atc.G03AA16]>1)*1
apu6<-(tmp.atc.lkm[,tmp.atc.G03AB03]>1)*1
apu7<-(tmp.atc.lkm[,tmp.atc.G03HB01]>1)*1
# COC E2
apu8<-(tmp.atc.lkm[,tmp.atc.G03AA14]>1)*1
apu9<-(tmp.atc.lkm[,tmp.atc.G03AB08]>1)*1
# POC
apu10<-(tmp.atc.lkm[,tmp.atc.G03AC01]>1)*1
apu11<-(tmp.atc.lkm[,tmp.atc.G03AC09]>1)*1
apu12<-(tmp.atc.lkm[,tmp.atc.G03AC10]>1)*1
# ring
apu13<-(tmp.atc.lkm[,tmp.atc.G02BB01]>1)*1
# patch
apu14<-(tmp.atc.lkm[,tmp.atc.G03AA13]>1)*1
# IUD
apu15<-(tmp.atc.lkm[,tmp.atc.G02BA03]>0)*1 # yksi riittää
# implant
apu16<-(tmp.atc.lkm[,tmp.atc.G03AC08]>0)*1 # yksi riittää
apu17<-(tmp.atc.lkm[,tmp.atc.G03AC03]>0)*1 # yksi riittää


apu.1.17<-cbind(apu1,apu2,apu3,apu4,apu5,apu6,apu7,apu8,apu9,apu10,apu11,apu12,apu13,apu14,apu15,apu16,apu17)
apu.1.17.max<-apply(apu.1.17,1,which.max)
apu.HC.3<-ifelse(apply(apu.1.17,1,sum)==0,0,apu.1.17.max)
apu.HC.3<-factor(apu.HC.3,levels = 0:17,
                 labels = c("no", "G03AA07", "G03AA09", "G03AA10", "G03AA12", "G03AA16", "G03AB03", "G03HB01", "G03AA14", 
                            "G03AB08", "G03AC01", "G03AC09", "G03AC10", "G02BB01", "G03AA13", "G02BA03", "G03AC08", "G03AC03"))

table(apu.HC.3,apu.1.17.max)

# apu<-(apply(tmp.atc.lkm[,tmp.atc.COC_EE],1,sum)>1)*1
# anemia.dt$COC_EE<-apu
# 
# apu<-(apply(tmp.atc.lkm[,tmp.atc.COC_E2],1,sum)>1)*1
# anemia.dt$COC_E2<-apu
# 
# apu<-(apply(tmp.atc.lkm[,tmp.atc.POC],1,sum)>1)*1
# anemia.dt$POC<-apu
# 
# apu<-(apply(tmp.atc.lkm[,tmp.atc.implant],1,sum)>0)*1
# anemia.dt$implant<-apu

anemia.dt$HormonalContraception<-apu.HC.3

table(anemia.dt$HormonalContraception,anemia.dt$HC.group.7)

save(anemia.dt,anemia.flow,file=paste(out,"anemiaFinalV20240527.RData",sep=''))# Korjatut arvo



###########################################################################################
########### HOW MANY USE DUPLICATE TYPES OF HC DURING THE 1 YEAR TIME WINDOW? #############
###########################################################################################

# here we check number of users who have multiple HC

# create a list of HC column names that end with .i.2 (marker of having the required number of HC purchases)
multiple_HC <- grep("\\.i.2$", names(anemia.dt),value=TRUE)

# add a column that sums the number of of HC use per individual
anemia.dt$totalHCuse <- rowSums(anemia.dt[multiple_HC])

table(anemia.dt$totalHCuse,anemia.dt$HormonalContraception)

# count number of individuals with more than one type of HC in use
totalHCuse<-sum(anemia.dt$totalHCuse > 1) # 225 individuals with multiple types of HC

multHC <- anemia.dt %>%
  group_by(Set) %>% 
  filter(totalHCuse > 1) %>% 
  count(Set) # the 225 individuals come from 217 different NCC sets 

# make a list of sets to be removed (can't remove just the individual, as this will break the NCC)
filtered_anemia.dt <- anemia.dt[anemia.dt$totalHCuse > 1, ] # filter rows with more than 1 type of HC 
tmp.set <- unique(filtered_anemia.dt$Set) # get the unique Set numbers of the filtered rows
anemia_multHCrm.dt <- anemia.dt[!anemia.dt$Set %in% tmp.set, ] # remove 1198 women with the identified Set numbers, 16243 remain

table(anemia_multHCrm.dt$HormonalContraception,anemia_multHCrm.dt$HC.group.7)

percentageNew <- anemia_multHCrm.dt %>% 
  group_by(HormonalContraception) %>% 
  count %>% 
  summarise(n = (n)) %>% 
  mutate(PercentageHC = formattable::percent(n / sum(n))) %>% 
  arrange(desc(PercentageHC))

percentageOld <- anemia.dt %>% 
  group_by(HormonalContraception) %>% 
  count %>% 
  summarise(n = (n)) %>% 
  mutate(PercentageHC = formattable::percent(n / sum(n))) %>% 
  arrange(desc(PercentageHC))


                            # AFTER REMOVAL                  BEFORE REMOVAL
# HormonalContraception     n             %                 n             % 
# no                        11136         68.56%            11744         67.34%  
# G03AA12                   1279          7.87%             1446          8.29%
# G03AC09                   1217          7.49%             1323          7.59%
# G03AA10                   555           3.42%             621           3.56%
# G03HB01                   366           2.25%             406           2.33%
# G03AA09                   328           2.02%             372           2.13%
# G02BB01                   312           1.92%             341           1.96%
# G02BA03                   219           1.35%             231           1.32%
# G03AB08                   198           1.22%             223           1.28%
# G03AA14                   178           1.10%             222           1.27%
# G03AA16                   104           0.64%             118           0.68%
# G03AC01                   101           0.62%             112           0.64%
# G03AC03                   74            0.46%             80            0.46%
# G03AA07                   67            0.41%             84            0.48%
# G03AC08                   51            0.31%             55            0.32%
# G03AA13                   41            0.25%             41            0.24%
# G03AC10                   14            0.09%             18            0.10%
# G03AB03                   3             0.02%             4             0.02%
 
# this solved most of the issues. but there is still some overlap, particularly in "no" for HormonalContraception there are some EE, E2 and POC users. This is due to us making the HC.group.7 in an early stage of the code (before removing individuals with multiple HC types), whereby someone using several different HC types was put in a group depending randomly based on the different HC types. Let's instead make a new variable based on the current filtering.

anemia_multHCrm.dt <- anemia_multHCrm.dt %>% 
  mutate(Classification = case_when(
    HormonalContraception == "no" ~ "no",
    HormonalContraception == "G03AA07" ~ "COC_EE",
    HormonalContraception == "G03AA09" ~ "COC_EE",
    HormonalContraception == "G03AA10" ~ "COC_EE",
    HormonalContraception == "G03AA12" ~ "COC_EE",
    HormonalContraception == "G03AA16" ~ "COC_EE",
    HormonalContraception == "G03AB03" ~ "COC_EE",
    HormonalContraception == "G03HB01" ~ "COC_EE",
    HormonalContraception == "G03AA14" ~ "COC_E2",
    HormonalContraception == "G03AB08" ~ "COC_E2",
    HormonalContraception == "G03AC01" ~ "POC",
    HormonalContraception == "G03AC09" ~ "POC",
    HormonalContraception == "G03AC10" ~ "POC",
    HormonalContraception == "G02BB01" ~ "ring",
    HormonalContraception == "G03AA13" ~ "patch",
    HormonalContraception == "G02BA03" ~ "IUD",
    HormonalContraception == "G03AC08" ~ "implant",
    HormonalContraception == "G03AC03" ~ "implant"
     )) 

anemia_multHCrm.dt <- anemia_multHCrm.dt %>% 
  mutate(Classification = factor(Classification, levels = c("no","implant","IUD","COC_EE","COC_E2","POC","ring","patch")
   ))
  
    
table(anemia_multHCrm.dt$HormonalContraception,anemia_multHCrm.dt$Classification)

anemia.dt <- anemia_multHCrm.dt
save(anemia.dt,anemia.flow,file=paste(out,"anemiaFinalV20240527.RData",sep=''))# Tämä versio käytössä 21082024

