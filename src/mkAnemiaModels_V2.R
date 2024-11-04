#------------------------------------------------------------------------------
# Data for HC Anemia Study (HCA)
# 
# Created: 21-05-2024, original code by JKH
# Updated:  29-05-2024 (SE), 03-06-2024 (SE), 11-06-2024 (SE), 12-06-2024 (SE), 13-08-2024 (SE, manuscript V1)
# Authors: JKH & SE
#
#------------------------------------------------------------------------------


#------------------------
# 
# Load data
# 
#------------------------
# load(file="anemiaFinal.RData")
#load(file="anemiaFinalV20240527.RData")
load(file="W:/Veripalvelu/mikko/anemia_and_hc/out20240821/anemiaFinalV20240527.RData")


# Some functions
# source('W:/multipleNCC/includefunctionsNCC.R')

mk.ncc<-function(data,formula,tabvar,out.df=FALSE,level=.95,small.N=5){
  require(Epi)
  require(survival)
  lv.failvar<-as.character(formula)[2]
  lv.1<-clogit(formula=formula,data=data,iter.max = 10000) # original 200
  lv.2<-exp(cbind(coefficients(lv.1),confint(lv.1,level = level)))
  if(!missing(tabvar)){
    lv.2<-lv.2[grepl(tabvar,rownames(lv.2),fixed = TRUE),]
    eval(parse(text=paste0("lv.3<-with(",deparse(substitute(data)),",table(",tabvar,",",lv.failvar,"))")))
    lv.3<-ifelse(lv.3<small.N,999999,lv.3)
    lv.ulos<-cbind(lv.3,rbind(1,lv.2))
    if(out.df){
      lv.ulos<-data.frame(Var=rownames(lv.ulos),lv.ulos,stringsAsFactors = FALSE)
      names(lv.ulos)<-c(tabvar,"N.control","N.case","OR","OR.lo","OR.hi")
      return(lv.ulos)
    }
    lv.ulos<-cbind(rownames(lv.ulos),lv.ulos)
    return(lv.ulos)
  }
  lv.2
}


#------------------------
# NCC modeling
#------------------------
library(survival)
library(codeCollection)
library(Epi)

# dput(names(anemia.dt))
tmp.expovars<-c(
                #"HC.group.7", 
                "HormonalContraception", 
                "Classification"
                # "G02BA03.i.2", 
                # "G02BB01.i.2", 
                # "G03AA07.i.2", 
                # "G03AA09.i.2",
                # "G03AA10.i.2", 
                # "G03AA12.i.2", 
                # "G03AA13.i.2", 
                # "G03AA14.i.2", 
                # "G03AA16.i.2",
                # "G03AB03.i.2", 
                # "G03AB08.i.2", 
                # "G03AC01.i.2", 
                # "G03AC03.i.2", 
                # "G03AC08.i.2", 
                # "G03AC09.i.2", 
                # "G03AC10.i.2", 
                # "G03HB01.i.2"
            )#,
            

# Only exposure variables

# mk.ncc(data=anemia.dt,formula=Fail~G02BA03.i.2+strata(Set), tabvar='G02BA03.i.2',out.df=TRUE,small.N = 5)

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt,formula=Fail~",x,
                 "+strata(Set), tabvar='",x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m0<-do.call(rbind,apu1)

# SOSECO and history of pregnancy added

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt,formula=Fail~",x,
                 "+strata(Set)+sivsr+soseR+ututkuR+N.deliveries.1, tabvar='",x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m1<-do.call(rbind,apu1)

# Chronic diseases (special reimbursement rights) added
# Kela reimbursment codes include the following diagnoses:
# *  SK103: Diabetes mellitus; E10-E14, E89.1 
# *  SK104: Hypothyreosis; C73, E03, E89.0  
# *  SK115: Breast cancer; C50, D05.1   
# *  SK117: Leukemia, lymphoma and other malignancy related to blood or bone marrow; C81-C85, C88, C90-C96, D45-D47, D72.1, D75 
# *  SK128: Gynecologic cancer; D39, C51-C58 
# *  SK130: Other malignancies; C00-C26, C30-C34, C37-C41, C43-C49, C64-C80, C97 (list excluded malignancies of male reproductive organs)
# *  SK202: Connective tissue disease; A04.6, A39.8, A50.5, D76, H20, H30, I33.0, I40.8, J84, K50.9, K51.9, K73.2, K74.3, K75.4, K83.0, L40.5, M02, M05, M06, M08, M13, M30-M35, M45, M46.1, M46.9, M86.6, M94.1, N03, N04, Q44.2  
# *  SK208: IBD; K50, K51  
# *  SK215: Diabetes mellitus, other treatment; E10-E14, E89.1
# Available ICD codes:
#   * I80 Phlebitis and thrombophlebitis
# * I82 Other venous embolism and thrombosis
# * I26 Pulmonary embolism
# * I21 Myocardial infaction
# * I63 Cerebral infarction


apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt,formula=Fail~",x,
                 "+strata(Set)+sivsr+soseR+ututkuR+N.deliveries.1+DM+SK104+SK115+SK117+SK128+SK130+SK202+SK208+I82+I26+I21+I63+E66, tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m2<-do.call(rbind,apu1)


# HMB added

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt,formula=Fail~",x,
                 "+strata(Set)+sivsr+soseR+ututkuR+N.deliveries.1+DM+SK104+SK115+SK117+SK128+SK130+SK202+SK208+I80+I81+I82+I26+I21+I63+E66+D25+N80+N84+N85+N92+N93, tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m3<-do.call(rbind,apu1)

  


sink(file="AnemiaNCCResults.html")
cat("<br><pre>Updated:",format(Sys.time()),"</pre><hr><br>")
knitr::kable(anemia.m0,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             No adjustments. ")
cat("<hr><br>")
knitr::kable(anemia.m1,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Adjusted by sivsr+soseR+ututkuR+N.deliveries.1.")
cat("<hr><br>")
knitr::kable(anemia.m2,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Adjusted by sivsr+soseR+ututkuR+N.deliveries.1+DM+SK104+SK115+SK117+SK128+SK130+SK202+SK208+I82+I26+I21+I63+E66.")
cat("<hr><br>")
knitr::kable(anemia.m3,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Adjusted by sivsr+soseR+ututkuR+N.deliveries.1+DM+SK104+SK115+SK117+SK128+SK130+SK202+SK208+I80+I81+I82+I26+I21+I63+E66+D25+N80+N84+N85+N92+N93.")
sink()


save(anemia.m0, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/models/anemia.m0.Rdata")
save(anemia.m1, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/models/anemia.m1.Rdata")
save(anemia.m2, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/models/anemia.m2.Rdata")
save(anemia.m3, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/models/anemia.m3.Rdata")

####################################################################################################################################
###################################### ADDITIONAL MODELS FOR TESTING AND SENSITIVITY ANALYSIS ######################################
####################################################################################################################################


####### MODEL m4: number of pregnancy instead of history of pregnancy ###########


apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt,formula=Fail~",x,
                 "+strata(Set)+sivsr+soseR+ututkuR+N.deliveries, tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m4<-do.call(rbind,apu1)


sink(file="AnemiaNCCResults.html")
cat("<hr><br>")
knitr::kable(anemia.m4,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Adjusted by sivsr+soseR+ututkuR+Number of previous pregnancies.")
sink()


save(anemia.m4, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/sensitivity_analysis/models_sensitivity/anemia.m4.Rdata")



####### MODEL m5: exclude age <25 years (possible free HC) ###########

anemia.dt_25 <- anemia.dt %>% 
  filter(age>25) # we are matching based on age (one year bins), so we don't have to worry about the Set as everyone in the set is the same age

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt_25,formula=Fail~",x,
                 "+strata(Set), tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m5<-do.call(rbind,apu1)


sink(file="AnemiaNCCResults.html")
cat("<hr><br>")
knitr::kable(anemia.m5,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Non-adjusted, women <25 years old excluded")
sink()


save(anemia.m5, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/sensitivity_analysis/models_sensitivity/anemia.m5.Rdata")




####### MODEL m6: exclude pensioners (possible underlying pathology causing anemia) ###########

# count number of individuals who are pensioners
pensioner<-sum(anemia.dt$soseR == "Pensioners") # 360 individuals 

multHC <- anemia.dt %>%
  group_by(Set) %>% 
  filter(soseR == "Pensioners") %>% 
  count(Set) # the 360 individuals come from 342 different NCC sets 

# make a list of sets to be removed (can't remove just the individual, as this will break the NCC)
filtered_anemia.dt <- anemia.dt[anemia.dt$soseR == "Pensioners", ] # filter pensioners
tmp.set <- unique(filtered_anemia.dt$Set) # get the unique Set numbers of the filtered rows
anemia.dt_pen <- anemia.dt[!anemia.dt$Set %in% tmp.set, ] # remove 4178 women with the identified Set numbers, 14356 remain

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt_pen,formula=Fail~",x,
                 "+strata(Set)+sivsr+soseR+ututkuR+N.deliveries.1+DM+SK104+SK115+SK117+SK128+SK130+SK202+SK208+I80+I81+I82+I26+I21+I63+E66+D25+N80+N84+N85+N92+N93, tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m6<-do.call(rbind,apu1)


sink(file="AnemiaNCCResults.html")
cat("<hr><br>")
knitr::kable(anemia.m6,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Pensioners excluded, fully adjusted (compare to m3)")
sink()


save(anemia.m6, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/sensitivity_analysis/models_sensitivity/anemia.m6.Rdata")



####### Model m7: remove <35 year-olds (examine if difference in POC and COCs is caused by age difference, or the product) #########


anemia.dt_35 <- anemia.dt %>% 
  filter(age>=35) # we are matching based on age (one year bins), so we don't have to worry about the Set as everyone in the set is the same age

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt_35,formula=Fail~",x,
                 "+strata(Set), tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m7<-do.call(rbind,apu1)


sink(file="AnemiaNCCResults.html")
cat("<hr><br>")
knitr::kable(anemia.m7,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Non-adjusted, women <35 years old excluded")
sink()


save(anemia.m7, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/sensitivity_analysis/models_sensitivity/anemia.m7.Rdata")



####### Model m8: only Helsinki #########


anemia.dt_hki <- anemia.dt %>% 
  filter(seutukunta=="011") # we are matching based on municipality, so we don't have to worry about the Set

apu1<-lapply(tmp.expovars,function(x){
  lv.txt<-paste0("lv.m<-mk.ncc(data=anemia.dt_hki,formula=Fail~",x,
                 "+strata(Set), tabvar='",
                 x,"',out.df=TRUE)")
  eval(parse(text=lv.txt))
  lv.m<-cbind(mkcodelabel(gsub(x=names(lv.m)[1],pattern = ".i.2",replacement = "",fixed = TRUE),
                          ATCKoodit,as.factor.B = FALSE),lv.m)
  names(lv.m)[1:2]<-c("Var","Level")
  lv.m
})

anemia.m8<-do.call(rbind,apu1)


sink(file="AnemiaNCCResults.html")
cat("<hr><br>")
knitr::kable(anemia.m8,format="html",row.names = FALSE,digits=c(3),
             caption="NCC odds ratios based on conditional logistic regression. 
             Non-adjusted, only Helsinki")
sink()


save(anemia.m8, file = "W:/Veripalvelu/sofie/anemia_and_hc/results/sensitivity_analysis/models_sensitivity/anemia.m8.Rdata")




#------------------------
# Models taking into account HMB
# Updated: 05-07-2024
# Author: Jari Haukka
#------------------------

# Exposure var: HC.group.7

# Add HMB 
# Load updated HILMO and Avo-HILMO data
load("W:/UpdateData/THL0224hilmo.RData")
load("W:/UpdateData/THL0224avohilmo.RData")

thl2023_5525_hilmo$TUPVA.pvm<-Epi::cal.yr(substr(thl2023_5525_hilmo$TUPVA,1,10),format="%d.%m.%Y")
thl2023_5525_avohilmo$KAYNTI_ALKOI.pvm<-Epi::cal.yr(substr(thl2023_5525_avohilmo$KAYNTI_ALKOI,1,10),format="%d.%m.%Y")

head(thl2023_5525_hilmo_icd10,2)

apu.hilmo<-subset(thl2023_5525_hilmo,shnro%in%anemia.dt$shnro)
apu.hilmo.icd<-subset(thl2023_5525_hilmo_icd10,(HILMO_ID%in%apu.hilmo$HILMO_ID)&grepl("D25|N80|N84|N85|N92|N93",KOODI))

apu.avohilmo<-subset(thl2023_5525_avohilmo,shnro%in%anemia.dt$shnro)
apu.avohilmo.icd<-subset(thl2023_5525_ah_icd10,(AVOHILMO_ID%in%apu.avohilmo$AVOHILMO_ID)&grepl("D25|N80|N84|N85|N92|N93",ICD10))

apu.i.1<-(apu.hilmo$HILMO_ID%in%apu.hilmo.icd$HILMO_ID)
apu.i.2<-(apu.avohilmo$AVOHILMO_ID%in%apu.avohilmo.icd$AVOHILMO_ID)

apu.HMB<-data.frame(shnro=c(apu.hilmo$shnro[apu.i.1],apu.avohilmo$shnro[apu.i.2]),
                    pvm=c(apu.hilmo$TUPVA.pvm[apu.i.1],apu.avohilmo$KAYNTI_ALKOI.pvm[apu.i.2])
                    )
apu.HMB$Time<-anemia.dt$Time[match(apu.HMB$shnro,anemia.dt$shnro)]


anemia.dt$any.HMB<-1*(anemia.dt$shnro%in%apu.HMB$shnro)
anemia.dt$before.HMB<-1*(anemia.dt$shnro%in%apu.HMB$shnro[apu.HMB$pvm<apu.HMB$Time])

anemia.dt.any.HMB<-subset(anemia.dt,any.HMB==1)
anemia.dt.before.HMB<-subset(anemia.dt,before.HMB==1)

  
# tmp.atc.COC_EE <- c("G03AA07","G03AA09","G03AA10","G03AA12","G03AA16","G03AB03","G03HB01")
# tmp.atc.COC_E2 <-c("G03AA14","G03AB08")
# tmp.atc.ring <-c("G02BB01")#
# tmp.atc.patch <- c("G03AA13")
# tmp.atc.POC <-c("G03AC01","G03AC09","G03AC10")
# tmp.atc.implant <- c("G03AC08","G03AC03")# yksi
# tmp.atc.IUD <- c("G02BA03")# yksi

apu.i.<-(anemia.dt.any.HMB$HC.group.7%in%c("IUD","COC_EE"))
anemia.dt.any.HMB.1<-subset(anemia.dt.any.HMB,apu.i.)
anemia.dt.any.HMB.1$HC.group.7<-droplevels(anemia.dt.any.HMB.1$HC.group.7)


apu.i.<-(anemia.dt.before.HMB$HC.group.7%in%c("IUD","COC_EE"))
anemia.dt.before.HMB.1<-subset(anemia.dt.before.HMB,apu.i.)
anemia.dt.before.HMB.1$HC.group.7<-droplevels(anemia.dt.before.HMB.1$HC.group.7)


mk.ncc(data=anemia.dt.any.HMB.1,formula=Fail~HC.group.7+strata(Set),tabvar='HC.group.7',out.df=TRUE)
mk.ncc(data=anemia.dt.before.HMB.1,formula=Fail~HC.group.7+strata(Set),tabvar='HC.group.7',out.df=TRUE)


