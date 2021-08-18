setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
library(asreml)
library(reshape)
dfall<- read.csv('predicted values table.csv', row.names=1)
bluesYld<- dfall[which(dfall$studyName %in% c('Adv_Car_18','Adv_Neo_18','Adv_StJ_18','Adv_Urb_18',
                                              'Adv_Car_20','Adv_Neo_20','Adv_StJ_20','Adv_Urb_20', 'Adv_Scb_20',"AdvHY_Urb_20",
                                              'Adv_Car_19','Adv_Neo_19','Adv_Stj_19','Adv_Urb_19',"AdvHY_Urb_19",
                                              'Pr[0-9]_Car_19','Pr[0-9]_Neo_19','Pr[0-9]_Stj_19','Pr[0-9]_Scb_19','Pr[0-9]_Urb_19',
                                              'Pr[0-9]_Car_18','Pr[0-9]_Neo_18','Pr[0-9]_Stj_18','Pr[0-9]_Scb_18','Pr[0-9]_Urb_18',
                                              'Adv_Scb_19', 'Adv_Scb_20', "Adv_Sbmv_20")),]
bluesYld$name[grep('Pio',bluesYld$name)]
bluesYld$name<- gsub("Pio ", "Pio", bluesYld$name)
##Functions to use later
#function to check model convergence and update until converged (tolerate a 1.5% change in components)
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#find out what lines overlap across trials
yr<- matrix(unlist(strsplit(as.character(bluesYld$studyName), split="_")), nrow=3)[3,]
bluesYld<- data.frame(bluesYld, yr)
tab<- cast(bluesYld, name~yr, value='predicted.value')
cks<- tab[intersect(intersect(which(tab[,2]!=0), which(tab[,3]!=0)), which(tab[,4]!=0)),'name']
head(bluesYld)
isck<- rep(0, nrow(bluesYld))
isck[which(bluesYld$name %in% cks)]<- 1
bluesYld<- data.frame(bluesYld, isck)
bluesYld$isck<- as.factor(as.character(bluesYld$isck))


#add rows for the missing values
a<- cast(bluesYld, name+studyName+isck~trait, value='predicted.value')
b<- cast(bluesYld, name+studyName+isck~trait, value='standard.error')
#a[,-c(1:3)]<- scale(a[,-c(1:3)], scale=FALSE)
mta<- melt(a)
mtb<- melt(b)


#estimate weights based on the standard error 
wt<- 1/(mtb$value^2)

#set missing weights to zero
wt[which(is.na(wt))]<- 0

#add weight to the data-set
df<- data.frame(mta, wt)
df<- droplevels.data.frame(df)

####fit model across locations for agronomic traits
Loc<- matrix(unlist(strsplit(as.character(df$studyName), split="_")), nrow=3)[2,]
Loc<- gsub('StJ', "Stj", Loc)
Loc2<- Loc
Loc2<- gsub("Urb", "North", Loc2)
Loc2<- gsub("Neo", "North", Loc2)
Loc2<- gsub("Stj", "South", Loc2)
Loc2<- gsub("Car", "South", Loc2)
df<- data.frame(df, Loc, Loc2)

#subset the traits
trts<- c('Grain.test.weight...lbs.bu', 'Grain.yield...bu.ac','Heading.time...Julian.date..JD..','Plant.height.inches',
         'FHB.DON.content...ppm.','FHB.incidence.....','FHB.severity.....','FHB.grain.incidence.....', 'SBMV')
df<- df[which(df$trait %in% trts),]

#consider yld in Urb/Neo, Stj, and Car as different traits
ixyld<- intersect(which(df$trait=='Grain.yield...bu.ac'), which(!is.na(df$value)))
df$trait<- as.character(df$trait)
df$trait[ixyld]<- paste(df$trait[ixyld], df$Loc2[ixyld])
df$trait<- as.factor(df$trait)

#fit the model
df$name<- as.factor(df$name)
mod<- asreml(fixed=value~1+trait, random=~at(isck, "1"):studyName+us(trait):name, 
             data=df, weights=wt, family = asr_gaussian(dispersion = 1), maxiter=500)
mod<- mkConv(mod)
p<- predict(mod, classify = 'name:trait', pworkspace=64e7)
allblup<- na.omit(p$pvals)

#summarize in wide format
wide<- cast(allblup, name~trait, value='predicted.value')
wide0<- wide
for(i in 2:ncol(wide)){
  wide[,i]<- scale(wide[, i], scale=FALSE)
  if(NA %in% wide[,i]){
    wide[which(is.na(wide[,i])),i]<- mean(wide[,i], na.rm=T)
  }
}

cr<- cor(wide[,-1])
colnames(cr)<- colnames(wide)[-1]
rownames(cr)<- colnames(wide)[-1]

##calculate indies
balancedIX<-wide$Grain.test.weight...lbs.bu * .5 +
  wide$`Grain.yield...bu.ac North` * 3.5 +
  wide$`Grain.yield...bu.ac South` * 7 +
  wide$Heading.time...Julian.date..JD.. * -7 +
  wide$FHB.grain.incidence.....* -0.1  +
  wide$FHB.incidence.....* -.1  +
  wide$FHB.severity.....* -.5  +
  wide$FHB.DON.content...ppm. * -7+ wide$Plant.height.inches * -4 + wide$SBMV * -9

#add index to the table
wide2<- cbind(wide0, balancedIX)
wide3<- data.frame(wide2, ranking=rank(-balancedIX))
wide3<- wide3[order(wide3$ranking),]

expectedGain<- colMeans(wide3[which(wide3$ranking<=30),-1])- colMeans(wide3[which(wide3$ranking>30),-1])
expectedGain

#indicate number of years of testing
sets<- c()
for(i in 1:nrow(wide3)){
  a<- strsplit(as.character(unique(bluesYld[which(wide3[i,1]==bluesYld$name),'studyName'])), split="_")
  b<- t(matrix(unlist(a), nrow=3))[,c(1,3)]
  b<- unique(b)
  b<- as.matrix(b)
  b[,1]<- gsub("HY", "", b[,1])
  b<- unique(b)
  if(ncol(b)>1){
    b<- t(b)
  }
  set<- unique(paste(b[1,], b[2,], sep=""))
  set<- paste(set, collapse=",")
  sets<- append(sets, set)
}
wide3<- data.frame(wide3, sets, multiyr='false', stringsAsFactors = FALSE)
wide3[grep(",", sets),'multiyr']<- 'true'


#re-arrange columns
cols<-c("name","Grain.yield...bu.ac.North" , "Grain.yield...bu.ac.South","Grain.test.weight...lbs.bu", "Heading.time...Julian.date..JD..","Plant.height.inches","FHB.DON.content...ppm.",
        "FHB.grain.incidence.....","FHB.incidence.....","FHB.severity.....","SBMV", 
        "Grain.yield...bu.ac.North" , "Grain.yield...bu.ac.South", "balancedIX","ranking","multiyr", "sets")         
wide4<- wide3[,c(cols)]


#add column for licensing status
lic<- c('12-26004', '13-1910', '13-1960','14-11830',
        '14-28468', '15-17909', '15-18976', '15-23771','15-2639', '15-27270', '15-4957',
        '16-11226', '16-1922', '16-28753', '16-8048', '16-8737')
      
licenced<- wide3[,1] %in% lic
wide4<- data.frame(wide4, licenced)

#plot correlation between yields and index
library(GGally)
ggpairs(wide4[,c(12:14)])



write.csv(wide4, file='all 2018-2020 rslts summary.csv')

