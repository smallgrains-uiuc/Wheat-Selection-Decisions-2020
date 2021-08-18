setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
library(asreml)
library(reshape)
dfall<- read.csv('predicted values table.csv', row.names=1)
bluesYld<- dfall[which(dfall$studyName %in% c('Adv_Car_20','Adv_Neo_20','Adv_StJ_20','Adv_Urb_20', 'Adv_Scb_20',"AdvHY_Urb_20",
                                              'Adv_Car_19','Adv_Neo_19','Adv_Stj_19','Adv_Urb_19',"AdvHY_Urb_19",
                                              'Pr[0-9]_Car_19','Pr[0-9]_Neo_19','Pr[0-9]_Stj_19','Pr[0-9]_Scb_19','Pr[0-9]_Urb_19',
                                              'Adv_Scb_19', 'Adv_Scb_20', "Adv_Sbmv_20")),]

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


#add rows for the missing values
a<- cast(bluesYld, name+studyName~trait, value='predicted.value')
b<- cast(bluesYld, name+studyName~trait, value='standard.error')
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

#subset the yield trials and agronomic traits
trts<- c('Grain.test.weight...lbs.bu', 'Grain.yield...bu.ac','Heading.time...Julian.date..JD..','Plant.height.inches',
         'FHB.DON.content...ppm.','FHB.incidence.....','FHB.severity.....','FHB.grain.incidence.....', 'SBMV')
df<- df[which(df$trait %in% trts),]

#fit the model
mod<- asreml(fixed=value~1+trait, random=~studyName+studyName:trait+us(trait):name, 
             data=df, weights=wt, family = asr_gaussian(dispersion = 1),maxiter=500)
mod<- mkConv(mod)
p<- predict(mod, classify = 'name:trait')
allblup<- p$pvals

#summarize in wide format
wide<- cast(allblup, name~trait, value='predicted.value')
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
balancedIX<-wide$Grain.test.weight...lbs.bu * 0.5 +
  wide$Grain.yield...bu.ac * 5 +
  wide$Heading.time...Julian.date..JD.. * -6 +
  wide$FHB.grain.incidence.....* -0.1  +
  wide$FHB.incidence.....* -.5  +
  wide$FHB.severity.....* -1.5  +
  wide$FHB.DON.content...ppm. * -1
  wide$Plant.height.inches * -12+
  wide$SBMV * 2

#add index to the table
wide2<- cbind(wide, balancedIX)

#subset lines that have been in advanced for only one year
ixadv19<- which(wide$name %in% unique(df[which(df$studyName == 'Adv_Urb_19'),'name']))
ixadv20<- which(wide$name %in% unique(df[which(df$studyName == 'Adv_Urb_20'),'name']))
lines<- wide2[ixadv20[-which(ixadv20 %in% ixadv19)],1]
lines<- setdiff(lines, 'Pio25R74')
wide2<- wide2[match(lines, wide2[,1]),]

#matrix
wide3<- data.frame(wide2, ranking=rank(-wide2$balancedIX))

#estimate expected gain
expectedGain<- colMeans(wide3[which(wide3$ranking<=40),-1])- colMeans(wide3[which(wide3$ranking>40),-1])
as.matrix(round(expectedGain,2))

#get the selected
selected<- wide3[which(wide3$ranking<=40),]
selected<- selected[1:40,]
write.csv(selected, file="selected from advanced 40.csv")

#get the next best for the Nifa North (NN)
nonselected<- wide3[which(wide3$ranking<=70),]
nonselected<- nonselected[order(nonselected$ranking),]
nonselected<- nonselected[41:70,1]
library(gtools)
nonselected<- nonselected[mixedorder(gsub("-", "A",nonselected))]
write.csv(nonselected, file="selected from advanced 30 for NN trial.csv")

#create entry list for advanced 2021
selected1<- read.csv("selected from prelim 130.csv", row.names=1)
cks<- c('07-4415', '02-18228', '07-19334', 'Pio 25R74', 'Kaskaskia')
advlist<- c(cks, as.character(selected[,1]), as.character(selected1[,1]))

library(gtools)
advlist<- c(advlist[c(1:5)], advlist[-c(1:5)][mixedorder(gsub("-", "A", advlist[-c(1:5)]))])
source<- c(rep("Checks",5), rep('Selected from Advanced 2019-2020', 170))
source[which(advlist %in% selected1[,1])]<- "Selected from Prelim 2019-2020"
advlist<- data.frame(advlist, source)
write.csv(advlist, file='advanced list 2020.csv')

getwd()


