setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
library(reshape)
library(asreml)
library(asremlPlus)
dfall<- read.csv('predicted values table.csv', row.names=1)
bluesYld<- dfall[which(dfall$studyName %in% c("VE_Urb_20", "VE_Scb_20","VE_Urb_19", "VE_Scb_19")),]

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
mod<- asreml(fixed=value~1+trait+studyName+studyName:trait, random=~ us(trait):name, 
             data=df, weights=wt, family = asr_gaussian(dispersion = 1),maxiter=500)
mod<- mkConv(mod)
pscb<- predict(mod, classify = 'name:trait', present=list(c("studyName"),prwts=c(1, 0)))
scbblup<- pscb$pvals
scbblup[which(scbblup$trait=='Heading.time...Julian.date..JD..'),'predicted.value']<- NA
pyld<- predict(mod, classify = 'name:trait', present=list(c("studyName"),prwts=c(0, 1)))
yldblup<- pyld$pvals
allblup<- na.omit(rbind(scbblup, yldblup))

#summarize in wide format
wide<- cast(allblup, name~trait, value='predicted.value', fun.aggregate = 'mean')
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
  wide$Grain.yield...bu.ac * 8 +
  wide$Heading.time...Julian.date..JD.. * -6 +
  wide$FHB.grain.incidence.....* -0.3  +
  wide$FHB.incidence.....* -1  +
  wide$FHB.severity.....* -1  +
  wide$Plant.height.inches * -10

#add index to the table
wide2<- cbind(wide, balancedIX)

#matrix
wide3<- data.frame(wide2, ranking=rank(-wide2$balancedIX))

#selected
vec<- c('17-29544', 'US16-IL-064-005', 'US16-IL-063-221', 
  'US16-IL-063-123', '17-8930', '17-23868', '17-8832')

unique(dfall[which(dfall[,1] %in% vec),'studyName'])

write.csv(wide3[match(c('17-29544', 'US16-IL-064-005', 'US16-IL-063-221', 
                                   'US16-IL-063-123', '17-8930', '17-23868', '17-8832'), wide3[,1]),], file='best from VE test 2020.csv')
