setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
dfall<- read.csv('predicted values table.csv', row.names=1)
bluesYld<- dfall[which(dfall$studyName %in% c('Pr_Car_20','Pr_Neo_20','Pr_Stj_20','Pr_Urb_20', 'Pr_Scb_20','Pr_Sbmv_20')),]

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

#add weight to the dataset
df<- data.frame(mta, wt)
df<- droplevels.data.frame(df)

#fit model across locations
mod<- asreml(fixed=value~1+studyName+trait+studyName:trait, random=~us(trait):name, 
             data=df, weights=wt, family = asr_gaussian(dispersion = 1))
mod<- mkConv(mod)


#get the blups + means across locations for yield and test weight
p<- predict(mod, classify = 'name:trait', levels= list(trait=c('Grain.yield...bu.ac','Grain.test.weight...lbs.bu')), 
            pworkspace=64e6,
            present=list(c("studyName"),prwts=c(.25,.25,0,0,.25,.25)))
blup<- p$pvals

#get the blups + means across locations for height
p<- predict(mod, classify = 'name:trait', pworkspace=64e6, 
            present=list(c("studyName"),
            prwts=c(0,0.5,0,0,0, 0)))
blupht<- p$pvals
blupht<- blupht[which(blupht$trait=='Plant.height.inches'),]

#get the blups + means across locations for julain
p<- predict(mod, classify = 'name:trait', pworkspace=64e6, 
            present=list(c("studyName"),
                         prwts=c(0,0,0,0,0,1)))
blupjul<- p$pvals
blupjul<- blupjul[which(blupjul$trait=="Heading.time...Julian.date..JD.."),]

#get the blups + means across locations for scab
p<- predict(mod, classify = 'name:trait', pworkspace=64e6, 
            present=list(c("studyName"),
                         prwts=c(0,0,0,1,0,0)))
blupscb<- na.omit(p$pvals)
blupscb<- blupscb[-which(blupscb[,2]=='Heading.time...Julian.date..JD..'),]

#get the blups + means across locations for sbmv
p<- predict(mod, classify = 'name:trait', pworkspace=64e6, 
            present=list(c("studyName"),
                         prwts=c(0,0,1,0,0,0)))
blupsbmv<- na.omit(p$pvals)

#combine results
allblup<- rbind(blup, blupht, blupjul, blupscb, blupsbmv)
wide<- cast(allblup, name~trait, value='predicted.value')
for(i in 2:ncol(wide)){
  wide[,i]<- scale(wide[, i], scale=FALSE)
  if(NA %in% wide[,i]){
    wide[which(is.na(wide[,i])),i]<- mean(wide[,i], na.rm=T)
  }
}

##calculate indies
balancedIX<- wide$Grain.test.weight...lbs.bu * 0.5 +
  wide$Grain.yield...bu.ac * 2.5 +
  wide$Heading.time...Julian.date..JD.. * -3.7 +
  wide$FHB.grain.incidence.....* -0.06  +
  wide$FHB.incidence.....* -0.03  +
  wide$FHB.severity.....* -0.02  +
  wide$Plant.height.inches * -3.2+
  wide$SBMV * .1

#add index to the table
wide2<- cbind(wide, balancedIX)

#check correlation matrix
as.matrix(cor(wide2[,-1])[,'balancedIX'])

#matrix
wide3<- data.frame(wide2, ranking=rank(-wide2$balancedIX))
rmv<- c('18-14132', '18-14647','18-15357','18-1050','18-1435','18-370','18-2872',
         '18-16216','18-16368','18-4772','18-1108','18-4022','18-8149',
         '18-3314','18-8907','18-10148','18-14308','18-12636','18-15271','18-26107','18-870',
        'Pio25R74', 'Pio 25R74', '07-4415','02-18228')
wide3<- wide3[-which(wide3[,1] %in% rmv),]
selected<- wide3[which(wide3$ranking<=139),]

#calculate expected gain
expectedGain<- colMeans(selected[,-1])- colMeans(wide3[,-1])
as.matrix(round(expectedGain,2))
nrow(selected)
write.csv(selected, file="selected from prelim 130.csv")

#select from the non-selected for the Nifa North
nonselected<- wide3[which(wide3$ranking>139),]
regionalAdvanced1list<- nonselected[order(nonselected$ranking),][1:60,1]
library(gtools)
regionalAdvanced1list<- regionalAdvanced1list[mixedorder(gsub("-", "A",regionalAdvanced1list))]
write.csv(regionalAdvanced1list, file="selected from prelim 60 for NN trial.csv")

#sample the 111 for don testing
ents<- wide3[-which(wide3[,1] %in% c('Pio25R74', 'Pio 25R74', '07-4415','02-18228', '07-20728','Kaskaskia','Pio25R47')),1]

#prelim_for_don<- sample(sample(ents), 111)
df<- scb[which(scb$replicate==2),c('observationUnitName', 'studyName', 'plotNumber', 'germplasmName')]
df<- data.frame(df, send_for_DON_testing='No', stringsAsFactors = FALSE)
df[which(df$germplasmName %in% prelim_for_don),'send_for_DON_testing']<- 'Yes'
write.csv(df, file='Pr_Scb_20 rep2 plots for DON testing.csv')

