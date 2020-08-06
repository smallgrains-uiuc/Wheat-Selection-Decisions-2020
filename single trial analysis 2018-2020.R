setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
library(asreml)
library(asremlPlus)
library(reshape)
data<- read.csv('All 2018-2020 data as of July 31 2020.csv', as.is=TRUE)
row.names(data)<- data$observationUnitName

#exclude the Y3N and Y4N for now
rowixY3Y4<- grep(c("Y3N|Y4N"), data$studyName)
data<- data[-rowixY3Y4,]

#Exclude data-points declared as outliers
data['A5S_Urb_20-plot216',"Heading.time...Julian.date..JD..CO_321.0001233"]<- NA
data['A6S_Urb_19-plot105',"Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot1048',"Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot1049', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot1096', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2025', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_Neo_18-plot2114', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Neo_20-plot518', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Neo_20-plot559', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_Scb_20-plot3019', "FHB.grain.incidence.....CO_321.0001155"] <- NA
data['Adv_Stj_18-plot2009', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_StJ_20-plot317', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Adv_StJ_20-plot322', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_StJ_20-plot520', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Adv_StJ_20-plot539', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['AdvHY_Urb_20-plot30163', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['ORG_Urb_20-plot112', "Septoria.nodorum.leaf.blotch.severity...0.9.percentage.scale.CO_321.0501145"] <- NA
data['Pr_Sbmv_20-plot216', "Soil.borne.mosaic.plant.response...0.9.Response.Scale.CO_321.0501140"] <- NA
data['Pr_Urb_20-plot763', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr2_Neo_18-plot201', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr3_Urb_18-plot214', "Grain.yield...kg.ha.CO_321.0001218"] <- NA
data['Pr5_Urb_19-plot162', "Grain.test.weight...g.l.CO_321.0001210"] <- NA
data['Pr6_Urb_19-plot189', "Grain.test.weight...g.l.CO_321.0001210"] <- NA


#exclude the Aug for now because it requires high memory
rowixAug<- grep(c("Aug"), data$studyName)
data<- data[-rowixAug,]

#sort plot no within trial
data<- data[order(data$plotNumber),]
data<- data[order(data$studyName),]

#make vector of ids that were submitted to cooperatives
coopgids<- c("14-28444","14-28468","14-28307",
  "15-23803","07-19334","07-4415","02-18228",
  "15-2639","15-4957","15-17909","15-30529",
  "14-DC-64-95-118","15-27270","14-11830",
  "14-28462","12-17257", "13-20616")

############################
##Data corrections
############################

pltsGy<- c('AdvHY_Urb_20-plot20162', 'AdvHY_Urb_20-plot30106', 'Adv_Neo_18-plot1095',
'Adv_Neo_18-plot1097', 'Adv_Neo_18-plot1120', 'Adv_Neo_20-plot267',
'Adv_Neo_20-plot340', 'Adv_StJ_20-plot118', 'Adv_StJ_20-plot139', 'Adv_StJ_20-plot511')
data[which(data$observationUnitName %in% pltsGy), "Grain.yield...kg.ha.CO_321.0001218"]<- NA

pltsTw<- c('Adv_Neo_18-plot1024', 'Adv_Neo_18-plot1073', 'Adv_Neo_18-plot1120', 'Adv_Neo_18-plot2001',
'Adv_Neo_18-plot2021', 'Adv_Neo_18-plot2045', 'Adv_Neo_20-plot565', 'Pr_Car_20-plot135')
data[which(data$observationUnitName %in% pltsTw), "Grain.test.weight...g.l.CO_321.0001210"]<- NA

############################
##Functions to be used later
############################

#mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#function to check model convergence and update until converged (tolerate a 1.5% change in components)
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#create vector of traits recorded in the trial
selectTraitcols<- function(trl, trtnms, thresh=0.2){
  ppres<- c()
  for(j in 1:length(trtnms)){
    vec<-trl[,trtnms[j]]
    ppres<- append(ppres, length(na.omit(vec))/length(vec))
  }
  return(trtnms[which(ppres>thresh)])
}

#function to add means, MSE, LSD, and CV
addRows<- function(df, varNm, label, vec){
  lenvec<- length(vec)
  df[nrow(df)+1,]<-df[nrow(df),]
  df[nrow(df), c(1:varNm)]<- rep("", varNm)
  df[nrow(df), varNm]<- label
  df[nrow(df), -c(1:varNm)]<- round(vec, 5)
  return(df)
}

#convert yld to bu/acre
convYld<- function(y){
  x<- y/c(60 * 0.453592 * 2.47105)
  return(x)
}

#convert tw to lbs/bu
convTwt<- function(y){
  x<- y/1000 *2.2046 *35.2391
  return(x)
}  

############################
## Data curation
############################

#treat prelims from before 2020 as one trial with separate blocks in the un-replicated sites
studgrpExp<- c("Pr[0-9]_Car_19","Pr[0-9]_Stj_19","Pr[0-9]_Rid_18","Pr[0-9]_Stj_18",
               "Pr[0-9]_Urb_18","Pr[0-9]_Urb_19", "Pr[0-9]_Scb_19")
for(x in 1:length(studgrpExp)){
  
  #get the study names and row positions of the study group
  studGrp<- unique(data[grep(c(studgrpExp[x]), data$studyName, fixed=FALSE),'studyName'])
  ixgrp<- which(data$studyName %in% studGrp)
  stdNms<- as.character(data[ixgrp,'studyName'])
  
  #make block the prelim number
  block<- matrix(unlist(strsplit(stdNms, split="_")), nrow=3)[1,]
  data[ixgrp,'blockNumber']<- block
  data[ixgrp,'studyName']<- studgrpExp[x]
  
  #determine checks and change entryType info
  mtck<- unique(data.frame(nm=data[ixgrp,'germplasmName'], block))
  cks<- names(which(table(mtck$nm)==length(unique(block))))
  ixgrpCk<- which(data[ixgrp, 'germplasmName'] %in% cks)
  data[ixgrp,'entryType'][ixgrpCk]<- 'check'
  
  #edit design info
  data[ixgrp,'studyDesign']<- 'Augmented RCBD'
}


#edit entry type and design info for the Augmented trials
ixAug<- grep('Aug', data$studyName)
ixAugCk<- which(data[ixAug,'germplasmName'] %in% c("Kaskaskia", "07-4415"))
data[ixAug,][ixAugCk, 'entryType']<- 'check'
data[ixAug,'studyDesign']<- 'Augmented RCBD'

#convert yield and test weight to common units
data[,'Grain.yield...kg.ha.CO_321.0001218']<- convYld(data[,'Grain.yield...kg.ha.CO_321.0001218'])
data[,'Grain.test.weight...g.l.CO_321.0001210']<- convTwt(data[,'Grain.test.weight...g.l.CO_321.0001210'])
data[,"Plant.height...cm.CO_321.0001301"]<- data[,"Plant.height...cm.CO_321.0001301"] *0.393701
colnames(data)<- gsub('kg.ha.CO_321.0001218', "bu.ac", colnames(data))
colnames(data)<- gsub('g.l.CO_321.0001210', "lbs.bu", colnames(data))
colnames(data)<- gsub('..cm.CO_321.0001301', "inches", colnames(data))

#shorten the trait names to avoid errors in model fitting
ixCOs<- grep('CO_', colnames(data))
colnames(data)[ixCOs]<- matrix(unlist(strsplit(colnames(data)[ixCOs], split="CO_")), nrow=2)[1,]
colnames(data)<- gsub("Soil.borne.mosaic.plant.response...0.9.Response.Scale." , "SBMV", colnames(data))

############################
##Subset each trial
############################

#get unique study names
stdnms<- unique(data$studyName)

#indicate which require a single trial analysis summary file
stdnmsCoop<- stdnms[grep(c("NU|PNU|A5S|P5S|A6S|P6S|SU|UE|VT|VE"), stdnms)]


for(i in 1:length(stdnms)){

  #subset single trial
  trl<- data[data$studyName==stdnms[i],]

  ############################
  ##Extract design information
  ############################
  
  #get vector of all possible traits
  trtnms<- colnames(data)[c(41:ncol(data)-1)]
  
  #exclude plots with low plant stands
  ps<- trl[,'Plant.stand...0.9.density.scale.']
  if(any(ps<5, na.rm=T)){
    trl<- trl[-which(ps<5),]
  }

  #exclude milling and baking traits, notes,and plant stand as response variables
  trtnms<-setdiff(trtnms, c("Flour.protein.content.....","Flour.yield.score.....",
          "Grain.hardness...skcs.index.","Lactic.Acid.SRC.score.....",
          "Softness.equivalent.score.....","Sucrose.SRC.score.....",
          'Plant.stand...0.9.density.scale.', 'notes')) 
  
  #vector of traits used in the selected trial
  ttrt<- selectTraitcols(trl, trtnms)
  
  #convert missing lodging scores to zero
  if("Lodging.incidence...0.9.percentage.scale." %in% ttrt){
    #convert missing lodging scores to zero
    vrldg<- var(trl$Lodging.incidence...0.9.percentage.scale., na.rm=T)
    if(vrldg>0){
      trl[which(is.na(trl$Lodging.incidence...0.9.percentage.scale.)),'Lodging.incidence...0.9.percentage.scale.']<- 0
    }
  }

  #check if there is blocking for the traits measured
  minBlkno<- length(unique(na.omit(trl[,c('blockNumber', ttrt)])$blockNumber))
  
  #single-trait or multitrait model?
  if(length(ttrt)==1){
    uvvmv<- "UV"
    clasfy<- 'germplasmDbId'
  }else{
    uvvmv<- "MV"
    clasfy<- 'germplasmDbId:trait'
  }
  
  #############################
  ## Create the fixed formula
  #############################
  if(uvvmv== "UV"){
    fxform<- paste(ttrt, "~1+germplasmDbId", sep="")
  }
  if(uvvmv== "MV"){
    fxform<- paste('cbind(', paste(ttrt, collapse=", "), ")~1+trait+us(trait):germplasmDbId", sep="")
  }
  #add the blocking factor if any
  if(minBlkno>1 & minBlkno<4){
    #for augmented designs use the checks to estimate the block effect
    if(trl$studyDesign[1]=='Augmented RCBD'){ 
      fxform<- paste(fxform, "at(entryType, 'check'):blockNumber", sep="+")
      if(length(unique(trl$replicate))>1){
        fxform<- paste(fxform, "blockNumber:replicate", sep="+")
      }
    }else{
      fxform<- paste(fxform, "blockNumber", sep="+")
    }
  }

  #convert to formula
  fxform<- as.formula(fxform)
  
  #############################
  ## Create the random formula
  #############################
  rform<-NA
  #add the blocking factor if any
  if(minBlkno>=4){
    rform<- "~blockNumber"
    #for augmented designs use the checks to estimate the block effect
    if(trl$studyDesign[1]=='Augmented RCBD'){ 
      rform<- "~at(entryType, 'check'):blockNumber"
      clasfy<- paste(clasfy, ":entryType", sep="")
      if(length(unique(trl$replicate))>1){
        rform<- paste(rform, "blockNumber:replicate", sep="+")
      }
    }
    rform<- as.formula(rform)
  }
  
  #############################
  ## Convert variables to factors
  #############################
  trl$germplasmDbId<- as.factor(as.character(trl$germplasmDbId))
  trl$blockNumber<- as.factor(as.character(trl$blockNumber))
  trl$entryType<- as.factor(as.character(trl$entryType))
  trl$replicate<- as.factor(as.character(trl$replicate))
  
  #################################
  ## Fit model and extract results
  #################################

  if(uvvmv=='MV'){
    if(class(rform)=='logical'){
      mod<- suppressWarnings(asreml(fixed=fxform, residual=~id(units):us(trait), data=trl, trace=FALSE, aom=T, workspace=64e6))
    }else{
      mod<- suppressWarnings(asreml(fixed=fxform, random=rform, residual=~id(units):us(trait), data=trl, trace=FALSE, aom=T, workspace=64e6))
    }
    mod<- mkConv(mod)
    p<- suppressWarnings(predictPlus(mod, classify = clasfy, meanLSD.type='factor.combination', LSDby = 'trait', pworkspace=64e7))
    
  }
  if(uvvmv=='UV'){
    if(class(rform)=='logical'){
      mod<- suppressWarnings(asreml(fixed=fxform, data=trl, trace=FALSE, aom=T,workspace=64e6))
    }else{
      mod<- suppressWarnings(asreml(fixed=fxform, random=rform, data=trl, trace=FALSE, aom=T,workspace=64e6))
    }
    mod<- mkConv(mod)
    p<- suppressWarnings(predictPlus(mod, classify = clasfy, pworkspace=64e6))

  }
  blues<- p$predictions
  
  if('entryType' %in% colnames(blues)){
    blues<- blues[,-match('entryType', colnames(blues))]
    blues<- unique(blues)
  }
  
  #compute the LSD
  LSDs<- p$LSD[,'meanLSD']
  names(LSDs)<- row.names(p$LSD)
  
  #add study name to the blues table
  df<- data.frame(studyName=stdnms[i], blues)
  
  if(uvvmv=='UV'){
    df<- data.frame(studyName=stdnms[i], trait=ttrt, blues)
    df<- df[,c(1,3,2, 4:ncol(df))]
  }
  if(uvvmv=='MV'){
    df<- data.frame(studyName=stdnms[i], blues)
  }
  
  #get residuals
  jpeg(file=paste(stdnms[i], "-residuals.jpeg", sep=""))
  resids<- resid(mod, type="stdCond")
  plot(resid(mod, type="stdCond"), main=stdnms[i])
  dev.off()
  
  #make potential outlier table
  mltTrl<- melt(trl, id.vars=c('observationUnitName','germplasmDbId'), measure.vars=ttrt)
  mltTrl<-mltTrl[order(mltTrl$variable),]
  mltTrl<-mltTrl[order(mltTrl$observationUnitName),]
  residsTab<- cbind(mltTrl, resids)
  outTab<- residsTab[which(sqrt(resids^2)>3),]
  
  if(stdnms[i] %in% stdnmsCoop){
    #################################
    ## create basic single trial analysis summary table 
    ## Not very useful for selection, but some people like to see it
    #################################
    #means
    smryA<- cast(df, studyName+germplasmDbId~trait, value='predicted.value')
    Means<- colMeans(smryA[,-c(1:2)], na.rm=TRUE)
    #standard errors
    smryB<- cast(df, studyName+germplasmDbId~trait, value='standard.error')
    SE<- colMeans(smryB[,-c(1:2)], na.rm=T)
    #trial metadata
    metaCols<- c("studyName", "studyDescription", "studyYear", "studyDesign", "locationName", "germplasmName","germplasmDbId")
    meta<- trl[match(smryA$germplasmDbId, trl$germplasmDbId),metaCols]
    #merge results with meta
    df2<- merge(meta, smryA[,-1], by='germplasmDbId')
    #convert factors to characters
    df2[,c(1:7)] <- lapply(df2[,c(1:7)], as.character)
    #add the IL to Illinois names
    df2$germplasmName[which(df2$germplasmName %in% coopgids)]<- paste("IL", df2$germplasmName[which(df2$germplasmName %in% coopgids)], sep="")
    trl$germplasmName[which(trl$germplasmName %in% coopgids)]<- paste("IL", trl$germplasmName[which(trl$germplasmName %in% coopgids)], sep="")
    #re-order the rows according to the raw data
    trl$replicate<- as.numeric(as.character(trl$replicate))
    trl1<- trl[which(trl$replicate==min(trl$replicate, na.rm=T)),]
    df2<- df2[match(trl1$germplasmName,df2$germplasmName),]

    #count number of rows
    nrowdf2<- nrow(df2)
    
    #add means, SE, LSD, and CV
    df2<- addRows(df2, 7, "MEAN", Means)
    df2<- addRows(df2, 7, "SE", SE)
    df2<- addRows(df2, 7, "LSD", LSDs)
    df2[,-c(1:7)]<- round(df2[,-c(1:7)],2)
  
    #add the number of replicates
    nrep<-c()
    for(k in 1:length(ttrt)){
      tb<- table(na.omit(trl[,c('germplasmDbId', ttrt[k])])[,1])
      tb[which(tb==0)]<- NA
      nrep<- append(nrep, mean(tb, na.rm=T))
    }
    df2<- addRows(df2, 7, "No. of Reps", nrep)
    df2[,-c(1:7)]<- round(df2[,-c(1:7)],2)
    
    #add CV
    df2<- addRows(df2, 7, "CV", round(sqrt(nrep)*SE / Means *100,1))
    
    #add rankings
    df2rnk<- df2[, ttrt]
    for(k in 1:ncol(df2rnk)){
      const<- 1
      if(ttrt[k] %in% c("Grain.test.weight...lbs.bu", "Grain.yield...bu.ac") ){
        const<- -1
      }
      rnk<- rank(df2rnk[1:nrowdf2,k]*const, ties.method='min')
      rnk[which(is.na(df2rnk[1:nrowdf2,k]))]<- NA
    df2rnk[,k]<- c(rnk, rep(NA, nrow(df2)-nrowdf2)) 
    }
    colnames(df2rnk)<- paste(colnames(df2rnk), "Rank", sep="-")
    df3<- cbind(df2, df2rnk)
    df3<- df3[,c(colnames(df3)[1:7], sort(colnames(df3[,-c(1:7)])))]

    #rename rank cols
    colnames(df3)[grep("Rank", colnames(df3))]<- 'rank'
    
    #write to an excel sheet
    excellist<- list(results=df3, rawdata=trl)
    WriteXLS::WriteXLS(excellist, ExcelFileName=paste(stdnms[i], ".xls", sep=""), 
                       SheetNames=c('results', 'rawdata'))
  }  

  #################################
  ## combine all analysis results into one table
  #################################    
  if(i==1){
    dfall<- df
    outTabs<- outTab
  }else{
    dfall<- rbind(dfall, df)
    outTabs<- rbind(outTabs, outTab)
  } 
}

name<- data[match(dfall$germplasmDbId, data$germplasmDbId),'germplasmName']
dfall<- data.frame(name, dfall)
write.csv(dfall, file='predicted values table.csv')
write.csv(outTabs, file='possible outliers.csv')



                     
                     