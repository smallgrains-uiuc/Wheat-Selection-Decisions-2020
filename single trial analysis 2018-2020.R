setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
library(asreml)
library(asremlPlus)
library(reshape)
data<- read.csv('All 2018-2020 data as of July 9 2020.csv', as.is=TRUE)

############################
##Functions to be used later
############################

#function to check model convergence and update until converged (tolerate a 1.5% change in components)
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >1.5)){
    mod<-suppressWarnings(update(mod))
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
addRows<- function(df, varNm, vec){
  lenvec<- length(vec)
  df[nrow(df)+1,]<-df[nrow(df),]
  df[nrow(df), c(1:(lenvec+1))]<- rep("", lenvec+1)
  df[nrow(df), lenvec+2]<- varNm
  df[nrow(df), -c(1:(lenvec+2))]<- round(vec, 5)
  return(df)
}

#function to get error co-variance matrix
getEcov<- function(mod){  
  vc<- summary(mod)$varcomp
  g<- vc[-c(1:which(row.names(vc)=='units:trait!R')),]
  g<- cbind(data.frame(t(matrix(unlist(strsplit(matrix(unlist(strsplit(row.names(g), split="trait_")), nrow=2)[2,], split=":")), nrow=2))), g)
  m1<- cast(g, X1~X2, value='component')[,-1]
  m1<- as.matrix(m1)
  m2<- cast(g, X2~X1, value='component')[,-1]
  m2<- as.matrix(m2)
  vcovE<- ifelse(is.na(m1), ifelse(is.na(m2), NA, m2), ifelse(is.na(m2), m1, m1 + m2))
  labs<- cast(g, X1~X2, value='component')[,1]
  colnames(vcovE)<- labs
  row.names(vcovE)<- labs
  return(vcovE) 
}

############################
## Data curation
############################

#treat prelims as one trial with separate blocks
#get study names for the prelims
studgrpExp<- c("Pr[0-9]_Urb_19","Pr[0-9]_Scb_19","Pr[0-9]_Car_19","Pr[0-9]_Stj_19")

studGrp<- unique(data[grep(c(studgrpExp[x]), data$studyName),'studyName'])
studGrp<-unique(data[grep(c(studgrpExp[x]), data$studyName),'studyName'])
studGrp<-unique(data[grep(c(studgrpExp[x]), data$studyName),'studyName'])
studGrp<-unique(data[grep(c(studgrpExp[x]), data$studyName),'studyName'])


#add column for check
#pheno<- data.frame(pheno, check="0", stringsAsFactors = FALSE)
#pheno[which(pheno$germplasmName=='07-4415'),'check']<- '1'
#pheno[which(pheno$germplasmName=='Kaskaskia'),'check']<- '1'

############################
##Subset each trial
############################

#get unique study names
stdnms<- unique(data$studyName) 

#indicate which require a single trial analysis summary file
stdnmsCoop<- stdnms[grep(c("A6S|P6S|SU|UE|VT"), stdnms)]

#indicate which have an augmented design
stdnmsAug<- stdnms[grep(c("Aug"), stdnms)]

for(i in 1:length(stdnms)){
  
  #subset single trial
  trl<- data[data$studyName==stdnms[i],]
  
  ############################
  ##Extract design information
  ############################
  
  #get vector of all possible traits
  trtnms<- colnames(data)[c(41:ncol(data)-1)]
  
  #vector of traits used in the selected trial
  ttrt<- selectTraitcols(trl, trtnms)
  
  #check if there is blocking for the traits measured
  minBlkno<- length(unique(na.omit(trl[,c('blockNumber', ttrt)])$blockNumber))
  
  #single-trait or multi-trait model?
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
    fxform<- paste(ttrt, "~1", sep="")
  }
  if(uvvmv== "MV"){
    fxform<- paste('cbind(', paste(ttrt, collapse=", "), ")~1+trait+us(trait):germplasmDbId", sep="")
  }
  #add the blocking factor if any
  if(minBlkno>1 & minBlkno<5){
    fxform<- paste(fxform, "block", sep="+")
  }
  #convert to formula
  fxform<- as.formula(fxform)
  
  #############################
  ## Create the random formula
  #############################
  rform<-NA
  #add the blocking factor if any
  if(minBlkno>=5){
    rform<- "~block"
    #for augmented designs use the checks to estimate the block effect
    if(stdnms[i] %in% stdnmsAug){ 
      rform<- "~at(check, '1'):block"
    }
    rform<- as.formula(rform)
  }
  
  #############################
  ## Convert variables to factors
  #############################
  trl$germplasmDbId<- as.factor(as.character(trl$germplasmDbId))
  trl$block<- as.factor(as.character(trl$block))
  
  #################################
  ## Fit model and extract results
  #################################
  if(uvvmv=='MV'){
    mod<- suppressWarnings(asreml(fixed=fxform, residual=~id(units):us(trait), data=trl, trace=FALSE))
  }
  if(uvvmv=='UV'){
    mod<- suppressWarnings(asreml(fixed=fxform, data=trl, trace=FALSE))
  }
  mod<- mkConv(mod)
  p<- suppressWarnings(predictPlus(mod, classify = clasfy, meanLSD.type='factor.combination', LSDby = 'trait'))
  blues<- p$predictions

  #compute the LSD
  LSDs<- p$LSD[,'meanLSD']
  names(LSDs)<- row.names(p$LSD)
  
  #add study name to the blues table
  df<- data.frame(studyName=stdnms[i], blues)
  
  if(stdnms[i] %in% stdnmsCoop){
    #################################
    ## create basic single trial analysis summary table 
    #################################
    #means
    smryA<- cast(df, studyName+germplasmDbId~trait, value='predicted.value')
    Means<- colMeans(smryA[,-c(1:2)])
    #standard errors
    smryB<- cast(df, studyName+germplasmDbId~trait, value='standard.error')
    MSE<- colMeans(smryB[,-c(1:2)])
    #trial metadata
    metaCols<- c("studyName", "studyDescription", "studyYear", "studyDesign", "locationName", "germplasmName","germplasmDbId")
    meta<- trl[match(smryA$germplasmDbId, trl$germplasmDbId),metaCols]
    #merge results with meta
    df2<- merge(meta, smryA[,-1], by='germplasmDbId')
    #convert factors to characters
    df2[,c(1:7)] <- lapply(df2[,c(1:7)], as.character)
    
    #add means, MSE, LSD, and CV
    df2<- addRows(df2, "MEAN", Means)
    df2<- addRows(df2, "MSE", MSE)
    df2<- addRows(df2, "LSD", LSDs)
    df2<- addRows(df2, "CV", MSE/Means *100)
    df2[,-c(1:7)]<- round(df2[,-c(1:7)],3)
    
    #write to an excel sheet
    excellist<- list(results=df2, rawdata=trl)
    WriteXLS::WriteXLS(excellist, ExcelFileName=paste(stdnms[i], ".xls", sep=""), SheetNames=c('results', 'rawdata'))
  }  

  #################################
  ## combine all analysis results into one table
  #################################    
  if(i==1){
    dfall<- df
  }else{
    dfall<- rbind(dfall, df)
  } 
}



    
                     
                     