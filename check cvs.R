
cv<- dfall$standard.error/dfall$predicted.value *100
feedback<- data.frame(study=dfall$studyName, trait=dfall$trait, CV=cv)
wide<- cast(feedback, study~trait, value='CV', fun.aggregate = mean)
wide20<- wide[c(grep('20', wide[,1])), ]
wide19<- wide[c(grep('19', wide[,1])), ]
wide18<- wide[c(grep('18', wide[,1])), ]

a<- data.frame( trait=colnames(wide20[,-c(1)]), avgCV= colMeans(wide20[,-c(1)], na.rm=T), year='20')
b<- data.frame( trait=colnames(wide19[,-c(1)]), avgCV= colMeans(wide19[,-c(1)], na.rm=T), year='19')
c<- data.frame( trait=colnames(wide18[,-c(1)]), avgCV= colMeans(wide18[,-c(1)], na.rm=T), year='18')
all<- rbind(a, b, c)


write.csv(cast(all, trait~year, value='avgCV'), file='cvSummary.csv')
