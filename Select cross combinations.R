#read in analysis results
setwd("~/Documents/GitHub/Wheat-Selection-Decisons-2020")
blups<- read.csv('all 2018-2020 rslts summary.csv')

#read in lines already selected as parents 
setwd("~/Documents/Wheat/2021")
prnts<- read.csv('crossing parents 2021.csv')

#starting price of wheat and soybean, five year average based on macrotrends.net
wheat_price0<- mean(c(5.4621, 4.9414, 4.9757, 4.4014, 4.3945))
soybean_price<- mean(c(9.3785, 8.9298, 9.3456, 9.7820, 9.8753))

#wheat price fcn
wheatPrice<- function(fdk, don, twt, wheat_price0){
  donDiscount<- sqrt(don)*-0.2
  fdkDiscount<- sqrt(fdk)*-0.04
  twtDiscount<- c(58-twt)*-.2
  twtDiscount[which(twtDiscount>0)]<- 0
  wheat_price<- wheat_price0+donDiscount+fdkDiscount+twtDiscount
  return(wheat_price)
}

#net merit function
netMerit<- function(headings, yields, dons, fdks, twt, wheat_price0, soybean_price){
  wheat_price1<- wheatPrice(fdks, dons, twt, wheat_price0)
  soy_yld_gain<- 0.5* (max(headings)-headings)
  soy_profit_gain<- soy_yld_gain*soybean_price
  wheat_profit<- yields*wheat_price1
  total_profit<- wheat_profit + soy_profit_gain
  return(total_profit)
}

###compute net merit of all possible crosses
nms<- as.character(prnts[1:39,2])
combos<- combn(nms,2)
nmvec<-c()
for(i in 1:ncol(combos)){
  pair<- combos[,i]
  ix<- c(match(pair[1], blups[,2]), match(pair[2], blups[,2]))
  F1<- colMeans(blups[ix,3:13])
  nm<- netMerit(F1['Heading.time...Julian.date..JD..'], 
           F1['Grain.yield...bu.ac.South'], F1['FHB.DON.content...ppm.'], 
           F1['FHB.grain.incidence.....'], F1['Grain.test.weight...lbs.bu'],
           wheat_price0, soybean_price)
  if(i==1){
    F1s<- F1
    pairs<-  paste(pair, collapse="/")
  }else{
    F1s<- rbind(F1s, F1)
    pairs<- append(pairs, paste(pair, collapse="/"))
  }
  nmvec<- append(nmvec, nm)
}

df<- data.frame(pairs, F1s)

#make desired crosses table
rankorder<- order(nmvec, decreasing=TRUE)
combos<- rbind(combos, nmvec)
combos<- combos[,rankorder]
table<- t(combos)
table<- as.data.frame(table)
row.names(table)<- c(1:nrow(table))
table<- data.frame(table, rank=c(row.names(table)))

#add the assigned entry numbers to the table
V1no<- prnts[match(table$V1, prnts$names),'Entry.number']
V2no<- prnts[match(table$V2, prnts$names),'Entry.number']
table<- data.frame(V1no, V2no, table)

#eliminate crosses with rank 351 or higher
table[which(as.numeric(as.character(table$rank))>500),'nmvec']<- NA

#sort by entry number
table<- table[order(table$V1no),]
table<- table[order(table$V2no),]

#make cross list by entry number
table$V1no<- as.numeric(as.character(table$V1no))
table$V2no<- as.numeric(as.character(table$V2no))
tball<- matrix(ncol=6)
colnames(tball)<- c('V1no','V2no','V1','V2','nmvec','rank')
for(i in 1:39){
  a<- table[c(which(table$V1no==i)),]
  b<- table[c(which(table$V2no==i)),]
  b2<- b
  b2$V1no<- b$V2no
  b2$V2no<- b$V1no
  b2$V1<- b$V2
  b2$V2<- b$V1
  tb<- rbind(a, b2)
  tb<- tb[order(as.numeric(as.character(tb$rank))),]
  tball<- rbind(tball, tb)
  
}
tball<- tball[-1,]
tball$nmvec<- round(as.numeric(as.character(tball$nmvec)),2)
write.csv(tball, file='selected cross combinations 2021.csv')

tball<- data.frame(pairs=paste(tball$V1, tball$V2, sep="/"), tball)
tball<- tball[order(as.numeric(tball$rank)),]
infotb<- merge(tball, df, by='pairs', sort=FALSE)
infotb[,6]<- as.numeric(infotb[,6])
infotb[,7]<- as.numeric(infotb[,7])
infotb<- rbind(infotb[1:25,], c(rep(NA, 5), colMeans(infotb[,-c(1:5)])))
write.csv(infotb, file='top25.csv')

#make colorized matrix
tball$V1no<- factor(as.character(tball$V1no), levels=as.character(c(1:39)))
tball$V2no<- factor(as.character(tball$V2no), levels=as.character(c(1:39)))
library(wesanderson)
library(ggplot2)
pal <- wes_palette("Zissou1", 100, type = "continuous")
p <- ggplot(tball, aes(V1no, V2no)) + geom_tile(aes(fill = nmvec), colour = "black") + 
  scale_fill_gradientn(colours =pal, na.value="white")
base_size <- 32
p <- p + theme_grey(base_size=base_size)
ggsave("crossingMatrix.jpeg", units="in", width=40, height=27, dpi=300)

p <- p + labs(x="", y="")
p <- p + scale_x_discrete(expand = c(0,0))
p <- p + theme(legend.position="none", axis.ticks=element_blank(), 
               axis.text.x=element_text(size=base_size*0.8, angle=330, hjust = 0, colour="grey50"))
####
q<- read.csv('intercrossparents2021X1.csv')
q[,1]<- gsub("d", "-", q[,1])
q[,1]<- gsub("2021X", "2021X_", q[,1])
write.csv(q, file='v5intercrossparents2021X1.csv',fileEncoding = "UTF-8", row.names=FALSE)

q<- read.csv('intercrosswishlist2021X1.csv')
q[,1]<- gsub("d", "-", q[,1])
q[,1]<- gsub("2021X", "2021X_", q[,1])
q[,2]<- gsub("d", "-", q[,2])
q[,2]<- gsub("2021X", "2021X_", q[,2])
write.csv(q, file='v4intercrosswishlist2021X1.csv', row.names=FALSE,fileEncoding = "UTF-8")


