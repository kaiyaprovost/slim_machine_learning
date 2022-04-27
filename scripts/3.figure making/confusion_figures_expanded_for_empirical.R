dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package,  ")"))))
}
packages = c("igraph", "statnet", "ggnetwork", "ggplot2","intergraph","RColorBrewer")
for (p in packages) {
  dynamic_require(p)
}

## making confusion figures from the probability data for scutulatus
{
scut_age="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-LOCIMISSING-COMBO.STATS_AGE_PROB.txt"
scut_dem="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-LOCIMISSING-COMBO.STATS_DEMOG_PROB.txt"
scut_ibd="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-LOCIMISSING-COMBO.STATS_IBD_PROB.txt"

age=read.table(scut_age,header=F,sep="\t"); colnames(age) = c('1,000,000','120,000','21,000','6,000')
dem=read.table(scut_dem,header=F,sep="\t"); colnames(dem) = c('G','I','P','S'); dem = dem[,c("P","I","G","S")]
ibd=read.table(scut_ibd,header=F,sep="\t"); colnames(ibd) = c('NO','YES')

barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)

age = age[order(age$`1,000,000`,age$`120,000`,age$`21,000`,age$`6,000`),]
dem = dem[order(-dem$P,-dem$I,-dem$G,-dem$S),]
ibd = ibd[order(ibd$NO,ibd$YES),]

barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)

scut_ful="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-LOCIMISSING-COMBO.STATS_PREDICTED.txt"
ful=read.table(scut_ful,header=F,sep="\t"); colnames(ful)=c("DEM","IBD","AGE")
table(ful$DEM)/nrow(ful)
table(ful$IBD)/nrow(ful)
table(ful$AGE)/nrow(ful)
library(dplyr)
ful %>% group_by_all() %>% summarise(COUNT = n())

pdf("age-dem-ibd_cumulative-maximum_scut.pdf",height=8,width=6.5)
par(mfrow=c(3,2),mar=c(1,4,3,0))
barplot(as.matrix(colSums(age)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="AGE",main="By Cumulative Probability")
barplot(as.matrix(table(ful$AGE)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="AGE",main="By Maximum Probability",yaxt="n")

barplot(as.matrix(colSums(dem)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="DEM")
barplot(as.matrix(table(ful$DEM)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="DEM",yaxt="n")

barplot(as.matrix(colSums(ibd)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="IBD")
barplot(as.matrix(table(ful$IBD)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="IBD",yaxt="n")
dev.off()
}

## that but for cardinalis
{
  card_age="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CARDINALIS.COMBO.STATS_AGE_PROB.txt"
  card_dem="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CARDINALIS.COMBO.STATS_DEMOG_PROB.txt"
  card_ibd="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CARDINALIS.COMBO.STATS_IBD_PROB.txt"
  
  age=read.table(card_age,header=F,sep="\t"); colnames(age) = c('1,000,000','120,000','21,000','6,000')
  dem=read.table(card_dem,header=F,sep="\t"); colnames(dem) = c('G','I','P','S'); dem = dem[,c("P","I","G","S")]
  ibd=read.table(card_ibd,header=F,sep="\t"); colnames(ibd) = c('NO','YES')
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  age = age[order(age$`1,000,000`,age$`120,000`,age$`21,000`,age$`6,000`),]
  dem = dem[order(-dem$P,-dem$I,-dem$G,-dem$S),]
  ibd = ibd[order(ibd$NO,ibd$YES),]
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  card_ful="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CARDINALIS.COMBO.STATS_PREDICTED.txt"
  ful=read.table(card_ful,header=F,sep="\t"); colnames(ful)=c("DEM","IBD","AGE")
  table(ful$DEM)/nrow(ful)
  table(ful$IBD)/nrow(ful)
  table(ful$AGE)/nrow(ful)
  library(dplyr)
  ful %>% group_by_all() %>% summarise(COUNT = n())
  
  pdf("age-dem-ibd_cumulative-maximum_card.pdf",height=8,width=6.5)
  par(mfrow=c(3,2),mar=c(1,4,3,0))
  barplot(as.matrix(colSums(age)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="AGE",main="By Cumulative Probability")
  barplot(as.matrix(table(ful$AGE)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="AGE",main="By Maximum Probability")
  
  barplot(as.matrix(colSums(dem)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="DEM")
  barplot(as.matrix(table(ful$DEM)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="DEM")
  
  barplot(as.matrix(colSums(ibd)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="IBD")
  barplot(as.matrix(table(ful$IBD)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="IBD")
  dev.off()
}


## doing the other three species now tat I have them
## ATROX
{
  atrl_age="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-ATROX-LOCI-COMBO.STATS_AGE_PROB.txt"
  atrl_dem="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-ATROX-LOCI-COMBO.STATS_DEMOG_PROB.txt"
  atrl_ibd="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-ATROX-LOCI-COMBO.STATS_IBD_PROB.txt"
  
  age=read.table(atrl_age,header=F,sep="\t"); colnames(age) = c('1,000,000','120,000','21,000','6,000')
  dem=read.table(atrl_dem,header=F,sep="\t"); colnames(dem) = c('G','I','P','S'); dem = dem[,c("P","I","G","S")]
  ibd=read.table(atrl_ibd,header=F,sep="\t"); colnames(ibd) = c('NO','YES')
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  age = age[order(age$`1,000,000`,age$`120,000`,age$`21,000`,age$`6,000`),]
  dem = dem[order(-dem$P,-dem$I,-dem$G,-dem$S),]
  ibd = ibd[order(ibd$NO,ibd$YES),]
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  atrl_ful="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-ATROX-LOCI-COMBO.STATS_PREDICTED.txt"
  ful=read.table(atrl_ful,header=F,sep="\t"); colnames(ful)=c("DEM","IBD","AGE")
  table(ful$DEM)/nrow(ful)
  table(ful$IBD)/nrow(ful)
  table(ful$AGE)/nrow(ful)
  library(dplyr)
  #ful %>% group_by_all() %>% summarise(COUNT = n())
  
  colMeans(age)
  colMeans(dem)
  colMeans(ibd)
  
  pdf("age-dem-ibd_cumulative-maximum_atrl.pdf",height=8,width=6.5)
  par(mfrow=c(3,2),mar=c(1,4,3,0))
  barplot(as.matrix(colSums(age)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="AGE",main="By Cumulative Probability")
  barplot(as.matrix(table(ful$AGE)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="AGE",main="By Maximum Probability")
  
  barplot(as.matrix(colSums(dem)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="DEM")
  barplot(as.matrix(table(ful$DEM)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="DEM")
  
  barplot(as.matrix(colSums(ibd)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="IBD")
  barplot(as.matrix(table(ful$IBD)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="IBD")
  dev.off()
}
## CATENIFER
{
  catl_age="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-CATENIFER-LOCI-COMBO.STATS_AGE_PROB.txt"
  catl_dem="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-CATENIFER-LOCI-COMBO.STATS_DEMOG_PROB.txt"
  catl_ibd="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-CATENIFER-LOCI-COMBO.STATS_IBD_PROB.txt"
  
  age=read.table(catl_age,header=F,sep="\t"); colnames(age) = c('1,000,000','120,000','21,000','6,000')
  dem=read.table(catl_dem,header=F,sep="\t"); colnames(dem) = c('G','I','P','S'); dem = dem[,c("P","I","G","S")]
  ibd=read.table(catl_ibd,header=F,sep="\t"); colnames(ibd) = c('NO','YES')
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  age = age[order(age$`1,000,000`,age$`120,000`,age$`21,000`,age$`6,000`),]
  dem = dem[order(-dem$P,-dem$I,-dem$G,-dem$S),]
  ibd = ibd[order(ibd$NO,ibd$YES),]
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  catl_ful="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-CATENIFER-LOCI-COMBO.STATS_PREDICTED.txt"
  ful=read.table(catl_ful,header=F,sep="\t"); colnames(ful)=c("DEM","IBD","AGE")
  table(ful$DEM)/nrow(ful)
  table(ful$IBD)/nrow(ful)
  table(ful$AGE)/nrow(ful)
  library(dplyr)
  #ful %>% group_by_all() %>% summarise(COUNT = n())
  
  colMeans(age)
  colMeans(dem)
  colMeans(ibd)
  
  pdf("age-dem-ibd_cumulative-maximum_catl.pdf",height=8,width=6.5)
  par(mfrow=c(3,2),mar=c(1,4,3,0))
  barplot(as.matrix(colSums(age)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="AGE",main="By Cumulative Probability")
  barplot(as.matrix(table(ful$AGE)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="AGE",main="By Maximum Probability")
  
  barplot(as.matrix(colSums(dem)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="DEM")
  barplot(as.matrix(table(ful$DEM)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="DEM")
  
  barplot(as.matrix(colSums(ibd)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="IBD")
  barplot(as.matrix(table(ful$IBD)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="IBD")
  dev.off()
}
## GETULA
{
  getl_age="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-GETULA-LOCI-COMBO.STATS_AGE_PROB.txt"
  getl_dem="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-GETULA-LOCI-COMBO.STATS_DEMOG_PROB.txt"
  getl_ibd="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-GETULA-LOCI-COMBO.STATS_IBD_PROB.txt"
  
  age=read.table(getl_age,header=F,sep="\t"); colnames(age) = c('1,000,000','120,000','21,000','6,000')
  dem=read.table(getl_dem,header=F,sep="\t"); colnames(dem) = c('G','I','P','S'); dem = dem[,c("P","I","G","S")]
  ibd=read.table(getl_ibd,header=F,sep="\t"); colnames(ibd) = c('NO','YES')
  
  colMeans(age)
  colMeans(dem)
  colMeans(ibd)
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  age = age[order(age$`1,000,000`,age$`120,000`,age$`21,000`,age$`6,000`),]
  dem = dem[order(-dem$P,-dem$I,-dem$G,-dem$S),]
  ibd = ibd[order(ibd$NO,ibd$YES),]
  
  #barplot(t(as.matrix(age)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(dem)),col=c("grey","red","goldenrod","cornflowerblue"),border=NA)
  #barplot(t(as.matrix(ibd)),col=c("grey","red"),border=NA)
  
  getl_ful="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-GETULA-LOCI-COMBO.STATS_PREDICTED.txt"
  ful=read.table(getl_ful,header=F,sep="\t"); colnames(ful)=c("DEM","IBD","AGE")
  table(ful$DEM)/nrow(ful)
  table(ful$IBD)/nrow(ful)
  table(ful$AGE)/nrow(ful)
  library(dplyr)
  #ful %>% group_by_all() %>% summarise(COUNT = n())
  
  pdf("age-dem-ibd_cumulative-maximum_getl.pdf",height=8,width=6.5)
  par(mfrow=c(3,2),mar=c(1,4,3,0))
  barplot(as.matrix(colSums(age)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="AGE",main="By Cumulative Probability")
  barplot(as.matrix(table(ful$AGE)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="AGE",main="By Maximum Probability")
  
  barplot(as.matrix(colSums(dem)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="DEM")
  barplot(as.matrix(table(ful$DEM)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="DEM")
  
  barplot(as.matrix(colSums(ibd)),col=c("grey","red","goldenrod","cornflowerblue"),ylab="IBD")
  barplot(as.matrix(table(ful$IBD)),col=c("grey","red","goldenrod","cornflowerblue"),xlab="IBD")
  dev.off()
}





## can we do them both at the same time?
age_c=read.table(card_age,header=F,sep="\t"); colnames(age_c) = c('1,000,000','120,000','21,000','6,000')
dem_c=read.table(card_dem,header=F,sep="\t"); colnames(dem_c) = c('G','I','P','S')
ibd_c=read.table(card_ibd,header=F,sep="\t"); colnames(ibd_c) = c('NO','YES')
age_s=read.table(scut_age,header=F,sep="\t"); colnames(age_s) = c('1000K','120K','21K','6K')
dem_s=read.table(scut_dem,header=F,sep="\t"); colnames(dem_s) = c('G','I','P','S')
ibd_s=read.table(scut_ibd,header=F,sep="\t"); colnames(ibd_s) = c('NO','YES')
age_a=read.table(atrl_age,header=F,sep="\t"); colnames(age_a) = c('1,000,000','120,000','21,000','6,000')
dem_a=read.table(atrl_dem,header=F,sep="\t"); colnames(dem_a) = c('G','I','P','S')
ibd_a=read.table(atrl_ibd,header=F,sep="\t"); colnames(ibd_a) = c('NO','YES')
age_g=read.table(getl_age,header=F,sep="\t"); colnames(age_g) = c('1,000,000','120,000','21,000','6,000')
dem_g=read.table(getl_dem,header=F,sep="\t"); colnames(dem_g) = c('G','I','P','S')
ibd_g=read.table(getl_ibd,header=F,sep="\t"); colnames(ibd_g) = c('NO','YES')
age_p=read.table(catl_age,header=F,sep="\t"); colnames(age_p) = c('1,000,000','120,000','21,000','6,000')
dem_p=read.table(catl_dem,header=F,sep="\t"); colnames(dem_p) = c('G','I','P','S')
ibd_p=read.table(catl_ibd,header=F,sep="\t"); colnames(ibd_p) = c('NO','YES')



ful_c = read.table(card_ful,header=F,sep="\t"); colnames(ful_c)=c("DEM","IBD","AGE")
ful_s = read.table(scut_ful,header=F,sep="\t"); colnames(ful_s)=c("DEM","IBD","AGE")
ful_a = read.table(atrl_ful,header=F,sep="\t"); colnames(ful_a)=c("DEM","IBD","AGE")
ful_g = read.table(getl_ful,header=F,sep="\t"); colnames(ful_g)=c("DEM","IBD","AGE")
ful_p = read.table(catl_ful,header=F,sep="\t"); colnames(ful_p)=c("DEM","IBD","AGE")



age_tab = list((table(ful_c$AGE)),(table(ful_s$AGE)),
                (table(ful_a$AGE)),(table(ful_g$AGE)),
                (table(ful_p$AGE)))
age_tab[[3]][4] = 0; names(age_tab[[3]]) = rownames(age_tab[[1]])
age_tab[[4]][3:4] = 0; names(age_tab[[4]]) = rownames(age_tab[[1]])
age_tab[[5]][3:4] = 0; names(age_tab[[5]]) = rownames(age_tab[[1]])
age_tab=matrix(unlist(age_tab),4,5)
rownames(age_tab) = c('1000K','120K','21K','6K')

dem_tab = cbind(as.matrix(table(ful_c$DEM)),as.matrix(table(ful_s$DEM)),
                as.matrix(table(ful_a$DEM)),as.matrix(table(ful_g$DEM)),
                as.matrix(table(ful_p$DEM)))
ibd_tab = cbind(as.matrix(table(ful_c$IBD)),as.matrix(table(ful_s$IBD)),
                as.matrix(table(ful_a$IBD)),as.matrix(table(ful_g$IBD)),
                as.matrix(table(ful_p$IBD)))

library(RColorBrewer)
fullrange = brewer.pal(12,"Paired")
blues=colorRampPalette(fullrange[1:2]) ## pleist
greens=colorRampPalette(fullrange[3:4]) ## plio
reds=colorRampPalette(fullrange[5:6]) ## unclear
oranges=colorRampPalette(fullrange[7:8]) ## plio-pleist
purples=colorRampPalette(fullrange[9:10]) ## mio-plio
browns=colorRampPalette(fullrange[11:12]) ## mio 



pdf("loci_only_cumulative_maximum_prob_comparison_with_legend.pdf",height=8,width=6.5)
par(mfrow=c(3,3),mar=c(2,4,1,0))
barplot(t(rbind(colMeans(age_c),colMeans(age_s),colMeans(age_a),colMeans(age_g),colMeans(age_p))),
        ylab="Support for Predicted Generations Since Divergence",main="By Cumulative Probability",col=brewer.pal(4,"YlOrRd"))
barplot(t(t(age_tab)/colSums(age_tab)),main="By Maximum Probability",yaxt="n",col=brewer.pal(4,"YlOrRd"))
plot.new()
legend("center",rev(names(colMeans(age_c))),fill=rev(brewer.pal(4,"YlOrRd")),bty="n")

barplot(t(rbind(colMeans(dem_c),colMeans(dem_s),colMeans(dem_a),colMeans(dem_g),colMeans(dem_p))),
        ylab="Support for Predicted Demography",col=brewer.pal(4,"PuBuGn"))
barplot(t(t(dem_tab)/colSums(dem_tab)),yaxt="n",col=brewer.pal(4,"PuBuGn"))
plot.new()
legend("center",rev(names(colMeans(dem_c))),fill=rev(brewer.pal(4,"PuBuGn")),bty="n")

barplot(t(rbind(colMeans(ibd_c),colMeans(ibd_s),colMeans(ibd_a),colMeans(ibd_g),colMeans(ibd_p))),
        ylab="Support for Predicted Isolation-by-Distance",names=c("Cc","Cs","Ca","Lg","Pc"),col=brewer.pal(2,"Greys")[c(1,3)])
barplot(t(t(ibd_tab)/colSums(ibd_tab)),names=c("Cc","Cs","Ca","Lg","Pc"),yaxt="n",col=brewer.pal(2,"Greys")[c(1,3)])
plot.new()
legend("center",rev(names(colMeans(ibd_c))),fill=rev(brewer.pal(2,"Greys")[c(1,3)]),bty="n")
dev.off()



## now for the other species
atr_age=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_ATROX.COMBO.STATS_AGE_PROB_MONOONLY.txt",sep="\t",header=F)
atr_dem=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_ATROX.COMBO.STATS_DEMOG_PROB_MONOONLY.txt",sep="\t",header=F)
atr_ibd=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_ATROX.COMBO.STATS_IBD_PROB_MONOONLY.txt",sep="\t",header=F)

get_age=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_GETULA.COMBO.STATS_AGE_PROB_MONOONLY.txt",sep="\t",header=F)
get_dem=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_GETULA.COMBO.STATS_DEMOG_PROB_MONOONLY.txt",sep="\t",header=F)
get_ibd=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_GETULA.COMBO.STATS_IBD_PROB_MONOONLY.txt",sep="\t",header=F)

cat_age=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CATENIFER.COMBO.STATS_AGE_PROB_MONOONLY.txt",sep="\t",header=F)
cat_dem=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CATENIFER.COMBO.STATS_DEMOG_PROB_MONOONLY.txt",sep="\t",header=F)
cat_ibd=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_DEM-IBD-AGE_empiricalcomparison_CATENIFER.COMBO.STATS_IBD_PROB_MONOONLY.txt",sep="\t",header=F)

scu_age=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-GENOME-COMBO.STATS_AGE_PROB.txt",sep="\t",header=F)
scu_dem=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-GENOME-COMBO.STATS_DEMOG_PROB.txt",sep="\t",header=F)
scu_ibd=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_MONO-SCUTULATUS-GENOME-COMBO.STATS_IBD_PROB.txt",sep="\t",header=F)

car_age=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_cardcard16.fasta.FULLSTATS.COMBO.STATS_AGE_PROB.txt",sep="\t",header=F)
car_dem=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_cardcard16.fasta.FULLSTATS.COMBO.STATS_DEMOG_PROB.txt",sep="\t",header=F)
car_ibd=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/DEM-IBD-AGEtrained_neural_network_**NN_expand_1604531680.joblib_X_train.**NN_expand_1604531680.temp_cardcard16.fasta.FULLSTATS.COMBO.STATS_IBD_PROB.txt",sep="\t",header=F)


ages=rbind(atr_age,get_age,cat_age,scu_age,car_age); rownames(ages) = c("ATROX","GETULA","CATENIFER","SCUTULATUS","CARDINALIS"); colnames(ages) = c("1,000,000","120,000","21,000","6,000")
dems=rbind(atr_dem,get_dem,cat_dem,scu_dem,car_dem); rownames(dems) = c("ATROX","GETULA","CATENIFER","SCUTULATUS","CARDINALIS"); colnames(dems) = c('G','I','P','S'); dems = dems[,c("P","I","G","S")]
ibds=rbind(atr_ibd,get_ibd,cat_ibd,scu_ibd,car_ibd); rownames(ibds) = c("ATROX","GETULA","CATENIFER","SCUTULATUS","CARDINALIS"); colnames(ibds) = c("NO","YES")

barplot(t(as.matrix(ages)),col=c("grey","red","goldenrod","cornflowerblue"))
barplot(t(as.matrix(dems)),col=c("grey","red","goldenrod","cornflowerblue"))
barplot(t(as.matrix(ibds)),col=c("grey","red","goldenrod","cornflowerblue"))

means_age = rbind(ages,colMeans(age_c),colMeans(age_s),colMeans(age_a),colMeans(age_g),colMeans(age_p)); rownames(means_age) = c("Ca.G","Lg.G","Pc.G","Cs.G","Cc.G", "Cc.L","Cs.L","Ca.L","Lg.L","Pc.L")
means_dem = rbind(dems,colMeans(dem_c),colMeans(dem_s),colMeans(dem_a),colMeans(dem_g),colMeans(dem_p)); rownames(means_dem) = c("Ca.G","Lg.G","Pc.G","Cs.G","Cc.G", "Cc.L","Cs.L","Ca.L","Lg.L","Pc.L")
means_ibd = rbind(ibds,colMeans(ibd_c),colMeans(ibd_s),colMeans(ibd_a),colMeans(ibd_g),colMeans(ibd_p)); rownames(means_ibd) = c("Ca.G","Lg.G","Pc.G","Cs.G","Cc.G", "Cc.L","Cs.L","Ca.L","Lg.L","Pc.L")

means_age = means_age[c("Ca.G","Ca.L","Cs.G","Cs.L","Lg.G","Lg.L","Pc.G","Pc.L", "Cc.G","Cc.L"),]
means_dem = means_dem[c("Ca.G","Ca.L","Cs.G","Cs.L","Lg.G","Lg.L","Pc.G","Pc.L", "Cc.G","Cc.L"),]
means_ibd = means_ibd[c("Ca.G","Ca.L","Cs.G","Cs.L","Lg.G","Lg.L","Pc.G","Pc.L", "Cc.G","Cc.L"),]

pdf("all_species_comparison_with_legend.pdf",height=8,width=7)
par(mfrow=c(3,2),mar=c(4,4,0.5,0))
layout(matrix(c(1,2,3,4,5,6),3,2,byrow=T),widths=c(3,1))
barplot(t(as.matrix(means_age)),ylab="Proportion Support",col=brewer.pal(4,"YlOrRd"),names.arg=rep("",10),las=1,space=c(0.5,0,0.5,0,0.5,0,0.5,0,0.5,0))
plot.new()
legend("center",rev(colnames(means_age)),fill=rev(brewer.pal(4,"YlOrRd")),bty="n",title="Predicted Age")
barplot(t(as.matrix(means_dem)),ylab="Proportion Support",col=brewer.pal(4,"PuBuGn"),names.arg=rep("",10),las=1,space=c(0.5,0,0.5,0,0.5,0,0.5,0,0.5,0))
plot.new()
legend("center",rev(colnames(means_dem)),fill=rev(brewer.pal(4,"PuBuGn")),bty="n",title="Predicted Demography")
barplot(t(as.matrix(means_ibd)),ylab="Proportion Support",col=brewer.pal(4,"Greys")[c(1,3)],las=1,space=c(0.5,0,0.5,0,0.5,0,0.5,0,0.5,0))
plot.new()
legend("center",rev(colnames(means_ibd)),fill=rev(brewer.pal(4,"Greys")[c(1,3)]),bty="n",title="Predicted IBD")
dev.off()



