library(igraph)
library(statnet)
library(ggnetwork)
library(ggplot2)
library(intergraph)
library(RColorBrewer)

## make a comparison between sims and empricail
dffile_all = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/alldf.txt" ## with or without _trim
df = read.csv(dffile_all,sep="\t",skip=1)
head(df)
dim(df)
sum(complete.cases(df))

dffile_trim = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/alldf_trim.txt" ## with or without _trim
dft = read.csv(dffile_trim,sep="\t",skip=0)
head(dft)
dim(dft)
sum(complete.cases(dft))

dfsub = df[,c(1:40)]
sum(complete.cases(dfsub))

dftsub = dft[,c(1:3,5,8:9,11:13,16:17,19:24,27:29,31:33,35,37:40)]
sum(complete.cases(dftsub))

## pcaplot of empirical vs not
dffull = dfsub[complete.cases(dfsub),]
dftfull = dftsub[complete.cases(dftsub),]

df_meta = dffull[,c(1:15)]
dft_meta = dftfull[,1:5]

dffull = sapply(dffull[,c(16:ncol(dffull))], as.numeric)
dftfull = sapply(dftfull[6:ncol(dftfull)], as.numeric)

## sort the data
dffull = dffull[,order(colnames(dffull))]
df_meta = df_meta[order(dffull[,1],dffull[,2],dffull[,3]),]
dffull = dffull[order(dffull[,1],dffull[,2],dffull[,3]),]
df_meta = df_meta[,order(colnames(df_meta))]


dftfull = dftfull[,order(colnames(dftfull))]
dft_meta = dft_meta[order(dftfull[,1],dftfull[,2],dftfull[,3]),]
dftfull = dftfull[order(dftfull[,1],dftfull[,2],dftfull[,3]),]
dft_meta = dft_meta[,order(colnames(dft_meta))]

all.equal(dffull,dftfull)
## dffull and dftfull are exactly the same in theory but not in pactice
pca = prcomp(dffull,center=T,scale=T)
summary(pca)
data = as.data.frame(pca$x)
withpca = cbind(df_meta,data)

empircal = withpca[withpca$DEMOG=="EMPIRICAL",]
sims = withpca[withpca$DEMOG!="EMPIRICAL",]
sims$TREATMENT = paste(sims$DEMOG,sims$IBD,sims$YEAR)
sims$YEAR = droplevels(sims$YEAR)
empircal$SPECIES = droplevels(empircal$SPECIES)
levels(empircal$SPECIES) = sort(unique(empircal$SPECIES))


par(mfrow=c(2,2),mar=c(0,0,0,0))
for(part in 1:4){
plot(withpca$PC1,withpca$PC2,type="n")
palette(c("purple","magenta","cyan","red","green"))
test = lapply(1:length(levels(sims$SIMDEFAULT.)),FUN=function(y) {
  year=levels(sims$SIMDEFAULT.)[y]
  k1000=sims[sims$SIMDEFAULT.==year,16:17]
  hpts=chull(k1000)
  hpts <- c(hpts, hpts[1])
  mean=colMeans(k1000)
  meanx=median(k1000[,1])
  meany=median(k1000[,2])
  
  lines(k1000[hpts,],
        col=y,lty=1,lwd=2)
  coords=k1000[hpts,]
  text(x=meanx,y=meany,labels=year,col=y)
  
  return(coords)
})
palette(c("blue","darkred","darkblue","goldenrod",
          "darkgreen","darkmagenta","cornflowerblue",
          "orange","darkorange"))
test2 = lapply(1+(10*(part-1)):min(10+(10*(part-1)),length(levels(empircal$CHROMOSOME))),FUN=function(y) {
  year=levels(empircal$CHROMOSOME)[y]
  k1000=empircal[empircal$CHROMOSOME==year,16:17]
  hpts=chull(k1000)
  hpts <- c(hpts, hpts[1])
  coords=k1000[hpts,]
  mean=colMeans(k1000)
  lines(k1000[hpts,],
          col=y,lty=2,lwd=1)
  text(x=mean[1],y=mean[2],labels=year,col=y,cex=0.5)
  return(coords)
})
}


par(mfrow=c(2,5),mar=c(0,0,0,0))
for(sppnum in 1:10) {
  print(sppnum)
  plot(withpcat$PC1,withpcat$PC2,type="n")
  palette(c("purple","magenta","cyan","red","green"))
  #points(sims$PC1,
  #       sims$PC2,
  #       col=as.numeric(as.factor(sims$YEAR)),
  #       pch=(as.numeric(as.factor(sims$YEAR))),
  #       cex=1)
  
  test = lapply(1:length(levels(sims$YEAR)),FUN=function(y) {
    year=levels(sims$YEAR)[y]
    k1000=sims[sims$YEAR==year,16:17]
    hpts=chull(k1000)
    hpts <- c(hpts, hpts[1])
    lines(k1000[hpts,],
          col=y,lty=1,lwd=2)
    coords=k1000[hpts,]
    return(coords)
  })
  names(test) = levels(sims$YEAR)
  palette(c("darkred","darkorange","goldenrod",
            "black","brown","darkblue","darkgreen","darkcyan",
            "darkmagenta","darkslategrey"))
  test2 = lapply(1:length(levels(empircal$SPECIES)),FUN=function(y) {
    
    year=levels(empircal$SPECIES)[y]
    k1000=empircal[empircal$SPECIES==year,16:17]
    hpts=chull(k1000)
    hpts <- c(hpts, hpts[1])
    coords=k1000[hpts,]
    if(y==sppnum){
      lines(k1000[hpts,],
            col="black",lty=2,lwd=3)
      text(x=10,y=-5,labels=year)
    }
    return(coords)
  })
  names(test2) = levels(empircal$SPECIES)
}


## just use dftfull

pcat = prcomp(dftfull,center=T,scale=T)
summary(pcat)
datat = as.data.frame(pcat$x)
withpcat = cbind(dft_meta,datat)

empircal = withpcat[withpcat$DEMOG=="EMPIRICAL",]
sims = withpcat[withpcat$DEMOG!="EMPIRICAL",]
sims$TREATMENT = paste(sims$DEMOG,sims$IBD,sims$YEAR)
sims$YEAR = droplevels(sims$YEAR)
empircal$SPECIES = droplevels(empircal$SPECIES)

## plot the pcas relative to the treatments 
par(mfrow=c(2,5),mar=c(0,0,0,0))
for(sppnum in 1:10) {
plot(withpcat$PC1,withpcat$PC2,type="n")
palette(c("purple","magenta","cyan","red","green"))
points(sims$PC1,
       sims$PC2,
       col=as.numeric(as.factor(sims$YEAR)),
       pch=(as.numeric(as.factor(sims$YEAR))),
       cex=1,type="n")

test = lapply(1:length(levels(sims$YEAR)),FUN=function(y) {
  year=levels(sims$YEAR)[y]
  k1000=sims[sims$YEAR==year,6:7]
  hpts=chull(k1000)
  hpts <- c(hpts, hpts[1])
  lines(k1000[hpts,],
        col=y,lty=1,lwd=2)
  coords=k1000[hpts,]
  return(coords)
})
names(test) = levels(sims$YEAR)
palette(c("darkred","darkorange","goldenrod",
                      "black","brown","darkblue","darkgreen","darkcyan",
                      "darkmagenta","darkslategrey"))
test2 = lapply(1:length(levels(empircal$SPECIES)),FUN=function(y) {
  
  year=levels(empircal$SPECIES)[y]
  k1000=empircal[empircal$SPECIES==year,6:7]
  hpts=chull(k1000)
  hpts <- c(hpts, hpts[1])
  coords=k1000[hpts,]
  if(y==sppnum){
  lines(k1000[hpts,],
        col="black",lty=2,lwd=3)
  text(x=10,y=-5,labels=year)
  }
  return(coords)
})
names(test2) = levels(empircal$SPECIES)
}


points_in_polys = lapply(1:nrow(empircal),FUN=function(x) {
  polyS = empircal[x,6:7]
  species=as.character(empircal$SPECIES[x])
  
  df = data.frame(matrix(nrow=nrow(polyS),ncol=length(test)+3))
  for(i in 1:length(test)){
    polyY = as.data.frame(test[i][1])
    pointtest = point.in.polygon(point.x=polyS$PC1,
                     point.y=polyS$PC2,
                     pol.x=polyY[,1],
                     pol.y=polyY[,2])
    df[,i] = c(pointtest)
  }
  df[,length(test)+1] = polyS$PC1
  df[,length(test)+2] = polyS$PC2
  df[,length(test)+3] = species
  colnames(df) = c(names(test),"PC1","PC2","SPECIES")
  return(df)
}
) 
points_in_polys = do.call("rbind", points_in_polys)

poi_in_pol = points_in_polys[,c(1:4,7)]
poi_in_pol %>% dplyr::group_by_all() %>% dplyr::summarise(COUNT = dplyr::n()) -> poi_in_pol
poi_in_pol = data.frame(poi_in_pol)

poi_in_pol=poi_in_pol[order(poi_in_pol[,5],-poi_in_pol[,6]),]




aggd = aggregate(points_in_polys[,1:4],by=list(points_in_polys[,7]),FUN=sum)


aggd = aggregate(points_in_polys[,1:4],by=list(points_in_polys[,7]),FUN=sum)
aggd = cbind(aggd,summary(empircal$SPECIES))
colnames(aggd)[ncol(aggd)] = "TOTAL"
aggd = aggd[,2:ncol(aggd)]

# for (i in 1:length(sort(unique(empircal$SPECIES)))) {
#   category = sort(unique(empircal$SPECIES))[i]
#   d = empircal[empircal$SPECIES==category,1:7]
#   colors=c("darkred","darkorange","goldenrod",
#            "black","brown","darkblue","darkgreen","darkcyan",
#            "darkmagenta","darkslategrey")
#   car::ellipse(center = colMeans(d[,6:7],na.rm = T), 
#                shape = cov(d[,6:7],use="pairwise.complete.obs"),
#                center.pch=i,center.cex=2,lwd=2,
#                radius = sqrt(qchisq(.5, df=2)),col = colors[i],
#                lty=3)
# }

## now plot pca with respect to empircal vs simulation

png("/Users/kprovost/Dropbox (AMNH)/Dissertation/simulations_vs_empirical.png",
    width=500,height=500)
plot(withpcat$PC1,withpcat$PC2,type="n")
points(sims$PC1,sims$PC2,col="lightgrey",pch=as.numeric(as.factor(sims$TREATMENT)))
points(empircal$PC1,empircal$PC2,col="darkgrey",pch=as.numeric(as.factor(empircal$SPECIES)))
for (i in 1:length(sort(unique(empircal$SPECIES)))) {
  category = sort(unique(empircal$SPECIES))[i]
  d = empircal[empircal$SPECIES==category,1:7]
  colors=c("darkred","darkorange","goldenrod",
           "black","brown","darkblue","darkgreen","darkcyan",
           "darkmagenta","darkslategrey")
  car::ellipse(center = colMeans(d[,6:7],na.rm = T), 
               shape = cov(d[,6:7],use="pairwise.complete.obs"),
               center.pch=i,center.cex=2,lwd=2,
               radius = sqrt(qchisq(.5, df=2)),col = colors[i],
               lty=3)
}
for (i in 1:length(sort(unique(sims$TREATMENT)))) {
  category = sort(unique(sims$TREATMENT))[i]
  d = sims[sims$TREATMENT==category,1:7]
  colors=c(rep("red",4),
           rep("orange",4),
           rep("yellow",4),
           rep("green",4),
           rep("blue",4),
           rep("cyan",4),
           rep("magenta",4),
           rep("purple",4))
  car::ellipse(center = colMeans(d[,6:7],na.rm = T), 
               shape = cov(d[,6:7],use="pairwise.complete.obs"),
               center.pch=ceiling(i/4),center.cex=2,lwd=2,
               radius = sqrt(qchisq(.5, df=2)),col = colors[i],
               lty=2)
}
dev.off()

## make a pca plot of the sumstats! 

{
  ## with everything at once
  file = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/TESTING/RECAPSTATS/MERGEDPOP_FINAL_extratrimme_notcollinear.txt",sep="\t")
  #file = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/Confusion_Matrix_DEM-IBD-AGE_fulltrain.csv",sep=",")
  
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/correlations_final_all11_26sept2019.pdf",
      width=15,height=5)
  par(xpd=T)
  par(mfrow=c(1,4))
  for (i in 6:ncol(file)) {
    boxplot(file[,i]~file$year,main=names(file)[i],notch=T)
    boxplot(file[,i]~file$demog,main=names(file)[i],notch=T)
    boxplot(file[,i]~file$ibd,main=names(file)[i],notch=T)
    boxplot(file[,i]~file$numrun,main=names(file)[i],notch=T)
  }
  dev.off()
  
  #k21 = file[file$y=="21000",]
  
  pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/correlations_graph_final_all11_26sept2019.pdf",
      width=7,height=7)
  for (j in 6:(ncol(file)-1)) {
    #if (!(j %in% c(32:37))) {
    #par(xpd=T)
    print(j)
    
    plot(file[,j],file[,j+1],
         #col=as.factor(file$numrun),
         col="lightgrey",
         ylab=names(file)[j+1],xlab=names(file)[j])
    for (i in 1:8) {
      category = i
      d = file[file$numrun==i,]
      
      colors=c("darkred","darkred","goldenrod","goldenrod",
               "black","black","blue","blue")
      
      car::ellipse(center = colMeans(file[file$numrun==category,c(j,j+1)],na.rm = T), 
                   shape = cov(file[file$numrun==category,c(j,j+1)],use="pairwise.complete.obs"),
                   center.pch=i,center.cex=2,lwd=2,
                   radius = sqrt(qchisq(.5, df=2)),col = colors[i],
                   lty=1+(i %% 2))
    }
    legend("top", inset=c(0,-0.1),legend=unique(file$numrun),horiz=T,
           col=colors,lty=c(1,2,1,2,1,2,1,2),lwd=1,xpd=T,bty="n",
           pch=1:length(unique(file$numrun)))
    
    #}
  }
  dev.off()
  
  temp = na.omit(file)
  vals = temp[,1:5]
  temp = (temp[,6:ncol(temp)])
  df = sapply(temp, as.numeric)
  pca = prcomp(df,center=T,scale=T)
  summary(pca)
  data = as.data.frame(pca$x)
  withpca = cbind(vals,data)
  
  plot(withpca$PC1,withpca$PC2,col=as.factor(withpca$year),
       pch=as.numeric(withpca$year))
  
  plot(withpca$PC3,withpca$PC2,col=as.factor(withpca$year),
       pch=as.numeric(withpca$year))
  
  aggregate(withpca$PC1 + withpca$PC2 + withpca$PC3 ~ withpca$year, FUN=mean)
  
  rgl::plot3d(withpca$PC1, withpca$PC2, withpca$PC3,col=as.numeric(as.factor(withpca$year)))
  
  boxplot(withpca$PC4~withpca$year)
  
  #withpca = withpca[withpca$year=="120000",]
  pdf("results_of_scaling.pdf")
  par(mfrow=c(2,2),mar=c(2.2,2.2,2.2,0.2))
  
  models=c("x","x","Isolation","Isolation IBD","Gene Flow","Gene Flow IBD")
  
  for (modelnum in 3:6) {
    
    plot(withpca$PC2,withpca$PC1,col=rgb(0.9,0.9,0.9),
         type="n",main=models[modelnum],
         ylab="PC1",xlab="PC2")
    points(withpca$PC2[withpca$numrun==modelnum],withpca$PC1[withpca$numrun==modelnum],col=rgb(0.9,0.9,0.9))
    for(scale in c(0.25,0.5,1,2,4)) {
      d = withpca[withpca$scale==scale,]
      i = which(c(0.25,0.5,1,2,4) == scale)
      
      colors=c("red","goldenrod","lightblue","blue",
               "black","grey","cornflowerblue","lightblue")
      
      car::ellipse(center = colMeans(d[d$numrun==modelnum,c("PC2","PC1")]), 
                   shape = cov(d[d$numrun==modelnum,c("PC2","PC1")]),
                   center.pch=i,center.cex=2,lwd=1,lty=2,
                   radius = sqrt(qchisq(.5, df=2)),col = colors[i])
    }
    legend("topleft",legend=c(0.25,0.5,1,2,4),
           col=colors,lwd=1,lty=2,
           pch=1:length(unique(withpca$scale)),
           bty="n")
  }
  
  dev.off()
  
  
  
  
  plot(withpca$PC2,withpca$PC1,col=as.factor(withpca$numrun),
       type="n")
  for (i in 3:6) {
    category = i
    d = withpca[withpca$numrun==i,]
    
    colors=c("red","pink","goldenrod","lightgoldenrod",
             "grey","lightgrey","cornflowerblue","lightblue")
    
    car::ellipse(center = colMeans(withpca[withpca$numrun==category,c("PC2","PC1")]), 
                 shape = cov(withpca[withpca$numrun==category,c("PC2","PC1")]),
                 center.pch=i,center.cex=2,lwd=1,
                 radius = sqrt(qchisq(.5, df=2)),col = colors[i],
                 lty=1+(i %% 2))
  }
  legend("topright",legend=unique(withpca$numrun),
         col=colors,lty=c(1,2,1,2,1,2,1,2),lwd=1,
         pch=1:length(unique(withpca$numrun)))
  
  
  model = lm(numrun ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + 
               PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14 + PC15, data = withpca)
  summary(model)
  
  
  gn6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-GENEFLOW-NOIBD-6K.POPGENOME.STATS",sep="	",header=F)
  gn21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-GENEFLOW-NOIBD-21K.POPGENOME.STATS",sep="	",header=F)
  gn120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-GENEFLOW-NOIBD-120K.POPGENOME.STATS",sep="	",header=F)
  gy6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-GENEFLOW-YESIBD-6K.POPGENOME.STATS",sep="	",header=F)
  gy21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-GENEFLOW-YESIBD-21K.POPGENOME.STATS",sep="	",header=F)
  gy120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-GENEFLOW-YESIBD-120K.POPGENOME.STATS",sep="	",header=F)
  in6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-ISOLATION-NOIBD-6K.POPGENOME.STATS",sep="	",header=F)
  in21= read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-ISOLATION-NOIBD-21K.POPGENOME.STATS",sep="	",header=F)
  in120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-ISOLATION-NOIBD-120K.POPGENOME.STATS",sep="	",header=F)
  iy6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-ISOLATION-YESIBD-6K.POPGENOME.STATS",sep="	",header=F)
  iy21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-ISOLATION-YESIBD-21K.POPGENOME.STATS",sep="	",header=F)
  iy120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-ISOLATION-YESIBD-120K.POPGENOME.STATS",sep="	",header=F)
  pn6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-PANMIXIA-NOIBD-6K.POPGENOME.STATS",sep="	",header=F)
  pn21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-PANMIXIA-NOIBD-21K.POPGENOME.STATS",sep="	",header=F)
  pn120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-PANMIXIA-NOIBD-120K.POPGENOME.STATS",sep="	",header=F)
  py6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-PANMIXIA-YESIBD-6K.POPGENOME.STATS",sep="	",header=F)
  py21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-PANMIXIA-YESIBD-21K.POPGENOME.STATS",sep="	",header=F)
  py120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-PANMIXIA-YESIBD-120K.POPGENOME.STATS",sep="	",header=F)
  sn6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-SECCON-NOIBD-6K.POPGENOME.STATS",sep="	",header=F)
  sn21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-SECCON-NOIBD-21K.POPGENOME.STATS",sep="	",header=F)
  sn120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-SECCON-NOIBD-120K.POPGENOME.STATS",sep="	",header=F)
  sy6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-SECCON-YESIBD-6K.POPGENOME.STATS",sep="	",header=F)
  sy21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-SECCON-YESIBD-21K.POPGENOME.STATS",sep="	",header=F)
  sy120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/DEM-IBD-AGE_COMBINED-SECCON-YESIBD-120K.POPGENOME.STATS",sep="	",header=F)
  
  namesrows = c(rep("gn6",nrow(gn6)),
                rep("gn21",nrow(gn21)),
                rep("gn120",nrow(gn120)),
                rep("gy6",nrow(gy6)),
                rep("gy21",nrow(gy21)),
                rep("gy120",nrow(gy120)),
                rep("in6",nrow(in6)),
                rep("in21",nrow(in21)),
                rep("in120",nrow(in120)),
                rep("iy6",nrow(iy6)),
                rep("iy21",nrow(iy21)),
                rep("iy120",nrow(iy120)),
                rep("pn6",nrow(pn6)),
                rep("pn21",nrow(pn21)),
                rep("pn120",nrow(pn120)),
                rep("py6",nrow(py6)),
                rep("py21",nrow(py21)),
                rep("py120",nrow(py120)),
                rep("sn6",nrow(sn6)),
                rep("sn21",nrow(sn21)),
                rep("sn120",nrow(sn120)),
                rep("sy6",nrow(sy6)),
                rep("sy21",nrow(sy21)),
                rep("sy120",nrow(sy120)))
  
  temp = rbind(gn6,gn21,gn120,gy6,gy21,gy120,in6,in21,in120,iy6,iy21,iy120,pn6,pn21,pn120,py6,py21,py120,sn6,sn21,sn120,sy6,sy21,sy120)
  temp$NAMES = namesrows
  
  temp1 = temp[,c(1:2,8:9,11:12,17:18,22,24:27,29:32,34:44,47)]
  temp1[temp1 == -100] = NA
  temp1 = na.omit(temp1)
  
  pca = prcomp(temp1[,-ncol(temp1)],center=T,scale=T)
  summary(pca)
  data = as.data.frame(pca$x)
  data$NAMES = temp1$NAMES
  
  pca6 = data[data$NAMES %in% c("pn6","py6","in6","iy6","gn6","gy6","sn6","sy6"),]
  pca21 = data[data$NAMES %in% c("pn21","py21","in21","iy21","gn21","gy21","sn21","sy21"),]
  withpca = data[data$NAMES %in% c("pn120","py120","in120","iy120","gn120","gy120","sn120","sy120"),]
  
  plot120 = withpca[,c("PC2","PC3","NAMES")]
  plot(plot120[,1],plot120[,2],col=as.factor(plot120[,3]),
       type="n",main="120000")
  for (i in 1:length(unique(plot120$NAMES))) {
    category = unique(plot120$NAMES)[i]
    d = plot120[plot120$NAMES==category,]
    
    colors=c("red","pink","goldenrod","lightgoldenrod",
             "grey","lightgrey","cornflowerblue","lightblue")
    
    car::ellipse(center = colMeans(plot120[plot120$NAMES==category,1:2]), 
                 shape = cov(plot120[plot120$NAMES==category,1:2]),
                 center.pch=i,center.cex=2,lwd=1,
                 radius = sqrt(qchisq(.95, df=2)),col = colors[i],
                 lty=1+(i %% 2))
  }
  legend("bottomleft",legend=unique(plot120$NAMES),
         col=colors,lty=c(1,2,1,2,1,2,1,2),lwd=1,
         pch=1:length(unique(plot120$NAMES)))
  
  
}


#####

{
  
  path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/"
  setwd(path)
  
  matlist = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_6K_11.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_6K_33.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_21K_11.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_21K_33.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_120K_11.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_120K_33.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_1000K_11.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_1000K_33.txt")
  
  #matlist = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/SLIMMATRIX_ALL_11.txt"
  
  names_mod = c("PAN","PAN_I","ISO","ISO_I","GF","GF_I","SEC","SEC_I")
  
  for (mat in matlist) {
    
    print(mat)  
    suffix = paste(strsplit(mat,"_")[[1]][4:5],collapse="_")
    
    mat_6k = as.matrix(read.csv(mat,
                                sep="\t",
                                header=T,
                                row.names = 1))
    colnames(mat_6k) = names_mod
    rownames(mat_6k) = names_mod
    
    vertcolors = c("lightgrey","grey","pink","red","lightgoldenrod","goldenrod","lightblue","cornflowerblue")
    reorder = c(2,6,1,5,7,3,8,4)
    
    mat_6k = mat_6k[reorder,reorder]
    
    vertcolors = vertcolors[reorder]
    
    pdfname = paste(path,"confusion_figure_",suffix,".pdf",sep="")
    pngname = paste(path,"confusion_figure_",suffix,".png",sep="")
    
    #pdf(pdfname,
    #    h=8,w=10)
    png(pngname,h=600,w=750)
    par(mar=c(0,2.5,1.5,2.5),bg=NA)
    
    net=igraph::graph.adjacency(mat_6k,mode="directed",weighted=TRUE,diag=TRUE) 
    igraph::plot.igraph(net,vertex.label=V(net)$name,main=suffix,
                #layout=layout_on_grid,
                #layout=x,
                layout=layout_in_circle,
                vertex.shape="circle",
                vertex.size=35,
                #vertex.size=5,
                vertex.size2=20,
                #vertex.color="white",
                vertex.color=vertcolors,
                #vertex.color=c("grey","red","goldenrod","cornflowerblue"),
                #vertex.label.color="black",
                edge.color="black",
                edge.width=E(net)$weight*2, 
                #edge.with=1,
                #edge.lty=2,
                edge.arrow.size=0.5,edge.curved=0.2)
    
    ## need to add a legend with weights 
    ## intro -- really broad barriers or sentence 1 CFB -- CFB easier to write but could be hard to be impactful if ot super well done 
    ## figures without IBD and one with? 
    ## igraph plot in ggplot to add image? 
    
    #legend("bottomleft",legend=c(seq(10,100,10)),col="black",lwd=c(seq(1,10)),
    #       bty="n")
    legend("bottomleft",legend=c("10-24%","25-49%","50-74%","75%+"),col="black",lwd=c(1,2,3,4)*2,
           bty="n")
    
    dev.off()
  }
  
  
  
  #/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/SLIMMATRIX_6K_DEMOG.txt
  
  ##DEMOGRAPHY
  
  matlist = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_6K_DEMOG.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_21K_DEMOG.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_120K_DEMOG.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_1000K_DEMOG.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_ALL_DEMOG.txt")
  
  names_mod = c("PAN","ISO","GF","SEC")
  
  for (mat in matlist) {
    
    print(mat)  
    suffix = paste(strsplit(mat,"_")[[1]][4:5],collapse="_")
    
    mat_6k = as.matrix(read.csv(mat,
                                sep="\t",
                                header=T,
                                row.names = 1))
    colnames(mat_6k) = names_mod
    rownames(mat_6k) = names_mod
    
    vertcolors = c("grey","red","goldenrod","cornflowerblue")
    
    reorder = c(1,3,4,2)
    
    mat_6k = mat_6k[reorder,reorder]
    vertcolors = vertcolors[reorder]
    
    pdfname = paste(path,"confusion_figure_",suffix,".pdf",sep="")
    pngname = paste(path,"confusion_figure_",suffix,".png",sep="")
    
    #pdf(pdfname,
    #    h=8,w=10)
    png(pngname,h=600,w=750)
    par(mar=c(0,2.5,1.5,2.5),bg=NA)
    
    net=graph.adjacency(mat_6k,mode="directed",weighted=TRUE,diag=TRUE) 
    plot.igraph(net,vertex.label=V(net)$name,main=suffix,
                #layout=layout_on_grid,
                #layout=x,
                layout=layout_in_circle,
                vertex.shape="circle",
                vertex.size=35,
                vertex.size2=20,
                vertex.color=vertcolors,
                edge.color="black",
                edge.width=E(net)$weight*2, 
                #edge.with=1,
                #edge.lty=2,
                edge.arrow.size=0.5,edge.curved=T)
    
    ## need to add a legend with weights 
    ## intro -- really broad barriers or sentence 1 CFB -- CFB easier to write but could be hard to be impactful if ot super well done 
    ## figures without IBD and one with? 
    ## igraph plot in ggplot to add image? 
    
    #legend("bottomleft",legend=c(seq(10,100,10)),col="black",lwd=c(seq(1,10)),
    #       bty="n")
    legend("bottomleft",legend=c("10-24%","25-49%","50-74%","75%+"),col="black",lwd=c(1,2,3,4)*2,
           bty="n")
    
    dev.off()
  }
  
  
  matlist = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_6K_IBD.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_21K_IBD.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_120K_IBD.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_1000K_IBD.txt",
              "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_ALL_IBD.txt")
  
  names_mod = c("NO","YES")
  
  for (mat in matlist) {
    
    print(mat)  
    suffix = paste(strsplit(mat,"_")[[1]][4:5],collapse="_")
    
    mat_6k = as.matrix(read.csv(mat,
                                sep="\t",
                                header=T,
                                row.names = 1))
    colnames(mat_6k) = names_mod
    rownames(mat_6k) = names_mod
    
    pdfname = paste(path,"confusion_figure_",suffix,".pdf",sep="")
    pngname = paste(path,"confusion_figure_",suffix,".png",sep="")
    
    #pdf(pdfname,
    #    h=8,w=10)
    png(pngname,h=600,w=750)
    par(mar=c(0,2.5,1.5,2.5),bg=NA)
    net=graph.adjacency(mat_6k,mode="directed",weighted=TRUE,diag=TRUE) 
    plot.igraph(net,vertex.label=V(net)$name,main=suffix,
                #layout=layout_on_grid,
                #layout=x,
                layout=layout_in_circle,
                vertex.shape="circle",
                vertex.size=35,
                vertex.size2=20,
                vertex.color=c("grey","red","goldenrod","cornflowerblue"),
                edge.color="black",
                edge.width=E(net)$weight*2, 
                #edge.with=1,
                #edge.lty=2,
                edge.arrow.size=0.5,edge.curved=T)
    
    ## need to add a legend with weights 
    ## intro -- really broad barriers or sentence 1 CFB -- CFB easier to write but could be hard to be impactful if ot super well done 
    ## figures without IBD and one with? 
    ## igraph plot in ggplot to add image? 
    
    #legend("bottomleft",legend=c(seq(10,100,10)),col="black",lwd=c(seq(1,10)),
    #       bty="n")
    legend("bottomleft",legend=c("10-24%","25-49%","50-74%","75%+"),col="black",lwd=c(1,2,3,4)*2,
           bty="n")
    
    dev.off()
  }
  
  matlist = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/images/datasets/SLIMMATRIX_ALL_AGE.txt")
  
  names_mod = c("1000K","120K","21K","6K")
  
  for (mat in matlist) {
    
    print(mat)  
    suffix = paste(strsplit(mat,"_")[[1]][4:5],collapse="_")
    
    mat_6k = as.matrix(read.csv(mat,
                                sep="\t",
                                header=T,
                                row.names = 1))
    colnames(mat_6k) = names_mod
    rownames(mat_6k) = names_mod
    
    pdfname = paste(path,"confusion_figure_",suffix,".pdf",sep="")
    pngname = paste(path,"confusion_figure_",suffix,".png",sep="")
    
    #pdf(pdfname,
    #    h=8,w=10)
    png(pngname,h=600,w=750)
    par(mar=c(0,2.5,1.5,2.5))
    net=graph.adjacency(mat_6k,mode="directed",weighted=TRUE,diag=TRUE) 
    plot.igraph(net,vertex.label=V(net)$name,main=suffix,
                #layout=layout_on_grid,
                #layout=x,
                layout=layout_in_circle,
                vertex.shape="circle",
                vertex.size=35,
                vertex.size2=20,
                vertex.color=c("grey","red","goldenrod","cornflowerblue"),
                edge.color="black",
                edge.width=E(net)$weight, 
                #edge.with=1,
                #edge.lty=2,
                edge.arrow.size=0.5,edge.curved=T)
    
    ## need to add a legend with weights 
    ## intro -- really broad barriers or sentence 1 CFB -- CFB easier to write but could be hard to be impactful if ot super well done 
    ## figures without IBD and one with? 
    ## igraph plot in ggplot to add image? 
    
    #legend("bottomleft",legend=c(seq(10,100,10)),col="black",lwd=c(seq(1,10)),
    #       bty="n")
    legend("bottomleft",legend=c("10-24%","25-49%","50-74%","75%+"),col="black",lwd=c(1,2,3,4)*2,
           bty="n")
    
    dev.off()
  }
  
  
  
  
  
}
#####

x = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/test_accuracy.txt",sep="\t")
names(x) = c("model","year","dem","ibd","modelnum","RIGHT.D","RIGHT.I","RIGHT.M","CONFD","CONFI","CONF")

#CONF=c(53.27102804,21.42857143,75,71.28712871,69.44444444,66.66666667,77.31958763,74.74747475,27.27272727,16.37931034,53.40909091,69.87951807,41.11111111,71.59090909,76.31578947,75.58139535,24.75247525,30.76923077,54.05405405,56.70103093,33.33333333,64.64646465,86.02150538,89.47368421)
#CONFD=c(34.57943925,3.571428571,68.75,62.37623762,61.11111111,57.47126437,72.16494845,68.68686869,17.27272727,0.862068966,46.59090909,61.44578313,36.66666667,64.77272727,67.10526316,70.93023256,4.95049505,12.82051282,38.73873874,41.2371134,28.88888889,52.52525253,81.72043011,87.36842105)
#CONFI=c(21.4953271,18.75,18.75,22.77227723,25.92592593,14.94252874,17.5257732,12.12121212,12.72727273,16.37931034,15.90909091,18.07228916,11.11111111,19.31818182,36.84210526,11.62790698,19.8019802,24.35897436,27.92792793,22.68041237,14.07407407,19.19191919,31.1827957,14.73684211)
#dem=c("P","P","I","I","G","G","S","S","P","P","I","I","G","G","S","S","P","P","I","I","G","G","S","S")
#ibd=c("N","Y","N","Y","N","Y","N","Y","N","Y","N","Y","N","Y","N","Y","N","Y","N","Y","N","Y","N","Y")
#model=c("M1","M2","M3","M4","M5","M6","M7","M8","M1","M2","M3","M4","M5","M6","M7","M8","M1","M2","M3","M4","M5","M6","M7","M8")
#year=c("6K","6K","6K","6K","6K","6K","6K","6K","21K","21K","21K","21K","21K","21K","21K","21K","120K","120K","120K","120K","120K","120K","120K","120K")
#modelnum=c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8)
#x = data.frame(CONF,CONFD,CONFI,dem,ibd,model,year,modelnum)

k6 = x[x$year=="6K",]
k21 = x[x$year=="21K",]
k120 = x[x$year=="120K",]
k1000 = x[x$year=="1000K",]
kALL = x[x$year=="ALLK",]

colors = RColorBrewer::brewer.pal(12, "Paired")

#png("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/model_misclassification_figures_noall.png",
#    height = 300,width = 900,pointsize=24)
pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/model_misclassification_figures_noall.pdf",
    height=2,width=7)
par(mfrow=c(1,3),mar=c(4,4,1,1))

plot(k6$modelnum,k6$CONF,ylim=c(0,100),col=colors[1],type="l",pch="a",lty=1,lwd=2,
     xlab=c(""),ylab="Misclassification",main="Overall Confusion",
     xaxt="n")
axis(1,at=1:8,labels=c("P","P.I","I","I.I","G","G.I","S","S.I"))
points(k6$modelnum,k21$CONF,ylim=c(0,100),col=colors[4],type="l",pch="b",lty=2,lwd=2)
points(k6$modelnum,k120$CONF,ylim=c(0,100),col=colors[8],type="l",pch="c",lty=3,lwd=2)
points(k6$modelnum,k1000$CONF,ylim=c(0,100),col=colors[6],type="l",pch="d",lty=4,lwd=2)
#points(k6$modelnum,kALL$CONF,ylim=c(0,100),col=colors[10],type="l",pch="e",lty=5,lwd=2)
#abline(h=100*(7/8),col="black",lty=6)

plot(k6$modelnum,k6$CONFD,ylim=c(0,100),col=colors[1],type="l",pch="a",lty=1,lwd=2,
     xlab=c("Model"),ylab="",main="Demography Confusion",
     xaxt="n")
axis(1,at=1:8,labels=c("P","P.I","I","I.I","G","G.I","S","S.I"))
points(k6$modelnum,k21$CONFD,ylim=c(0,100),col=colors[4],type="l",pch="b",lty=2,lwd=2)
points(k6$modelnum,k120$CONFD,ylim=c(0,100),col=colors[8],type="l",pch="c",lty=3,lwd=2)
points(k6$modelnum,k1000$CONFD,ylim=c(0,100),col=colors[6],type="l",pch="d",lty=4,lwd=2)
#points(k6$modelnum,kALL$CONFD,ylim=c(0,100),col=colors[10],type="l",pch="e",lty=5,lwd=2)
#abline(h=100*(3/4),col="black",lty=6)

plot(k6$modelnum,k6$CONFI,ylim=c(0,100),col=colors[1],type="l",pch="a",lty=1,lwd=2,
     xlab=c(""),ylab="",main="IBD Confusion",
     xaxt="n")
axis(1,at=1:8,labels=c("P","P.I","I","I.I","G","G.I","S","S.I"))
points(k6$modelnum,k21$CONFI,ylim=c(0,100),col=colors[4],type="l",pch="b",lty=2,lwd=2)
points(k6$modelnum,k120$CONFI,ylim=c(0,100),col=colors[8],type="l",pch="c",lty=3,lwd=2)
points(k6$modelnum,k1000$CONFI,ylim=c(0,100),col=colors[6],type="l",pch="d",lty=4,lwd=2)
#points(k6$modelnum,kALL$CONFI,ylim=c(0,100),col=colors[10],type="l",pch="e",lty=5,lwd=2)
#abline(h=100*(1/2),col="black",lty=6)

#legend("topright",legend=c("6K","21K","120K","1000K","All"),
#       col=colors[c(1,4,8,6,10)],lty=c(1,2,3,4,5),lwd=2)
legend("topright",legend=c("6K","21K","120K","1000K"),
       col=colors[c(1,4,8,6)],lty=c(1,2,3,4),lwd=2)

dev.off()

#####

k6 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/old/AGE-DEM-IBD_COMBINED-6K_ISOLATION-NOIBD.POPGENOME.STATS",
              sep="\t",header=F)
k21 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/old/AGE-DEM-IBD_COMBINED-21K_ISOLATION-NOIBD.POPGENOME.STATS",
               sep="\t",header=F)
k120 = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/combostats_8models/old/AGE-DEM-IBD_COMBINED-120K_ISOLATION-NOIBD.POPGENOME.STATS",
                sep="\t",header=F)

k6$year=6
k21$year=21
k120$year=120

forp = rbind(k6,k21,k120)
forp = sapply(forp,as.numeric)
forp = as.data.frame(forp)
forp=forp[complete.cases(forp),]


p = prcomp(forp)
d = as.data.frame(p$x)
plot(d$PC1,d$PC2,col=forp$year)
