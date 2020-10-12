library(gtools)
library(R.utils)

doPopgenome = FALSE
doSumstat = TRUE
doSmallMerge = FALSE 
doMerges = FALSE
doTrimmed = FALSE

## import popgenome
if (doPopgenome==TRUE) {
  
  pathlist = c(
    #"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/BELLII/"
    "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/scutulatus/stats/" 
    #"/Users/kprovost/Dropbox (AMNH)/cardinalis vcf/stats/" 
    )
  
  
  for (path in pathlist) {
    print(path)
    setwd(path)
    files = c(list.files(pattern = "popgenome.stats$", recursive = F,full.names = T),
              list.files(pattern = ".STATS$", recursive = F,full.names = T))
    #files = list.files(pattern = "STATS.*.txt", recursive = F,full.names = T)
    #files = sample(files)
    print(length(files))
    print("list")
    df_list <- lapply(files, FUN = function(x){
      csv = data.table::fread(file = x,sep="\t",na.strings = c("NA","NaN","NAN","-100","NULL"),fill=T)
      csv$FILE = paste(rownames(csv),basename(x),sep="-")
      return(csv)
    })
    print("merge")
    merged = plyr::rbind.fill(df_list)
    merged = unique(merged)
    
    print("OUTPUTTING")
    
    write.table(merged,paste("MERGED_basics",".txt",sep=""),append=T,row.names=F)
    
    lapply(files,FUN=function(x){gzip(x,skip=T)})
    
    ## BELOW HERE IS COMMENTED OUT DUE TO BEING OLD AND CLUNKY
    # blocksize = 300
    # blocks = ceiling(length(files)/blocksize)
    # 
    # for (i in 1:blocks) {
    #   
    #   first = ((i-1)*blocksize+1)
    #   last = min((i*blocksize),length(files))
    #   print(paste("block",i,"first",first,"last",last))
    #   ptm <- proc.time()
    #   
    #   
    #   for (f in first:last) {
    #     #for (f in 1:5) {
    #     #for (f in 1:length(files)) {
    #     if (f %% 100 == 0) {
    #       print(f)
    #       #print(colnames(merged))
    #     }
    #     filename = files[f]
    #     #print(filename)
    #     #csv = read.csv(filename,sep="\t",na.strings = c("NA","NaN","NAN","-100","NULL"))
    #     csv = data.table::fread(file = filename,sep="\t",na.strings = c("NA","NaN","NAN","-100","NULL"))
    #     
    #     csv$FILE = paste(rownames(csv),basename(filename),sep="-")
    #     
    #     if (nrow(csv) > 0) {
    #       
    #       if (f==1) {
    #         merged = csv
    #       } else if (f==first) {
    #         merged = plyr::rbind.fill(merged[1,],csv)
    #       } else {
    #         merged = plyr::rbind.fill(merged,csv)
    #       }
    #       
    #     }
    #     
    #     #gzip(filename,skip=T)
    #     
    #   }
    #   print(proc.time() - ptm)
      
      
      
    #}
    
  }
  
}

## import sumstats
if (doSumstat == TRUE) {
  
  pathlist=c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/scutulatus/stats/" )
  
  
  for (path in pathlist) {
    
    print(path)
    
    setwd(path)
    files = list.files(pattern = "*sumstats.stats$", recursive = TRUE,full.names=T)
        print(length(files))
    print("list")
    df_list <- lapply(files[1:41], FUN = function(x){
      csv = data.table::fread(file = x,sep=" ",na.strings = c("NA","NaN","NAN","-100","NULL"),fill=T)
      csv$FILE = paste(rownames(csv),basename(x),sep="-")
      return(csv)
    })
    print("merge")
    merged = plyr::rbind.fill(df_list)
    merged = unique(merged)
    
    print("OUTPUTTING")
    
    write.table(merged,paste("MERGED_basics_sumstats",".txt",sep=""),append=T,row.names=F)
    
    print("zips")
    lapply(files,FUN=function(x){gzip(x,skip=T)})
    
  }
  
}
## merge popgenome and sumstats
if (doSmallMerge == TRUE) {
  
  pathlist = c("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/TESTING/")
  
  mergeacross = FALSE
  
  for (path in pathlist) {
    
    #path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/TESTING/WITHRECAP/"
    
    
    
    setwd(path)
    files = list.files(pattern = "MERGEDPOP_YYY.txt", recursive = TRUE)
    for (i in 1:length(files)) {
      file = files[i]
      print(file)
      if (i==1) {
        #merged = read.table(file,header=T)
        merged = read.csv(file,header=T,sep="\t",stringsAsFactors=F,row.names = NULL)
      } else {
        csv = read.csv(file,header=T,sep="\t",stringsAsFactors=F,row.names = NULL)
        merged = smartbind(merged,csv)
      }
      print("UNIQ")
      merged = unique(merged)
      #print("AGG")
      #merged = aggregate(merged,by=list(merged$RUN),FUN=c)
      #merged = aggregate(merged,by=list(merged$RUN),
      #                   FUN=function(x){
      #                     c(unique(c(x)[!(is.na(c(x)))]),"")
      #})
    }
    write.table(merged,"MERGEDPOP_QQQ.txt",row.names = F)
    
    #xx = aggregate(merged,by=list(merged$RUN),FUN=function(x) {unique(x)})
    #write.table(xx,"MERGEDPOP_2A.txt",row.names = F)
    
    if (mergeacross == T) {
      
      
      start = Sys.time()
      for (i in c(22658)) {
        print(i)
        
        for (j in 8:ncol(merged)) {
          
          
          merged[i,j] = paste((unique(unlist(merged[i,j]))),collapse=", ")
          
        }
      }
      
      for (i in 1:ncol(merged)) {
        
        
        print(i)
        merged[,i] = unlist(merged[,i])
      }
      
      write.table(merged,"MERGEDPOP_2BBB.txt",row.names = F)
    }
    
    # path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/POPGENOME/"
    # setwd(path)
    # files = list.files(pattern = "MERGEDPOP_2.txt", recursive = TRUE)
    # for (i in 1:length(files)) {
    #   file = files[i]
    #   print(file)
    #   if (i==1) {
    #     merged = read.table(file,header=T)
    #   } else {
    #     csv = read.table(file,header=T)
    #     merged = smartbind(merged,csv)
    #   }
    # }
    # write.table(merged,"MERGEDSUMSTAT_3.txt",row.names = F)
    # 
    # 
    # 
    # path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/SUMSTAT/120000/"
    # setwd(path)
    # files = list.files(pattern = "MERGEDSUMSTAT.txt", recursive = TRUE)
    # for (i in 1:length(files)) {
    #   file = files[i]
    #   print(file)
    #   if (i==1) {
    #     merged = read.table(file,header=T)
    #   } else {
    #     csv = read.table(file,header=T)
    #     merged = smartbind(merged,csv)
    #   }
    # }
    # write.table(merged,"MERGEDSUMSTAT_2.txt",row.names = F)
    # 
    # path = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/SUMSTAT/"
    # setwd(path)
    # files = list.files(pattern = "MERGEDSUMSTAT_2.txt", recursive = TRUE)
    # for (i in 1:length(files)) {
    #   file = files[i]
    #   print(file)
    #   if (i==1) {
    #     merged = read.table(file,header=T)
    #   } else {
    #     csv = read.table(file,header=T)
    #     merged = smartbind(merged,csv)
    #   }
    # }
    # write.table(merged,"MERGEDSUMSTAT_3.txt",row.names = F)
    # 
  }
}
if (doMerges == TRUE) {
  
  print("reading!")
  
  path1 = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/SUMSTAT/"
  setwd(path1)
  files1 = list.files(pattern = "MERGEDSUMSTAT_3.txt", recursive = TRUE,full.names = F)
  
  path2 = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/POPGENOME/"
  setwd(path2)
  files2 = list.files(pattern = "MERGEDPOP_3.txt", recursive = TRUE,full.names = F)
  
  end1 = files1[1]
  end2 = files2[1]
  
  print("uploading")
  csv1 = unique(read.csv(paste(path1,end1,sep=""),sep=" "))
  csv2 = unique(read.csv(paste(path2,end2,sep=""),sep=" "))
  
  ## make sure for csv2 that you get rid of the ones with doubles
  
  print("merging")
  csv3 = merge(csv1,csv2)
  
  print("writing")
  write.table(csv3,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/MERGEDCOMBO_3.txt",append=T,row.names=F)
  
}
if (doTrimmed==TRUE) {

print("trimming")
#trimmed = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/MERGEDPOP_TOGETHER.txt",
#                   sep="\t")
trimmed = read.csv("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/TESTING/RECAPSTATS/MERGEDPOP_FINAL_extratrimme_notcollinear.txt",
                   sep="\t")
nrow(trimmed)

## need to convert MS.haplotype.diversity 
#b4 = trimmed$MS.haplotype.diversity


#corrm = cor(trimmed[,c(6:41,43:55)],use="pairwise.complete.obs")
#corrm[is.na(corrm)] = 0

#bin = corrm
#bin[corrm >=  0.8] = 1
#bin[corrm <   0.8] = 0
#bin[corrm <= -0.8] = -1

#corrplot::corrplot(bin,method="color",order="hclust",tl.cex=0.5)

#trimmed = trimmed[,order(colnames(trimmed))]
#names(trimmed)

#corrd = cor(trimmed[,c("n.biallelic.sites","n.segregating.sites","ss.")],
#            use="pairwise.complete.obs")
#corrplot::corrplot(corrd,method="number")

trimmed = na.omit(trimmed)
nrow(trimmed)

ibdyes = trimmed[trimmed$ibd==TRUE,]
ibdno = trimmed[trimmed$ibd==FALSE,]

panyes = ibdyes[ibdyes$demog=="panmixia",]
panno = ibdno[ibdno$demog=="panmixia",]
isoyes = ibdyes[ibdyes$demog=="isolation",]
isono = ibdno[ibdno$demog=="isolation",]
gfyes = ibdyes[ibdyes$demog=="geneflow",]
gfno = ibdno[ibdno$demog=="geneflow",]
secyes = ibdyes[ibdyes$demog=="seccon",]
secno = ibdno[ibdno$demog=="seccon",]

panyes.1000 = panyes[panyes$year=="1000k",]
panno.1000 = panno[panno$year=="1000k",]
isoyes.1000 = isoyes[isoyes$year=="1000k",]
isono.1000 = isono[isono$year=="1000k",]
gfyes.1000 = gfyes[gfyes$year=="1000k",]
gfno.1000 = gfno[gfno$year=="1000k",]
secyes.1000 = secyes[secyes$year=="1000k",]
secno.1000 = secno[secno$year=="1000k",]

panyes.120 = panyes[panyes$year=="120k",]
panno.120 = panno[panno$year=="120k",]
isoyes.120 = isoyes[isoyes$year=="120k",]
isono.120 = isono[isono$year=="120k",]
gfyes.120 = gfyes[gfyes$year=="120k",]
gfno.120 = gfno[gfno$year=="120k",]
secyes.120 = secyes[secyes$year=="120k",]
secno.120 = secno[secno$year=="120k",]

panyes.21 = panyes[panyes$year=="21k",]
panno.21 = panno[panno$year=="21k",]
isoyes.21 = isoyes[isoyes$year=="21k",]
isono.21 = isono[isono$year=="21k",]
gfyes.21 = gfyes[gfyes$year=="21k",]
gfno.21 = gfno[gfno$year=="21k",]
secyes.21 = secyes[secyes$year=="21k",]
secno.21 = secno[secno$year=="21k",]

panyes.6 = panyes[panyes$year=="6k",]
panno.6 = panno[panno$year=="6k",]
isoyes.6 = isoyes[isoyes$year=="6k",]
isono.6 = isono[isono$year=="6k",]
gfyes.6 = gfyes[gfyes$year=="6k",]
gfno.6 = gfno[gfno$year=="6k",]
secyes.6 = secyes[secyes$year=="6k",]
secno.6 = secno[secno$year=="6k",]


setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS_testing/TESTING/")

write.table(panyes.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-YES-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(panno.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-NO-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isoyes.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-YES-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isono.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-NO-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfyes.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-YES-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfno.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-NO-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secyes.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-YES-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secno.1000[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-NO-1000K.COMBO.STATS",col.names=F,row.names=F,append=T)


write.table(panyes.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-YES-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(panno.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-NO-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isoyes.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-YES-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isono.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-NO-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfyes.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-YES-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfno.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-NO-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secyes.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-YES-120K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secno.120[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-NO-120K.COMBO.STATS",col.names=F,row.names=F,append=T)

write.table(panyes.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-YES-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(panno.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-NO-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isoyes.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-YES-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isono.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-NO-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfyes.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-YES-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfno.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-NO-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secyes.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-YES-21K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secno.21[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-NO-21K.COMBO.STATS",col.names=F,row.names=F,append=T)

write.table(panyes.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-YES-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(panno.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_PAN-NO-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isoyes.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-YES-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(isono.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_ISO-NO-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfyes.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-YES-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(gfno.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_GF-NO-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secyes.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-YES-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
write.table(secno.6[,c(-1:-5)],"DEM-IBD-AGE_COMBINED_SEC-NO-6K.COMBO.STATS",col.names=F,row.names=F,append=T)
}
