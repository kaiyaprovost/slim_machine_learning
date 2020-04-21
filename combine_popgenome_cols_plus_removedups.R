setwd("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/")

f1="alldf_trim.temp.txt"
f2="alldf_trim.txt"
f3="alldf.txt"
f4="MERGED_empirical_1.txt"
f5="MERGED_empirical_BEL.txt"
f6="MERGED_empirical_BIL.txt"
f7="MERGED_empirical_BRU.txt"
f8="MERGED_empirical_CRI.txt"
f9="MERGED_empirical_CUR.txt"
f10="MERGED_empirical_FLA.txt"
f11="MERGED_empirical_FUS.txt"
f12="MERGED_empirical_MEL.txt"
f13="MERGED_empirical_NIT.txt"
f14="MERGED_empirical_sims.txt"
f15="MERGED_empirical_SIN.txt"
f16="merged_empirical_stats_same_cols_TRIMMED_extra.txt"
f17="merged_empirical_stats_same_cols_TRIMMED.txt copy"
f18="merged_empirical_stats_same_cols_TRIMMED.txt"
f19="merged_empirical_stats_same_cols.txt"
#f20="MERGED_empirical.txt" ## bad header
#f21="merged_window_stats_empirical.txt"
f22="MERGEDPOP_scale0.25.txt"
f23="MERGEDPOP_scale0.50.txt"
f24="MERGEDPOP_scale1.txt"
f25="MERGEDPOP_scale2.txt"
f26="MERGEDPOP_scale4.txt"
f27="MERGEDPOP_TOGETHER.txt"

df1=read.table(f1,header=T,sep="\t",stringsAsFactors = F)
df2=read.table(f2,header=T,sep="\t",stringsAsFactors = F)
df3=read.table(f3,header=T,sep="\t",stringsAsFactors = F,skip=1) ## weird?
df4=read.table(f4,header=T,stringsAsFactors = F)
df5=read.table(f5,header=T,stringsAsFactors = F)
df6=read.table(f6,header=T,stringsAsFactors = F)
df7=read.table(f7,header=T,stringsAsFactors = F)
df8=read.table(f8,header=T,stringsAsFactors = F)
df9=read.table(f9,header=T,stringsAsFactors = F)
df10=read.table(f10,header=T,stringsAsFactors = F)
df11=read.table(f11,header=T,stringsAsFactors = F)
df12=read.table(f12,header=T,stringsAsFactors = F)
df13=read.table(f13,header=T,stringsAsFactors = F)
df14=read.table(f14,header=T,stringsAsFactors = F)
df15=read.table(f15,header=T,stringsAsFactors = F)
df16=read.table(f16,header=T,sep="\t",stringsAsFactors = F)
df17=read.table(f17,header=T,sep="\t",stringsAsFactors = F)
df18=read.table(f18,header=T,sep="\t",stringsAsFactors = F)
df19=read.table(f19,header=T,sep="\t",stringsAsFactors = F)
#df20=read.table(f20,header=T,stringsAsFactors = F,sep="\t")
#df21=read.table(f21,header=T,sep="\t",stringsAsFactors = F) #file is wrong
df22=read.table(f22,header=T,stringsAsFactors = F)
df23=read.table(f23,header=T,sep="\t",stringsAsFactors = F)
df24=read.table(f24,header=T,sep="\t",stringsAsFactors = F)
df25=read.table(f25,header=T,sep="\t",stringsAsFactors = F)
df26=read.table(f26,header=T,sep="\t",stringsAsFactors = F)
df27=read.table(f27,header=T,sep="\t",stringsAsFactors = F)

colnames(df1) = toupper(colnames(df1))
colnames(df2) = toupper(colnames(df2))
colnames(df3) = toupper(colnames(df3))
colnames(df4) = toupper(colnames(df4))
colnames(df5) = toupper(colnames(df5))
colnames(df6) = toupper(colnames(df6))
colnames(df7) = toupper(colnames(df7))
colnames(df8) = toupper(colnames(df8))
colnames(df9) = toupper(colnames(df9))
colnames(df10) = toupper(colnames(df10))
colnames(df11) = toupper(colnames(df11))
colnames(df12) = toupper(colnames(df12))
colnames(df13) = toupper(colnames(df13))
colnames(df14) = toupper(colnames(df14))
colnames(df15) = toupper(colnames(df15))
colnames(df16) = toupper(colnames(df16))
colnames(df17) = toupper(colnames(df17))
colnames(df18) = toupper(colnames(df18))
colnames(df19) = toupper(colnames(df19))
#colnames(df20) = toupper(colnames(df20))
#colnames(df21) = toupper(colnames(df21))
colnames(df22) = toupper(colnames(df22))
colnames(df23) = toupper(colnames(df23))
colnames(df24) = toupper(colnames(df24))
colnames(df25) = toupper(colnames(df25))
colnames(df26) = toupper(colnames(df26))
colnames(df27) = toupper(colnames(df27))

mergedA=merge(df1,df2,all=T); mergedA=unique(mergedA)
#mergedB=merge(df3,df4,all=T); mergedB=unique(mergedB) ## choking
mergedB=merge(mergedA,df4,all=T)
mergedC=merge(df5,df6,all=T); mergedC=unique(mergedC)
mergedD=merge(df7,df8,all=T); mergedD=unique(mergedD)
mergedE=merge(df9,df10,all=T); mergedE=unique(mergedE)
mergedF=merge(df11,df12,all=T); mergedF=unique(mergedF)
mergedG=merge(df13,df14,all=T); mergedG=unique(mergedG) 
mergedH=merge(df15,df16,all=T); mergedH=unique(mergedH)
mergedI=merge(df17,df18,all=T); mergedI=unique(mergedI)
#mergedJ=merge(df19,df20,all=T); mergedJ=unique(mergedJ) 
mergedJ=merge(df19,df22,all=T); mergedJ=unique(mergedJ) 
mergedK=merge(df23,df24,all=T); mergedK=unique(mergedK)
mergedL=merge(df25,df26,all=T); mergedL=unique(mergedL)
mergedM=merge(mergedL,df27,all=T); mergedM=unique(mergedM)
mergedN = merge(mergedB,mergedC,all=T)
mergedO = merge(mergedD,mergedE,all=T)
mergedP = merge(mergedF,mergedG,all=T)
mergedQ = merge(mergedH,mergedI,all=T)
mergedR = merge(mergedJ,mergedK,all=T)
#mergedS = merge(mergedL,mergedM,all=T)
mergedS = merge(mergedR,mergedM,all=T)
mergedT = merge(mergedN,mergedO,all=T)
mergedU = merge(mergedP,mergedQ,all=T)
#mergedV = merge(mergedR,mergedS,all=T)
mergedV = merge(mergedU,mergedS,all=T)
#mergedW=merge(mergedT,mergedU,all=T)
#mergedX=merge(mergedW,mergedV,all=T)
mergedX=merge(mergedT,mergedV,all=T)

mergedX=unique(mergedX)

dim(mergedX)

#mergedY=gtools::smartbind(df3,mergedX)
mergedY=merge(mergedX,df3,all=T)

mergedY=unique(mergedY)
mergedZ=mergedY[,-198]
colnames(mergedZ)
mergedZ=unique(mergedZ)
dim(mergedZ)

write.table(mergedZ,"all_merged_8_april_2020.txt",sep="\t",row.names = F,quote=F)




tomerge = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/all_merged_8_april_2020.txt"
newdata = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MERGEDPOP_2_new_sims.txt"

df1 = read.table(tomerge,sep="\t",skip=0,header=T)
df2 = read.table(newdata,sep=" ",header=T)

x=merge(df1,df2,all=T)
x=unique(x)
write.table(x,"all_merged_8_april_2020_2.txt",sep="\t",row.names = F,quote=F)

x=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/all_merged_8_april_2020_2.txt",
             sep="\t",             stringsAsFactors = T,             header=T)

x=x[order(x$FILE,x$CHROMOSOME,x$WALL.Q.1),]

takeout=c()
for(colnum in 1:ncol(x)){
  thiscol=x[,colnum]
  print(colnum)
  nums=as.numeric(unique(thiscol))
  nona=nums[!is.na(nums)]
  nona=nona[!is.infinite(nona)]
  if(length(nona)<=1){
    takeout=c(takeout,colnum)
  }
}
print(takeout)
newx=x[,-takeout]
write.table(newx,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/all_merged_8_april_2020_2.txt",
            sep="\t",row.names = F,quote=F)
x=newx

x=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/all_merged_8_april_2020_2.txt",
             sep="\t",header=T,stringsAsFactors = F)
df1=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MERGEDPOP_2_new_sims.txt",
               sep="\t",header=T)
df2=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/MERGEDPOP_3_new_sims.txt",
               sep="\t",header=T)

newx=merge(x,df1,all=T)
x=merge(newx,df2,all=T)
newx=unique(x)
x=newx

x=x[order(x$NUC.DIVERSITY.WITHIN,x$TIMESTAMP,x$CHROMOSOME,x$WALL.Q.1),]
for(colnum in 1:ncol(x)){
  print(paste("COLNUM",colnum,"OF",ncol(x)))
  x=x[order(x[,colnum]),]
  findmatches1=c()
  for (i in 2:nrow(x)) {
    j=i-1
    if(i %% 100 == 0) {print(paste(i,"of",nrow(x)))}
    temp=t(x[c(i,j),])
    temp[temp[,1]==""]=NA
    temp[temp[,2]==""]=NA
    temp2=temp[complete.cases(temp),]
    ## how many of the cols match 
    matches=sum(temp2[,1] == temp2[,2])
    #print(paste(matches,nrow(temp2)))
    if (matches == nrow(temp2)) {
      ## they are a match
      # if(is.null(findmatches1)){
      #   findmatches1=data.frame(i=i,j=j)
      # } else {
      #   add=c(i,j)
      #   findmatches1 = rbind(findmatches1,add)
      # }
      
      ## get them matching
      #print("MATCH")
      for (r in 1:nrow(temp)) {
        row=temp[r,]
        replace=unique(row[complete.cases(row)])
        replace=paste(replace,sep=";",collapse = ";")
        temp[r,1]=replace
        temp[r,2]=replace
      }
      
      x[c(i,j),] = t(temp)
      
      ## find the ones that do not match
      
    }
  }
  dim(findmatches1)
  x=unique(x)
}
newx=unique(x)

write.table(newx,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/all_merged_8_april_2020_2.txt",
            sep="\t",row.names = F,quote=F)

x=newx



x=x[order(x$NUC.DIVERSITY.WITHIN,x$TIMESTAMP,x$CHROMOSOME,x$WALL.Q.1),]
findmatches=c()
maxk=1000
for (i in 1:nrow(x)) {
  for (j in 1:i) {
    if (i > j) {
      if(j %% 100 == 0) {print(paste(k,i,j))}
      temp=t(x[c(i,j),])
      temp=temp[complete.cases(temp),]
      ## how many of the cols match 
      matches=sum(temp[,1] == temp[,2]) ## dont name temp
      if (matches == nrow(temp)) {
        ## they are a match
        print("MATCH")
        if(is.null(findmatches)){
          findmatches=data.frame(i=i,j=j)
        } else {
          add=c(i,j)
          findmatches = rbind(findmatches,add)
        } 
        
        ## your code needs better variables
        for (r in 1:nrow(temp)) {
          row=temp[r,]
          replace=unique(row[complete.cases(row)])
          replace=paste(replace,sep=";",collapse = ";")
          temp[r,1]=replace
          temp[r,2]=replace
        }
        
        x[c(i,j),] = t(temp)
      }
      ## add progress marker 
      
    }
  }
  print(paste("DONE WITH",i,"OF",nrow(x),": UNIQUE"))
  x=unique(x)
}

dim(findmatches)
newx=unique(x)
x=newx

findmatches=c()
x=x[sample(nrow(x)),]
##intersect instead of for loop, or apply 
maxk=1000
for(k in 1:maxk){
  for (i in floor(1+((nrow(x)/maxk)*(k-1))):ceiling((nrow(x)/maxk)*k)) {
    for (j in i:ceiling((nrow(x)/maxk)*k)) {
      if (i < j) {
        if(j %% 100 == 0) {print(paste(k,i,j))}
        temp=t(x[c(i,j),])
        temp=temp[complete.cases(temp),]
        ## how many of the cols match 
        matches=sum(temp[,1] == temp[,2]) ## dont name temp
        if (matches == nrow(temp)) {
          ## they are a match
          if(is.null(findmatches)){
            findmatches=data.frame(i=i,j=j)
          } else {
            add=c(i,j)
            findmatches = rbind(findmatches,add)
          } 
            
          ## your code needs better variables
            for (r in 1:nrow(temp)) {
              row=temp[r,]
              replace=unique(row[complete.cases(row)])
              replace=paste(replace,sep=";",collapse = ";")
              temp[r,1]=replace
              temp[r,2]=replace
            }
            
            (x[c(i,j),]) = t(temp)
           ## add progress marker 
          }
        
      }
    }
    print("UNIQUE")
    x=unique(x)
  }
  print("UNIQUE")
  x=unique(x)
}
dim(findmatches)
newx=unique(x)
x=newx
