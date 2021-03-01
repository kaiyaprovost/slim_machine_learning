## note: this file assumes you have a folder called "simulated_statistics" which contains your statistics files

## load library
library(gtools)
library(R.utils)

pathlist = c(
  "~/simulated_statistics/"
)

for (path in pathlist) {
  print(path)
  setwd(path)
  files = c(list.files(pattern = "stats$", recursive = T,full.names = T),
            list.files(pattern = "STATS$", recursive = T,full.names = T))
  files = files[!(grepl("COMBO",files))]
  print(length(files))
  print("list")
  df_list <- lapply(files, FUN = function(x){
    csv = data.table::fread(file = x,sep="\t",na.strings = c("NA","NaN","NAN","-100","NULL"),fill=T)
    return(csv)
  })
  print("merge")
  merged = plyr::rbind.fill(df_list)
  merged = unique(merged)
  
  print("OUTPUTTING")
  
  write.table(merged,paste("MERGED.txt",sep=""),append=T,row.names=F)
  
  print("zipping")
  lapply(files,FUN=function(x){gzip(x,skip=T)})
  
}


