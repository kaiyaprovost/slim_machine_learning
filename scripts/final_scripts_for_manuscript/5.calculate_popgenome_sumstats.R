

dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package, ")"))))
}
packages = c("PopGenome", "moments", "R.utils", "gtools")
for (p in packages) {
  dynamic_require(p)
}

files=c()
path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/cardinalis vcf" ## 2957
setwd(path)
files = c(files,list.files(path = path,pattern = ".*split.*ms$",recursive = F,full.names = T))
#files = c(files,list.files(path = path,pattern = "*fulltemp$",recursive = TRUE,full.names = T))
#files = c(files,list.files(path = path,pattern = "*subsettemp$",recursive = TRUE,full.names = T))
#files = c(files,list.files(path = path,pattern = "*window.temp$",recursive = TRUE,full.names = T))
#files = c(files,list.files(path = path,pattern = "*vcf.temp$",recursive = TRUE,full.names = T))
#files = c(files,list.files(path = path,pattern = "*converted.temp$",recursive = TRUE,full.names = T))
#path="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/fulltemps/"
#setwd(path)
#files = c(files,list.files(path = path,pattern = "*fulltemp$",recursive = TRUE,full.names = T))
#files = c(files,list.files(path = path,pattern = "*subsettemp$",recursive = TRUE,full.names = T))
files=unique(files)
files = files[!grepl("finished", files)]
files = files[!grepl("DONE", files)]


x <- file.info(files)
x <- x[order(x$size), ]
files = rownames(x)
## linkage works with missingtrue monoeither but not missingfalse monoeither?
#myfiles = sample(files)
myfiles=files
length(myfiles)
do_vcf = F
overwrite = F
doWindow = T
verbose = T
fixHeaderMS = function(msfile) {
  lines = readLines(msfile)
  numlines = sum(lines == "//")
  segsites = as.numeric(strsplit(lines[2], ":")[[1]][2])
  
  if(is.na(segsites)){
    if(substr(lines[2],1,4)=="slim"){
      lines=lines[c(1,3:length(lines))]
      writeLines(lines, msfile)
      lines = readLines(msfile)
      header = lines[1]
      segsites = as.numeric(strsplit(lines[2], ":")[[1]][2])
    } else if (lines[2]=="//") {
      lines=lines[c(1,3:length(lines),2)]
      writeLines(lines, msfile)
      lines = readLines(msfile)
      header = lines[1]
      segsites = as.numeric(strsplit(lines[2], ":")[[1]][2])
    }
  }
  
  header = lines[1]
  if (sum(lines == paste(rep(".", segsites), collapse = "")) > 0) {
    lines = lines[lines != paste(rep(".", segsites), collapse = "")]
    writeLines(lines, msfile)
    lines = readLines(msfile)
    header = lines[1]
  }
  if (substr(header, 1, 8) == "slim 100") {
    header = "slim 20 1"
    lines[1] = header
    writeLines(lines, msfile)
    lines = readLines(msfile)
    header = lines[1]
  }
  if (substr(header, 1, 4) != "slim") {
    lines[1] = paste("slim", header, sep = " ")
    writeLines(lines, msfile)
    numreps = as.numeric(strsplit(header, " ")[[1]][2])
    lines = readLines(msfile)
    header = lines[1]
  } else {
    numreps = as.numeric(strsplit(header, " ")[[1]][3])
  }
  header = lines[1]
  if (header == "slim 1") {
    header = "slim 36 1"
    lines[1] = header
    writeLines(lines, msfile)
    lines = readLines(msfile)
    header = lines[1]
  }
  if (numlines == 1) {
    newheader = strsplit(header, " ")[[1]]
    newheader[length(newheader)] = 1
    lines[1] = paste(newheader, sep = " ", collapse = " ")
    writeLines(lines, msfile)
    numreps = 1
    lines = readLines(msfile)
    header = newheader
  }
  numinds = as.numeric(strsplit(header, " ")[[1]][2])
  expected_lines = 3 + (numinds * 2)
  non_blank_lines = sum(lines != "")
  if (is.na(expected_lines)) {
    numinds = (non_blank_lines - 3) / 2
    expected_lines = 3 + (numinds * 2)
    header = paste("slim ", numinds, " 1", sep = "")
    lines[1] = header
    writeLines(lines, msfile)
    lines = readLines(msfile)
    header = lines[1]
  }
  if (non_blank_lines != expected_lines) {
    new_numinds = (non_blank_lines - 3) / 2
    header = paste("slim", new_numinds, "1", sep = " ")
    lines[1] = header
    writeLines(lines, msfile)
    lines = readLines(msfile)
    header = lines[1]
  }
  can_run = (
    numreps >= numlines &
      lines[3] != "positions: " &
      non_blank_lines >= expected_lines &
      expected_lines >= 4 & segsites > 3)
  return(list(can_run = can_run, lines = lines))
}
runStatsMS = function(MS.class_temp=MS.class, window = F, sum = T, neut = T, link = T, recom = T, fst = T, div = T, mkt = T,
                      sweep = T, r2 = T, hap = T, min = T, nuc = div, detail = T,verbose=F) {
  skip = F
  if (sum == T) {
    MS.class_temp_outtable = get.sum.data(MS.class_temp) ## 1 row
  } else {
    MS.class_temp_outtable = NULL
  }
  if (neut == T) {
    if(verbose==T){print("neutrality stats")}
    MS.class_temp <- neutrality.stats(MS.class_temp)
    MS.class_temp_outtable = cbind(
      MS.class_temp_outtable, ## 1 line
      get.neutrality(MS.class_temp, theta = T, stats =
                       T)[[1]], get.neutrality(MS.class_temp, theta = T, stats =
                                                 T)[[2]])
  }
  if (link == T) {
    if(verbose==T) {print("linkage stats")}
    x = withTimeout(
      linkage.stats(MS.class_temp, detail = T)
      , timeout = 30, onTimeout = "silent") ## lengthy and failing
    if (is.null(x)) {
      if(verbose==T) {print("LINKAGE DID NOT WORK")}
      if(verbose==T) {print("WILL SKIP R2")}
      skip = T
    } else {
      MS.class_temp = x
      skip = F
      MS.class_temp_outtable = cbind(
        MS.class_temp_outtable, get.linkage(MS.class_temp)[[1]], get.linkage(MS.class_temp)[[2]]) ## 1 line
    }
  }
  if (recom == T) {
    if(verbose==T) {print("recomb stats")}
    x <- withTimeout(recomb.stats(MS.class_temp), timeout = 30, onTimeout = "silent") ## lengthy
    if (is.null(x)) {
      if(verbose==T) {print("RECOMB DID NOT WORK")}
    } else {
      MS.class_temp = x
      MS.class_temp_outtable = cbind(
        MS.class_temp_outtable, get.recomb(MS.class_temp)[[1]], get.recomb(MS.class_temp)[[2]]) ## 1 line
    }
  }
  if (fst == T) {
    if(verbose==T) {print("F_ST stats")}
    MS.class_temp <-
      F_ST.stats(MS.class_temp) ## not relevant with one pop?
    MS.class_temp_outtable = cbind(MS.class_temp_outtable, get.F_ST(MS.class_temp, mode = F, pairwise =
                                                                      F)) ## 1 line
  }
  if (div == T) {
    if(verbose==T) {print("diversity stats")}
    
    MS.class_temp <-
      diversity.stats(MS.class_temp, pi = T, keep.site.info = T,new.populations=F) ## doesn't work with pops specified? didn't have new.pop last time
    
    variable_names=colnames(get.diversity(MS.class_temp)[[1]])
    data=unlist(get.diversity(MS.class_temp))
    data = rbind(data)
    colnames(data) = rep(variable_names,ceiling(length(data)/length(variable_names)))[1:length(data)]
    
    MS.class_temp_outtable = cbind(MS.class_temp_outtable,data)
    
  }
  if (mkt == T) {
    if(verbose==T) {print("MKT stats")}
    MS.class_temp <-
      MKT(MS.class_temp) ## only works if multiple pops
    MS.class_temp_outtable = cbind(MS.class_temp_outtable, get.MKT(MS.class_temp))
  }
  if (detail == T) {
    if(verbose==T) {print("detail stats")}
    MS.class_temp <-
      detail.stats(MS.class_temp, biallelic.structure = T)
    MS.class_temp_outtable = cbind(
      MS.class_temp_outtable, get.detail(MS.class_temp)[[1]], get.detail(MS.class_temp)[[2]])
  }
  if (sweep == T) {
    if(verbose==T) {print("sweep stats")}
    if (detail == T) {
      MS.class_temp <- sweeps.stats(MS.class_temp)
      freq <- MS.class_temp@region.stats@minor.allele.freqs[[1]]
      freq.table = list()
      freq.table <- table(freq)
      MS.class_temp <-
        sweeps.stats(MS.class_temp, freq.table = freq.table)
    } else {
      MS.class_temp <- sweeps.stats(MS.class_temp)
    }
    CL = MS.class_temp@CL
    colnames(CL) = paste("CL", colnames(CL), sep = ".")
    CLR = MS.class_temp@CLR
    colnames(CLR) = paste("CLR", colnames(CLR), sep = ".")
    MS.class_temp_outtable = cbind(MS.class_temp_outtable, CL, CLR)
  }
  if (skip == F & r2 == T) {
    if(verbose==T) {print("calc r2")}
    MS.class_temp <- withTimeout(calc.R2(MS.class_temp), timeout = 30, onTimeout = "silent") ## lengthy
    ## I don't know how this works
    #linkage=MS.class_temp@region.stats@linkage.disequilibrium[[1]]
    #print("LINKAGE:")
    #print(linkage)
    #MS.class_temp_outtable = cbind(MS.class_temp_outtable, # linkage[[1]], # linkage[[2]])
  }
  colnames(MS.class_temp_outtable) = make.unique(colnames(MS.class_temp_outtable))
  MS.REGION.class = MS.class_temp@region.data
  MS.STATS.class = MS.class_temp@region.stats
  if (hap == T) {
    MS.haplotype.diversity = as.data.frame(cbind((MS.STATS.class@haplotype.diversity)))
    colnames(MS.haplotype.diversity) = "MS.haplotype.diversity" ## exactly identical to haplotype.diversity.within
    rownames(MS.haplotype.diversity) = rownames(MS.class_temp_outtable)
    if(verbose==T) {print("moments hap count")}
    momentsMS_haplotype.counts = t(sapply(1:length(MS.STATS.class@haplotype.counts), function(i) {
      data = unlist(MS.STATS.class@haplotype.counts[i])
      if (is.null(data)) {
        if(verbose==T) {print("failed momentsMS_haplotype.counts")}
        print(i)
        return(c(1,-100,-100,-100,-100))
      } else {
        return(all.moments(
          unlist(MS.STATS.class@haplotype.counts[i]), order.max = 4))
      }
    }))
    colnames(momentsMS_haplotype.counts) = c(
      "prob_haplotype.counts", "mean_haplotype.counts", "var_haplotype.counts", "skew_haplotype.counts", "kurt_haplotype.counts")
    rownames(momentsMS_haplotype.counts) = rownames(MS.class_temp_outtable)
    MS.class_temp_outtable = cbind(MS.class_temp_outtable, MS.haplotype.diversity, momentsMS_haplotype.counts)
  }
  if (min == T) {
    if(verbose==T) {print("moments minor")}
    momentsMS_minor.allele.freqs = t(sapply(1:length(MS.STATS.class@minor.allele.freqs), function(i) {
      data = unlist(MS.STATS.class@minor.allele.freqs[i])
      if (is.null(data)) {
        if(verbose==T) {print(i)
          print("failed momentsMS_minor.allele.freqs")}
        return(c(1,-100,-100,-100,-100))
      } else {
        return(all.moments(
          unlist(MS.STATS.class@minor.allele.freqs[i]), order.max = 4))
      }
    }))
    colnames(momentsMS_minor.allele.freqs) = c(
      "prob_minor.allele.freqs", "mean_minor.allele.freqs", "var_minor.allele.freqs", "skew_minor.allele.freqs", "kurt_minor.allele.freqs")
    rownames(momentsMS_minor.allele.freqs) = rownames(MS.class_temp_outtable)
    MS.class_temp_outtable = cbind(MS.class_temp_outtable, momentsMS_minor.allele.freqs)
  }
  if (nuc == T) {
    if(verbose==T) {print("moments nuc div within")}
    momentsMS_nuc.diversity.within = t(sapply(1:length(MS.STATS.class@nuc.diversity.within), function(i) {
      data = unlist(MS.STATS.class@nuc.diversity.within[i])
      if (is.null(data)) {
        if(verbose==T) {print("failed momentsMS_nuc.diversity.within")}
        print(i)
        return(c(1,-100,-100,-100,-100))
      } else {
        return(all.moments(
          unlist(MS.STATS.class@nuc.diversity.within[i]), order.max = 4))
      }
    }))
    colnames(momentsMS_nuc.diversity.within) = c(
      "prob_nuc.diversity.within", "mean_nuc.diversity.within", "var_nuc.diversity.within", "skew_nuc.diversity.within", "kurt_nuc.diversity.within")
    rownames(momentsMS_nuc.diversity.within) = rownames(MS.class_temp_outtable)
    MS.nuc_diversity = momentsMS_nuc.diversity.within[, 2]
    MS.div_per_site = cbind(MS.nuc_diversity / MS.class_temp@n.sites)
    colnames(MS.div_per_site) = "div_per_site"
    if (window == T) {
      MS.class_temp_outtable = cbind(MS.class_temp_outtable, momentsMS_nuc.diversity.within)
    } else {
      MS.class_temp_outtable = cbind(MS.class_temp_outtable, momentsMS_nuc.diversity.within, MS.div_per_site)
    }
  }
  if (window == F) {
    MS.class_temp_outtable[, which(is.na(MS.class_temp_outtable))] = -100
    MS.class_temp_outtable[, which(is.null(MS.class_temp_outtable))] = -100
  } else {
    MS.class_temp_outtable[is.na(MS.class_temp_outtable)] = NA
    MS.class_temp_outtable[MS.class_temp_outtable == "NULL"] = NA
    MS.class_temp_outtable = as.data.frame(MS.class_temp_outtable)
    for (i in 1:nrow(MS.class_temp_outtable)) {
      are100 = as.numeric(MS.class_temp_outtable[i, ]) == -100
      are100[is.na(are100)] = F
      if (sum(are100, na.rm = T) > 0) {
        MS.class_temp_outtable[i, ][are100] = NA
      }
    }
  }
  for (col in 1:ncol(MS.class_temp_outtable)) {
    MS.class_temp_outtable[, col] = as.character(MS.class_temp_outtable[, col])
  }
  if(verbose==T) {print("number of cols final:")
    print(ncol(MS.class_temp_outtable))}
  return(MS.class_temp_outtable)
}
momentsWindowStatsMS = function(MS.class_outtable_window) {
  moments_df = NULL
  for (colnumber in 1:ncol(MS.class_outtable_window)) {
    MS.class_outtable_window[, colnumber]
    number_filled = sum(!(is.na(
      as.numeric(MS.class_outtable_window[, colnumber]))))
    if (number_filled > 0) {
      colname = colnames(MS.class_outtable_window)[colnumber]
      data = unlist(as.numeric(MS.class_outtable_window[, colnumber]))
      moments = rbind(all.moments(data, order.max = 4, na.rm = T)[2:5])
      colnames(moments) = c("mean", "var", "skew", "kurt")
      colnames(moments) = paste(colnames(moments), paste(colname, "window", sep = "."), sep = "_")
      moments_df = cbind(moments_df, moments)
    }
  }
  return(moments_df)
}

## from popgenome package trying to play nice
readMS.custom = function (file, big.data = FALSE,verbose=FALSE) {
  if (!big.data) {
    out <- read.ms.output.custom(file.ms.output = file)
    gametes <- out$gametes
    dir.create("SwapMS",showWarnings = F)
    for (xx in 1:length(gametes)) {
      d <- gametes[[xx]]
      d <- list(matrix = d, positions = NaN)
      samplename <- paste("ms_sample_", xx, ".RD", sep = "")
      save(d, file = file.path("SwapMS", samplename))
    }
    test <- readData.custom(
      "SwapMS", SNP.DATA = TRUE, FAST = TRUE, format = "RData", big.data = big.data)
    unlink("SwapMS", recursive = TRUE)
    return(test)
  }
  if (big.data) {
    out <- read.big.ms.output(file)
    gametes <- out$gametes
    dir.create("SwapMS")
    for (xx in 1:length(gametes)) {
      open(gametes[[xx]])
      d <- gametes[[xx]][,]
      close(gametes[[xx]])
      d <- list(matrix = d, positions = NaN)
      samplename <- paste("ms_sample_", xx, ".RD", sep = "")
      save(d, file = file.path("SwapMS", samplename))
    }
    test <- readData.custom(
      "SwapMS", SNP.DATA = TRUE, FAST = TRUE, format = "RData", big.data = big.data)
    unlink("SwapMS", recursive = TRUE)
    return(test)
  }
}
read.ms.output.custom = function(txt = NA, file.ms.output = NA, MSMS = FALSE,verbose=FALSE) {
  if (!is.na(file.ms.output))
    txt <- scan(
      file = file.ms.output, what = character(0), sep = "\n", quiet = TRUE)
  if (is.na(txt[1])) {
    print("Usage: read.ms.output(txt), or read.ms.output(file=filename)")
    return()
  }
  if (MSMS[1] == FALSE) {
    nsam <- as.integer(strsplit(txt[1], split = " ")[[1]][2])
    ndraws <- as.integer(strsplit(txt[1], split = " ")[[1]][3])
  }
  #print(strsplit(txt[1], split=" "))
  h <- numeric()
  result <- list()
  gamlist <- list()
  positions <- list()
  #marker <- grep("prob",txt)
  #probs <- sapply(strsplit(txt[marker], split=":"), function(vec) as.numeric(vec[2]))
  #marker <- grep("time",txt)
  #times <- sapply(strsplit(txt[marker], split="\t"), function(vec){ as.numeric(vec[2:3])})
  times <- NaN
  probs <- NaN
  ## THE OUTPUT TEXT FOR EACH DRAW SHOULD CONTAIN THE WORD "segsites"
  marker <- grep("segsites", txt)
  if (MSMS[1] != FALSE) {
    ndraws <- length(marker)
    nsam <- MSMS$nsam
  } # MSMS
  stopifnot(length(marker) == ndraws)
  ## GET NUMBERS OF SEGREGATING SITES IN EACH DRAW
  segsites <-
    sapply(strsplit(txt[marker], split = " "), function(vec)
      as.integer(vec[2]))
  for (draw in seq(along = marker)) {
    # if(!(draw %% 100)) cat(draw, " ")
    if (segsites[draw] > 0) {
      tpos <- strsplit(txt[marker[draw] + 1], split = " ")
      positions[[draw]] <-
        as.numeric(tpos[[1]][2:(segsites[draw] + 1)])
      haplotypes <-
        txt[(marker[draw] + 2):(marker[draw] + 2 + nsam - 1)]
      haplotypes <- strsplit(haplotypes, split = "")
      h <- sapply(haplotypes, function(el)
        c(as.integer(el)))
      ## IF THERE'S 1 SEGREGATING SITE, THIS WON'T BE A MATRIX
      if (segsites[draw] == 1)
        h <- as.matrix(h)
      ## OTHERWISE, IT NEEDS TO BE TRANSPOSED
      else
        h <- t(h)
    }
    else {
      h <- matrix(nrow = nsam, ncol = 0)
      positions[[draw]] <- NA
    }
    gamlist[[draw]] <- h
    stopifnot(all(dim(h) == c(nsam, segsites[draw])))
  }
  list(
    segsites = segsites, gametes = gamlist, probs = probs, times = t(times), positions = positions, nsam = nsam, nreps = ndraws)
}
readData.custom <- function(path,populations=FALSE,outgroup=FALSE,include.unknown=FALSE,
                            gffpath=FALSE,format="fasta",parallized=FALSE,progress_bar_switch=FALSE,
                            FAST=FALSE,big.data=FALSE,SNP.DATA=FALSE,verbose=FALSE){
  
  
  ## CHECK INPUT
  if(!file.exists(path)){stop("Cannot find path !")}
  if(!file.info(path)[2]){stop("Put your file/files in a folder !")}
  
  if(gffpath[1]!=FALSE){
    if(!file.exists(gffpath)){stop("Cannot find gff path !")}
    if(!file.info(gffpath)[2]){stop("Put your file/files in a folder ! (GFF)")}
  }
  
  if(format=="HapMap"){SNP.DATA=TRUE;FAST=TRUE}
  if(format=="VCF")   {SNP.DATA=TRUE;FAST=TRUE}
  #if(format=="VCFhap"){SNP.DATA=TRUE;FAST=TRUE}
  
  if(SNP.DATA){big.data=TRUE}
  
  # Parallized version of readData
  
  if(parallized){
    
    #library(parallel)
    n.cores <- parallel::detectCores() #multicore:::detectCores()
    #n.cores <- n.cores - 1 
    #options(cores=n.cores)
    #getOption('cores')
    
    files       <- list.files(path)
    split_files <- split(files,sort(rep(1:n.cores,ceiling(length(files)/n.cores))))
    
    if(gffpath[1]!=FALSE){
      gff_files       <- list.files(gffpath)
      split_files_gff <- split(gff_files,sort(rep(1:n.cores,ceiling(length(gff_files)/n.cores))))
    }
    
    
    xxx <- NULL
    yyy <- NULL
    
    for (i in 1:n.cores){
      
      command <- paste("mkdir split",i,sep="")      
      system(command)
      filename      <- paste("split",i,sep="")
      filename_path <- file.path(getwd(),filename)
      sapply(split_files[[i]],function(x){
        command <- file.path(path,x)
        #command <- paste("mv",command,filename_path,sep=" ")
        command <- paste("cp",command,filename_path,sep=" ")
        system(command)
        
      })
      
      xxx <- c(xxx,filename_path)
      
      if(gffpath[1]!=FALSE){
        
        command <- paste("mkdir GFFRObjects_split",i,sep="")      
        system(command)
        filename      <- paste("GFFRObjects_split",i,sep="")
        filename_path <- file.path(getwd(),filename)
        sapply(split_files_gff[[i]],function(x){
          command <- file.path(gffpath,x)
          #command <- paste("mv",command,filename_path,sep=" ")
          command <- paste("cp",command,filename_path,sep=" ")
          system(command)
        })
        
        yyy <- c(yyy,filename_path)
        
      }
      
    }
    
    
    if(gffpath[1]==FALSE){
      
      #print(xxx)
      #print(format)
      #print(progress_bar_switch)
      #stop("Jo")
      
      rueck    	      <- parallel::mclapply(xxx,readData,
                                            format=format,parallized=FALSE,progress_bar_switch=FALSE,
                                            FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA,
                                            include.unknown=include.unknown, mc.cores = n.cores, mc.silent = TRUE)
      
    }else{
      
      if(verbose==T){cat("Calculation ... \n")}
      INPUT <- vector("list",n.cores)
      for(iii in 1:n.cores){
        INPUT[[iii]] <- c(xxx[iii],yyy[iii])
      }
      
      rueck    	       <- parallel::mclapply(INPUT,function(x){
        daten    <- x[1]
        gffdaten <- x[2]
        TT <- readData(path=daten,
                       gffpath=gffdaten,format=format,parallized=FALSE,progress_bar_switch=FALSE,
                       FAST=FAST,big.data=big.data,SNP.DATA=SNP.DATA,include.unknown=include.unknown)
        
        return(TT) 
      },mc.cores = n.cores, mc.silent = TRUE, mc.preschedule = TRUE)
      
    }
    
    
    genome   	       <- concatenate(rueck,n.cores)
    genome@basepath        <- file.path(path)  
    genome@project         <- file.path(path) 
    liste    	       <- list.files(path,full.names=TRUE)
    genome@genelength      <- length(liste)
    
    if(FAST){
      genome@Pop_Slide$calculated      <- TRUE # nur wegen Einschraenkungen, die auch im Slid Modus gelten
    }else{
      genome@Pop_Slide$calculated      <- FALSE
    }
    
    if(!is.list(populations)){genome@populations <- list(NULL)}else{genome@populations <- populations}
    genome@outgroup        <- outgroup
    
    ## Move splits back to original folder
    ## and Delete splits afterwards #TODO
    
    for (i in 1:n.cores){  
      #command <- paste("mv split",i,"/*"," ",path,sep="")
      #system(command)
      command <- paste("rm -r split",i,sep="")      
      system(command)
      if(gffpath[1]!=FALSE){
        #command <- paste("mv split",i,"/*"," ",gffpath,sep="")
        #system(command)
        command <- paste("rm -r GFFRObjects_split",i,sep="")      
        system(command)
      }
    }# end of delete
    if(verbose==T){cat("\n")}
    return(genome)
    
    # SNOWFALL
    
    #def.cpus <- 3
    #sfInit(parallel=TRUE,cpus=def.cpus,type="SOCK")
    #sfLibrary(ape)
    #sfExportAll() # notwendig !
    #sfSource("C:/Users/NULLEINS/Desktop/POPGEN/GENOME.R")
    #sfSource("C:/Users/NULLEINS/Desktop/POPGEN/GEN.R")
    #sfSource("C:/Users/NULLEINS/Desktop/POPGEN/DATA.R")
    
    # split the Data
    #files   <- list.files(path)
    #split_files <- split(files,1:def.cpus)
    
    #xxx <- NULL
    #for (i in 1:def.cpus){
    
    #    command <- paste("md split",i,sep="")      
    #    shell(command)
    #    filename      <- paste("split",i,sep="")
    #    filename_path <- file.path(getwd(),filename)
    #    sapply(split_files[[i]],function(x){
    #      command <- file.path(path,x)
    #      command <- paste("cp",command,filename_path,sep=" ")
    #      shell(command)
    #    })
    #    
    #    xxx <- c(xxx,filename_path) 
    #}
    
    #genome <- sfLapply(xxx,readData,parallized=FALSE)
    #sfStop()
    #return(genome)
    #genome <- sfLapply()
    
  }
  
  
  
  #### End of parallization
  
  methods <- "DATA"
  
  
  ########################
  # if (.Platform$OS.type == "unix") {
  # path_C_code      <- file.path(.path.package("PopGenome"),"libs","PopGenome.so")
  # }else{
  # ## muss je nach 32 64 Win geaendert werden !
  # path_C_code      <- file.path(.path.package("PopGenome"),"libs","PopGenome.dll")
  # }
  # dyn.load(path_C_code,PACKAGE="PopGenome")
  ########################
  
  
  npops    <- length(populations)
  popnames <- paste("pop",1:npops)
  
  # Get the alignments
  liste    <- list.files(path,full.names=TRUE)
  liste2   <- list.files(path)
  liste3   <- gsub("\\.[a-zA-Z]*","",liste2)
  
  
  # sort the list -------------
  ordered            <- as.numeric(gsub("\\D", "", liste))
  if(!any(is.na(ordered))){
    names(ordered)    <- 1:length(ordered)
    ordered           <- sort(ordered)
    ids               <- as.numeric(names(ordered))
    liste             <- liste [ids]
    liste2            <- liste2[ids]
    liste3            <- liste3[ids]  
  }
  # --------------------------- 
  
  
  gff_objects  <- vector("list",length(liste))
  SNP.GFF      <- FALSE
  
  
  GFF.BOOL        <- FALSE
  ## GET THE GFF FILES -----------------------------
  if(gffpath[1]!=FALSE){ 
    
    GFF.BOOL     <- TRUE
    gff_liste    <- list.files(gffpath,full.names=TRUE)
    gff_liste2   <- list.files(gffpath)
    gff_liste3   <- gsub("\\.[a-zA-Z]*","",gff_liste2)
    # print(gff_liste3)
    # print("-------------")
    # print(liste3)             
    #treffer      <- match(gff_liste3,liste3)
    treffer       <- match(liste3,gff_liste3)
    
    ### Print a warning when GFF files are missing 
    if(any(is.na(treffer))){
      if(verbose==T){cat("WARNING:: Could not find GFF files for:\n")
        cat(liste3[is.na(treffer)],"\n",sep=",")}
    }
    
    #print(liste3)
    #print(gff_liste3)
    
    #print("------------")
    #print(treffer)
    gff_liste    <- gff_liste[treffer]
    #print(liste)
    #print(gff_liste)
    #stop("")
    
  }
  
  
  # ---------------------------------------------
  
  
  ################# Get the poppairs names ##################################
  if(npops>1){
    if(outgroup[1]!=FALSE){
      poppairs <- choose(npops+1,2) # Outgroup is included !!
      pairs    <- combn(1:(npops+1),2)
    }else{
      poppairs <- choose(npops,2)   # Outgroup is not included !!
      pairs    <- combn(1:(npops),2)
    } 
    ###########################################################################
    
    
    #### --- Names of population pairs --- ####################################### 
    nn <- paste("pop",pairs[1,1],"/pop",pairs[2,1],sep="")
    if(dim(pairs)[2]>1){ # more than 2 Populations
      for(xx in 2:dim(pairs)[2]){
        m  <- paste("pop",pairs[1,xx],"/pop",pairs[2,xx],sep="")
        nn <- c(nn,m)
      } 
    }#END if
  } else{poppairs <- 1;nn <- "pop1"} 
  ##### ------------------------------ ######################################### 
  
  
  
  ### --- ------------------------------ #### 
  sizeliste              <- length(liste)
  genome                 <- new("GENOME")  
  genome@basepath        <- file.path(path)  
  genome@project         <- file.path(path) 
  genome@genelength      <- sizeliste
  #populationsX           <- as.matrix(populations)
  #rownames(populationsX) <- popnames
  #colnames(populationsX) <- "Number of Samples"
  if(!is.list(populations)){genome@populations <- list(NULL)}else{genome@populations <- populations}
  genome@poppairs        <- nn
  genome@outgroup        <- outgroup
  genome@region.names    <- liste2
  DATABOOL               <- is.element("DATA",methods)
  
  nsites                 <- integer(sizeliste)
  nsites2                <- integer(sizeliste)
  
  # if(ALL_METHODS){
  # genelist  <- vector("list",length(liste))      # list of objects (class GEN )
  # datalist  <- vector("list",length(liste))    # list of objects (class DATA)
  # } 
  
  # INIT
  # do this calculation every time
  region.data  <- new("region.data")
  region.stats <- new("region.stats")
  init         <- vector("list",sizeliste)
  init2        <- numeric(sizeliste)
  init3        <- rep(NaN,sizeliste)
  
  # region.data init
  populationsX      <- init
  populations2      <- init
  popmissing        <- init 
  outgroupX         <- init
  outgroup2         <- init
  
  CodingSNPS        <- init
  UTRSNPS           <- init
  IntronSNPS        <- init
  ExonSNPS          <- init
  GeneSNPS          <- init
  
  reading.frame     <- init
  rev.strand        <- init 
  Coding.matrix     <- init
  Coding.matrix2    <- init
  UTR.matrix        <- init
  Intron.matrix     <- init
  Exon.matrix       <- init
  Gene.matrix       <- init
  
  #Coding.info       <- init
  #UTR.info          <- init
  #Intron.info       <- init
  #Exon.info         <- init
  #Gene.info         <- init
  
  
  Coding.region     <- init2
  UTR.region        <- init2
  Intron.region     <- init2
  Exon.region       <- init2
  Gene.region       <- init2
  
  transitions         <- init  # matrix_sv  transition war eine 1
  biallelic.matrix    <- init  # matrix_pol
  biallelic.sites     <- init  # matrix_pos
  biallelic.sites2    <- init
  
  reference           <- init
  matrix_codonpos     <- init # codonpos. of biallelics
  synonymous          <- init # synnonsyn
  matrix_freq         <- init
  n.singletons        <- init # unic
  polyallelic.sites   <- init # mhitbp
  n.nucleotides       <- init # sum_sam
  biallelic.compositions  <- init # TCGA
  biallelic.substitutions <- init # subst
  minor.alleles       <- init # mutations
  codons              <- init
  sites.with.gaps     <- init # gaps
  sites.with.unknowns <- init
  
  # GENOME data  init
  n.valid.sites       <- init2
  n.gaps              <- init2
  n.unknowns          <- init2
  n.polyallelic.sites <- init2
  n.biallelic.sites   <- init2
  trans.transv.ratio  <- init3
  
  
  ## PROGRESS #########################
  if(progress_bar_switch){ ### because of parallized 
    progr <- progressBar.custom()
  }
  #####################################
  
  # -----------------------------#
  #if(!is.list(populations)){
  #allseq <- T}else{allseq <- F}
  # -----------------------------#
  
  for(xx in 1:sizeliste){ # 
    
    if(!progress_bar_switch){
      print(liste[xx])
    }	
    
    #print(liste[xx])
    CCC <- try(PopGenread.custom(liste[xx],format),silent=TRUE) 
    if(is.na(CCC$matrix[1])){next}
    
    gen <- CCC$matrix
    pos <- CCC$positions
    ref <- CCC$reference 
    
    # GFF stuff ----------------------------------------------------------------------
    
    gff.object.exists    <- FALSE
    FIT                  <- FALSE
    
    if(GFF.BOOL){
      
      gff.object.exists <- TRUE
      # if(length(grep("GFFRObjects",gffpath))!=0){
      
      if(SNP.DATA){
        if(!is.na(gff_liste[xx])){
          if(length(grep("GFFRObjects",gffpath))!=0){
            
            Robj                       <- load(gff_liste[xx])
            gff_object                 <- get(Robj[1])
          }else{
            gff_object                 <- gffRead(gff_liste[xx])        
          }
          
          gff_object_fit             <- fitting_gff_fast(pos,gff_object)              
          FIT                        <- TRUE
        }else{gff.object.exists <- FALSE} 
      }else{
        
        if(!is.na(gff_liste[xx])){
          gff_object                  <- gffRead(gff_liste[xx])
        }else{gff.object.exists <- FALSE}
      }   
      # Parsing GFF file
      if(gff.object.exists){
        
        gff_object                 <- parse_gff(gff_object,SNP.DATA=SNP.DATA)
        if(FIT){
          gff_object_fit             <- parse_gff(gff_object_fit,SNP.DATA=SNP.DATA)   
          GLOBAL.GFF$GFF             <- NULL                  
        }else{
          gff_object_fit             <- gff_object # in case of non SNP.DATA
        }
        
      }else{
        gff_object     <- FALSE
        gff_object_fit <- FALSE
      }
      
    }else{gff_object <- FALSE;gff_object_fit <- FALSE}
    
    # GFF stuff END -------------------------------------------------
    
    if(is.matrix(gen)){
      
      nsites[xx]  <- dim(gen)[2] 
      nsites2[xx] <- dim(gen)[2] # important for SNPDATA
      
      #------------------------------------------------------#
      result    <- popgen.custom(gen,Populations=populations,outgroup=outgroup,methods=methods,include.unknown=include.unknown,gff=gff_object_fit,FAST,SNP.DATA)
      #------------------------------------------------------#
      rm(gen) # delete gen
      
    }else{result <- NA;next}
    
    # PROGRESS #######################################################
    if(progress_bar_switch){ # wegen parallized
      progr <- progressBar.custom(xx,sizeliste, progr)
    }
    
    if(xx==ceiling(sizeliste/2)){gc()}
    ###################################################################
    
    if(is.list(result)){ # biallelic.sites exist
      
      # fill region.data
      populationsX[[xx]]        <- result$populations
      populations2[[xx]]        <- result$populations2
      popmissing[[xx]]          <- result$popmissing
      outgroupX[[xx]]           <- result$outgroup
      outgroup2[[xx]]           <- result$outgroup2
      
      datt                      <- result$data.sum
      
      CodingSNPS[[xx]]          <- datt$CodingSNPS
      UTRSNPS[[xx]]             <- datt$UTRSNPS
      IntronSNPS[[xx]]          <- datt$IntronSNPS
      ExonSNPS[[xx]]            <- datt$ExonSNPS
      GeneSNPS[[xx]]            <- datt$GeneSNPS	
      
      if(GFF.BOOL & !gff.object.exists){
        fillinit                  <- vector(,length(datt$biallelic.sites))
        CodingSNPS[[xx]]          <- fillinit
        UTRSNPS[[xx]]             <- fillinit
        IntronSNPS[[xx]]          <- fillinit
        ExonSNPS[[xx]]            <- fillinit
        GeneSNPS[[xx]]            <- fillinit
      }
      
      transitions[[xx]]         <- datt$transitions       # matrix_sv  transition war eine 1
      
      
      if(is.na(pos[1])){
        biallelic.sites[[xx]]     <- datt$biallelic.sites   # matrix_pos
        polyallelic.sites[[xx]]   <- datt$polyallelic.sites
        sites.with.gaps[[xx]]     <- datt$sites.with.gaps
        sites.with.unknowns[[xx]] <- datt$sites.with.unknowns
      }else{
        biallelic.sites2[[xx]]    <- datt$biallelic.sites
        biallelic.sites[[xx]]     <- pos[datt$biallelic.sites]
        nsites[xx]                <- biallelic.sites[[xx]][length(biallelic.sites[[xx]])] 
        polyallelic.sites[[xx]]   <- pos[datt$polyallelic.sites] # mhitbp
        sites.with.gaps[[xx]]     <- pos[datt$sites.with.gaps]
        sites.with.unknowns[[xx]] <- pos[datt$sites.with.unknowns]          
      }
      
      if(!big.data){
        
        biallelic.matrix[[xx]]           <- datt$biallelic.matrix  # matrix_pol
        colnames(biallelic.matrix[[xx]]) <- biallelic.sites[[xx]]
        
        if(GFF.BOOL & gff.object.exists){
          
          reading.frame[[xx]]              <- gff_object$reading.frame
          rev.strand[[xx]]                 <- gff_object$rev.strand
          Coding.matrix[[xx]]              <- gff_object$Coding
          UTR.matrix[[xx]]                 <- gff_object$UTR
          Intron.matrix[[xx]]              <- gff_object$Intron
          Exon.matrix[[xx]]                <- gff_object$Exon
          Gene.matrix[[xx]]                <- gff_object$Gene
          
          #
          #Coding.info[[xx]]                <- gff_object$Coding.info
          #UTR.info[[xx]]                   <- gff_object$UTR.info
          #Intron.info[[xx]]                <- gff_object$Intron.info
          #Exon.info[[xx]]                  <- gff_object$Exon.info
          #Gene.info[[xx]]                  <- gff_object$Gene.info
          
        }
        
      }else{ # BIG DATA FF
        
        biallelic.matrix[[xx]]           <- ff(datt$biallelic.matrix,dim=dim(datt$biallelic.matrix))  # matrix_pol
        close(biallelic.matrix[[xx]])
        dimnames(biallelic.matrix[[xx]]) <- list(rownames(datt$biallelic.matrix),biallelic.sites[[xx]])
        
        if(GFF.BOOL & gff.object.exists){
          
          if(dim(gff_object$Coding)[1]>0){
            
            fill                    <- as.matrix(gff_object$Coding)
            Coding.matrix[[xx]]     <- ff(fill,dim=dim(fill))
            close(Coding.matrix[[xx]])
            
            if(dim(gff_object_fit$Coding)[1]>0){
              reading.frame[[xx]]     <- gff_object_fit$reading.frame
              rev.strand[[xx]]        <- gff_object_fit$rev.strand
              fill                    <- as.matrix(gff_object_fit$Coding)
              Coding.matrix2[[xx]]    <- ff(fill,dim=dim(fill)) # wegen set.synnonsyn
              close(Coding.matrix2[[xx]])
            }
            #  Coding.info[[xx]]       <- ff(as.factor(gff_object$Coding.info))
          }
          if(dim(gff_object$UTR)[1]>0){
            fill                    <- as.matrix(gff_object$UTR)
            UTR.matrix[[xx]]        <- ff(fill,dim=dim(fill))
            close(UTR.matrix[[xx]])
            #  UTR.info[[xx]]          <- ff(as.factor(gff_object$UTR.info))
          }
          if(dim(gff_object$Intron)[1]>0){
            fill                    <- as.matrix(gff_object$Intron)   
            Intron.matrix[[xx]]     <- ff(fill,dim=dim(fill))
            close(Intron.matrix[[xx]])
            #  Intron.info[[xx]]       <- ff(as.factor(gff_object$Intron.info))
          }
          if(dim(gff_object$Exon)[1]>0){
            fill                    <- as.matrix(gff_object$Exon) 
            Exon.matrix[[xx]]       <- ff(fill,dim=dim(fill))
            close(Exon.matrix[[xx]])
            #  Exon.info[[xx]]         <- ff(as.factor(gff_object$Exon.info))
          } 
          if(dim(gff_object$Gene)[1]>0){
            fill                    <- as.matrix(gff_object$Gene)  
            Gene.matrix[[xx]]       <- ff(fill,dim=dim(fill))
            close(Gene.matrix[[xx]])
            #  Gene.info[[xx]]         <- ff(as.factor(gff_object$Gene.info))
          }    
        }
      }
      
      
      
      if(!is.na(ref[1])){
        reference[[xx]]           <- ref[datt$biallelic.sites]
      }
      
      matrix_codonpos[[xx]]     <- datt$matrix_codonpos   # codonpos. of biallelics
      synonymous[[xx]]          <- datt$synonymous        # synnonsyn
      matrix_freq[[xx]]         <- datt$matrix_freq
      n.singletons[[xx]]        <- datt$n.singletons      # unic
      
      n.nucleotides[[xx]]       <- datt$n.nucleotides     # sum_sam
      biallelic.compositions[[xx]]  <- datt$biallelic.compositions # TCGA
      biallelic.substitutions[[xx]] <- datt$biallelic.substitutions# subst
      minor.alleles[[xx]]       <- datt$minor.alleles # mutations
      codons[[xx]]              <- datt$codons
      
      
      # fill Genome data
      n.valid.sites[xx]       <- datt$n.valid.sites
      n.gaps[xx]              <- length(datt$sites.with.gaps)
      n.unknowns[xx]          <- length(datt$sites.with.unknowns)
      n.polyallelic.sites[xx] <- length(datt$polyallelic.sites)
      n.biallelic.sites[xx]   <- length(datt$biallelic.sites)
      trans.transv.ratio[xx]  <- datt$trans.transv.ratio
      
      Coding.region[xx]       <- datt$Coding_region_length
      UTR.region[xx]          <- datt$UTR_region_length
      Intron.region[xx]       <- datt$Intron_region_length
      Exon.region[xx]         <- datt$Exon_region_length
      Gene.region[xx]         <- datt$Gene_region_length    
      
    }# else{warnings("No biallelic position !")}
    
    
  }# End of For  
  
  # region.data
  region.data@populations      <- populationsX
  region.data@populations2     <- populations2
  region.data@popmissing       <- popmissing
  region.data@outgroup         <- outgroupX
  region.data@outgroup2        <- outgroup2
  
  
  region.data@CodingSNPS       <- CodingSNPS
  region.data@UTRSNPS          <- UTRSNPS
  region.data@IntronSNPS       <- IntronSNPS
  region.data@ExonSNPS         <- ExonSNPS
  region.data@GeneSNPS         <- GeneSNPS
  
  region.data@reading.frame    <- reading.frame
  region.data@rev.strand       <- rev.strand
  region.data@Coding.matrix    <- Coding.matrix
  region.data@Coding.matrix2   <- Coding.matrix2
  region.data@UTR.matrix       <- UTR.matrix
  region.data@Intron.matrix    <- Intron.matrix
  region.data@Exon.matrix      <- Exon.matrix
  region.data@Gene.matrix      <- Gene.matrix
  
  #region.data@Coding.info      <- Coding.info
  #region.data@UTR.info         <- UTR.info
  #region.data@Intron.info      <- Intron.info
  #region.data@Exon.info        <- Exon.info
  #region.data@Gene.info        <- Gene.info
  
  region.data@transitions      <- transitions  # matrix_sv  transition war eine 1
  region.data@biallelic.matrix <- biallelic.matrix# matrix_pol
  region.data@biallelic.sites  <- biallelic.sites # matrix_pos
  region.data@biallelic.sites2 <- biallelic.sites2
  region.data@reference        <- reference
  region.data@matrix_codonpos  <- matrix_codonpos # codonpos. of biallelics
  region.data@synonymous       <- synonymous# synnonsyn
  region.data@matrix_freq      <- matrix_freq
  region.data@n.singletons     <- n.singletons # unic
  region.data@polyallelic.sites <- polyallelic.sites # mhitbp
  region.data@n.nucleotides    <- n.nucleotides # sum_sam
  region.data@biallelic.compositions  <- biallelic.compositions  # TCGA
  region.data@biallelic.substitutions <- biallelic.substitutions # subst
  region.data@minor.alleles    <- minor.alleles # mutations
  region.data@codons           <- codons
  region.data@sites.with.gaps  <- sites.with.gaps # gaps
  region.data@sites.with.unknowns <- sites.with.unknowns
  
  region.stats@nucleotide.diversity   <- init
  region.stats@haplotype.diversity    <- init
  region.stats@haplotype.counts       <- init       # sfreqh
  region.stats@minor.allele.freqs     <- init       # JFD
  region.stats@biallelic.structure    <- init       # SXX
  region.stats@linkage.disequilibrium <- init
  
  # GENOME data
  genome@big.data            <- big.data
  names(nsites)              <- liste2
  genome@n.sites             <- nsites
  genome@n.sites2            <- nsites2
  genome@n.valid.sites       <- n.valid.sites    
  genome@n.gaps              <- n.gaps     
  genome@n.unknowns          <- n.unknowns   
  genome@n.polyallelic.sites <- n.polyallelic.sites
  genome@n.biallelic.sites   <- n.biallelic.sites
  genome@trans.transv.ratio  <- trans.transv.ratio
  
  genome@Coding.region     <- Coding.region
  genome@UTR.region        <- UTR.region
  genome@Intron.region     <- Intron.region
  genome@Exon.region       <- Exon.region
  genome@Gene.region       <- Gene.region
  
  
  genome@region.data  <- region.data
  genome@region.stats <- region.stats
  
  
  genome@Pop_Neutrality$calculated <- FALSE
  genome@Pop_FSTN$calculated       <- FALSE
  genome@Pop_FSTH$calculated       <- FALSE
  genome@Pop_MK$calculated         <- FALSE
  genome@Pop_Linkage$calculated    <- FALSE
  genome@Pop_Recomb$calculated     <- FALSE
  genome@Pop_Slide$calculated      <- FALSE
  genome@Pop_Detail$calculated     <- FALSE
  
  if(FAST){genome@Pop_Slide$calculated <- TRUE}
  if(GFF.BOOL){genome@gff.info<-TRUE}else{genome@gff.info <- FALSE}
  
  genome@snp.data <- SNP.DATA
  
  if(verbose==T){cat("\n")}
  return(genome)
  
}# End of Function
progressBar.custom <- function (i, total, symb_drawn = 0, symb = "=", nl = "\n",verbose=FALSE) {
  
  # i:          loop variable
  # total:      maximum value of i
  # symb_drawn: number of printet symbs
  # symb:       printed symbol
  # nl:         newline character
  
  if (missing(i) && missing(total)) {
    init = TRUE
  }
  else {
    init = FALSE
  }
  
  if (init == TRUE) {
    # print a timeline
    if(verbose==T){cat("|            :            |            :            | 100 %", nl, "|", sep = "")}
  }
  else {
    
    if (symb_drawn < 51) {
      
      progress <- round((i/total) * 52, digits = 0) - symb_drawn
      symb_drawn <- progress + symb_drawn
      
      while (progress > 0) {
        if(verbose==T){
          cat(symb);
        }
        progress <- progress - 1
      }
      
      #return(symb_drawn)
    }
    else {
      if (i == total) {
        if(verbose==T){cat("| ;-) \n")}
      }
    }
  }
  
  return(symb_drawn)
}
PopGenread.custom <- function(filepath,format) {
  
  if(format=="VCF"){
    res  <- myReadVCF(filepath)
    mat  <- res$matrix
    ref  <- res$reference
    pos  <- res$positions
    return(list(matrix=mat,reference=ref,positions=pos)) 
  }
  
  ############################################################
  #if(format=="VCF"){
  # res  <- readVCFchunk(filepath)
  # mat  <- res$matrix
  # ref  <- res$reference
  # pos  <- res$positions
  # return(list(matrix=mat,reference=ref,positions=pos)) 
  #}
  
  #if(format=="VCFhap"){
  # res  <- readVCFchunkHap(filepath)
  # mat  <- res$matrix
  # ref  <- res$reference
  # pos  <- res$positions
  # return(list(matrix=mat,reference=ref,positions=pos)) 
  #}
  
  #if(format=="VCFtri"){
  # res  <- readVCFchunk_tri(filepath)
  # mat  <- res$matrix
  # ref  <- res$reference
  # pos  <- res$positions
  # return(list(matrix=mat,reference=ref,positions=pos)) 
  #}
  
  #if(format=="VCFtet"){
  # res  <- readVCFchunk_tet(filepath)
  # mat  <- res$matrix
  # ref  <- res$reference
  # pos  <- res$positions
  # return(list(matrix=mat,reference=ref,positions=pos)) 
  #}
  ###############################################################
  
  if(format=="HapMap"){
    res  <- parse_HapMap(filepath)
    mat  <- res$matrix
    ref  <- res$reference
    pos  <- res$positions
    return(list(matrix=mat,reference=ref,positions=pos)) 
  }
  
  if(format=="RData"){
    XXX   <- load(filepath)
    mmm   <- get(XXX[1])
    return(list(matrix=mmm$matrix,reference=NaN,positions=mmm$positions))
  }
  
  
  if(format=="nexus"){
    
    matrix     <- my_read.nexus(filepath)
    nn         <- rownames(matrix)
    
    number     <- c(1,1,1,1,2,2,3,3,4,4,5,5,5,6)
    nuc        <- c("T","t","U","u","C","c","G","g","A","a","N","n","?","-")
    
    
    matrix <- apply(matrix,1,function(x){return(as.integer(number[match(x,nuc)]))})
    matrix <- t(matrix)
    
    matrix[is.na(matrix)] <- 5
    rownames(matrix)      <- nn
    attr(matrix,"path")   <- filepath
    return(list(matrix=matrix,reference=NaN,positions=NaN))
  }
  
  
  gen            <-  .Call("readdna",filepath)
  #gen            <-  .Call("my_read_fasta",filepath)
  rownames(gen)  <-  gsub(" ","",rownames(gen))
  rownames(gen)  <-  gsub("\r","",rownames(gen))
  if(is.null(gen)){gen <- NaN}  
  
  #----ape
  #  gen     <- read.dna(filepath,format,as.character=TRUE) # ape package
  #  r.names <- rownames(gen)
  #  gen     <- .Call("code_nucs",gen)
  #  rownames(gen) <- r.names
  #----ape  
  
  return(list(matrix=gen,reference=NaN,positions=NaN))
  
  #############################################################
  
  
}
popgen.custom <- function(Code_matrix,Populations=FALSE,outgroup=FALSE,methods=FALSE,include.unknown=TRUE,gff=FALSE,FAST,SNP.DATA){
  
  
  ## KONVERTING ##########################################
  # Check if an allel doesnt exists in the alignement ####
  # ----------------------------------------------### ####
  # If it doesnt exist delete it #########################
  
  ALLROWS <- FALSE
  
  if(!is.list(Populations)){                     #### wenn keine Population definiert nimm alle !
    Populations  <- list(1:dim(Code_matrix)[1])
    Populations2 <- list(rownames(Code_matrix))
    npops        <- 1
    ALLROWS      <- TRUE    
    
  }else{
    
    npops        <- length(Populations)
    Populations2 <- vector("list",npops) # beinhaltet die Namen der Sequenzen
    
    for(xx in 1:npops){
      
      if(is.character(Populations[[xx]])){
        
        namesX  <- Populations[[xx]]
        isit    <- which(is.na(match(namesX,rownames(Code_matrix))))
        
        if(length(isit)>=1){
          # Sequence doesnt exist
          if(length(attr(Code_matrix,"path"))!=0){
            warning("GEN",attr(Code_matrix,"path"))
          }
          warning("The allele ----->: ",namesX[isit], " <---- doesnt exist")
          namesX <- namesX[-isit]
        }
        
        # unbekannte geloescht ! 
        Populations[[xx]]   <- match(namesX,rownames(Code_matrix))
        Populations2[[xx]]  <- namesX
      } # End Population is character
      
      if(is.numeric(Populations[[xx]])){
        Populations2[[xx]] <- rownames(Code_matrix)[Populations[[xx]]]
      }
      
    }# End of for npops
  }# else
  
  # --------------------------------------------------------------------------------#
  
  ########################################################
  # Check if Population == NULL # Verkleiner sonst die Population
  ########################################################
  
  if(is.list(Populations)){
    res          <- delNULLpop(Populations)
    Populations  <- res$Populations
    popmissing   <- res$popmissing
    res          <- delNULLpop(Populations2)
    Populations2 <- res$Populations
    npops        <- npops - length(popmissing)
  }
  #########################################################
  
  
  
  #---------------------------------------------------------------------------------#
  # Outgroup #################### --------------------------------------------------#
  
  outgroup2 <- vector()        # outgroup Namen
  
  if(is.character(outgroup)){
    
    namesX  <- outgroup
    isit    <- which(is.na(match(namesX,rownames(Code_matrix))))
    
    if(length(isit)>=1){
      
      if(length(attr(Code_matrix,"path"))!=0){
        warning("region.stats",attr(Code_matrix,"path"))
      }
      warning("The outgroup allele ----->: ",namesX[isit], "<---- doesnt exist")
      namesX   <- namesX[-isit] # unbekannte geloescht !
      
    } 
    outgroup  <- match(namesX,rownames(Code_matrix))
    outgroup2 <- rownames(Code_matrix)[outgroup]
    
  } # end of outgroup character
  
  if(is.numeric(outgroup)){
    outgroup2 <- rownames(Code_matrix)[outgroup]
  }
  
  
  #####################################
  #### NAMES ##########################
  #####################################
  
  NAME    <- paste("pop",1:npops)
  if(outgroup[1]){
    NAME <- c(NAME,"outgroup")
  }
  
  
  #dat   <- new("DATA") # create a new class of type DATA
  
  # Save genename # ----------------------------------
  if(length(attr(Code_matrix,"path"))==0){genename <- "unknown"}else{
    genename <- attr(Code_matrix,"path")}
  # ---------------------------------------------------
  
  # GET INFORMATION FROM THE CODING MATRIX    <---------------------------------------------- take only the defined population alleles
  if(outgroup[1]==FALSE){outgr <- NULL}else{outgr <- outgroup}
  #
  if(!ALLROWS){Code_matrix <- Code_matrix[unique(c(unlist(Populations),outgr)),,drop=FALSE]}
  #
  if(dim(Code_matrix)[1]<=1){return(NA)} 
  # --------------------------------------
  
  # Aktualisiere Populationen auf geschrumpfte Code_matrix
  
  if(!ALLROWS){          # Wenn nicht alle Sequencen (beziehungweise keine Population definiert wurde)
    for(xx in 1:npops){
      Populations[[xx]] <- match(Populations2[[xx]],rownames(Code_matrix))
    } 
  }
  
  if(outgroup[1]!=FALSE){outgroup <- match(outgroup2,rownames(Code_matrix))}
  if(length(outgroup)==0){outgroup <- FALSE} # Outgroup existiert garnicht !
  ###### GET DATA #############################
  
  
  obj <- get_data(Code_matrix,include.unknown,gff=gff,FAST,SNP.DATA)
  
  
  ## Exception
  if(length(obj)==1){ # No statistics calculated
    return(NA) # no biallelic sites
  }
  ############
  
  if(length(rownames(Code_matrix)) > 0){ 
    genes       <- rownames(Code_matrix)
  }else{
    genes       <- 1:dim(Code_matrix)[1]
  }
  
  #dat@MATRIX      <- Code_matrix
  
  #if(is.list(Populations)){
  populations  <- Populations
  populations2 <- Populations2
  popmissing   <- popmissing
  #}
  
  outgroup                       <- as.matrix(outgroup)   # bezieht sich auf das komplette Alignment
  outgroup2                      <- as.matrix(outgroup2)
  rownames(obj$biallelic.matrix) <- genes
  
  #if(!is.list(Populations)){return(dat)}
  
  
  return(list(data.sum=obj,genename=genename,genes=genes,populations=populations,populations2=populations2,popmissing=popmissing,outgroup=outgroup,outgroup2=outgroup2))
  
}
delNULLpop <- function(Populations){
  
  
  npops        <- length(Populations)
  PPopulations <- list()
  popmissing   <- NULL
  
  yy <- 1
  for(xx in 1:npops){
    if(length(Populations[[xx]])>0){
      PPopulations[[yy]] <- Populations[[xx]]
      yy <- yy + 1
    }else{
      popmissing <- c(popmissing,xx)
    }
  }
  
  Populations <- PPopulations
  if(length(popmissing)==0){
    popmissing <- integer(0)
  }
  
  return(list(Populations=Populations,popmissing=popmissing))
}
get_data <- function(matr,include.unknown=FALSE,gff=FALSE,FAST,SNP.DATA){
  
  reverse.codons <- FALSE
  
  # poly <- get.polymorph(matr)
  # matr <- matr[,poly,drop=FALSE]
  # cat("Suche nach polymorphen Stellen ist fertig !")
  # GLOBAL is an evironment
  
  if(FAST){
    
    
    if(include.unknown){
      bialpos <- .Call("polyCinclude",matr)
    }else{  
      bialpos <- .Call("polyC",matr)
    }
    
    
    #bialpos <- as.logical(bialpos)
    #RETURNLISTE <- list(n_site=n_site,transitions=matrix_sv, biallelic.matrix=matrix_pol,
    #biallelic.sites=matrix_pos,matrix_codonpos=matrix_codonpos,n.singletons=unic,totalmis=NaN,s_sites=NaN,mvbp=mvbp
    #,trans.transv.ratio=(transitions/transversions),n.valid.sites=algsites,n.biallelic.sites=bial_sites,
    #polyallelic.sites=mhitbp,n.nucleotides=sum_sam,biallelic.compositions=TCGA,ROUGH=erg,matrix_freq=matrix_freq,syn=syn,
    #nonsyn=nonsyn,synonymous=!as.logical(synnonsyn),biallelic.substitutions=subst,minor.alleles=mutations,codons=Codons,
    #CodingSNPS=CodingSNPS,UTRSNPS=UTRSNPS,IntronSNPS=IntronSNPS,ExonSNPS=ExonSNPS,GeneSNPS=GeneSNPS,Coding_region_length=Coding_region_length,
    #UTR_region_length=UTR_region_length,Intron_region_length=Intron_region_length,Exon_region_length=Exon_region_length, #Gene_region_length=Gene_region_length, sites.with.gaps=gaps,sites.with.unknowns=unknown) 
    
    
    biallelic.sites    <- which(bialpos==1)
    polyallelic.sites  <- which(bialpos==4)
    
    #FIXME
    gaps            <- which(bialpos==2)
    unknowns        <- which(bialpos==3)
    
    
    if(length(biallelic.sites)==0){
      # print("No biallelic positions !")
      return(NA)
    }
    
    
    SUBMAT          <- matr[,biallelic.sites,drop=FALSE]
    
    ## very fast ---------------
    if(include.unknown){
      res          <- .Call("makeBialMatrixinclude",SUBMAT)
      Bial.Mat     <- res[[1]]
      Bial.Mat[Bial.Mat==-1] <- NaN
    }else{
      res           <- .Call("makeBialMatrix",SUBMAT)
      Bial.Mat      <- res[[1]]
    }
    
    
    
    transitions             <- as.logical(res[[2]])
    biallelic.substitutions <- res[[3]]
    rownames(biallelic.substitutions) <- c("minor","major")
    n.transitions           <- sum(transitions)
    n.transversions         <- length(biallelic.sites) - n.transitions
    tt.ratio                <- n.transitions/n.transversions 
    
    
    
    ## GFF-File
    # GFF
    #if(is.list(gff)){
    #features <- parse_gff(gff)
    #}
    
    ## GFF
    if(is.list(gff)){
      
      features <- gff
      
      if(length(features$Gene)>0){
        Gene_region        <- as.vector(unlist(apply(features$Gene,1,function(x){return(x[1]:x[2])})))
        GeneSNPS           <- is.element(biallelic.sites,Gene_region)
        Gene_region_length <- length(Gene_region)
        #IntronSNPS    <- matrix_pos[IntronSNPS]
      }else{GeneSNPS<-NaN;Gene_region_length <- 0} 
      
      if(length(features$Intron)>0){
        Intron_region <- as.vector(unlist(apply(features$Intron,1,function(x){return(x[1]:x[2])})))
        IntronSNPS    <- is.element(biallelic.sites,Intron_region)
        Intron_region_length <- length(Intron_region)
        #IntronSNPS    <- matrix_pos[IntronSNPS]
      }else{IntronSNPS<-NaN;Intron_region_length <- 0}
      
      if(length(features$UTR)>0){	
        UTR_region        <- as.vector(unlist(apply(features$UTR,1,function(x){return(x[1]:x[2])})))
        UTRSNPS           <- is.element(biallelic.sites,UTR_region)
        UTR_region_length <- length(UTR_region) 
        #UTRSNPS       <- matrix_pos[UTRSNPS]
      }else{UTRSNPS<-NaN;UTR_region_length <- 0}
      
      if(length(features$Exon)>0){
        Exon_region        <- as.vector(unlist(apply(features$Exon,1,function(x){return(x[1]:x[2])})))
        ExonSNPS           <- is.element(biallelic.sites,Exon_region)
        Exon_region_length <- length(Exon_region)
      }else{ExonSNPS<-NaN;Exon_region_length<-0}
      
      if(length(features$Coding)>0){ 
        Coding_region  <- as.vector(unlist(apply(features$Coding,1,function(x){return(x[1]:x[2])})))
        CodingSNPS     <- is.element(biallelic.sites,Coding_region)
        Coding_region_length <- length(Coding_region)
        # CodingSNPS     <- biallelic.sites[CodingSNPS]
        # size           <- length(CodingSNPS)
      }else{CodingSNPS<-NaN;Coding_region_length <- 0}
      
    }else{CodingSNPS<- NaN;IntronSNPS <- NaN; UTRSNPS <- NaN;ExonSNPS <- NaN;GeneSNPS<-NaN;
    Coding_region_length<-NaN;Intron_region_length<-NaN;UTR_region_length<-NaN;Exon_region_length<-NaN;Gene_region_length<-NaN}
    
    
    return(list(biallelic.matrix=Bial.Mat,biallelic.sites=biallelic.sites,polyallelic.sites=polyallelic.sites,
                sites.with.gaps=gaps,sites.with.unknowns=unknowns,
                transitions=transitions,n.valid.sites=NaN,synonymous=rep(NaN,length(biallelic.sites)),
                trans.transv.ratio=tt.ratio,biallelic.substitutions=biallelic.substitutions,CodingSNPS=CodingSNPS,UTRSNPS=UTRSNPS,IntronSNPS=IntronSNPS,
                ExonSNPS=ExonSNPS,GeneSNPS=GeneSNPS,Coding_region_length=Coding_region_length,
                Gene_region_length=Gene_region_length, Exon_region_length=Exon_region_length, Intron_region_length=Intron_region_length,UTR_region_length=UTR_region_length))
    ## end of very fast
    
    #############################################################################
    
    
    ######################################################################
    num.rows        <- dim(matr)[1]       
    
    
    pyrimid         <- c(1,2)
    purin           <- c(3,4)
    XXX             <- new.env()
    XXX$transitions <- vector(,length(biallelic.sites))
    XXX$count       <- 1
    ret.vek         <- integer(num.rows)
    
    Bial.Mat <- 
      apply(SUBMAT,2,function(x){
        nucs     <- unique(x)
        nuc1     <- x[1]
        nuc2     <- x[2]
        trans1   <- setequal(nucs,pyrimid) # pyrimid
        trans2   <- setequal(nucs,purin)   # purin 
        if(trans1 | trans2){XXX$transitions[XXX$count]<-TRUE}
        XXX$count <- XXX$count + 1 
        nuc1.id   <- nuc1==x
        nuc2.id   <- nuc2==x
        num.nuc1  <- sum(nuc1.id)
        num.nuc2  <- sum(nuc2.id)
        if(num.nuc1<=num.nuc2){
          ret.vek[nuc1.id] <- 1
        }else{
          ret.vek[nuc2.id] <- 1
        }
        return(ret.vek)
      })
    
    # Transitions
    transitions        <- XXX$transitions
    n.transitions      <- sum(transitions)
    n.transversions    <- length(biallelic.sites) - n.transitions
    tt.ratio           <- n.transitions/n.transversions 
    
    return(list(biallelic.matrix=Bial.Mat,biallelic.sites=biallelic.sites,
                transitions=transitions,n.valid.sites=NaN,synonymous=rep(NaN,length(biallelic.sites)),
                trans.transv.ratio=tt.ratio,Coding_region_length=NaN,Gene_region_length=NaN,Exon_region_length=NaN,
                Intron_region_length=NaN,UTR_region_length=NaN))
  }
  
  #############################################################################
  # END OF FAST
  
  # GFF
  # if(is.list(gff)){
  # features <- gff
  #}
  
  
  
  #### INIT
  n_site  <- dim(matr)[2]
  nsamtot <- dim(matr)[1]
  
  matrix_pol <- integer(nsamtot)
  
  
  
  ### Apply var
  #matrix_pos      <- FALSE
  #mhitbp          <- NaN
  #mis             <- 0
  #matrix_freq     <- NaN
  #matrix_sv       <- NA
  #algsites        <- 0
  #mvbp            <- 0
  #s_sites         <- 0
  #bial_sites      <- 0
  #mis             <- 0
  #gap             <- 0 
  #############
  
  pyrimid         <- c(1,2)
  purin           <- c(3,4)
  err             <- c(5,6)
  
  GLOBAL <- new.env()
  GLOBAL$GLOBALmatrix_pol <- NULL         ### BIALLELIC MATRIX
  RETURN                  <-  integer(9)
  
  #count <<- 0
  
  # names(RETURN)    <-  c("matrix_sv","mhitbp","mis","mvbp","algsites","bial_sites")
  
  # 1 : matrix_sv
  # 2 : mhitbp
  # 3 : mis
  # 4 : mvbp
  # 5 : algsites
  # 6 : bial_sites
  
  
  
  ### APPLY: ITERATION OVER ALL COLUMNS
  erg <- apply(matr,2,function(check){
    
    
    #cat(count,"\n")
    #count <<- count + 1
    
    fuenfsechs   <- is.element(err,check) # err = c(5,6)
    fuenf        <- fuenfsechs[1]
    sechs        <- fuenfsechs[2]
    
    # keine 5 oder 6
    if(!fuenf & !sechs){
      
      test <- unique(check)
      size <- length(test)
      
      if(size==1){
        # mono     <- 1
        return(RETURN)
      }
      
      if(size>2){
        RETURN[5] <- 1 # algsites
        RETURN[2] <- 1 # mhitbp
        return(RETURN)
      }
      
      if(size==2){
        
        RETURN[6] <- 1 # bial_sites
        
        nuc1       <- test[1]
        nuc2       <- test[2]
        nuc12      <- c(nuc1,nuc2)
        # transition/transversion
        trans1        <- setequal(nuc12,pyrimid) # pyrimid
        trans2        <- setequal(nuc12,purin)   # purin        
        ifelse(trans1 | trans2,RETURN[1] <- 1,RETURN[1] <- 0) # matrix_sv
        
        countnuc1  <- sum(check == nuc1)
        countnuc2  <- nsamtot - countnuc1  
        
        if(countnuc1 <= countnuc2){
          RETURN[7]                 <- nuc1
          matrix_pol[check==nuc1]   <- 1 
          GLOBAL$GLOBALmatrix_pol   <- cbind(GLOBAL$GLOBALmatrix_pol,matrix_pol)
          
        }else{
          RETURN[7]                 <- nuc2
          matrix_pol[check==nuc2]   <- 1 
          GLOBAL$GLOBALmatrix_pol   <- cbind(GLOBAL$GLOBALmatrix_pol,matrix_pol)
        }
        
        return(RETURN)    
      } # end size==2
      
    }# end of !fuenfsechs
    
    
    if(sechs){   
      RETURN[8]   <- 1 # gap
      RETURN[5]   <- 1 # algsites
      
      return(RETURN)
    } # return
    
    if(fuenf && include.unknown){
      
      mis        <- 1  
      gapids     <- check==5
      check2     <- check[!gapids]
      hh         <- unique(check2)
      size       <- length(hh)
      RETURN[9]  <- 1       
      
      if(size>2){
        
        RETURN[5]   <- 1 # algsites
        RETURN[2]   <- 1 # third nucleotide mismatch # mhitbp  
        
        return(RETURN)
        
      }
      
      if(size==2){ # biallelic
        
        RETURN[4]             <- 1 # mvbp
        RETURN[6]             <- 1 # bial_sites
        
        nuc1       <- hh[1]
        nuc2       <- hh[2]
        
        nuc12         <- c(nuc1,nuc2)
        # hier transition/transversion
        trans1        <- setequal(nuc12,pyrimid) # pyrimid
        trans2        <- setequal(nuc12,purin)   # purin        
        ifelse(trans1 | trans2,RETURN[1] <- 1,RETURN[1] <- 0) # matrix_sv
        # -----------------------------------------------------------
        
        countnuc1  <- sum(check2 == nuc1)
        countnuc2  <- sum(check2 == nuc2)
        
        if(countnuc1 <= countnuc2){
          RETURN[7]                 <- nuc1
          matrix_pol[check==nuc1]   <- 1 
        }else{
          RETURN[7]                 <- nuc2
          matrix_pol[check==nuc2]   <- 1 
        }
        
        matrix_pol[gapids] <- NaN # GAPS
        GLOBAL$GLOBALmatrix_pol   <- cbind(GLOBAL$GLOBALmatrix_pol,matrix_pol)
        
        # biallelic
      } # end size == 2
      
      return(RETURN)
      
    }# end 
    
    if(fuenf && !include.unknown){
      RETURN[5]          <- 1 # algsites
      RETURN[9]          <- 1
      return(RETURN)
    }
    
    return(RETURN)
    
  }) # End of APPLY
  
  #cat("APPLY is durch \n")
  
  ## Function
  #####################################
  #erg <- as.matrix(erg)
  #####################################
  
  # row.names(erg) <- c("matrix_sv","mhitbp","mis","mvbp","algsites","bial_sites")
  
  ## BIALLELIC INDICES ######################
  matrix_pos <- which(erg[6, ]==1)
  bial_sites <- length(matrix_pos)
  ###########################################
  
  #cat("bial_sites")
  
  ### ----------------- Exception ----------------- ####
  ### if there are no biallelic positions return NA ####
  ### --------------------------------------------- ####
  
  if(bial_sites==0){
    # print("No biallelic positions !")
    return(NA)
  }
  #########################################################
  
  ### gaps ################################################
  gaps      <- which(erg[8,]==1)
  #########################################################
  
  #cat("gaps")
  
  ## unknown #############################################
  unknown  <- which(erg[9,]==1)
  ########################################################
  
  #cat("unknown")
  
  ## SUM_SAM ##############################################
  algsites  <- which(erg[5,]==1)
  if(length(algsites)==0){sum_sam  <- rowSums(matr!=5) }else{sum_sam <- rowSums(matr[,-algsites,drop=FALSE]!=5)}
  #########################################################
  
  #cat("sum_sam")
  
  ### MATRIX_POL (BIALLELIC MATRIX)
  matrix_pol           <- GLOBAL$GLOBALmatrix_pol
  rm(GLOBAL) # remove environment !
  colnames(matrix_pol) <- matrix_pos
  #--------------------------------
  #################################
  
  #cat("matrix_pol")
  
  ###########################################
  ### Generate codonpositions from matrix_pos
  ###----------------------------------------
  
  #if(!is.list(gff)){
  
  if(!SNP.DATA){ 
    
    y <- 1
    matrix_codonpos <- vector(,3*bial_sites)
    
    for(xx in 1:bial_sites){
      x  <- matrix_pos[xx]
      if(x%%3==0){matrix_codonpos[y:(y+2)] <- c(x-2,x-1,x);y <- y+3;next;}
      if(x%%3==1){matrix_codonpos[y:(y+2)] <- c(x,x+1,x+2);y <- y+3;next;}
      if(x%%3==2){matrix_codonpos[y:(y+2)] <- c(x-1,x,x+1);y <- y+3;next;}
    } 
    
    # matrix_codonpos <- unique(matrix_codonpos) important change !
    
    
    if(matrix_codonpos[length(matrix_codonpos)]>dim(matr)[2]){
      ende            <- length(matrix_codonpos)-3
      matrix_codonpos <- matrix_codonpos[1:ende]
    } # Schmeisse das letzte codon raus ! Kommt nur vor wenn die Original Matrix # #  
    
  }else{matrix_codonpos <- NaN}
  
  
  #cat("matrix_codon_pos")
  ##############################################
  
  ## GFF
  if(is.list(gff)){
    
    features <- gff
    
    if(length(features$Gene)>0){
      Gene_region        <- as.vector(unlist(apply(features$Gene,1,function(x){return(x[1]:x[2])})))
      GeneSNPS           <- is.element(matrix_pos,Gene_region)
      Gene_region_length <- length(Gene_region)
      #IntronSNPS    <- matrix_pos[IntronSNPS]
    }else{GeneSNPS<-NaN;Gene_region_length <- 0} 
    
    if(length(features$Intron)>0){
      Intron_region <- as.vector(unlist(apply(features$Intron,1,function(x){return(x[1]:x[2])})))
      IntronSNPS    <- is.element(matrix_pos,Intron_region)
      Intron_region_length <- length(Intron_region)
      #IntronSNPS    <- matrix_pos[IntronSNPS]
    }else{IntronSNPS<-NaN;Intron_region_length <- 0}
    
    if(length(features$UTR)>0){	
      UTR_region        <- as.vector(unlist(apply(features$UTR,1,function(x){return(x[1]:x[2])})))
      UTRSNPS           <- is.element(matrix_pos,UTR_region)
      UTR_region_length <- length(UTR_region) 
      #UTRSNPS       <- matrix_pos[UTRSNPS]
    }else{UTRSNPS<-NaN;UTR_region_length <- 0}
    
    if(length(features$Exon)>0){
      Exon_region        <- as.vector(unlist(apply(features$Exon,1,function(x){return(x[1]:x[2])})))
      ExonSNPS           <- is.element(matrix_pos,Exon_region)
      Exon_region_length <- length(Exon_region)
    }else{ExonSNPS<-NaN;Exon_region_length<-0}
    
    if(length(features$Coding)>0){ 
      Coding_region  <- as.vector(unlist(apply(features$Coding,1,function(x){return(x[1]:x[2])})))
      start.pos      <- as.vector(unlist(apply(features$Coding,1,function(x){return(rep(x[1],x[2]-x[1]+1))})))
      end.pos        <- as.vector(unlist(apply(features$Coding,1,function(x){return(rep(x[2],x[2]-x[1]+1))})))
      # pump reverse strand #FIXME Performance
      JJJ            <- new.env()
      JJJ$count      <- 1
      open(features$rev.strand)
      rev.strand     <-  as.vector(unlist(apply(features$Coding,1,function(x){back <- rep(features$rev.strand[JJJ$count],x[2]-x[1]+1);JJJ$count <- JJJ$count+1;return(back)}))) 
      close(features$rev.strand) 
      #---------------
      
      ids            <- match(matrix_pos, Coding_region)
      ids            <- ids[!is.na(ids)]
      
      # reverse strand
      rev.strand     <- rev.strand[ids]
      
      
      start.pos      <- start.pos[ids] 
      end.pos        <- end.pos[ids]
      
      
      CodingSNPS     <- is.element(matrix_pos,Coding_region)
      Coding_region_length <- length(Coding_region)
      CodingSNPS2     <- matrix_pos[CodingSNPS]
      size            <- length(CodingSNPS2)
    }else{CodingSNPS<-NaN;size<-0;Coding_region_length <- 0}
    
    if(size>0 & !SNP.DATA){  # wenn SNPS in den codierenden regionen existieren
      y <- 1 
      matrix_codonpos <- vector(,3*size)
      reverse.codons  <- vector(,3*size)
      rvt             <- c(T,T,T)
      
      for(xx in 1:size){
        x  <- CodingSNPS2[xx]
        if(rev.strand[xx]){# reverse strand
          
          if((end.pos[xx]-x)%%3==0){reverse.codons[y:(y+2)] <- rvt;matrix_codonpos[y:(y+2)] <- c(x,x-1,x-2);y <- y+3;next;}
          if((end.pos[xx]-x)%%3==1){reverse.codons[y:(y+2)] <- rvt;matrix_codonpos[y:(y+2)] <- c(x+1,x,x-1);y <- y+3;next;}
          if((end.pos[xx]-x)%%3==2){reverse.codons[y:(y+2)] <- rvt;matrix_codonpos[y:(y+2)] <- c(x+2,x+1,x);y <- y+3;next;}
          
        }else{# non-reverse strand
          
          if((x-start.pos[xx])%%3==0){matrix_codonpos[y:(y+2)] <- c(x,x+1,x+2);y <- y+3;next;}
          if((x-start.pos[xx])%%3==1){matrix_codonpos[y:(y+2)] <- c(x-1,x,x+1);y <- y+3;next;}
          if((x-start.pos[xx])%%3==2){matrix_codonpos[y:(y+2)] <- c(x-2,x-1,x);y <- y+3;next;}
          
        }
        
      } 
      
      
      # matrix_codonpos <-  unique(matrix_codonpos) IMPORTANT ! CHECK !
      
      if(matrix_codonpos[length(matrix_codonpos)]>dim(matr)[2]){
        ende            <- length(matrix_codonpos)-3
        matrix_codonpos <- matrix_codonpos[1:ende]
        
      } # Schmeisse das letzte codon raus ! Kommt nur vor wenn die Original Matrix # #  
      
    }else{matrix_codonpos <- NaN;CodingSNPS<- NaN;Coding_region_length <- NaN}# if size >0 
  }else{CodingSNPS<- NaN;IntronSNPS <- NaN; UTRSNPS <- NaN;ExonSNPS <- NaN;GeneSNPS<-NaN;
  Coding_region_length<-NaN;Intron_region_length<-NaN;UTR_region_length<-NaN;Exon_region_length<-NaN;Gene_region_length<-NaN}
  
  # if (bial_sites!=0 & is.list(gff))
  # END GFF
  
  
  ### TCGA #################  COUNT THE NUCLEOTIDES PER ROW
  
  nucmat <- matr[,matrix_pos,drop=FALSE]
  
  ### Get Nucleotide Substitutions and Mutations #
  
  subst            <- apply(nucmat,2,function(x){return(unique(x[x!=5]))})
  colnames(subst)  <- matrix_pos
  mutations        <- erg[7,matrix_pos]
  
  #cat("subst")
  ### -----------------------------#
  
  tcga1  <- rowSums(nucmat==1)
  tcga2  <- rowSums(nucmat==2)
  tcga3  <- rowSums(nucmat==3)
  tcga4  <- rowSums(nucmat==4)
  tcga5  <- rowSums(nucmat==5)
  tcga6  <- rowSums(nucmat==6)
  
  TCGA           <- cbind(tcga1,tcga2,tcga3,tcga4,tcga5,tcga6)
  colnames(TCGA) <- c("T","C","G","A","unknown","GAP")
  ###############################
  
  #cat("TCGA")
  
  # Syn nonSyn Sites #################
  
  synnonsyn            <- rep(NaN,bial_sites)
  
  if(!SNP.DATA){
    
    if(!is.na(matrix_codonpos[1]) ){
      
      # GFF 
      testmatrix           <- matr[,matrix_codonpos,drop=FALSE]
      
      #print(testmatrix)
      
      ### change rows for reverse strands
      if(any(reverse.codons)){
        komplement         <- c(4,3,2,1,5,6)   
        testmatrix.reverse <- testmatrix[,reverse.codons,drop=FALSE]
        testmatrix.reverse <- apply(testmatrix.reverse,2,function(x){return(komplement[x])})
        testmatrix[,reverse.codons] <- testmatrix.reverse
        rm(testmatrix.reverse)
      }
      
      rm(matr)          # -------------------------------------------------------> remove original matrix
      
      
      colnames(testmatrix) <- matrix_codonpos
      synnonsynL           <- getsyn(testmatrix)
      Codons               <- synnonsynL$Codons
      syn                  <- synnonsynL$synid
      nonsyn               <- synnonsynL$nonsynid
      synid                <- match(syn,matrix_pos)
      synid                <- synid[!is.na(synid)]
      nonsynid             <- match(nonsyn,matrix_pos)
      nonsynid             <- nonsynid[!is.na(nonsynid)]
      synnonsyn[synid]     <- 0 
      synnonsyn[nonsynid]  <- 1 
      
    }else{syn  <- NaN;nonsyn <- NaN; Codons <- as.list(NaN)}
  }else{syn  <- NaN;nonsyn <- NaN; Codons <- as.list(NaN)}
  #####################################
  
  #cat("synonsyn")
  
  # ALGSITES (SITES WITHOUT VALUES >5 OR MHIT) or ==5 and include.unknown==0)
  ############################################################################
  
  algsites <- n_site - length(algsites)
  
  ############################################################################
  
  
  ##### MHIT (IF THERE IS A >2 NUCLEOTIDE)
  mhitbp <- which(erg[2,]==1)
  mhit   <- length(mhitbp)
  ########################
  #cat("mhitp")
  
  #### MIS (SUM OF ALL MHITS) # if 5 and includeunknown==1
  # totalmis <- sum(erg["mis",])
  ##############
  
  
  ## MVBP: SITE POSITION -GAP INCLUDED IN THE BIALLELIC MATRIX
  mvbp <- which(erg[4,]==1)
  ###########################
  
  #cat("mvbp")
  
  #### BIALLELIC SUBMATRIX FOR ANALYSIS
  erg2 <- erg[ ,matrix_pos,drop=FALSE]
  #####################################
  
  ## Delete erg ###########################################
  rm(erg)
  erg <- as.matrix(NaN)
  #########################################################
  
  
  #-------------------------------------#
  # VALUES FOR THE BIALLELIC MATRIX !!! #
  #-------------------------------------#
  
  # Transition/Transversions ############
  # 0: transversions
  # 1: transitions
  #######################################
  matrix_sv      <- erg2[1,]
  transitions    <- sum(matrix_sv)
  transversions  <- bial_sites - transitions
  matrix_sv      <- as.logical(matrix_sv)
  #######################################
  
  #cat("transitions")
  
  rm(erg2) ########################### delete erg 2
  
  ## matrix_freq #############################
  matrix_freq <- colSums(matrix_pol,na.rm=TRUE)
  ############################################
  
  ## UNIC ###################################
  unicids      <- which(matrix_freq==1)
  unicmatrix   <- matrix_pol[,unicids,drop=FALSE]
  unic         <- rowSums(unicmatrix,na.rm=TRUE)
  ############################################
  
  #cat("singletons")
  #
  
  #print(IntronSNPS)
  
  RETURNLISTE <- list(n_site=n_site,transitions=matrix_sv, biallelic.matrix=matrix_pol,
                      biallelic.sites=matrix_pos,matrix_codonpos=matrix_codonpos,n.singletons=unic,totalmis=NaN,s_sites=NaN,mvbp=mvbp
                      ,trans.transv.ratio=(transitions/transversions),n.valid.sites=algsites,n.biallelic.sites=bial_sites,
                      polyallelic.sites=mhitbp,n.nucleotides=sum_sam,biallelic.compositions=TCGA,ROUGH=erg,matrix_freq=matrix_freq,syn=syn,
                      nonsyn=nonsyn,synonymous=!as.logical(synnonsyn),biallelic.substitutions=subst,minor.alleles=mutations,codons=Codons,
                      CodingSNPS=CodingSNPS,UTRSNPS=UTRSNPS,IntronSNPS=IntronSNPS,ExonSNPS=ExonSNPS,GeneSNPS=GeneSNPS,Coding_region_length=Coding_region_length,
                      UTR_region_length=UTR_region_length,Intron_region_length=Intron_region_length,Exon_region_length=Exon_region_length, Gene_region_length=Gene_region_length, sites.with.gaps=gaps,sites.with.unknowns=unknown)
  
  
  return(RETURNLISTE)
  
}

j=1
if (do_vcf == F) {
  for (i in (j:length(myfiles))) {
    msfile = myfiles[i]
    
    if(file.exists(msfile)) {
      
      if(verbose==T){print(paste("Starting ", i, " of ", length(myfiles), sep = ""))}
      can_run_lines = fixHeaderMS(msfile)
      can_run = can_run_lines[[1]]
      lines = can_run_lines[[2]]
      if (can_run == TRUE) {
        if(verbose==T){print("passes can run check")}
        folder = dirname(msfile)
        file = basename(msfile)
        prefix = basename(msfile)
        outfile_text = paste(folder, "/", prefix, ".popgenome.stats", sep = "")
        setwd(folder)
        if (overwrite != T && file.exists(outfile_text)) {
          if(verbose==T){print("DO NOT RUN, ALREADY DONE")}
          #gzip(msfile, overwrite = T)
        } else {
          #MS.class = readMS.custom(file)
          if(verbose==T){print("passes overwrite check")}
          MS.class = try(readMS.custom(file),TRUE)
          
          if(isTRUE(class(MS.class)=="try-error")) { 
            if(verbose==T){print("error reading file")}
            next
            } else { 
            
            #try(readMS(file),silent=F)
            if (length(get.individuals(MS.class)) != 0) {
              if(verbose==T){print("file not empty")}
              numinds = length(get.individuals(MS.class)[[1]])
              numpopstodo = 2
              perpop = numinds / numpopstodo
              pop1 = (1:((numinds) / 2))
              pop2 = (max(pop1) + 1):numinds
              MS.class <- set.populations(MS.class, list(pop1, pop2))
              
              ## don't do diversity if less than 4 inds
              if(numinds>=4){do_div = T} else {do_div = F}
              
              #MS.class_outtable = runStatsMS(MS.class,div=do_div)
              MS.class_outtable = try(runStatsMS(MS.class,div=do_div),TRUE)
              if(isTRUE(class(MS.class_outtable)=="try-error")) { 
                if(verbose==T){print("error on stats")}
                next } else { 
                
                if(verbose==T){print("outputting 1")}
                MS.class_outtable= MS.class_outtable[ , order(colnames(MS.class_outtable))]
                MS.class_outtable$file = basename(msfile)
                write.table(
                  as.matrix(MS.class_outtable), file = outfile_text, quote = F, row.names = F, sep = "\t")
                
                if (doWindow == T) {
                  ## automatically calclulate the windows by dividing by 10
                  n.sites = get.sum.data(MS.class)[1]
                  if (n.sites > 100000) {
                    width = 100000
                    jump = 10000
                  } else if (n.sites > 10000) {
                    width = 10000
                    jump = 1000
                  } else if (n.sites > 1000) {
                    width = 1000
                    jump = 100
                  } else if (n.sites > 100) {
                    width = 100
                    jump = 10
                  } else if (n.sites > 10) {
                    width = 10
                    jump = 1
                  } else {
                    width = 2
                    jump = 1
                  }
                  if (n.sites > 1) {
                    if(verbose==T){print(paste("SITES:", n.sites, "WIDTH:", width, "JUMP:", jump))}
                    
                    MS.class.win = try(sliding.window.transform(MS.class, width = width, jump = jump, type = 2, whole.data = FALSE),TRUE)
                    if(isTRUE(class(MS.class.win)=="try-error")) { 
                      if(verbose==T){print("sliding window conversion error")}
                      next
                      } else { 
                      
                      #MS.class.win <- sliding.window.transform(MS.class, width = width, jump = jump, type = 2, whole.data = FALSE)
                      MS.class.win = set.populations(MS.class.win, list(pop1, pop2))
                      
                      MS.class_outtable_window = try(runStatsMS(MS.class.win, sum = F, window = T,div=F),TRUE)
                      if(isTRUE(class(MS.class_outtable_window)=="try-error")) { 
                        if(verbose==T){print("sliding window stats error")}
                        next } else { 
                        
                        #MS.class_outtable_window = runStatsMS(MS.class.win, sum = F, window = T,div=F)
                        
                        window_moments = try(momentsWindowStatsMS(MS.class_outtable_window),TRUE)
                        
                        if(isTRUE(class(window_moments)=="try-error")) { next } else { 
                          
                          #window_moments = momentsWindowStatsMS(MS.class_outtable_window)
                          MS.class_outtable = cbind(MS.class_outtable, window_moments)
                        }
                      }
                    }
                  }
                }
                if(verbose==T){print("outputting 2")}
                MS.class_outtable= MS.class_outtable[ , order(colnames(MS.class_outtable))]
                write.table(
                  as.matrix(MS.class_outtable), file = outfile_text, quote = F, row.names = F, sep = "\t")
                rm(MS.class_outtable)
              }
            }
            rm(MS.class)
          }
        }
        gzip(msfile, overwrite = T)
      } else {  gzip(msfile,overwrite=T)  }
      
    }
  }
}
# for(i in 1:max(j,i)){
#   if(file.exists(myfiles[i])){
#     gzip(myfiles[i],overwrite=T)
#   }
# }
# j=i+1; print(j)



