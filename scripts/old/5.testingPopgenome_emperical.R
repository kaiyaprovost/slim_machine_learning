dynamic_require <- function(package,lib=NULL) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  if(is.null(lib)) {
    install.packages(package,repos='http://cran.us.r-project.org')
  } else {
    install.packages(package,lib=lib,repos='http://cran.us.r-project.org')
  }
  return(eval(parse(text = paste(
    "require(", package,  ")"))))
}


packages = c("PopGenome", "moments", "R.utils")

install.packages("moments")

for (p in packages) {
  dynamic_require(p,lib=.libPaths()[2])
}

#library(PopGenome)
#library(moments)
#library(R.utils)
#options(expressions = 10000)

## TODO:
## get this so that it can input either a VCF or a MS file and then choose which to output?
## for now, just do the MS files

## make it so that you can give it a list of the txt files to output

path="/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"



# locs = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/all_done/NE1000/LOCS/SUBSET/model3_isolation_6k-1558189269-53.subsetlocs"
# subset = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/all_done/NE1000/TEMPS/SUBSET/DONE/model3_isolation_6k-1558189269-53.withheader.subsettemp"
# locs1 = read.csv(locs,sep=" ")
#
# while(is.na(locs1[1,3])) {
#
# locs2 = readLines(locs)
# locs3 = paste(locs2,collapse = "\n",sep="")
# locsgrep = gsub("0\\.", " 0.",locs3)
# locsgrep = gsub("\n ", "\n",locsgrep)
# writeLines(locsgrep,locs,sep="")
# locs1 = read.csv(locs,sep=" ")
#
# }
#
# msfile = subset

## change the path to other folders and include the .txt files 1:4

setwd(path)

#files = list.files(pattern = "Amphispiza-bilineata-called.geno.NW_005087129.1.fulltemp", recursive = TRUE)
#files = list.files(pattern = "subsettemp$", recursive = FALSE)

files = list.files(pattern=".ms$",recursive=T)

myfiles = files

#myfiles = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TEMPS/SUBSET/model1_panmixia_6k-1558965538-1.withheader.subsettemp"

#### 21k = 962:1075 ## 1174, 416, 879, 910
#### 120 pan = 1045:1450 ##

length(myfiles)
#myfiles = files[14055:length(files)]

#myfiles = sample(myfiles)


## make this a for loop

do_vcf = F
overwrite = F

## /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/SINUATUS/WINDOWS/10
## Cardinalis-sinuatus-called.geno.PseudoNC_011474.1_Tgut_10.fixedchroms.converted.vcf_w100000_o20000_0.window.ms

if (do_vcf == F) {
  ##### with the MS data itself
  
  #msfolder  = "/Users/kprovost/Dropbox (AMNH)/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/pogenometest/ms/"
  msfolder = path
  
  ## for this to read, your header has to say "slim 20 100" or similar
  #msfile = "/Users/kprovost/Dropbox (AMNH)/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/pogenometest/ms/stable-neut-nospat-1pop-isol-1542226966.txt"
  #i = 1
  
  for (i in 1:length(myfiles)) {
    msfile = myfiles[i]
    print(paste("Starting ", i, " of ", length(myfiles), sep = ""))
    print(msfile)
    
    
    ## check if the file is temp or text
    
    ## check the header
    lines = readLines(msfile)
    
    numlines = 1
    
    header = lines[1]
    
    if (header == "//") {
      totinds = ((length(lines) - 3) / 2)
      header = paste("callgeno", totinds, 1,"\n//", sep = " ")
      lines[1] = header
      writeLines(lines, msfile)
      lines = readLines(msfile)
      header = lines[1]
    }
    
    split = strsplit(msfile, split = "/")[[1]]
    folder = paste(split[1:length(split) - 1], sep = "/", collapse = "/")
    file = split[length(split)]
    prefix = strsplit(file, split = "[.]")[[1]][1]
    
    if (overwrite != T &&
        file.exists(paste(path, folder, "/", prefix, ".popgenome.stats", sep = ""))) {
      print("DO NOT RUN, ALREADY DONE")
      #tomove = (paste("/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/DONE/",msfile,sep=""))
      #currently = (paste("/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/",msfile,sep=""))
      #print("moving")
      #file.rename(currently, tomove)
    } else {
      ## read in MS file
      MS.class = readMS(msfile,big.data=T)
      
      
      ## assumes you have read a locs file?
      ## honestly you don't need a locs file you can just assign them
      
      ## you need a populations file here though 
      
      
      numinds = length(get.individuals(MS.class)[[1]])
      numpopstodo = 2
      perpop = numinds / numpopstodo
      pop1 = (1:((numinds) / 2))
      pop2 = (((numinds / 2) + 1):numinds)
      
      MS.class <- set.populations(MS.class, list(pop1, pop2))
      
      
      
      #show.slots(MS.class)
      
      ## generate stats for MS
      
      print("neutrality stats")
      MS.class <- neutrality.stats(MS.class)
      print("linkage stats")
      
      x = withTimeout(
        linkage.stats(MS.class, detail = T)
        ,
        timeout = 30,
        onTimeout = "silent"
      ) ## lengthy
      if (is.null(x)) {
        print("LINKAGE DID NOT WORK")
        print("WILL SKIP R2")
        skip = T
      } else {
        MS.class = x
        skip = F
      }
      
      
      print("recomb stats")
      x <- withTimeout(recomb.stats(MS.class),
                       timeout = 30,
                       onTimeout = "silent") ## lengthy
      
      if (is.null(x)) {
        print("RECOMB DID NOT WORK")
      } else {
        MS.class = x
      }
      
      print("F_ST stats")
      
      MS.class <-
        F_ST.stats(MS.class) ## not relevant with one pop?
      #MS.class <- F_ST.stats.2(MS.class,snn=TRUE,Phi_ST=FALSE)
      print("diversity stats")
      MS.class <-
        diversity.stats(MS.class, pi = T, keep.site.info = T) ## doesn't work with pops specified?
      print("sweep stats")
      MS.class <- sweeps.stats(MS.class)
      print("MKT stats")
      MS.class <- MKT(MS.class) ## only works if multiple pops
      print("detail stats")
      MS.class <- detail.stats(MS.class, biallelic.structure = T)
      
      if (skip == F) {
        print("calc r2")
        MS.class <- withTimeout(calc.R2(MS.class),
                                timeout = 30,
                                onTimeout = "silent") ## lengthy
        
      } else {
        print("SKIPPING CALC R2")
      }
      
      ## get stats for MS
      # print(head(get.sum.data(MS.class)))
      # print(head(get.neutrality(MS.class)[[1]]))
      # print(head(get.linkage(MS.class)[[1]]))
      # print(head(get.recomb(MS.class)[[1]]))
      # head(get.F_ST(MS.class)) ## never relevant with one pop
      # print(head(get.diversity(MS.class)[[1]]))
      # head(get.sweeps(MS.class)) ## not working
      # head(get.MKT(MS.class))  ## only works if mult pops
      # head(get.detail(MS.class,biallelic.structure = T))
      
      ## THIS WORKS
      
      if (skip == F) {
        MS.class_outtable = cbind(
          get.sum.data(MS.class),
          get.neutrality(MS.class, theta = T, stats = T)[[1]],
          get.neutrality(MS.class, theta = T, stats = T)[[2]],
          get.F_ST(MS.class, mode = F, pairwise = F),
          get.MKT(MS.class),
          ## think only works if have actual ACTG values
          get.linkage(MS.class)[[1]],
          ##
          get.linkage(MS.class)[[2]],
          ##
          get.recomb(MS.class)[[1]],
          ##
          get.recomb(MS.class)[[2]],
          ##
          get.diversity(MS.class)[[1]],
          get.diversity(MS.class)[[2]]
        )
        
        colnames(MS.class_outtable) = make.unique(colnames(MS.class_outtable))
        
      } else {
        MS.class_outtable = cbind(
          get.sum.data(MS.class),
          get.neutrality(MS.class, theta = T, stats = T)[[1]],
          get.neutrality(MS.class, theta = T, stats = T)[[2]],
          c(
            "Wall.B" = NA,
            "Wall.Q" = NA,
            "Rozas.ZA" = NA,
            "Rozas.ZZ" = NA,
            "Kelly.Z_ns" = NA
          ),
          c("Hudson.Kaplan.RM" = NA),
          get.diversity(MS.class)[[1]]
          
          
        )
        colnames(MS.class_outtable) = make.unique(colnames(MS.class_outtable))
        
      }
      
      print("num cols A:")
      print(ncol(MS.class_outtable))
      
      #plot(MS.class@nuc.diversity.within)
      
      ## don't need to do sliding windows here
      
      
      ## get region.data as own class
      
      MS.REGION.class = MS.class@region.data
      
      #MS.REGION.class@populations ## list of individuals by pop
      #MS.REGION.class@populations2 ## list of individuals by pop?
      #MS.REGION.class@outgroup ## whether there is outgroup or not
      #MS.REGION.class@transitions ## transitions/transversions, not applicable
      #MS.REGION.class@biallelic.matrix[[1]] ## ?
      #MS.REGION.class@n.singletons ## blank
      #MS.REGION.class@biallelic.sites ## large
      #MS.REGION.class@reference ## blank
      #MS.REGION.class@n.nucleotides # #blank
      #MS.REGION.class@biallelic.compositions ## blank
      #MS.REGION.class@synonymous ## all NaN
      #MS.REGION.class@biallelic.substitutions ## large and probably not applicable
      #MS.REGION.class@polyallelic.sites ## blank
      #MS.REGION.class@sites.with.gaps  ## blank
      #MS.REGION.class@sites.with.unknowns ## blank
      #MS.REGION.class@minor.alleles ## blank
      #MS.REGION.class@codons ## blank
      #MS.REGION.class@IntronSNPS ## blank
      #MS.REGION.class@UTRSNPS ## blank
      #MS.REGION.class@CodingSNPS ## blank
      #MS.REGION.class@ExonSNPS ## blank
      #MS.REGION.class@GeneSNPS ## blank
      #MS.REGION.class@included ## blank
      
      #####
      
      #####
      #
      # msfile = "/Users/kprovost/Dropbox (AMNH)/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/slim-stable-neut-nospat-1pop-isol-1542226966/stable-neut-nospat-1pop-isol-1542226966-1.withheader.temp"
      #
      #
      #
      # #ms = readMS("/Users/kprovost/Dropbox (AMNH)/Classes/Machine Learning/SLiMTreeSeqPub-master/p1_ms_exponential.txt")
      #
      #
      #
      #
      # ## to get this to work you need header to read "slim number number"
      #
      #
      #
      #
      #
      #
      #
      #
      #
      # ## this code from pipemaster
      # ms = neutrality.stats(ms,FAST=T)
      # ms = diversity.stats(ms)
      # ms = F_ST.stats(ms,FAST=T)
      #
      # s.sites<-ms@n.segregating.sites
      # div.within = ms@nuc.diversity.within
      # #pi.witiin<-ms@nuc.diversity.within/as.numeric(model$loci[u,2])
      # Hap.div<-ms@hap.diversity.within
      # Taj.D<-ms@Tajima.D
      # Fu.Li.D<-ms@Fu.Li.D
      # Fu.Li.F<-ms@Fu.Li.F
      # Hap.Fst<-ms@haplotype.F_ST
      # nuc.Fst<-ms@nucleotide.F_ST
      # #OA.ms<-cbind(s.sites,pi.within,Hap.div,Taj.D,Fu.Li.D,Fu.Li.F)
      # OA.ms<-cbind(s.sites,div.within,Hap.div,Taj.D,Fu.Li.D,Fu.Li.F,Hap.Fst,nuc.Fst)
      #
      # SS<-NULL
      # kur<-NULL
      # vari<-NULL
      # skew<-NULL
      # for(jj in 1:sim.block.size){
      #   x<-NULL
      #   for(ii in 1:nrow(model$loci)){
      #     x<-rbind(x,OA.ss[[ii]][jj,])
      #   }
      #   SS<-rbind(SS,colMeans(x, na.rm=T))
      #   if(get.moments==T){
      #     vari<-rbind(vari,diag(var(x, na.rm=T)))
      #     kk<-NULL
      #     sk<-NULL
      #     for(uu in 1:ncol(x)){
      #       s<-skewness(x[,uu],na.rm=T)
      #       sk<-c(sk,s)
      #       k<-kurtosis(x[,uu],na.rm=T)
      #       kk<-c(kk,k)
      #     }
      #     kur<-rbind(kur,kk)
      #     skew<-rbind(skew,sk)
      
      
      
      ## get region.stats as own class
      
      MS.STATS.class = MS.class@region.stats
      
      #MS.STATS.class@nucleotide.diversity[[1]] ## blank
      MS.haplotype.diversity = as.data.frame(cbind((
        MS.STATS.class@haplotype.diversity
      )))
      colnames(MS.haplotype.diversity) = "MS.haplotype.diversity" ## exactly identical to haplotype.diversity.within
      rownames(MS.haplotype.diversity) = rownames(MS.class_outtable)
      
      #print(head(MS.STATS.class@haplotype.counts)) ## large
      ## probability, mean, variance, skewness, kurtosis
      print("moments hap count")
      momentsMS_haplotype.counts = t(sapply(1:length(MS.STATS.class@haplotype.counts), function(i) {
        #print(i)
        data = unlist(MS.STATS.class@haplotype.counts[i])
        
        if (is.null(data)) {
          print("failed momentsMS_haplotype.counts")
          print(i)
          return(c(1,-100,-100,-100,-100))
          
        } else {
          return(all.moments(
            unlist(MS.STATS.class@haplotype.counts[i]),
            order.max = 4
          ))
          
        }
        
      }))
      colnames(momentsMS_haplotype.counts) = c(
        "prob_haplotype.counts",
        "mean_haplotype.counts",
        "var_haplotype.counts",
        "skew_haplotype.counts",
        "kurt_haplotype.counts"
      )
      rownames(momentsMS_haplotype.counts) = rownames(MS.class_outtable)
      
      
      #print(head(MS.STATS.class@minor.allele.freqs)) ## large
      print("moments minor")
      momentsMS_minor.allele.freqs = t(sapply(1:length(MS.STATS.class@minor.allele.freqs), function(i) {
        data = unlist(MS.STATS.class@minor.allele.freqs[i])
        
        if (is.null(data)) {
          print(i)
          
          print("failed momentsMS_minor.allele.freqs")
          return(c(1,-100,-100,-100,-100))
          
        } else {
          return(all.moments(
            unlist(MS.STATS.class@minor.allele.freqs[i]),
            order.max = 4
          ))
          
        }
      }))
      colnames(momentsMS_minor.allele.freqs) = c(
        "prob_minor.allele.freqs",
        "mean_minor.allele.freqs",
        "var_minor.allele.freqs",
        "skew_minor.allele.freqs",
        "kurt_minor.allele.freqs"
      )
      rownames(momentsMS_minor.allele.freqs) = rownames(MS.class_outtable)
      
      
      #print(head(MS.STATS.class@linkage.disequilibrium))  ## HUGE amount of stuff - not sure how to get
      #MS.STATS.class@biallelic.structure[[1]] ## blank
      #MS.STATS.class@site.FST ## blank
      #MS.STATS.class@missing.freqs # #blank
      #MS.STATS.class@nuc.diversity.between ## blank
      #print(head(MS.STATS.class@nuc.diversity.within)) ## large
      print("moments nuc div within")
      momentsMS_nuc.diversity.within = t(sapply(1:length(MS.STATS.class@nuc.diversity.within), function(i) {
        data = unlist(MS.STATS.class@nuc.diversity.within[i])
        
        if (is.null(data)) {
          print("failed momentsMS_nuc.diversity.within")
          print(i)
          return(c(1,-100,-100,-100,-100))
          
        } else {
          return(all.moments(
            unlist(MS.STATS.class@nuc.diversity.within[i]),
            order.max = 4
          ))
        }
      }))
      colnames(momentsMS_nuc.diversity.within) = c(
        "prob_nuc.diversity.within",
        "mean_nuc.diversity.within",
        "var_nuc.diversity.within",
        "skew_nuc.diversity.within",
        "kurt_nuc.diversity.within"
      )
      rownames(momentsMS_nuc.diversity.within) = rownames(MS.class_outtable)
      
      
      MS.nuc_diversity = momentsMS_nuc.diversity.within[, 2]
      MS.div_per_site = cbind(MS.nuc_diversity / MS.class@n.sites)
      colnames(MS.div_per_site) = "div_per_site"
      #MS.STATS.class@D ## blank
      #MS.STATS.class@BDF ## blank
      
      ## THIS WORKS
      MS.STATS.class_outtable = cbind(
        MS.haplotype.diversity,
        momentsMS_haplotype.counts,
        ## remove the 0th moment, probability
        momentsMS_minor.allele.freqs,
        momentsMS_nuc.diversity.within,
        MS.div_per_site
      )
      
      print("Run num:")
      print(i)
      print("num cols B:")
      print(ncol(MS.STATS.class_outtable))
      
      ## THIS WORKS
      MS_full_outtable = cbind(MS.class_outtable, MS.STATS.class_outtable)
      sort(colnames(MS_full_outtable))
      
      ## this is probably the problem
      MS_trim_outtable = MS_full_outtable
      #MS_trim_outtable = MS_full_outtable[,!apply(is.na(MS_full_outtable), 2, all)]
      #MS_trim_outtable = (MS_trim_outtable[rowSums(is.na(MS_trim_outtable[, ])) == 0,])
      
      MS_trim_outtable[, which(is.na(MS_trim_outtable))] = -100
      MS_trim_outtable[, which(is.null(MS_trim_outtable))] = -100
      
      for (col in 1:ncol(MS_trim_outtable)) {
        #print(col)
        
        MS_trim_outtable[, col] = as.character(MS_trim_outtable[, col])
        #MS_trim_outtable[, col] = as.numeric(MS_trim_outtable[, col])
      }
      
      print("number of cols final:")
      print(ncol(MS_trim_outtable))
      
      print("outputting")
      write.table(
        as.matrix(MS_trim_outtable),
        file = paste(path, folder, "/", prefix, ".popgenome.stats", sep = ""),
        quote = F,
        row.names = F,
        sep = "\t"
      )
      
      #tomove = (paste("/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/DONE/",msfile,sep=""))
      #currently = (paste("/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/",msfile,sep=""))
      #print("moving")
      #file.rename(currently, tomove)
      
      print("#####")
      
    }
  } 
  
  
  print(paste("Done with ", i, " of ", length(myfiles), sep = ""))
  
  
  
}
