library(PopGenome)
library(moments)
#options(expressions = 10000)

## TODO:
## get this so that it can input either a VCF or a MS file and then choose which to output?
## for now, just do the MS files

## make it so that you can give it a list of the txt files to output

#path = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/"
#path = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/selection/"
#path = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/spatial/"
#path = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/demography/"
path = "/Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/runs/"

## change the path to other folders and include the .txt files 1:4

setwd(path)

files = list.files(pattern = "model.*txt$", recursive = TRUE)
myfiles = files
#myfiles = files[3]


## make this a for loop

do_vcf = F

if (do_vcf == T) {
  vcffolder = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/pogenometest/vcf/"
  ## the only thing in there must be the vcfs?
  
  vcffile = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/pogenometest/vcf/stable-neut-nospat-1pop-isol-1542226966-1-overlaid.vcf"
  GENOME.class = readData(vcffolder, format = "VCF")
  
  "
  1                  readData            Reading data           get.sum.data
  2          neutrality.stats        Neutrality tests         get.neutrality
  3             linkage.stats  Linkage disequilibrium            get.linkage
  4              recomb.stats           Recombination             get.recomb
  5                F_ST.stats          Fixation index get.F_ST,get.diversity
  6           diversity.stats             Diversities          get.diversity
  7              sweeps.stats        Selective sweeps             get.sweeps
  8                       MKT  McDonald-Kreitman test                get.MKT
  9              detail.stats        Mixed statistics             get.detail
  "
  
  
  ## generate all the data for main VCF file
  
  show.slots(GENOME.class)
  
  GENOME.class <- neutrality.stats(GENOME.class)
  GENOME.class <- linkage.stats(GENOME.class, detail = T)
  GENOME.class <- recomb.stats(GENOME.class)
  GENOME.class <- F_ST.stats(GENOME.class)
  GENOME.class <- diversity.stats(GENOME.class)
  GENOME.class <- sweeps.stats(GENOME.class)
  GENOME.class <- MKT(GENOME.class)
  GENOME.class <-
    detail.stats(GENOME.class, biallelic.structure = T)
  GENOME.class <- calc.R2(GENOME.class)
  GENOME.class <-
    introgression.stats(GENOME.class) ## only works if outgroup
  
  ## get all the data for main VCF file
  
  
  
  get.sum.data(GENOME.class)
  get.neutrality(GENOME.class)[[1]]
  get.linkage(GENOME.class)[[1]]
  get.recomb(GENOME.class)[[1]]
  #get.F_ST(GENOME.class) ## not if only one pop
  get.diversity(GENOME.class)[[1]]
  #get.sweeps(GENOME.class) ## not working
  #get.MKT(GENOME.class) ## doesn't work with single pop
  #get.detail(GENOME.class,biallelic.structure = T) ## all zeroes
  
  GENOME.class@n.sites
  GENOME.class@region.names
  
  GENOME.class_outtable = cbind(
    get.sum.data(GENOME.class),
    get.neutrality(GENOME.class)[[1]],
    get.linkage(GENOME.class)[[1]],
    get.recomb(GENOME.class)[[1]],
    get.diversity(GENOME.class)[[1]]
  )
  
  ## get region.data as own class
  
  REGION.class = GENOME.class@region.data
  
  #REGION.class@populations ## list of individuals by pop
  #REGION.class@populations2 ## list of individuals by pop?
  #REGION.class@outgroup ## whether there is outgroup or not
  #REGION.class@transitions ## transitions/transversions, not applicable
  #REGION.class@biallelic.matrix[[1]] ## ?
  #REGION.class@n.singletons ## blank
  #REGION.class@biallelic.sites ## large
  #REGION.class@reference ## blank
  #REGION.class@n.nucleotides # #blank
  #REGION.class@biallelic.compositions ## blank
  #REGION.class@synonymous ## all NaN
  #REGION.class@biallelic.substitutions ## large and probably not applicable
  #REGION.class@polyallelic.sites ## blank
  #REGION.class@sites.with.gaps  ## blank
  #REGION.class@sites.with.unknowns ## blank
  #REGION.class@minor.alleles ## blank
  #REGION.class@codons ## blank
  #REGION.class@IntronSNPS ## blank
  #REGION.class@UTRSNPS ## blank
  #REGION.class@CodingSNPS ## blank
  #REGION.class@ExonSNPS ## blank
  #REGION.class@GeneSNPS ## blank
  #REGION.class@included ## blank
  
  ## get region.stats as own class
  
  STATS.class = GENOME.class@region.stats
  
  #STATS.class@nucleotide.diversity[[1]] ## blank
  STATS.class@haplotype.diversity[[1]]
  STATS.class@haplotype.counts ## large
  ## probability, mean, variance, skewness, kurtosis
  moments_haplotype.counts = rbind(all.moments(unlist(STATS.class@haplotype.counts), order.max =
                                                 4))
  colnames(moments_haplotype.counts) = c(
    "prob_haplotype.counts",
    "mean_haplotype.counts",
    "var_haplotype.counts",
    "skew_haplotype.counts",
    "kurt_haplotype.counts"
  )
  
  
  STATS.class@minor.allele.freqs ## large
  moments_minor.allele.freqs = rbind(all.moments(unlist(STATS.class@minor.allele.freqs), order.max =
                                                   4))
  colnames(moments_minor.allele.freqs) = c(
    "prob_minor.allele.freqs",
    "mean_minor.allele.freqs",
    "var_minor.allele.freqs",
    "skew_minor.allele.freqs",
    "kurt_minor.allele.freqs"
  )
  
  
  STATS.class@linkage.disequilibrium  ## HUGE amount of stuff - not sure how to get
  #STATS.class@biallelic.structure[[1]] ## blank
  #STATS.class@site.FST ## blank
  #STATS.class@missing.freqs # #blank
  #STATS.class@nuc.diversity.between ## blank
  STATS.class@nuc.diversity.within ## large
  moments_nuc.diversity.within = rbind(all.moments(unlist(STATS.class@nuc.diversity.within), order.max =
                                                     4))
  colnames(moments_nuc.diversity.within) = c(
    "prob_nuc.diversity.within",
    "mean_nuc.diversity.within",
    "var_nuc.diversity.within",
    "skew_nuc.diversity.within",
    "kurt_nuc.diversity.within"
  )
  nuc_diversity = moments_nuc.diversity.within[2]
  div_per_site = nuc_diversity / GENOME.class@n.sites
  names(div_per_site) = "div_per_site"
  #STATS.class@D ## blank
  #STATS.class@BDF ## blank
  
  STATS.class_outtable = cbind(
    haplotype.diversity = STATS.class@haplotype.diversity[[1]][1],
    moments_haplotype.counts,
    ## remove the 0th moment, probability
    moments_minor.allele.freqs,
    moments_nuc.diversity.within,
    div_per_site
  )
  
  ## generate sliding windows
  
  GENOME.class.slide = sliding.window.transform(
    GENOME.class,
    width = GENOME.class@n.sites2,
    ## this works because it makes the whole MS thing one region
    jump = GENOME.class@n.sites2,
    ## this works because it makes the whole MS thing one region
    type = 1,
    whole.data = F
  )
  GENOME.class.slide@region.names
  
  show.slots(GENOME.class.slide)
  
  ## generate data for windows
  
  GENOME.class.slide <- neutrality.stats(GENOME.class.slide)
  GENOME.class.slide <-
    linkage.stats(GENOME.class.slide, detail = T) ## lengthy
  GENOME.class.slide <-
    recomb.stats(GENOME.class.slide) ## lengthy
  GENOME.class.slide <-
    F_ST.stats(GENOME.class.slide) ## never relevant with one pop
  GENOME.class.slide <- diversity.stats(GENOME.class.slide)
  GENOME.class.slide <- sweeps.stats(GENOME.class.slide)
  GENOME.class.slide <-
    MKT(GENOME.class.slide) ## only works if mult pops
  GENOME.class.slide <-
    detail.stats(GENOME.class.slide, biallelic.structure = T)
  GENOME.class.slide <- calc.R2(GENOME.class.slide) ## lengthy
  
  ## get data for windows
  
  #get.sum.data(GENOME.class.slide) ## not relevant
  get.neutrality(GENOME.class.slide)[[1]]
  get.linkage(GENOME.class.slide)[[1]]
  get.recomb(GENOME.class.slide)[[1]]
  #get.F_ST(GENOME.class.slide) ## never relevant with one pop
  get.diversity(GENOME.class.slide)[[1]]
  #get.sweeps(GENOME.class.slide) ## not working
  #get.MKT(GENOME.class.slide)  ## only works if mult pops
  #get.detail(GENOME.class.slide,biallelic.structure = T)
  
  GENOME.class.slide_outtable = cbind(
    ## this is identical to non-slide and therefore redundant
    get.neutrality(GENOME.class.slide)[[1]],
    get.linkage(GENOME.class.slide)[[1]],
    get.recomb(GENOME.class.slide)[[1]],
    get.diversity(GENOME.class.slide)[[1]]
  )
  
  
  plot(GENOME.class.slide@nuc.diversity.within)
  
  
  
}

if (do_vcf == F) {
  ##### with the MS data itself
  
  #msfolder  = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/pogenometest/ms/"
  msfolder = path
  
  ## for this to read, your header has to say "slim 20 100" or similar
  #msfile = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/pogenometest/ms/stable-neut-nospat-1pop-isol-1542226966.txt"
  #i = 1
  
  for (i in 1:length(myfiles)) {
    msfile = myfiles[i]
    print(msfile)
    
    ## check if the
    
    ## check the header
    lines = readLines(msfile)
    
    numlines = sum(lines == "//")
    
    header = lines[1]
    
    
    if (substr(header, 1, 4) != "slim") {
      #print("oops")
      lines[1] = paste("slim", header, sep = " ")
      writeLines(lines, msfile)
      numreps = as.numeric(strsplit(header, " ")[[1]][2])
    } else {
      numreps = as.numeric(strsplit(header, " ")[[1]][3])
    }
    
    if (numreps == numlines) {
      split = strsplit(msfile, split = "/")[[1]]
      folder = split[1]
      file = split[2]
      prefix = strsplit(file, split = "[.]")[[1]][1]
      
      if (file.exists(paste(path, folder, "/", prefix, ".popgenome.stats", sep = ""))) {
        print("DO NOT RUN, ALREADY DONE")
      } else {
        ## read in MS file
        MS.class = readMS(msfile)
        
        #show.slots(MS.class)
        
        ## generate stats for MS
        
        print("neutrality stats")
        MS.class <- neutrality.stats(MS.class)
        print("linkage stats")
        MS.class <- linkage.stats(MS.class, detail = T) ## lengthy
        print("recomb stats")
        MS.class <- recomb.stats(MS.class) ## lengthy
        print("F_ST stats")
        
        MS.class <-
          F_ST.stats(MS.class) ## not relevant with one pop?
        print("diversity stats")
        MS.class <- diversity.stats(MS.class)
        print("sweep stats")
        MS.class <- sweeps.stats(MS.class)
        print("MKT stats")
        #MS.class <- MKT(MS.class) ## only works if multiple pops
        print("detail stats")
        MS.class <- detail.stats(MS.class, biallelic.structure = T)
        print("calc r2")
        MS.class <- calc.R2(MS.class) ## lengthy
        
        ## get stats for MS
        #print(head(get.sum.data(MS.class)))
        #print(head(get.neutrality(MS.class)[[1]]))
        #print(head(get.linkage(MS.class)[[1]]))
        #print(head(get.recomb(MS.class)[[1]]))
        #head(get.F_ST(MS.class)) ## never relevant with one pop
        #print(head(get.diversity(MS.class)[[1]]))
        #head(get.sweeps(MS.class)) ## not working
        #head(get.MKT(MS.class))  ## only works if mult pops
        #head(get.detail(MS.class,biallelic.structure = T))
        
        MS.class_outtable = cbind(
          get.sum.data(MS.class),
          get.neutrality(MS.class)[[1]],
          get.linkage(MS.class)[[1]],
          get.recomb(MS.class)[[1]],
          get.diversity(MS.class)[[1]]
        )
        
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
            return(c(1, -100, -100, -100, -100))
            
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
            return(c(1, -100, -100, -100, -100))
            
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
            return(c(1, -100, -100, -100, -100))
            
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
        
        MS.STATS.class_outtable = cbind(
          MS.haplotype.diversity,
          momentsMS_haplotype.counts,
          ## remove the 0th moment, probability
          momentsMS_minor.allele.freqs,
          momentsMS_nuc.diversity.within,
          MS.div_per_site
        )
        
        print("num cols B:")
        print(ncol(MS.STATS.class_outtable))
        
        MS_full_outtable = cbind(MS.class_outtable, MS.STATS.class_outtable)
        sort(colnames(MS_full_outtable))
        
        ## this is probably the problem
        MS_trim_outtable = MS_full_outtable
        #MS_trim_outtable = MS_full_outtable[,!apply(is.na(MS_full_outtable), 2, all)]
        MS_trim_outtable = (MS_trim_outtable[rowSums(is.na(MS_trim_outtable[, ])) == 0,])
        
        for (col in 1:ncol(MS_trim_outtable)) {
          print(col)
          MS_trim_outtable[, col] = as.numeric(MS_trim_outtable[, col])
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
        print("#####")
        
      }
    } else {
      print("CANNOT RUN THIS FILE, NOT COMPLETE")
    }
  }
}



#####

#####
#
# msfile = "/Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/slim-stable-neut-nospat-1pop-isol-1542226966/stable-neut-nospat-1pop-isol-1542226966-1.withheader.temp"
#
#
#
# #ms = readMS("/Users/kprovost/Documents/Classes/Machine Learning/SLiMTreeSeqPub-master/p1_ms_exponential.txt")
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
