rm(list = ls())
folder = "./"
newworking= T
sims=T
overwrite = F
doWindows = F
doMultPop = F
specieslist=""

#install.packages("PopGenome_2.7.5.tar.gz",repos=NULL,type="source",lib="localR/")

dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package, ")"
  ))))
}
packages = c("PopGenome", "moments", "R.utils", "gtools","progress")
for (p in packages) {
  dynamic_require(p)
}

## TODO:

#do_vcf = T

## the way to get around needing a directory is probably to create a temporary directory, copy the file into it, run analyses, and then empty the directory
path="/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo"
#path=getwd()
for (path in path) {
  print(path)
  setwd(path)
  #files = list.files(pattern = "vcf$")
  files = list.files(pattern = "Lampropeltis.*vcf$")
  files = files[!grepl("generated",files)]
  
  x <- file.info(files)
  files <- files[order(x$size)]
  
  if(newworking==T) {
    if (!(file.exists("stats"))){
      dir.create(file.path("stats"))
      print("making stats")
    }
  }
  
  split = strsplit(path,"/")[[1]]
  split = split[length(split)]
  print(split)  
  
  for (filename in files) {
    
    print(filename)
    setwd(path)
    
    ## check if it has 7+ lines
    temp = readLines(filename, n = 10)
    if(length(temp) <= 6) {
      print("-----NODATA")
      gzip(filename)
    } else {
      
      if (file.exists(paste(path, "/",filename,"_STATS.txt", sep = "")) && overwrite == F ) {
        print("-----SKIPPING")
        gzip(filename)
      } else {
        if(newworking==T) {
          print("moving")
          file.copy(from=filename, to="stats")
        }
        start_time <- Sys.time()
        print(paste("READING IN", length(files), "FILES", sep = " "))
        
        ## you must make sure that the first X lines number matches the last X lines
        
        if(newworking==T) {
          pathtoreadin = paste(path,"/stats",sep="")
        } else {
          pathtoreadin = paste(path,sep="")
        }
        print("start read")
        GENOME.class = readData(
          pathtoreadin,
          populations = F,        outgroup = F,
          include.unknown = T,        gffpath = F,
          format = "VCF",        parallized = F,
          progress_bar_switch = T,        FAST = T,
          big.data = T,        SNP.DATA = T
        )
        print("done read")
        GENOME.class_outtable = cbind(get.sum.data(GENOME.class))
        
        numind = length(get.individuals(GENOME.class)[[1]]) ## this will be 2N
        print(paste(numind, "individuals", sep = " "))
        print(get.individuals(GENOME.class)[[1]])
        
        ## FIRST WE DO ALL OF THE SINGLE POP STUFF
        print("GET FULL STATS PANMIXIA")
        {
          print("neutrality stats")
          GENOME.class <-
            neutrality.stats(GENOME.class, detail = T, do.R2 = T)
          pantemp = get.neutrality(GENOME.class, theta = T, stats = T)[[1]]
          GENOME.class_outtable = cbind(GENOME.class_outtable, pantemp)
          print("linkage stats")
          x = try(withTimeout(
            linkage.stats(
              GENOME.class,
              detail = T, ## F makes it work on the empirical
              do.ZnS = F, ## f makes it work on the empirical -- having both these be T fails, see if you can get them to play nice
              do.WALL = T 
            ),
            timeout = 120,
            silent = T
          ),
          silent = F)
          ## lengthy
          if (class(x) != class(GENOME.class)) {
            print("LINKAGE DID NOT WORK")
            print("WILL SKIP R2")
            skip = T
          } else {
            GENOME.class = x
            skip = F
            pantemp = get.linkage(GENOME.class)[[1]]
            GENOME.class_outtable = cbind(GENOME.class_outtable, pantemp)
          }
          
          print("recomb stats")
          x <-
            try(withTimeout(recomb.stats(GENOME.class),
                            timeout = 120,
                            silent = T))
          ## lengthy
          
          if (class(x) != class(GENOME.class)) {
            print("RECOMB DID NOT WORK")
          } else {
            GENOME.class = x
            pantemp = get.recomb(GENOME.class)[[1]]
            GENOME.class_outtable = cbind(GENOME.class_outtable, pantemp)
          }
          
          print("F_ST stats")
          GENOME.class <-
            F_ST.stats.2(GENOME.class, snn = T, Phi_ST = T)
          GENOME.class <- F_ST.stats(GENOME.class, detail = T)
          pantemp = cbind(
            get.F_ST(GENOME.class, mode = F, pairwise = F),
            GENOME.class@Hudson.Snn,
            GENOME.class@Phi_ST,
            get.diversity(GENOME.class)[[1]]
          )
          GENOME.class_outtable = cbind(GENOME.class_outtable, pantemp)
          print("diversity stats")
          GENOME.class <-
            diversity.stats(GENOME.class, pi = T, keep.site.info = T) ## doesn't work with pops specified?
          panmtemp = get.diversity(GENOME.class)[[1]]
          GENOME.class_outtable = cbind(GENOME.class_outtable, pantemp)
          print("sweep stats")
          GENOME.class <-
            sweeps.stats(GENOME.class, freq.table = T, FST = T)
          print("detail stats")
          GENOME.class <-
            detail.stats(
              GENOME.class,
              biallelic.structure = T,
              mismatch.distribution = T,
              site.spectrum = T,
              site.FST = T
            )
          panmtemp = get.detail(GENOME.class)[[1]]
          GENOME.class_outtable = cbind(GENOME.class_outtable, pantemp)
          if (skip == F) {
            print("calc r2")
            x <-
              try(withTimeout(calc.R2(GENOME.class),
                              timeout = 120,
                              silent = T))
            
            
            if (class(x) != class(GENOME.class)) {
              print("CALC r2 DID NOT WORK")
            } else {
              GENOME.class = x
              linkage_matrix = as.data.frame(t(GENOME.class@region.stats@linkage.disequilibrium[[1]][[1]]))
              colnames(GENOME.class_outtable) = make.unique(colnames(GENOME.class_outtable))
            }
            ## lengthy
            
            
          } else {
            print("SKIPPING CALC R2")
          }
          
        }
        print("num cols panmixia:")
        print(ncol(GENOME.class_outtable))
        ## DO THE STATS AND REGIONS WITH ONE POP
        print("STATS AND REGION PANMIXIA")
        {
          ## get region.data as own class
          REGION.class = GENOME.class@region.data
          ## get region.stats as own class
          STATS.class = GENOME.class@region.stats
          STATS.class@nucleotide.diversity[[1]] ## NAN WITH BIL
          STATS.class@haplotype.diversity[[1]] ## NAN WITH BIL
          STATS.class@haplotype.counts ## POSSIBLE USEFUL
          STATS.class@minor.allele.freqs ## USEFUL
          STATS.class@linkage.disequilibrium[[1]][[1]] ## SAME AS LINKAGE_MARTIX
          STATS.class@biallelic.structure [[1]][[1]] ## NAN WITH BIL
          STATS.class@site.FST[[1]] ## NAN WITH BIL
          STATS.class@missing.freqs ## BLANK WITH BIL
          STATS.class@nuc.diversity.between ## BLANK WITH BIL
          STATS.class@nuc.diversity.within[[1]] ##USEFUL
          STATS.class@D ## BLANK WITH BIL
          #STATS.class@BDF ## BLANK WITH BIL
          GENOME.haplotype.diversity = as.data.frame(cbind((STATS.class@haplotype.diversity)))
          colnames(GENOME.haplotype.diversity) = "GENOME.haplotype.diversity" ## exactly identical to haplotype.diversity.within
          rownames(GENOME.haplotype.diversity) = rownames(GENOME.class_outtable)
          print("moments hap count")
          momentsGENOME_haplotype.counts = t(sapply(1:length(STATS.class@haplotype.counts), function(i) {
            #print(i)
            data = unlist(STATS.class@haplotype.counts[i])
            if (is.null(data)) {
              print("failed momentsGENOME_haplotype.counts")
              print(i)
              return(c(1,-100,-100,-100,-100))
            } else {
              return(all.moments(unlist(
                STATS.class@haplotype.counts[i]
              ),
              order.max = 4))
            }
          }))
          colnames(momentsGENOME_haplotype.counts) = c(
            "prob_haplotype.counts",
            "mean_haplotype.counts",
            "var_haplotype.counts",
            "skew_haplotype.counts",
            "kurt_haplotype.counts"
          )
          rownames(momentsGENOME_haplotype.counts) = rownames(GENOME.class_outtable)
          print("moments minor")
          momentsGENOME_minor.allele.freqs = t(sapply(1:length(STATS.class@minor.allele.freqs), function(i) {
            data = unlist(STATS.class@minor.allele.freqs[i])
            if (is.null(data)) {
              print(i)
              print("failed momentsGENOME_minor.allele.freqs")
              return(c(1,-100,-100,-100,-100))
            } else {
              return(all.moments(
                unlist(STATS.class@minor.allele.freqs[i]),
                order.max = 4
              ))
            }
          }))
          colnames(momentsGENOME_minor.allele.freqs) = c(
            "prob_minor.allele.freqs",
            "mean_minor.allele.freqs",
            "var_minor.allele.freqs",
            "skew_minor.allele.freqs",
            "kurt_minor.allele.freqs"
          )
          rownames(momentsGENOME_minor.allele.freqs) = rownames(GENOME.class_outtable)
          print("moments nuc div within")
          momentsGENOME_nuc.diversity.within = t(sapply(1:length(STATS.class@nuc.diversity.within), function(i) {
            data = unlist(STATS.class@nuc.diversity.within[i])
            if (is.null(data)) {
              print("failed momentsGENOME_nuc.diversity.within")
              print(i)
              return(c(1,-100,-100,-100,-100))
            } else {
              return(all.moments(
                unlist(STATS.class@nuc.diversity.within[i]),
                order.max = 4
              ))
            }
          }))
          colnames(momentsGENOME_nuc.diversity.within) = c(
            "prob_nuc.diversity.within",
            "mean_nuc.diversity.within",
            "var_nuc.diversity.within",
            "skew_nuc.diversity.within",
            "kurt_nuc.diversity.within"
          )
          rownames(momentsGENOME_nuc.diversity.within) = rownames(GENOME.class_outtable)
          GENOME.nuc_diversity = momentsGENOME_nuc.diversity.within[, 2]
          GENOME.div_per_site = cbind(GENOME.nuc_diversity / GENOME.class@n.sites)
          colnames(GENOME.div_per_site) = "div_per_site"
          STATS.class_outtable = cbind(
            GENOME.haplotype.diversity,
            momentsGENOME_haplotype.counts,
            ## remove the 0th moment, probability
            momentsGENOME_minor.allele.freqs,
            momentsGENOME_nuc.diversity.within,
            GENOME.div_per_site
          )
          GENOME.class_outtable = cbind(GENOME.class_outtable, STATS.class_outtable)
          colnames(GENOME.class_outtable) = make.unique(colnames(GENOME.class_outtable))
        }
        print("num cols panmixia 2:")
        print(ncol(GENOME.class_outtable))
        ## MAKE REGIONS WITH ONEPOP
        ## THEN WE DO ALL OF THE MULTI POP STUFF
        
        if(doMultPop == T) {
          
          print("GET FULL STATS 2POP")
          {
            ## this is for crissale below
            # pop1 = c("indiv1","indiv2","indiv6","indiv9","indiv12",
            #          "indiv13","indiv15","indiv17")
            # pop2=c("indiv0","indiv3", "indiv4", "indiv5",
            #        "indiv7","indiv8","indiv10",
            #        "indiv11","indiv14","indiv16")
            # change to bilineata -- not real just garbage
            pop1 = c(
              "indiv1",
              "indiv2",
              "indiv6",
              "indiv9",
              "indiv12",
              "indiv13",
              "indiv15",
              "indiv17",
              "indiv18",
              "indiv19"
            )
            pop2 = c(
              "indiv0",
              "indiv3",
              "indiv4",
              "indiv5",
              "indiv7",
              "indiv8",
              "indiv10",
              "indiv11",
              "indiv14",
              "indiv16"
            )
            GENOME.class <-
              set.populations(GENOME.class, list(pop1, pop2)) ## don't set this because will screw up stats below
            bimat = get.biallelic.matrix(GENOME.class, 1)
            ## generate all the data for main VCF file
            print("neutrality stats")
            GENOME.class <-
              neutrality.stats(
                GENOME.class,
                detail = T,
                do.R2 = T,
                new.populations = list(pop1, pop2)
              )
            pop1temp = get.neutrality(GENOME.class, theta = T, stats = T)[[1]]
            pop2temp = get.neutrality(GENOME.class, theta = T, stats = T)[[2]]
            GENOME.class_outtable = cbind(GENOME.class_outtable, pop1temp, pop2temp)
            
            
            
            print("linkage stats")
            x = try(withTimeout(
              linkage.stats(
                GENOME.class,
                detail = T,
                do.ZnS = T,
                do.WALL = T,
                new.populations = list(pop1, pop2)
              ),
              timeout = 120,
              silent = T
            ))
            ## lengthy
            if (class(x) != class(GENOME.class)) {
              print("LINKAGE DID NOT WORK")
              print("WILL SKIP R2")
              skip = T
            } else {
              GENOME.class = x
              skip = F
              pop1temp = get.linkage(GENOME.class)[[1]]
              pop2temp = get.linkage(GENOME.class)[[2]]
              GENOME.class_outtable = cbind(GENOME.class_outtable, pop1temp, pop2temp)
              
            }
            
            print("recomb stats")
            x <-
              try(withTimeout(
                recomb.stats(GENOME.class, new.populations = list(pop1, pop2)),
                timeout = 120,
                silent = T
              ))
            ## lengthy
            
            if (class(x) != class(GENOME.class)) {
              print("RECOMB DID NOT WORK")
            } else {
              GENOME.class = x
              pop1temp = get.recomb(GENOME.class)[[1]]
              pop2temp = get.recomb(GENOME.class)[[2]]
              GENOME.class_outtable = cbind(GENOME.class_outtable, pop1temp, pop2temp)
            }
            
            
            print("F_ST stats")
            GENOME.class <-
              F_ST.stats.2(
                GENOME.class,
                snn = T,
                Phi_ST = T,
                new.populations = list(pop1, pop2)
              )
            GENOME.class <-
              F_ST.stats(GENOME.class,
                         detail = T,
                         new.populations = list(pop1, pop2))
            popstemp = cbind(
              get.F_ST(GENOME.class, mode = F, pairwise = F),
              GENOME.class@Hudson.Snn,
              GENOME.class@Phi_ST,
              get.diversity(GENOME.class)[[1]],
              get.diversity(GENOME.class)[[2]]
            )
            GENOME.class_outtable = cbind(GENOME.class_outtable, popstemp)
            print("MKT stats")
            GENOME.class <-
              MKT(
                GENOME.class,
                do.fisher.test = T,
                new.populations = list(pop1, pop2)
              ) ## only works if multiple pops
            GENOME.class_outtable = cbind(GENOME.class_outtable, get.MKT(GENOME.class))
            print("detail stats")
            GENOME.class <-
              detail.stats(
                GENOME.class,
                biallelic.structure = T,
                mismatch.distribution = T,
                site.spectrum = T,
                site.FST = T,
                new.populations = list(pop1, pop2)
              )
            pop1temp = get.detail(GENOME.class, biallelic.structure = F)[[1]]
            pop2temp = get.detail(GENOME.class, biallelic.structure = F)[[2]]
            GENOME.class_outtable = cbind(GENOME.class_outtable, pop1temp, pop2temp)
          }
          print("num cols 2POP:")
          print(ncol(GENOME.class_outtable))
          print("STATS AND REGIONS 2POP")
          {
            GENOME.class <-
              set.populations(GENOME.class, list(pop1, pop2)) ## don't set this because will screw up stats below
            ## get region.data as own class
            REGION.class = GENOME.class@region.data
            ## get region.stats as own class
            STATS.class = GENOME.class@region.stats
            STATS.class@nucleotide.diversity[[1]]
            STATS.class@haplotype.diversity[[1]] ## NAN WITH BIL
            STATS.class@haplotype.counts[[1]] ## NAN WITH BIL
            STATS.class@minor.allele.freqs[[1]] ## USEFUL
            STATS.class@linkage.disequilibrium[[1]][[2]] ## NAN WITH BIL
            STATS.class@biallelic.structure[[1]]$POP ## MAYBE USEFUL
            STATS.class@site.FST[[1]] ## NAN WITH BIL
            STATS.class@missing.freqs ## BLANK WITH BIL
            STATS.class@nuc.diversity.between ## BLANK WITH BIL
            STATS.class@nuc.diversity.within[[1]] ##USEFUL
            STATS.class@D ## BLANK WITH BIL
            #STATS.class@BDF ## BLANK WITH BIL
            GENOME.haplotype.diversity = as.data.frame(cbind((STATS.class@haplotype.diversity)))
            colnames(GENOME.haplotype.diversity) = "GENOME.haplotype.diversity" ## exactly identical to haplotype.diversity.within
            rownames(GENOME.haplotype.diversity) = rownames(GENOME.class_outtable)
            print("moments hap count")
            momentsGENOME_haplotype.counts = t(sapply(1:length(STATS.class@haplotype.counts), function(i) {
              #print(i)
              data = unlist(STATS.class@haplotype.counts[i])
              if (is.null(data)) {
                print("failed momentsGENOME_haplotype.counts")
                print(i)
                return(c(1,-100,-100,-100,-100))
              } else {
                return(all.moments(unlist(
                  STATS.class@haplotype.counts[i]
                ),
                order.max = 4))
              }
            }))
            colnames(momentsGENOME_haplotype.counts) = c(
              "prob_haplotype.counts",
              "mean_haplotype.counts",
              "var_haplotype.counts",
              "skew_haplotype.counts",
              "kurt_haplotype.counts"
            )
            rownames(momentsGENOME_haplotype.counts) = rownames(GENOME.class_outtable)
            print("moments minor")
            momentsGENOME_minor.allele.freqs = t(sapply(1:length(STATS.class@minor.allele.freqs), function(i) {
              data = unlist(STATS.class@minor.allele.freqs[i])
              if (is.null(data)) {
                print(i)
                print("failed momentsGENOME_minor.allele.freqs")
                return(c(1,-100,-100,-100,-100))
              } else {
                return(all.moments(
                  unlist(STATS.class@minor.allele.freqs[i]),
                  order.max = 4
                ))
              }
            }))
            colnames(momentsGENOME_minor.allele.freqs) = c(
              "prob_minor.allele.freqs",
              "mean_minor.allele.freqs",
              "var_minor.allele.freqs",
              "skew_minor.allele.freqs",
              "kurt_minor.allele.freqs"
            )
            rownames(momentsGENOME_minor.allele.freqs) = rownames(GENOME.class_outtable)
            print("moments nuc div within")
            momentsGENOME_nuc.diversity.within = t(sapply(1:length(STATS.class@nuc.diversity.within), function(i) {
              data = unlist(STATS.class@nuc.diversity.within[i])
              if (is.null(data)) {
                print("failed momentsGENOME_nuc.diversity.within")
                print(i)
                return(c(1,-100,-100,-100,-100))
              } else {
                return(all.moments(
                  unlist(STATS.class@nuc.diversity.within[i]),
                  order.max = 4
                ))
              }
            }))
            colnames(momentsGENOME_nuc.diversity.within) = c(
              "prob_nuc.diversity.within",
              "mean_nuc.diversity.within",
              "var_nuc.diversity.within",
              "skew_nuc.diversity.within",
              "kurt_nuc.diversity.within"
            )
            rownames(momentsGENOME_nuc.diversity.within) = rownames(GENOME.class_outtable)
            GENOME.nuc_diversity = momentsGENOME_nuc.diversity.within[, 2]
            GENOME.div_per_site = cbind(GENOME.nuc_diversity / GENOME.class@n.sites)
            colnames(GENOME.div_per_site) = "div_per_site"
            STATS.class_outtable = cbind(
              GENOME.haplotype.diversity,
              momentsGENOME_haplotype.counts,
              ## remove the 0th moment, probability
              momentsGENOME_minor.allele.freqs,
              momentsGENOME_nuc.diversity.within,
              GENOME.div_per_site
            )
            GENOME.class_outtable = cbind(GENOME.class_outtable, STATS.class_outtable)
            colnames(GENOME.class_outtable) = make.unique(colnames(GENOME.class_outtable))
          }
          print("num cols 2POP 2:")
          print(ncol(GENOME.class_outtable))
          
        }
        
        
        if (doWindows == T) {
          print("SLIDING WINDOWS -- NOT WORKING?")
          ## generate sliding windows
          #GENOME.class.slide = sliding.window.transform(GENOME.class,
          #                                              width=GENOME.class@n.sites2, ## this works because it makes the whole MS thing one region
          #                                              jump=GENOME.class@n.sites2, ## this works because it makes the whole MS thing one region
          #                                              type=2,
          #                                              whole.data=T)
          GENOME.class.slide = sliding.window.transform(
            GENOME.class,
            width = 100,
            ## this works because it makes the whole MS thing one region
            jump = 50,
            ## this works because it makes the whole MS thing one region
            type = 2,
            whole.data = T
          )
          GENOME.class.slide@region.names
          show.slots(GENOME.class.slide)
          GENOME.class.slide <-
            neutrality.stats(
              GENOME.class.slide,
              detail = T,
              do.R2 = T,
              FAST = F
            )
          GENOME.class.slide <-
            linkage.stats(
              GENOME.class.slide,
              detail = T,
              do.ZnS = T,
              do.WALL = T
            ) ## lengthy
          GENOME.class.slide <-
            recomb.stats(GENOME.class.slide) ## lengthy
          GENOME.class.slide <-
            F_ST.stats.2(GENOME.class.slide, snn = T, Phi_ST = T) ## never relevant with one pop
          GENOME.class.slide <-
            F_ST.stats(
              GENOME.class.slide,
              detail = T,
              only.haplotype.counts = F,
              FAST = F
            ) ## never relevant with one pop
          GENOME.class.slide <-
            diversity.stats(GENOME.class.slide,
                            pi = T,
                            keep.site.info = T)
          GENOME.class.slide <-
            sweeps.stats(GENOME.class.slide,
                         freq.table = T,
                         FST = T)
          GENOME.class.slide <-
            MKT(
              GENOME.class.slide,
              new.populations = list(pop1, pop2),
              do.fisher.test = T,
              fixed.threshold.fst = T
            ) ## only works if mult pops
          GENOME.class.slide <-
            detail.stats(
              GENOME.class.slide,
              biallelic.structure = T,
              mismatch.distribution = T,
              site.spectrum = T,
              site.FST = T
            )
          GENOME.class.slide <- calc.R2(GENOME.class.slide) ## lengthy
          ## get data for windows
          get.neutrality(GENOME.class.slide, theta = T, stats = T)[[1]]
          get.neutrality(GENOME.class.slide, theta = T, stats = T)[[2]]
          get.linkage(GENOME.class.slide)[[1]]
          get.linkage(GENOME.class.slide)[[2]]
          get.recomb(GENOME.class.slide)[[1]]
          get.recomb(GENOME.class.slide)[[2]]
          get.F_ST(GENOME.class.slide) ## never relevant with one pop
          get.diversity(GENOME.class.slide)[[1]]
          #get.sweeps(GENOME.class.slide) ## not working
          #get.MKT(GENOME.class.slide)  ## only works if mult pops
          #get.detail(GENOME.class.slide,biallelic.structure = T)
          GENOME.class.slide_outtable = cbind(
            #get.sum.data(GENOME.class.slide),
            get.neutrality(GENOME.class.slide, theta = T, stats = T)[[1]],
            get.neutrality(GENOME.class.slide, theta = T, stats = T)[[2]],
            get.F_ST(
              GENOME.class.slide,
              mode = F,
              pairwise = F
            ),
            get.MKT(GENOME.class.slide),
            ## think only works if have actual ACTG values
            get.linkage(GENOME.class.slide)[[1]],
            ##
            get.linkage(GENOME.class.slide)[[2]],
            ##
            get.recomb(GENOME.class.slide)[[1]],
            ##
            get.recomb(GENOME.class.slide)[[2]],
            ##
            get.diversity(GENOME.class.slide)[[1]],
            get.diversity(GENOME.class.slide)[[2]]
          )
          colnames(GENOME.class.slide_outtable) = make.unique(colnames(GENOME.class.slide_outtable))
          View(GENOME.class.slide_outtable)
          plot(GENOME.class.slide@nuc.diversity.within)
          GENOME.class_outtable = cbind(GENOME.class_outtable, GENOME.class.slide_outtable)
          colnames(GENOME.class_outtable) = make.unique(colnames(GENOME.class_outtable))
        } else {
          print("SKIPPING WINDOWS")
        }
        print("num cols windows:")
        print(ncol(GENOME.class_outtable))
        
        head(GENOME.class_outtable)
        
        write.table(
          as.matrix(GENOME.class_outtable),
          paste(path, "/",filename,"_STATS.txt", sep = ""),
          sep = "\t",
          row.names = T,
          col.names = T,
          append=F
        )
        
        
        
        end_time <- Sys.time()
        print(end_time - start_time)
        
        if(newworking==T){
          setwd(paste(path,"/working",sep=""))
        } else {
          setwd(path) 
        }
        
        file.remove(filename)
        setwd(path)
        
        rm(GENOME.class)
        rm(STATS.class)
        rm(GENOME.class_outtable)
        rm(STATS.class_outtable)
        
        
      }
    }
  }
  
}


