rm(list = ls())
## Rscript "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/5.testingPopgenome_VCF.R"
## Rscript "/home/kprovost/nas2/convert_vcf_to_temp/5.testingPopgenome_VCF.R"
#folder = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/VCFS/"
folder = "/home/kprovost/nas2/Analysis_SLiM/VCFS/"
newworking= T
sims=T
overwrite = F
doWindows = F
doMultPop = F
#specieslist=c(  "BELLII",  "BILINEATA",  "BRUNNEICAPILLUS",  "CRISSALE",  "CURVIROSTRE",  "FLAVICEPS",  "FUSCA",  "MELANURA",  "NITENS","SINUATUS")
specieslist=""

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

if(newworking==T) {
pathlist = paste(folder,
  specieslist,"/",sep = "")
} else {
  if (sims==F) {
  pathlist = paste(folder,"working/",sep = "")
  } else {
    pathlist=folder
  }
}
#do_vcf = T


## the way to get around needing a directory is probably to create a temporary directory, copy the file into it, run analyses, and then empty the directory

for (path in pathlist) {
  ## this will work on some but not all chromosiomes?
  print(path)
  setwd(path)
  files = list.files(pattern = "vcf$")
  x <- file.info(files)
  files <- files[order(x$size)]
  
  if(newworking==T) {
  if (!(file.exists("working"))){
    dir.create(file.path("working"))
  }
  }
  
  split = strsplit(path,"/")[[1]]
  split = split[length(split)]
  
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
      file.copy(from=filename, to="working")
      }
      start_time <- Sys.time()
      #print(paste("READING IN", length(files), "FILES", sep = " "))
      
      if(newworking==T) {
      pathtoreadin = paste(path,"/working",sep="")
      } else {
        pathtoreadin = paste(path,sep="")
      }
      GENOME.class = readData(
        pathtoreadin,
        populations = F,        outgroup = F,
        include.unknown = T,        gffpath = F,
        format = "VCF",        parallized = F,
        progress_bar_switch = T,        FAST = T,
        big.data = T,        SNP.DATA = T
      )
      GENOME.class_outtable = cbind(get.sum.data(GENOME.class))
      
      numind = length(get.individuals(GENOME.class)[[1]]) ## this will be 2N
      #print(paste(numind, "individuals", sep = " "))
      #print(get.individuals(GENOME.class)[[1]])
      
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



## custom readData to debug
myReadVCF <- function(vcffile){
  
  vcf       <- .Call("myReadVCFC",vcffile)
  matrix    <- vcf[[1]]
  positions <- vcf[[2]]
  rownames(matrix) <- vcf[[3]]
  
  rm(vcf)
  gc()
  
  o_b_j_sub    <- list(matrix=matrix,reference=NaN,positions=positions)
  
  return(o_b_j_sub)
}
progressBar <- function (i, total, symb_drawn = 0, symb = "=", nl = "\n") {
  
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
    cat("|            :            |            :            | 100 %", nl, "|", sep = "")
  }
  else {
    
    if (symb_drawn < 51) {
      
      progress <- round((i/total) * 52, digits = 0) - symb_drawn
      symb_drawn <- progress + symb_drawn
      
      while (progress > 0) {
        cat(symb);
        progress <- progress - 1
      }
      
      #return(symb_drawn)
    }
    else {
      if (i == total) {
        cat("| ;-) \n")
      }
    }
  }
  
  return(symb_drawn)
}
PopGenread <- function(filepath,format) {
  
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
readData_2 <- function(path, populations = FALSE, outgroup = FALSE, include.unknown = FALSE, 
                       gffpath = FALSE, format = "fasta", parallized = FALSE, progress_bar_switch = TRUE, 
                       FAST = FALSE, big.data = FALSE, SNP.DATA = FALSE) 
{
  if (!file.exists(path)) {
    stop("Cannot find path !")
  }
  if (!file.info(path)[2]) {
    stop("Put your file/files in a folder !")
  }
  if (gffpath[1] != FALSE) {
    if (!file.exists(gffpath)) {
      stop("Cannot find gff path !")
    }
    if (!file.info(gffpath)[2]) {
      stop("Put your file/files in a folder ! (GFF)")
    }
  }
  if (format == "HapMap") {
    SNP.DATA = TRUE
    FAST = TRUE
  }
  if (format == "VCF") {
    SNP.DATA = TRUE
    FAST = TRUE
  }
  if (SNP.DATA) {
    big.data = TRUE
  }
  if (parallized) {
    n.cores <- parallel::detectCores()
    files <- list.files(path)
    split_files <- split(files, sort(rep(1:n.cores, ceiling(length(files)/n.cores))))
    if (gffpath[1] != FALSE) {
      gff_files <- list.files(gffpath)
      split_files_gff <- split(gff_files, sort(rep(1:n.cores, 
                                                   ceiling(length(gff_files)/n.cores))))
    }
    xxx <- NULL
    yyy <- NULL
    for (i in 1:n.cores) {
      command <- paste("mkdir split", i, sep = "")
      system(command)
      filename <- paste("split", i, sep = "")
      filename_path <- file.path(getwd(), filename)
      sapply(split_files[[i]], function(x) {
        command <- file.path(path, x)
        command <- paste("cp", command, filename_path, 
                         sep = " ")
        system(command)
      })
      xxx <- c(xxx, filename_path)
      if (gffpath[1] != FALSE) {
        command <- paste("mkdir GFFRObjects_split", i, 
                         sep = "")
        system(command)
        filename <- paste("GFFRObjects_split", i, sep = "")
        filename_path <- file.path(getwd(), filename)
        sapply(split_files_gff[[i]], function(x) {
          command <- file.path(gffpath, x)
          command <- paste("cp", command, filename_path, 
                           sep = " ")
          system(command)
        })
        yyy <- c(yyy, filename_path)
      }
    }
    if (gffpath[1] == FALSE) {
      rueck <- parallel::mclapply(xxx, readData, format = format, 
                                  parallized = FALSE, progress_bar_switch = FALSE, 
                                  FAST = FAST, big.data = big.data, SNP.DATA = SNP.DATA, 
                                  include.unknown = include.unknown, mc.cores = n.cores, 
                                  mc.silent = TRUE)
    }
    else {
      cat("Calculation ... \n")
      INPUT <- vector("list", n.cores)
      for (iii in 1:n.cores) {
        INPUT[[iii]] <- c(xxx[iii], yyy[iii])
      }
      rueck <- parallel::mclapply(INPUT, function(x) {
        daten <- x[1]
        gffdaten <- x[2]
        TT <- readData(path = daten, gffpath = gffdaten, 
                       format = format, parallized = FALSE, progress_bar_switch = FALSE, 
                       FAST = FAST, big.data = big.data, SNP.DATA = SNP.DATA, 
                       include.unknown = include.unknown)
        return(TT)
      }, mc.cores = n.cores, mc.silent = TRUE, mc.preschedule = TRUE)
    }
    genome <- concatenate(rueck, n.cores)
    genome@basepath <- file.path(path)
    genome@project <- file.path(path)
    liste <- list.files(path, full.names = TRUE)
    genome@genelength <- length(liste)
    if (FAST) {
      genome@Pop_Slide$calculated <- TRUE
    }
    else {
      genome@Pop_Slide$calculated <- FALSE
    }
    if (!is.list(populations)) {
      genome@populations <- list(NULL)
    }
    else {
      genome@populations <- populations
    }
    genome@outgroup <- outgroup
    for (i in 1:n.cores) {
      command <- paste("rm -r split", i, sep = "")
      system(command)
      if (gffpath[1] != FALSE) {
        command <- paste("rm -r GFFRObjects_split", i, 
                         sep = "")
        system(command)
      }
    }
    cat("\n")
    return(genome)
  }
  methods <- "DATA"
  npops <- length(populations)
  popnames <- paste("pop", 1:npops)
  liste <- list.files(path, full.names = TRUE)
  liste2 <- list.files(path)
  liste3 <- gsub("\\.[a-zA-Z]*", "", liste2)
  ordered <- as.numeric(gsub("\\D", "", liste))
  if (!any(is.na(ordered))) {
    names(ordered) <- 1:length(ordered)
    ordered <- sort(ordered)
    ids <- as.numeric(names(ordered))
    liste <- liste[ids]
    liste2 <- liste2[ids]
    liste3 <- liste3[ids]
  }
  gff_objects <- vector("list", length(liste))
  SNP.GFF <- FALSE
  GFF.BOOL <- FALSE
  if (gffpath[1] != FALSE) {
    GFF.BOOL <- TRUE
    gff_liste <- list.files(gffpath, full.names = TRUE)
    gff_liste2 <- list.files(gffpath)
    gff_liste3 <- gsub("\\.[a-zA-Z]*", "", gff_liste2)
    treffer <- match(liste3, gff_liste3)
    if (any(is.na(treffer))) {
      cat("WARNING:: Could not find GFF files for:\n")
      cat(liste3[is.na(treffer)], "\n", sep = ",")
    }
    gff_liste <- gff_liste[treffer]
  }
  if (npops > 1) {
    if (outgroup[1] != FALSE) {
      poppairs <- choose(npops + 1, 2)
      pairs <- combn(1:(npops + 1), 2)
    }
    else {
      poppairs <- choose(npops, 2)
      pairs <- combn(1:(npops), 2)
    }
    nn <- paste("pop", pairs[1, 1], "/pop", pairs[2, 1], 
                sep = "")
    if (dim(pairs)[2] > 1) {
      for (xx in 2:dim(pairs)[2]) {
        m <- paste("pop", pairs[1, xx], "/pop", pairs[2, 
                                                      xx], sep = "")
        nn <- c(nn, m)
      }
    }
  }
  else {
    poppairs <- 1
    nn <- "pop1"
  }
  sizeliste <- length(liste)
  genome <- new("GENOME")
  genome@basepath <- file.path(path)
  genome@project <- file.path(path)
  genome@genelength <- sizeliste
  if (!is.list(populations)) {
    genome@populations <- list(NULL)
  }
  else {
    genome@populations <- populations
  }
  genome@poppairs <- nn
  genome@outgroup <- outgroup
  genome@region.names <- liste2
  DATABOOL <- is.element("DATA", methods)
  nsites <- integer(sizeliste)
  nsites2 <- integer(sizeliste)
  region.data <- new("region.data")
  region.stats <- new("region.stats")
  init <- vector("list", sizeliste)
  init2 <- numeric(sizeliste)
  init3 <- rep(NaN, sizeliste)
  populationsX <- init
  populations2 <- init
  popmissing <- init
  outgroupX <- init
  outgroup2 <- init
  CodingSNPS <- init
  UTRSNPS <- init
  IntronSNPS <- init
  ExonSNPS <- init
  GeneSNPS <- init
  reading.frame <- init
  rev.strand <- init
  Coding.matrix <- init
  Coding.matrix2 <- init
  UTR.matrix <- init
  Intron.matrix <- init
  Exon.matrix <- init
  Gene.matrix <- init
  Coding.region <- init2
  UTR.region <- init2
  Intron.region <- init2
  Exon.region <- init2
  Gene.region <- init2
  transitions <- init
  biallelic.matrix <- init
  biallelic.sites <- init
  biallelic.sites2 <- init
  reference <- init
  matrix_codonpos <- init
  synonymous <- init
  matrix_freq <- init
  n.singletons <- init
  polyallelic.sites <- init
  n.nucleotides <- init
  biallelic.compositions <- init
  biallelic.substitutions <- init
  minor.alleles <- init
  codons <- init
  sites.with.gaps <- init
  sites.with.unknowns <- init
  n.valid.sites <- init2
  n.gaps <- init2
  n.unknowns <- init2
  n.polyallelic.sites <- init2
  n.biallelic.sites <- init2
  trans.transv.ratio <- init3
  print("CHECKPOINT 1")
  if (progress_bar_switch) {
    progr <- progressBar()
  }
  print("CHECKPOINT 2")
  #print(sizeliste)
  #print(head(liste))
  for (xx in 1:sizeliste) {
    if (!progress_bar_switch) {
      print(liste[xx])
    }
    print("CHECKPOINT 3")
    CCC <- try(PopGenread(liste[xx], format), silent = TRUE)
    print("CHECKPOINT 4")
    #print(head(CCC))
    print("CHECKPOINT 4.1")
    #print(head(CCC$matrix))
    print("CHECKPOINT 4.2")
    if (is.na(CCC$matrix[1])) {
      next
    }
    print("CHECKPOINT 5")
    gen <- CCC$matrix
    pos <- CCC$positions
    ref <- CCC$reference
    gff.object.exists <- FALSE
    FIT <- FALSE
    print("CHECKPOINT 6")
    if (GFF.BOOL) {
      gff.object.exists <- TRUE
      if (SNP.DATA) {
        if (!is.na(gff_liste[xx])) {
          if (length(grep("GFFRObjects", gffpath)) != 
              0) {
            Robj <- load(gff_liste[xx])
            gff_object <- get(Robj[1])
          }
          else {
            gff_object <- gffRead(gff_liste[xx])
          }
          gff_object_fit <- fitting_gff_fast(pos, gff_object)
          FIT <- TRUE
        }
        else {
          gff.object.exists <- FALSE
        }
      }
      else {
        if (!is.na(gff_liste[xx])) {
          gff_object <- gffRead(gff_liste[xx])
        }
        else {
          gff.object.exists <- FALSE
        }
      }
      if (gff.object.exists) {
        gff_object <- parse_gff(gff_object, SNP.DATA = SNP.DATA)
        if (FIT) {
          gff_object_fit <- parse_gff(gff_object_fit, 
                                      SNP.DATA = SNP.DATA)
          GLOBAL.GFF$GFF <- NULL
        }
        else {
          gff_object_fit <- gff_object
        }
      }
      else {
        gff_object <- FALSE
        gff_object_fit <- FALSE
      }
    }
    else {
      gff_object <- FALSE
      gff_object_fit <- FALSE
    }
    if (is.matrix(gen)) {
      nsites[xx] <- dim(gen)[2]
      nsites2[xx] <- dim(gen)[2]
      result <- popgen(gen, Populations = populations, 
                       outgroup = outgroup, methods = methods, include.unknown = include.unknown, 
                       gff = gff_object_fit, FAST, SNP.DATA)
      rm(gen)
    }
    else {
      result <- NA
      next
    }
    if (progress_bar_switch) {
      progr <- progressBar(xx, sizeliste, progr)
    }
    if (xx == ceiling(sizeliste/2)) {
      gc()
    }
    if (is.list(result)) {
      populationsX[[xx]] <- result$populations
      populations2[[xx]] <- result$populations2
      popmissing[[xx]] <- result$popmissing
      outgroupX[[xx]] <- result$outgroup
      outgroup2[[xx]] <- result$outgroup2
      datt <- result$data.sum
      CodingSNPS[[xx]] <- datt$CodingSNPS
      UTRSNPS[[xx]] <- datt$UTRSNPS
      IntronSNPS[[xx]] <- datt$IntronSNPS
      ExonSNPS[[xx]] <- datt$ExonSNPS
      GeneSNPS[[xx]] <- datt$GeneSNPS
      if (GFF.BOOL & !gff.object.exists) {
        fillinit <- vector(, length(datt$biallelic.sites))
        CodingSNPS[[xx]] <- fillinit
        UTRSNPS[[xx]] <- fillinit
        IntronSNPS[[xx]] <- fillinit
        ExonSNPS[[xx]] <- fillinit
        GeneSNPS[[xx]] <- fillinit
      }
      transitions[[xx]] <- datt$transitions
      if (is.na(pos[1])) {
        biallelic.sites[[xx]] <- datt$biallelic.sites
        polyallelic.sites[[xx]] <- datt$polyallelic.sites
        sites.with.gaps[[xx]] <- datt$sites.with.gaps
        sites.with.unknowns[[xx]] <- datt$sites.with.unknowns
      }
      else {
        biallelic.sites2[[xx]] <- datt$biallelic.sites
        biallelic.sites[[xx]] <- pos[datt$biallelic.sites]
        nsites[xx] <- biallelic.sites[[xx]][length(biallelic.sites[[xx]])]
        polyallelic.sites[[xx]] <- pos[datt$polyallelic.sites]
        sites.with.gaps[[xx]] <- pos[datt$sites.with.gaps]
        sites.with.unknowns[[xx]] <- pos[datt$sites.with.unknowns]
      }
      if (!big.data) {
        biallelic.matrix[[xx]] <- datt$biallelic.matrix
        colnames(biallelic.matrix[[xx]]) <- biallelic.sites[[xx]]
        if (GFF.BOOL & gff.object.exists) {
          reading.frame[[xx]] <- gff_object$reading.frame
          rev.strand[[xx]] <- gff_object$rev.strand
          Coding.matrix[[xx]] <- gff_object$Coding
          UTR.matrix[[xx]] <- gff_object$UTR
          Intron.matrix[[xx]] <- gff_object$Intron
          Exon.matrix[[xx]] <- gff_object$Exon
          Gene.matrix[[xx]] <- gff_object$Gene
        }
      }
      else {
        biallelic.matrix[[xx]] <- ff(datt$biallelic.matrix, 
                                     dim = dim(datt$biallelic.matrix))
        close(biallelic.matrix[[xx]])
        dimnames(biallelic.matrix[[xx]]) <- list(rownames(datt$biallelic.matrix), 
                                                 biallelic.sites[[xx]])
        if (GFF.BOOL & gff.object.exists) {
          if (dim(gff_object$Coding)[1] > 0) {
            fill <- as.matrix(gff_object$Coding)
            Coding.matrix[[xx]] <- ff(fill, dim = dim(fill))
            close(Coding.matrix[[xx]])
            if (dim(gff_object_fit$Coding)[1] > 0) {
              reading.frame[[xx]] <- gff_object_fit$reading.frame
              rev.strand[[xx]] <- gff_object_fit$rev.strand
              fill <- as.matrix(gff_object_fit$Coding)
              Coding.matrix2[[xx]] <- ff(fill, dim = dim(fill))
              close(Coding.matrix2[[xx]])
            }
          }
          if (dim(gff_object$UTR)[1] > 0) {
            fill <- as.matrix(gff_object$UTR)
            UTR.matrix[[xx]] <- ff(fill, dim = dim(fill))
            close(UTR.matrix[[xx]])
          }
          if (dim(gff_object$Intron)[1] > 0) {
            fill <- as.matrix(gff_object$Intron)
            Intron.matrix[[xx]] <- ff(fill, dim = dim(fill))
            close(Intron.matrix[[xx]])
          }
          if (dim(gff_object$Exon)[1] > 0) {
            fill <- as.matrix(gff_object$Exon)
            Exon.matrix[[xx]] <- ff(fill, dim = dim(fill))
            close(Exon.matrix[[xx]])
          }
          if (dim(gff_object$Gene)[1] > 0) {
            fill <- as.matrix(gff_object$Gene)
            Gene.matrix[[xx]] <- ff(fill, dim = dim(fill))
            close(Gene.matrix[[xx]])
          }
        }
      }
      if (!is.na(ref[1])) {
        reference[[xx]] <- ref[datt$biallelic.sites]
      }
      matrix_codonpos[[xx]] <- datt$matrix_codonpos
      synonymous[[xx]] <- datt$synonymous
      matrix_freq[[xx]] <- datt$matrix_freq
      n.singletons[[xx]] <- datt$n.singletons
      n.nucleotides[[xx]] <- datt$n.nucleotides
      biallelic.compositions[[xx]] <- datt$biallelic.compositions
      biallelic.substitutions[[xx]] <- datt$biallelic.substitutions
      minor.alleles[[xx]] <- datt$minor.alleles
      codons[[xx]] <- datt$codons
      n.valid.sites[xx] <- datt$n.valid.sites
      n.gaps[xx] <- length(datt$sites.with.gaps)
      n.unknowns[xx] <- length(datt$sites.with.unknowns)
      n.polyallelic.sites[xx] <- length(datt$polyallelic.sites)
      n.biallelic.sites[xx] <- length(datt$biallelic.sites)
      trans.transv.ratio[xx] <- datt$trans.transv.ratio
      Coding.region[xx] <- datt$Coding_region_length
      UTR.region[xx] <- datt$UTR_region_length
      Intron.region[xx] <- datt$Intron_region_length
      Exon.region[xx] <- datt$Exon_region_length
      Gene.region[xx] <- datt$Gene_region_length
    }
  }
  print("CHECKPOINT 10")
  region.data@populations <- populationsX
  region.data@populations2 <- populations2
  region.data@popmissing <- popmissing
  region.data@outgroup <- outgroupX
  region.data@outgroup2 <- outgroup2
  region.data@CodingSNPS <- CodingSNPS
  region.data@UTRSNPS <- UTRSNPS
  region.data@IntronSNPS <- IntronSNPS
  region.data@ExonSNPS <- ExonSNPS
  region.data@GeneSNPS <- GeneSNPS
  region.data@reading.frame <- reading.frame
  region.data@rev.strand <- rev.strand
  region.data@Coding.matrix <- Coding.matrix
  region.data@Coding.matrix2 <- Coding.matrix2
  region.data@UTR.matrix <- UTR.matrix
  region.data@Intron.matrix <- Intron.matrix
  region.data@Exon.matrix <- Exon.matrix
  region.data@Gene.matrix <- Gene.matrix
  region.data@transitions <- transitions
  region.data@biallelic.matrix <- biallelic.matrix
  region.data@biallelic.sites <- biallelic.sites
  region.data@biallelic.sites2 <- biallelic.sites2
  region.data@reference <- reference
  region.data@matrix_codonpos <- matrix_codonpos
  region.data@synonymous <- synonymous
  region.data@matrix_freq <- matrix_freq
  region.data@n.singletons <- n.singletons
  region.data@polyallelic.sites <- polyallelic.sites
  region.data@n.nucleotides <- n.nucleotides
  region.data@biallelic.compositions <- biallelic.compositions
  region.data@biallelic.substitutions <- biallelic.substitutions
  region.data@minor.alleles <- minor.alleles
  region.data@codons <- codons
  region.data@sites.with.gaps <- sites.with.gaps
  region.data@sites.with.unknowns <- sites.with.unknowns
  region.stats@nucleotide.diversity <- init
  region.stats@haplotype.diversity <- init
  region.stats@haplotype.counts <- init
  region.stats@minor.allele.freqs <- init
  region.stats@biallelic.structure <- init
  region.stats@linkage.disequilibrium <- init
  genome@big.data <- big.data
  names(nsites) <- liste2
  genome@n.sites <- nsites
  genome@n.sites2 <- nsites2
  genome@n.valid.sites <- n.valid.sites
  genome@n.gaps <- n.gaps
  genome@n.unknowns <- n.unknowns
  genome@n.polyallelic.sites <- n.polyallelic.sites
  genome@n.biallelic.sites <- n.biallelic.sites
  genome@trans.transv.ratio <- trans.transv.ratio
  genome@Coding.region <- Coding.region
  genome@UTR.region <- UTR.region
  genome@Intron.region <- Intron.region
  genome@Exon.region <- Exon.region
  genome@Gene.region <- Gene.region
  genome@region.data <- region.data
  genome@region.stats <- region.stats
  genome@Pop_Neutrality$calculated <- FALSE
  genome@Pop_FSTN$calculated <- FALSE
  genome@Pop_FSTH$calculated <- FALSE
  genome@Pop_MK$calculated <- FALSE
  genome@Pop_Linkage$calculated <- FALSE
  genome@Pop_Recomb$calculated <- FALSE
  genome@Pop_Slide$calculated <- FALSE
  genome@Pop_Detail$calculated <- FALSE
  if (FAST) {
    genome@Pop_Slide$calculated <- TRUE
  }
  if (GFF.BOOL) {
    genome@gff.info <- TRUE
  }
  else {
    genome@gff.info <- FALSE
  }
  genome@snp.data <- SNP.DATA
  cat("\n")
  return(genome)
}