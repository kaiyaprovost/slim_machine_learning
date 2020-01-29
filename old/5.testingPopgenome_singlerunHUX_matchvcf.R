dynamic_require <- function(package) {
  if (eval(parse(text = paste("require(", package, ")"))))
    return(TRUE)
  install.packages(package)
  return(eval(parse(text = paste(
    "require(", package, ")"
  ))))
}
packages = c("PopGenome", "moments", "R.utils")
for (p in packages) {
  dynamic_require(p)
}
path = "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/"
setwd(path)
files = list.files(pattern = "subsettemp$", recursive = FALSE)
myfiles = files
length(myfiles)
#myfiles = sample(myfiles)

do_vcf = F
overwrite = F
if (do_vcf == F) {
  ##### with the MS data itself
  msfolder = path
  for (i in 1:length(myfiles)) {
    msfile = myfiles[i]
    print(paste("Starting ", i, " of ", length(myfiles), sep = ""))
    print(msfile)
    ## check the header
    lines = readLines(msfile)
    numlines = sum(lines == "//")
    header = lines[1]
    if (substr(header, 1, 8) == "slim 100") {
      print("fixing number individuals")
      header = "slim 20 1"
      lines[1] = header
      writeLines(lines, msfile)
      lines = readLines(msfile)
      header = lines[1]
    }
    if (substr(header, 1, 4) != "slim") {
      print("oops")
      lines[1] = paste("slim", header, sep = " ")
      writeLines(lines, msfile)
      numreps = as.numeric(strsplit(header, " ")[[1]][2])
      lines = readLines(msfile)
      header = lines[1]
    } else {
      numreps = as.numeric(strsplit(header, " ")[[1]][3])
    }
    if (numlines == 1) {
      print("singleline")
      print(header)
      newheader  = strsplit(header, " ")[[1]]
      newheader[length(newheader)] = 1
      lines[1] = paste(newheader, sep = " ", collapse = " ")
      print(newheader)
      writeLines(lines, msfile)
      numreps = 1
      lines = readLines(msfile)
      header = newheader
    }
    print(numreps)
    print(numlines)
    if (numreps == numlines) {
      split = strsplit(msfile, split = "/")[[1]]
      folder = paste(split[1:length(split) - 1], sep = "/", collapse = "/")
      file = split[length(split)]
      prefix = strsplit(file, split = "[.]")[[1]][1]
      if (overwrite != T &&
          file.exists(paste(path, folder, "/", prefix, ".popgenome.stats", sep = ""))) {
        print("DO NOT RUN, ALREADY DONE")
        tomove = (
          paste(
            "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/DONE/",
            msfile,
            sep = ""
          )
        )
        currently = (
          paste(
            "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/",
            msfile,
            sep = ""
          )
        )
        print("moving")
        file.rename(currently, tomove)
      } else {
        ## read in MS file
        GENOME.class = readMS(msfile)
        ## assumes you have read a locs file?
        ## honestly you don't need a locs file you can just assign them
        numinds = length(get.individuals(GENOME.class)[[1]])
        numpopstodo = 2
        perpop = numinds / numpopstodo
        pop1 = (1:((numinds) / 2))
        pop2 = (((numinds / 2) + 1):numinds)
        GENOME.class_outtable = cbind(get.sum.data(GENOME.class))
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
              detail = T,
              do.ZnS = T,
              do.WALL = T
            ),
            timeout = 30,
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
                            timeout = 30,
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
          GENOME.class <- F_ST.stats.2(GENOME.class, snn = T, Phi_ST = T)
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
            GENOME.class <-
              try(withTimeout(calc.R2(GENOME.class),
                              timeout = 30,
                              silent = T))
            ## lengthy
            linkage_matrix = as.data.frame(t(GENOME.class@region.stats@linkage.disequilibrium[[1]][[1]]))
            colnames(GENOME.class_outtable) = make.unique(colnames(GENOME.class_outtable))
            
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
          STATS.class@BDF ## BLANK WITH BIL
          GENOME.haplotype.diversity = as.data.frame(cbind((
            STATS.class@haplotype.diversity
          )))
          colnames(GENOME.haplotype.diversity) = "GENOME.haplotype.diversity" ## exactly identical to haplotype.diversity.within
          rownames(GENOME.haplotype.diversity) = rownames(GENOME.class_outtable)
          print("moments hap count")
          momentsGENOME_haplotype.counts = t(sapply(1:length(STATS.class@haplotype.counts), function(i) {
            #print(i)
            data = unlist(STATS.class@haplotype.counts[i])
            if (is.null(data)) {
              print("failed momentsGENOME_haplotype.counts")
              print(i)
              return(c(1, -100, -100, -100, -100))
            } else {
              return(all.moments(
                unlist(STATS.class@haplotype.counts[i]),
                order.max = 4
              ))
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
              return(c(1, -100, -100, -100, -100))
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
              return(c(1, -100, -100, -100, -100))
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
        print("GET FULL STATS 2POP")
        {
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
          x = try(with(withTimeout(
            linkage.stats(
              GENOME.class,
              detail = T,
              do.ZnS = T,
              do.WALL = T,
              new.populations = list(pop1, pop2)
            ),
            timeout = 30,
            silent = T
          )))
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
              timeout = 30,
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
          STATS.class@BDF ## BLANK WITH BIL
          GENOME.haplotype.diversity = as.data.frame(cbind((
            STATS.class@haplotype.diversity
          )))
          colnames(GENOME.haplotype.diversity) = "GENOME.haplotype.diversity" ## exactly identical to haplotype.diversity.within
          rownames(GENOME.haplotype.diversity) = rownames(GENOME.class_outtable)
          print("moments hap count")
          momentsGENOME_haplotype.counts = t(sapply(1:length(STATS.class@haplotype.counts), function(i) {
            #print(i)
            data = unlist(STATS.class@haplotype.counts[i])
            if (is.null(data)) {
              print("failed momentsGENOME_haplotype.counts")
              print(i)
              return(c(1, -100, -100, -100, -100))
            } else {
              return(all.moments(
                unlist(STATS.class@haplotype.counts[i]),
                order.max = 4
              ))
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
              return(c(1, -100, -100, -100, -100))
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
              return(c(1, -100, -100, -100, -100))
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
        
        print("num cols windows:")
        print(ncol(GENOME.class_outtable))
        
        print("outputting")
        write.table(
          as.matrix(GENOME.class_outtable),
          paste(path, folder, "/", prefix, ".popgenome.stats", sep = ""),
          sep = "\t",
          row.names = T,
          col.names = T
        )
        
        
        
        tomove = (
          paste(
            "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/DONE/",
            msfile,
            sep = ""
          )
        )
        currently = (
          paste(
            "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TEMPS/SUBSET/",
            msfile,
            sep = ""
          )
        )
        print("moving")
        file.rename(currently, tomove)
        
        print("#####")
        
      }
    } else {
      print("CANNOT RUN THIS FILE, NOT COMPLETE")
    }
  }
  
  print(paste("Done with ", i, " of ", length(myfiles), sep = ""))
  
  
  
}
