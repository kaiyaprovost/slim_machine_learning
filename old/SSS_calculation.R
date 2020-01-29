#R script to calculate my set of SSS on the results of fastsimcoal
#Usage: Rscript SSS_calculation.R GeneticsOutputs general_dataset/deme_coords.csv 4
#Diego Alvarado-S., 22-agu-2014

args <- commandArgs(trailingOnly = F)	#pass commands from the shell toward the script and reads script path
scriptPath <- normalizePath(dirname(sub("^--file=", "", args[grep("^--file=", args)])))
workPath = args[6]      #defines the folder where the arlequin and snp files are
coordfile = sprintf('%s/%s', scriptPath, args[7])     #defines the .csv file with the demes names and coordinates
nproc = as.numeric(args[8])         #defines the number of processors to use
setwd(paste(scriptPath,workPath, sep='/'))
if (is.na(args[9])){
  startgrp = ''
} else{
  startgrp = sprintf('StartGroup%s', args[9])          #defines which startgroup to run when running in parallel    
}


### DEFINING FUNCTIONS TO RUN IN PARALLEL
#################################################################################################################################
nss_independent_SSS = function(index0, index1, scriptPath, workPath, coordfile, suffix=''){
  #AFS-based spatial summary statistics
  #index0: intial index of list of snpfiles to be run
  #index1: intial index of list of snpfiles to be run
  #scriptPath: directory where script is
  #workPath: directory where the files are found
  #xcoordfile: name of file with coordinates of the samples the demes being analyzed
  #suffix: if an output file is desired, what suffix should it have
  
  #source function to calculate SNPs sum stats
  source(paste(scriptPath,'spatial_sumstats_functions.R',sep='/'))
  #getting demes coordinates and distances
  xy.samples = read.csv(coordfile, header=T)
  npops = nrow(xy.samples)
  nsamples = mean(xy.samples[,4])     #note that this assumes all populations have the same sample size
  
  #function to catch errors
  tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
  }
  
  #getting spatial sumstats
  if (index0!=1){
    index0 = as.numeric(index0) + 1
  }
  
  snp.files = list.files(paste(scriptPath,workPath,sep='/'), full.names=T, pattern='.snp')
  snp.files = snp.files[grep(startgrp, snp.files)]
  snp.files = snp.files[as.numeric(index0):as.numeric(index1)]
  num.AFS.SSS = 4 + (((npops*2)+((((npops^2)-npops)/2)*2)) + 2)  #4 for the AFS-based stats and the rest for the sPCA-based SSS (with the last 2 corresponding to the significant values of the global and local tests)
  AFS.SSS = matrix(rep(NA,length(snp.files)*(1+num.AFS.SSS)),ncol=(1+num.AFS.SSS))
  
  #opening a connection to a file to write line by line
  if (suffix!=''){
    file.create(paste(scriptPath, '/', workPath, '/', 'NSS-independent_SSS-', suffix, '.tsv', sep=''))
    fileConn = paste(scriptPath, '/', workPath, '/', 'NSS-independent_SSS-', suffix, '.tsv', sep='')
  } else {
    file.create(paste(scriptPath, '/', workPath, '/', 'NSS-independent_SSS', '.tsv', sep=''))
    fileConn = paste(scriptPath, '/', workPath, '/', 'NSS-independent_SSS', '.tsv', sep='')
  }
  
  for (i in 1:length(snp.files)){
    file.name = tail(strsplit(snp.files[i],'/')[[1]], n=1)
    var.name = paste('S', sub('.snp','',sub(paste('COR-', workPath, '_', sep=''),'',file.name)), sep='')
    table  = read.table(snp.files[i],sep='\t',header=F)
    table = table[,-2]      #removes the column of haplotype frequencies
    popnums = NULL
    for (j in 1:length(table[,1])){
      popnums = as.numeric(c(popnums,as.numeric(strsplit(as.character(table[j,1]),'_')[[1]][1])))
    }
    assign(var.name, cbind(popnums, table))
    cat (sub('_GeneSamples', '', file.name))
    AFS.SSS[i, 1] =  sub('.snp', '', sub('_GeneSamples', '', file.name))
    
    AFS.stats = tryCatch.W.E(stats.AFS(get(var.name), 'sss', nsamples, xy.samples[,2:3]))
    if (!is.null(AFS.stats$message)){
      cat(paste('AFS.stats calculation failed:\n\t', AFS.stats$warning))
      AFS.SSS[i,2:5] = rep(NaN,4)
    } else{ if (is.vector(AFS.stats$value)){
      AFS.SSS[i,2:5] = AFS.stats$value
    } else{
      cat(paste('AFS.stats calculation failed:\n\t', AFS.stats$value))
      AFS.SSS[i,2:5] = rep(NaN,4)
    }}
    
    sPCA.stats = tryCatch.W.E(sPCA.dist(get(var.name), xy.samples[,2:3], nsamples, cpos=2, cneg=2, plot=F))
    if (!is.null(sPCA.stats$message)){
      cat(paste('sPCA.stats calculation failed:\n\t', sPCA.stats$warning))
      AFS.SSS[i, 6:ncol(AFS.SSS)] = rep(NaN,(ncol(AFS.SSS)-5))
    } else{ if (is.vector(sPCA.stats$value)){
      AFS.SSS[i, 6:ncol(AFS.SSS)] = sPCA.stats$value
    } else{
      cat(paste('sPCA.stats calculation failed:\n\t', sPCA.stats$value))
      AFS.SSS[i, 6:ncol(AFS.SSS)] = rep(NaN,(ncol(AFS.SSS)-5))
    }}
    if (i == 1 & suffix!=''){
      write(paste(c('SimN', names(AFS.stats$value), names(sPCA.stats$value)), collapse='\t'), file=fileConn, append=T)
    }
    if (suffix!=''){
      write(paste(AFS.SSS[i,], collapse='\t'), file=fileConn, append=T)
    }}
  colnames(AFS.SSS) = c('SimN', names(AFS.stats$value), names(sPCA.stats$value))
  return(AFS.SSS)
}
#################################################################################################################################




#################################################################################################################################
nss_dependent_SSS = function(index0, index1, scriptPath, workPath, coordfile, suffix=''){
  #SSS based on NSS
  #index0: intial index of list of nssfile to be run
  #index1: intial index of list of nssfile to be run
  #scriptPath: directory where scripts are
  #workPath: directory where the NSS file is found
  #coordfile: name of the file with coordinates of the samples the demes being analyzed
  #suffix: if an output file is desired, what suffix should it have
  
  #source function to calculate SNPs sum stats
  source(paste(scriptPath,'spatial_sumstats_functions.R',sep='/'))
  
  #getting demes coordinates and distances
  xy.samples = read.csv(coordfile, header=T)
  nsamples = nrow(xy.samples)    
  
  #reading the arlequin output file
  if (index0!=1){
    index0 = as.numeric(index0) + 1
  }
  NSS.files = list.files(paste(scriptPath,workPath,sep='/'), full.names=T, pattern='outSumStats')
  NSS.files = NSS.files[grep(startgrp,NSS.files)]
  NSS.files = NSS.files[as.numeric(index0):as.numeric(index1)]
  
  #opening a connection to a file to write line by line
  if (suffix!=''){
    file.create(paste(scriptPath, '/', workPath, '/', 'NSS-dependent_SSS-', suffix, '.tsv', sep=''))
    fileConn = paste(scriptPath, '/', workPath, '/', 'NSS-dependent_SSS-', suffix, '.tsv', sep='')
  } else {
    file.create(paste(scriptPath, '/', workPath, '/', 'NSS-dependent_SSS', '.tsv', sep=''))
    fileConn = paste(scriptPath, '/', workPath, '/', 'NSS-dependent_SSS', '.tsv', sep='')
  }
  
  for (N in 1:length(NSS.files)){
    
    #just to get the file headers    
    if (N == 1){
      con <- file(NSS.files[N], "r")
      stats.names = strsplit(readLines(con, n=1), split='\t')[[1]]
      close(con)
    }
    
    READ1 = parse(text=paste("read.table(sep='\t', header=F, file=pipe(\"sed -n -e'", (index0+1), ",", (index1+1), "p' ", NSS.files[N], "\"))", sep=''))
    sumstats.file = eval(READ1)
    cat(paste0('\n',tail(strsplit(NSS.files[N],'/')[[1]],n=1), ': lines ', (index0+1), '-', index1+1))
    names(sumstats.file) = stats.names
    names(sumstats.file)[(ncol(sumstats.file)-1)] = 'SimN'
    names(sumstats.file)[ncol(sumstats.file)] = 'X'
    
    #Spatial SumStats based on NSS
    popwise.stats = c('K','H','S')   #change to define which POPWISE stats to use
    pairwise.stats = c('FST','PI')  #change to define which PAIRWISE stats to use
    stats.used = c(popwise.stats,pairwise.stats)
    combs = combn(seq(nrow(xy.samples)),2)
    sss.groups = c('MOR', 'VAR', 'MON', 'MAN')
    num.SSS = (length(popwise.stats)*1) + (length(popwise.stats)*4) + (length(pairwise.stats)*8)  + (length(pairwise.stats)*4)      #1 stat per popwise stat for Moran,  4 stats per popwise stat for Variogram, 4 stats per pairwise stat for Monmonier, and 4 stats per pairwise stat for Mantel correlogram
    SSS.combined = matrix(rep(NA,length(NSS.files)*(1+num.SSS)),ncol=(1+num.SSS))
    for (l in 1:nrow(sumstats.file)){
      SSS.combined[l,1] = paste(sumstats.file[l,(ncol(sumstats.file)-1)], l, sep='_')
      SSS = NULL; sss.names = NULL; initial.col = 1
      for (y in 1:length(stats.used)){
        if (y<=length(popwise.stats)){
          data = paste(stats.used[y],seq(nrow(xy.samples)), sep='_')
          gen = NULL
          for (datum in data){
            gen = c(gen, sumstats.file[l,datum])
          }
          MOR = moran_SSS(gen, xy.samples[,c('X','Y')])
          VAR = semivar_SSS(gen, xy.samples[,c('X','Y')])
          for (i in 1:2){
            for (j in 1:length(names(get(sss.groups[i])))){
              SSS = c(SSS, get(sss.groups[i])[[j]])
              sss.names = c(sss.names, names(get(sss.groups[i]))[j])
            }}
          names(SSS)[initial.col:length(SSS)] = paste(stats.used[y], sss.names, sep='.')
          initial.col = length(SSS) + 1; sss.names = NULL
        } else { if (y>length(popwise.stats)){
          comps = NULL
          for (m1 in 1:ncol(combs)){
            comps = c(comps,paste(stats.used[y],combs[2,m1],combs[1,m1], sep='_'))
          }
          gen = matrix(rep(NA,nrow(xy.samples)^2),ncol=nrow(xy.samples))
          for (comp in comps){
            pops.comp = as.numeric(strsplit(comp,'_')[[1]][2:3])
            gen[pops.comp[1],pops.comp[2]] = sumstats.file[l,comp]
          }
          gen = as.dist(gen)
          MON = monmonier_SSS(gen, xy.samples[,c('X','Y')], floor(nrow(xy.samples)/2))
          MAN = mantel_SSS(gen, xy.samples[,c('X','Y')])
          for (i in 3:4){
            for (j in 1:length(names(get(sss.groups[i])))){
              SSS = c(SSS, get(sss.groups[i])[[j]])
              sss.names = c(sss.names, names(get(sss.groups[i]))[j])
            }}
          names(SSS)[initial.col:length(SSS)] = paste(stats.used[y], sss.names, sep='.')
          initial.col = length(SSS) + 1; sss.names = NULL
        }}}
      SSS.combined[l,2:ncol(SSS.combined)] = SSS
      if (N==1 &l == 1 & suffix!=''){
        write(paste(c(names(sumstats.file[-(ncol(sumstats.file)-1)]), 'SimN', names(SSS)), collapse='\t'), file=fileConn, append=T)    
      }
      if (suffix!=''){
        write(paste(paste(sumstats.file[l,-(ncol(sumstats.file)-1)], collapse='\t'), paste(SSS.combined[l,], collapse='\t'), sep='\t'), file=fileConn, append=T)
      }}
  }        
  colnames(SSS.combined) = c('SimN', names(SSS))
  SSS.combined = merge(sumstats.file, as.data.frame(SSS.combined, by='SimN'))
  return(SSS.combined)
}
#################################################################################################################################



### RUNNING IN PARALLEL
require(parallel)

#SSS-set1
SNP.files = list.files(paste(scriptPath,workPath,sep='/'), full.names=T, pattern='.snp')
if (nproc == 1){
  groups = c(1, length(SNP.files))
} else {
  groups = seq(1,length(SNP.files), by=length(SNP.files)%/%nproc)
  groups = c(groups[1:(length(groups)-1)], length(SNP.files))   
}

tasks1 = list()
jobs = NULL
grouping = matrix(rep(NA,(length(groups)-1)*3), ncol=3)
colnames(grouping) = c('Group', 'StartSimNum', 'StopSimNum')
for (J in 1:(length(groups)-1)){
  grouping[J,] = c(J, groups[J], groups[J+1])
  func = parse(text=paste('function() nss_independent_SSS(', groups[J], ', ', groups[J+1], ", '", scriptPath, "', '", workPath, "', '", coordfile, "', '", paste(startgrp,J,sep='-'), "')", sep=''))
  tasks1 = c(tasks1, eval(func))
  jobs = c(jobs, paste('afs.SSS.Group',J,sep=''))
}
#write.csv(grouping, paste(scriptPath, workPath, 'Simulations_split_record.csv',sep='/'), row.names=F)
names(tasks1) = jobs
cat('\nCalculation SSS-set1 started\n')
out.SSS1 = mclapply(
  tasks1,
  function(f) f(),
  mc.cores = length(tasks1)
)
cat('\nCalculation SSS-set1 finished\n')


#SSS-set2
arl.files = list.files(paste(scriptPath,workPath,sep='/'), full.names=T, pattern='outSumStats')
if (nproc == 1){
  groups = c(1, length(arl.files))
} else {
  groups = seq(1,length(arl.files), by=length(arl.files)%/%nproc)
  groups = c(groups[1:(length(groups)-1)], length(arl.files))   
}

tasks2 = list()
jobs = NULL
for (J in 1:(length(groups)-1)){
  func = parse(text=paste('function() nss_dependent_SSS(', groups[J], ', ', groups[J+1], ", '", scriptPath, "', '", workPath, "', '", coordfile, "', '", paste(startgrp,J,sep='-'), "')", sep='')) 
  tasks2 = c(tasks2, eval(func))
  jobs = c(jobs, paste('nss.SSS.Group',J,sep=''))
}
names(tasks2) = jobs
cat('\nCalculation SSS-set2 started\n')
out.SSS2 = mclapply(
  tasks2,
  function(f) f(),
  mc.cores = length(tasks2)
)
cat('\nCalculation SSS-set2 finished\n\n')

#combining output files
combine = F
if (combine){
  for (type in c('dependent', 'independent')){
    cmd = c('cat')
    for (nss in list.files(paste(scriptPath,workPath,sep='/'), full.names = T, pattern=sprintf('NSS-%s_SSS', type))){
      cmd = c(cmd, nss)
    }
    cmd = paste(c(cmd, sprintf('> %s/%s/ALL_NSS-%s_%s_1_SSS-temp.tsv', scriptPath, workPath, type, startgrp)), collapse=' ')
    system(paste(cmd, collapse=' '))
    system(sprintf('%s/rm_dup_headers.py %s/%s/ALL_NSS-%s_%s_1_SSS-temp.tsv %s/%s/ALL_NSS-%s_%s_2_SSS-temp.tsv', scriptPath, scriptPath, workPath, type, startgrp, scriptPath, workPath, type))
  }
  
  SS_set1 = read.table(paste(scriptPath, workPath, sprintf('ALL_NSS-dependent_%s_2_SSS-temp.tsv', startgrp) ,sep='/'), sep='t', header=T)
  SS_set2 = read.table(paste(scriptPath, workPath, sprintf('ALL_NSS-independent_%s_2_SSS-temp.tsv', startgrp) ,sep='/'), sep='t', header=T)
  FINAL.SS = merge(SS_set1, SS_set2, by='SimN', all=T)
  write.table(FINAL.SS, 'ALL_sum_stats.tsv', quote=F, row.names=F, sep='\t')
}