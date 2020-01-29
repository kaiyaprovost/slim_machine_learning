library(PopGenome)
library(moments)

args <- commandArgs(TRUE)
msfile <- as.character(args[1])

print(msfile)

## check the header
lines = readLines(msfile)

numlines = sum(lines=="//")

header = lines[1]

if (substr(header, 1, 4) != "slim") {
  #print("oops")
  lines[1] = paste("slim", header, sep = " ")
  writeLines(lines, msfile)
  numreps = as.numeric(strsplit(header," ")[[1]][2])
} else {
  numreps = as.numeric(strsplit(header," ")[[1]][3])
}

if (numreps == numlines) {

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

MS.class <- F_ST.stats(MS.class) ## not relevant with one pop?
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
MS.haplotype.diversity = as.data.frame(cbind((MS.STATS.class@haplotype.diversity)))
colnames(MS.haplotype.diversity) = "MS.haplotype.diversity" ## exactly identical to haplotype.diversity.within
rownames(MS.haplotype.diversity) = rownames(MS.class_outtable)

#print(head(MS.STATS.class@haplotype.counts)) ## large
## probability, mean, variance, skewness, kurtosis
print("moments hap count")
momentsMS_haplotype.counts = t(sapply(1:length(MS.STATS.class@haplotype.counts), function(i) {
  #print(i)
  data = unlist(MS.STATS.class@haplotype.counts[i])
  
  if(is.null(data)) {
	
	print("failed momentsMS_haplotype.counts")
	print(i)
	return(c(1,-100,-100,-100,-100))
	
  } else {
  
   return(all.moments(unlist(MS.STATS.class@haplotype.counts[i]), order.max = 4))
	
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
  
  if(is.null(data)) {
	print(i)

	print("failed momentsMS_minor.allele.freqs")
	return(c(1,-100,-100,-100,-100))
	
  } else {
  
  
   return(all.moments(unlist(MS.STATS.class@minor.allele.freqs[i]),
			  order.max = 4))
	
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
  
  if(is.null(data)) {
	print("failed momentsMS_nuc.diversity.within")
	print(i)
	return(c(1,-100,-100,-100,-100))
	
  } else {
  
  return(all.moments(unlist(MS.STATS.class@nuc.diversity.within[i]),
			  order.max = 4))
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
MS_trim_outtable = (MS_trim_outtable[rowSums(is.na(MS_trim_outtable[,])) == 0, ])

for (col in 1:ncol(MS_trim_outtable)){
  print(col)
  MS_trim_outtable[,col] = as.numeric(MS_trim_outtable[,col])
}

print("number of cols final:")
print(ncol(MS_trim_outtable))

print("outputting")
write.table(
  as.matrix(MS_trim_outtable),
  file = paste(msfile, ".popgenome.stats", sep = ""),
  quote = F,
  row.names = F,
  sep = "\t"
)
print("#####")

} else {
  print("CANNOT RUN THIS FILE, NOT COMPLETE")
}
}

