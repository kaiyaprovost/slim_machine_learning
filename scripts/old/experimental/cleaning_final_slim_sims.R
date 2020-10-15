
read_original_file=T

if(read_original_file == T){
  file=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020.txt",
                  sep="\t",header=T,
                  na.strings=c("NA","NAN","-100",""),
                  fill=T)
  df=file[ , order(names(file))]
} else {
  df=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020_CLEANED.txt",
                sep="\t",header=T,
                na.strings=c("NA","NAN","-100","","NaN","NULL"),
                fill=T)
  df=df[ , order(names(df))]
}
dim(df)





## ALREADY DONE TO CLEAN:
df=df[,(colSums(is.na(df) | df==-100))!=nrow(df)]

df$SPATIAL = 1
df$SPATIAL[grepl("NOSPAT",df$FILE) & !is.na(df$FILE)] = 0
df$SPATIAL[df$YEAR=="EMPIRICAL" & !(is.na(df$YEAR))] = "EMPIRICAL"
df$SIMDEFAULT.[df$SPATIAL==0] = "NO"


df$CHROM[is.na(df$CHROM)] = df$CHROMOSOME[is.na(df$CHROM)]
df$CHROMOSOME[is.na(df$CHROMOSOME)] = df$CHROM[is.na(df$CHROMOSOME)]
if(sum(colnames(df) %in% "CHROM",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("CHROM"))] }## CANNOT RUN TWICE

df$DEMOG[df$DEMOG=="GF"] = "GENEFLOW"
df$DEMOG[df$DEMOG=="ISO"] = "ISOLATION"
df$DEMOG[df$DEMOG=="PAN"] = "PANMIXIA"
df$DEMOG[df$DEMOG=="SEC"] = "SECCON"

df$DEMOG[df$CHROMOSOME=="SIMULATION" & is.na(df$DEMOG)] = "SIMULATION"
df$DEMOG[df$CHROMOSOME!="" & is.na(df$DEMOG)] = "EMPIRICAL"

df$DEMOG[grepl("locus",df$FILE) & is.na(df$DEMOG)] = "EMPIRICAL"
df$DEMOG[grepl("Crotalus",df$FILE) & is.na(df$DEMOG)] = "EMPIRICAL"
df$DEMOG[grepl("Lampropeltis",df$FILE) & is.na(df$DEMOG)] = "EMPIRICAL"
df$DEMOG[grepl("Pituophis",df$FILE) & is.na(df$DEMOG)] = "EMPIRICAL"


df$DEMOG[grepl("CALLED.GENO",df$FILE) & is.na(df$DEMOG)] = "EMPIRICAL"
df$DEMOG[grepl("POP-1",df$FILE) & is.na(df$DEMOG)] = "PANMIXIA"
df$DEMOG[grepl("POP-2",df$FILE) & grepl("MIGRATE-0.0-",df$FILE) & is.na(df$DEMOG)] = "ISOLATION"
df$DEMOG[grepl("MIGRATE-0.001-",df$FILE) & grepl("SECCON-0-POP-2",df$FILE) & is.na(df$DEMOG)] = "GENEFLOW"
df$DEMOG[grepl("MIGRATE-0.01-",df$FILE) & grepl("SECCON-0-POP-2",df$FILE) & is.na(df$DEMOG)] = "GENEFLOW"
df$DEMOG[grepl("MIGRATE-0.1-",df$FILE) & grepl("SECCON-0-POP-2",df$FILE) & is.na(df$DEMOG)] = "GENEFLOW"

df$DEMOG[grepl("MIGRATE-0.1-",df$FILE) & grepl("SECCON-1-POP-2",df$FILE) & is.na(df$DEMOG)] = "SECCON"
df$DEMOG[grepl("MIGRATE-0.01-",df$FILE) & grepl("SECCON-1-POP-2",df$FILE) & is.na(df$DEMOG)] = "SECCON"
df$DEMOG[grepl("MIGRATE-0.001-",df$FILE) & grepl("SECCON-1-POP-2",df$FILE) & is.na(df$DEMOG)] = "SECCON"

df$SPECIES[df$DEMOG %in% c("ISOLATION","PANMIXIA","GENEFLOW","SECCON")] = "SIMULATION"
df$SPECIES=stringr::str_replace(df$SPECIES,"-CALLED","")


df$WINDOWSIZE[is.na(df$WINDOWSIZE)] = df$WINDOW.SIZE[is.na(df$WINDOWSIZE)]
df$WINDOW.SIZE[is.na(df$WINDOW.SIZE)] = df$WINDOWSIZE[is.na(df$WINDOW.SIZE)]

if(sum(colnames(df) %in% "WINDOWSIZE",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("WINDOWSIZE"))] }## CANNOT RUN TWICE


df$WINDOW[is.na(df$WINDOW)] = df$WINDOWNUM[is.na(df$WINDOW)]
df$WINDOWNUM[is.na(df$WINDOWNUM)] = df$WINDOW[is.na(df$WINDOWNUM)]

if(sum(colnames(df) %in% "WINDOWNUM",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("WINDOWNUM"))] }## CANNOT RUN TWICE

unique(df[is.na(df$CHROMOSOME),c("CHROMOSOME","DEMOG")])

df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="GENEFLOW"] = "SIMULATION"
df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="PANMIXIA"] = "SIMULATION"
df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="ISOLATION"] = "SIMULATION"
df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="SECCON"] = "SIMULATION"
df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="SCALE1"] = "SIMULATION"
df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="SIMULATION"] = "SIMULATION"

df$IBD[df$IBD=="TRUE"] = 1
df$IBD[df$IBD=="FALSE"] = 0

table(colSums(is.na(df)))

df$IBD[is.na(df$IBD) & grepl("IBD-0-",df$FILE)] = 0
df$IBD[is.na(df$IBD) & grepl("IBD-1-",df$FILE)] = 1
df$IBD[is.na(df$IBD) & df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"


df$EFFPOPSIZE[grepl("POPSIZE-100000-",df$FILE) & is.na(df$EFFPOPSIZE)] = 100000
df$EFFPOPSIZE[grepl("POPSIZE-200000-",df$FILE) & is.na(df$EFFPOPSIZE)] = 200000
df$EFFPOPSIZE[df$EFFPOPSIZE=="400000" & !(is.na(df$EFFPOPSIZE))] = 400000
df$EFFPOPSIZE[df$EFFPOPSIZE=="200000" & !(is.na(df$EFFPOPSIZE))] = 200000
df$EFFPOPSIZE[df$EFFPOPSIZE=="100000" & !(is.na(df$EFFPOPSIZE))] = 100000
df$EFFPOPSIZE[df$EFFPOPSIZE=="4E+05" & !(is.na(df$EFFPOPSIZE))] = 400000
df$EFFPOPSIZE[df$EFFPOPSIZE=="2E+05" & !(is.na(df$EFFPOPSIZE))] = 200000
df$EFFPOPSIZE[df$EFFPOPSIZE=="1E+05" & !(is.na(df$EFFPOPSIZE))] = 100000


unique(df[,c("FILE","RUN")])

df$DEMOG[grepl("ISOLATION",df$RUN) & is.na(df$DEMOG)] = "ISOLATION"
df$DEMOG[grepl("PANMIXIA",df$RUN) & is.na(df$DEMOG)] = "PANMIXIA"
df$DEMOG[grepl("SECCON",df$RUN) & is.na(df$DEMOG)] = "SECCON"
df$DEMOG[grepl("GENEFLOW",df$RUN) & is.na(df$DEMOG)] = "GENEFLOW"


df$OVERLAP[is.na(df$OVERLAP)] = df$OVERLAP.SIZE[is.na(df$OVERLAP)]
df$OVERLAP.SIZE[is.na(df$OVERLAP.SIZE)] = df$OVERLAP[is.na(df$OVERLAP.SIZE)]

if(sum(colnames(df) %in% "OVERLAP",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("OVERLAP"))] }## CANNOT RUN TWICE


df$YEAR[df$YEAR=="6K"] = 6000
df$YEAR[df$YEAR=="21K"] = 21000
df$YEAR[df$YEAR=="120K"] = 120000
df$YEAR[df$YEAR=="1000K"] = 1000000
df$YEAR[df$YEAR=="2000K"] = 2000000
df$YEAR[df$YEAR=="4000K"] = 4000000
df$YEAR[df$YEAR=="6.00E+03"] = 6000

df$YEAR[grepl("GEN-6000-",df$FILE) & is.na(df$YEAR)] = 6000
df$YEAR[grepl("GEN-21000-",df$FILE) & is.na(df$YEAR)] = 21000
df$YEAR[grepl("GEN-120000-",df$FILE) & is.na(df$YEAR)] = 120000
df$YEAR[grepl("GEN-1000000-",df$FILE) & is.na(df$YEAR)] = 1000000
df$YEAR[grepl("GEN-2000000-",df$FILE) & is.na(df$YEAR)] = 2000000
df$YEAR[grepl("GEN-4000000-",df$FILE) & is.na(df$YEAR)] = 4000000

df$DEMOG[grepl("SCALE1",df$RUN) & is.na(df$DEMOG)] = "SCALE1"
df$DEMOG[grepl("EMPIRICAL",df$RUN) & is.na(df$DEMOG)] = "EMPIRICAL"

df$FILE[df$FILE=="FILE"] = NA
df$MUTATIONPOWER[grepl("MUT-2.21E-10-",df$FILE) & is.na(df$MUTATIONPOWER)] = -10
df$MUTATIONPOWER[grepl("MUT-2.21E-9-",df$FILE) & is.na(df$MUTATIONPOWER)] = -9
df$MUTATIONPOWER[grepl("MUT-2.21E-8-",df$FILE) & is.na(df$MUTATIONPOWER)] = -8

df$MIGRATION.RATE[grepl("MIGRATE-0.1-",df$FILE) & is.na(df$MIGRATION.RATE)] = 0.1
df$MIGRATION.RATE[grepl("MIGRATE-0.0-",df$FILE) & is.na(df$MIGRATION.RATE)] = 0
df$MIGRATION.RATE[grepl("MIGRATE-0.001-",df$FILE) & is.na(df$MIGRATION.RATE)] = 0.001
df$MIGRATION.RATE[grepl("MIGRATE-0.01-",df$FILE) & is.na(df$MIGRATION.RATE)] = 0.01


df$NAME[is.na(df$NAME)] = df$FILE[is.na(df$NAME)]
df$FILE[is.na(df$FILE)] = df$NAME[is.na(df$FILE)]

df$RUN[is.na(df$RUN)] = df$FILE[is.na(df$RUN)]
df$FILE[is.na(df$FILE)] = df$RUN[is.na(df$FILE)]

df$NAME[is.na(df$NAME)] = df$RUN[is.na(df$NAME)]
df$RUN[is.na(df$RUN)] = df$NAME[is.na(df$RUN)]

#df$NAME[is.na(df$NAME)]=""
#df$FILE[is.na(df$FILE)]=""
#df$RUN[is.na(df$RUN)]=""

if(sum(colnames(df) %in% "NAME",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NAME"))] }## CANNOT RUN TWICE


df$NUM.POP[grepl("-POP-1-",df$FILE) & is.na(df$NUM.POP)] = 1
df$NUM.POP[grepl("-POP-2-",df$FILE) & is.na(df$NUM.POP)] = 2 

df$NUMRUN[df$NUMRUN=="1.00E+00" & !(is.na(df$NUMRUN))] = 1
df$NUMRUN[df$NUMRUN=="2.00E+00" & !(is.na(df$NUMRUN))] = 2
df$NUMRUN[df$NUMRUN=="3.00E+00" & !(is.na(df$NUMRUN))] = 3
df$NUMRUN[df$NUMRUN=="4.00E+00" & !(is.na(df$NUMRUN))] = 4
df$NUMRUN[df$NUMRUN=="5.00E+00" & !(is.na(df$NUMRUN))] = 5
df$NUMRUN[df$NUMRUN=="6.00E+00" & !(is.na(df$NUMRUN))] = 6
df$NUMRUN[df$NUMRUN=="7.00E+00" & !(is.na(df$NUMRUN))] = 7
df$NUMRUN[df$NUMRUN=="8.00E+00" & !(is.na(df$NUMRUN))] = 8

df$NUMRUN[df$IBD=="0" & is.na(df$NUMRUN) & df$DEMOG=="PANMIXIA"] = 1
df$NUMRUN[df$IBD=="0" & is.na(df$NUMRUN) & df$DEMOG=="ISOLATION"] = 3
df$NUMRUN[df$IBD=="0" & is.na(df$NUMRUN) & df$DEMOG=="GENEFLOW"] = 5
df$NUMRUN[df$IBD=="0" & is.na(df$NUMRUN) & df$DEMOG=="SECCON"] =7 
df$NUMRUN[df$IBD=="1" & is.na(df$NUMRUN) & df$DEMOG=="PANMIXIA"] =2 
df$NUMRUN[df$IBD=="1" & is.na(df$NUMRUN) & df$DEMOG=="ISOLATION"] = 4
df$NUMRUN[df$IBD=="1" & is.na(df$NUMRUN) & df$DEMOG=="GENEFLOW"] = 6
df$NUMRUN[df$IBD=="1" & is.na(df$NUMRUN) & df$DEMOG=="SECCON"] = 8

df$NUMRUN[df$IBD=="0" & df$NUMRUN==0 & df$DEMOG=="PANMIXIA"] = 1
df$NUMRUN[df$IBD=="0" & df$NUMRUN==0 & df$DEMOG=="ISOLATION"] = 3
df$NUMRUN[df$IBD=="0" & df$NUMRUN==0 & df$DEMOG=="GENEFLOW"] = 5
df$NUMRUN[df$IBD=="0" & df$NUMRUN==0 & df$DEMOG=="SECCON"] =7 
df$NUMRUN[df$IBD=="1" & df$NUMRUN==0 & df$DEMOG=="PANMIXIA"] =2 
df$NUMRUN[df$IBD=="1" & df$NUMRUN==0 & df$DEMOG=="ISOLATION"] = 4
df$NUMRUN[df$IBD=="1" & df$NUMRUN==0 & df$DEMOG=="GENEFLOW"] = 6
df$NUMRUN[df$IBD=="1" & df$NUMRUN==0 & df$DEMOG=="SECCON"] = 8

df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==1 & !is.na(df$DEMOG)] = "PANMIXIA"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==2 & !is.na(df$DEMOG)] = "PANMIXIA"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==3 & !is.na(df$DEMOG)] = "ISOLATION"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==4 & !is.na(df$DEMOG)] = "ISOLATION"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==5 & !is.na(df$DEMOG)] = "GENEFLOW"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==6 & !is.na(df$DEMOG)] = "GENEFLOW"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==7 & !is.na(df$DEMOG)] = "SECCON"
df$DEMOG[df$DEMOG=="SIMULATION" & df$NUMRUN==8 & !is.na(df$DEMOG)] = "SECCON"


df$DEMOG[df$NUMRUN==1 & is.na(df$DEMOG)] = "PANMIXIA"
df$DEMOG[df$NUMRUN==2 & is.na(df$DEMOG)] = "PANMIXIA"
df$DEMOG[df$NUMRUN==3 & is.na(df$DEMOG)] = "ISOLATION"
df$DEMOG[df$NUMRUN==4 & is.na(df$DEMOG)] = "ISOLATION"
df$DEMOG[df$NUMRUN==5 & is.na(df$DEMOG)] = "GENEFLOW"
df$DEMOG[df$NUMRUN==6 & is.na(df$DEMOG)] = "GENEFLOW"
df$DEMOG[df$NUMRUN==7 & is.na(df$DEMOG)] = "SECCON"
df$DEMOG[df$NUMRUN==8 & is.na(df$DEMOG)] = "SECCON"

df$DEMOG[is.na(df$DEMOG) & df$IBD == "EMPIRICAL" & !(is.na(df$IBD))] = "EMPIRICAL"
df$IBD[is.na(df$IBD) & df$DEMOG == "EMPIRICAL" & !(is.na(df$DEMOG))] = "EMPIRICAL"
df$NUMRUN[is.na(df$NUMRUN) & df$IBD == "EMPIRICAL" & !(is.na(df$IBD))] = "EMPIRICAL"
df$NUMRUN[is.na(df$NUMRUN) & df$DEMOG == "EMPIRICAL" & !(is.na(df$DEMOG))] = "EMPRICAL"

df$NUMRUN[is.na(df$NUMRUN) & df$DEMOG == "SIMULATION" & !(is.na(df$DEMOG))] = "SIMULATION"

df$NUMRUN[!(is.na(df$NUMRUN)) & df$NUMRUN=="EMPRIICAL"] = "EMPIRICAL"
df$IBD[!(is.na(df$IBD)) & df$IBD=="EMPRIICAL"] = "EMPIRICAL"
df$DEMOG[!(is.na(df$DEMOG)) & df$DEMOG=="EMPRIICAL"] = "EMPIRICAL"

df$NUMRUN[!(is.na(df$NUMRUN)) & df$NUMRUN=="EMPRICAL"] = "EMPIRICAL"
df$IBD[!(is.na(df$IBD)) & df$IBD=="EMPRICAL"] = "EMPIRICAL"
df$DEMOG[!(is.na(df$DEMOG)) & df$DEMOG=="EMPRICAL"] = "EMPIRICAL"

df$NUMRUN[is.na(df$NUMRUN) & df$DEMOG == "SCALE1" & !(is.na(df$DEMOG))] = "SIMULATION"
df$IBD[is.na(df$IBD) & df$DEMOG == "SCALE1" & !(is.na(df$DEMOG))] = "SIMULATION"


df$NUMRUN[is.na(df$NUMRUN) & df$IBD == "1" & !(is.na(df$IBD))] = "SIMULATION"
df$NUMRUN[is.na(df$NUMRUN) & df$IBD == "0" & !(is.na(df$IBD))] = "SIMULATION"

df$DEMOG[is.na(df$DEMOG) & df$IBD == "1" & !(is.na(df$IBD))] = "SIMULATION"
df$DEMOG[is.na(df$DEMOG) & df$IBD == "0" & !(is.na(df$IBD))] = "SIMULATION"

df$TIMESTAMP[is.na(df$TIMESTAMP)]=as.character(sapply(df$RUN[is.na(df$TIMESTAMP)],FUN=function(x){strsplit(x,"-")[[1]][2]}))

df$TIMESTAMP[is.na(df$TIMESTAMP) & df$RUN == "EMPIRICAL" & !(is.na(df$RUN))] = "EMPIRICAL"
df$TIMESTAMP[is.na(df$TIMESTAMP) & df$RUN == "MERGED_EMPIRICAL_1" & !(is.na(df$RUN))] = "EMPIRICAL"
df$TIMESTAMP[is.na(df$TIMESTAMP) & df$RUN == "SIMULATION" & !(is.na(df$RUN))] = "SIMULATION"
df$TIMESTAMP[is.na(df$TIMESTAMP) & df$RUN == "MERGEDPOP_SCALE1" & !(is.na(df$RUN))] = "SIMULATION"
df$TIMESTAMP[is.na(df$TIMESTAMP) & df$RUN == "1" & !(is.na(df$RUN))] = "1"

df$SCALE[df$SCALE=="1.00E+00" & !(is.na(df$SCALE))] = 1
df$SCALE[df$SCALE=="5.00E-01" & !(is.na(df$SCALE))] = 0.5

df$SCALE[grepl("-0.02.VCF",df$FILE) & is.na(df$SCALE)] = 0.02

df$SECCON[is.na(df$SECCON) & df$RUN == "EMPIRICAL" & !(is.na(df$RUN))] = "EMPIRICAL"
df$SECCON[is.na(df$SECCON) & df$DEMOG == "SECCON" & !(is.na(df$DEMOG))] = "1"
df$SECCON[is.na(df$SECCON) & df$DEMOG == "PANMIXIA" & !(is.na(df$DEMOG))] = "0"
df$SECCON[is.na(df$SECCON) & df$DEMOG == "ISOLATION" & !(is.na(df$DEMOG))] = "0"
df$SECCON[is.na(df$SECCON) & df$DEMOG == "GENEFLOW" & !(is.na(df$DEMOG))] = "0"

df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MUTATIONPOWER == "-10" & !(is.na(df$MUTATIONPOWER))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MUTATIONPOWER == "-8" & !(is.na(df$MUTATIONPOWER))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MUTATIONPOWER == "EMPIRICAL" & !(is.na(df$MUTATIONPOWER))] = "EMPIRICAL"

df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$RECOMPOWER == "-9" & !(is.na(df$RECOMPOWER))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$RECOMPOWER == "-7" & !(is.na(df$RECOMPOWER))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$RECOMPOWER == "EMPIRICAL" & !(is.na(df$RECOMPOWER))] = "EMPIRICAL"

df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$YEAR == 1 & !(is.na(df$YEAR))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$YEAR == 2e+06 & !(is.na(df$YEAR))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$YEAR == 4e+06 & !(is.na(df$YEAR))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$YEAR == "EMPIRICAL" & !(is.na(df$YEAR))] = "EMPIRICAL"

df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$EFFPOPSIZE == 1e+05 & !(is.na(df$EFFPOPSIZE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$EFFPOPSIZE == 2e+05 & !(is.na(df$EFFPOPSIZE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$EFFPOPSIZE == "EMPIRICAL" & !(is.na(df$EFFPOPSIZE))] = "EMPIRICAL"

df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MIGRATION.RATE == 0.01 & !(is.na(df$MIGRATION.RATE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MIGRATION.RATE == 0.001 & !(is.na(df$MIGRATION.RATE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MIGRATION.RATE == "EMPIRICAL" & !(is.na(df$MIGRATION.RATE))] = "EMPIRICAL"

df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE == 0.25 & !(is.na(df$SCALE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE == 0.5 & !(is.na(df$SCALE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE == 1 & !(is.na(df$SCALE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE == 2 & !(is.na(df$SCALE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE == 4 & !(is.na(df$SCALE))] = "NO"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE == "EMPIRICAL" & !(is.na(df$SCALE))] = "EMPIRICAL"

df$YEAR[df$YEAR=="UNK" & !(is.na(df$YEAR))] = NA
df$EFFPOPSIZE[df$EFFPOPSIZE=="UNK" & !(is.na(df$EFFPOPSIZE))] = NA


df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE==0.02 & df$MUTATIONPOWER==-9 & df$RECOMPOWER==-8 & df$EFFPOPSIZE==4e+05] = "YES"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SCALE==0.02] = "YES"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$MUTATIONPOWER==-9] = "YES"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$RECOMPOWER==-8] = "YES"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$YEAR %in% c(6000,21000,120000,1000000)] = "YES"


df$MISSING = 0
df$MISSING[grepl("_MISSING0.5.",df$FILE)] = 0.5
df$MISSING[df$MISSING==0 & df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"

df$BASES = 0
df$BASES[grepl("_BASES.",df$FILE)] = 1
df$BASES[df$BASES==0 & df$DEMOG=="EMPIRICAL"] = "EMPIRICAL"



## need to do migration rate empirical
## need to do recom power empirical
# need to do scale empirical

df$MIGRATION.RATE[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$MIGRATION.RATE)] = "EMPIRICAL"
df$RECOMPOWER[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$RECOMPOWER)] = "EMPIRICAL"
df$SCALE[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$SCALE)] = "EMPIRICAL"
df$SIMDEFAULT.[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$SIMDEFAULT.)] = "EMPIRICAL"
df$MUTATIONPOWER[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$MUTATIONPOWER)] = "EMPIRICAL"
df$MISSING[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$MISSING)] = "EMPIRICAL"
df$BASES[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$BASES)] = "EMPIRICAL"
df$EFFPOPSIZE[df$CHROMOSOME!="SIMULATION" & !(is.na(df$CHROMOSOME)) & is.na(df$EFFPOPSIZE)] = "EMPIRICAL"

## chromosome, demog, effpopsize, file, ibd, mutationpower, migration.rate, num.pop, 
## numrun, run, scale, seccon, simdefault., species, timestamp, window, window.size, year

## now make sure emp files okay
#AMPHISPIZA-BILINEATA-CALLED.GENO.PSEUDONC_007897.1_TGUT_MTDNA.VCF.FIXEDCHROMS.CONVERTED_W100000_O100000_0.WINDOW.VCF

rev(unique(df$FILE[is.na(df$SPECIES)]))

df$SPECIES[df$DEMOG!="EMPIRICAL" & !(is.na(df$DEMOG)) & is.na(df$SPECIES)] = "SIMULATION"

df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("MELOZONE-FUSCA",df$FILE)] = "MELOZONE-FUSCA"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("POLIOPTILA-MELANURA",df$FILE)] = "POLIOPTILA-MELANURA"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("CAMPYLORHYNCHUS-BRUNNEICAPILLUS",df$FILE)] = "CAMPYLORHYNCHUS-BRUNNEICAPILLUS"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("TOXOSTOMA-CRISSALE",df$FILE)] = "TOXOSTOMA-CRISSALE"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("TOXOSTOMA-CURVIROSTRE",df$FILE)] = "TOXOSTOMA-CURVIROSTRE"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("VIREO-BELLII",df$FILE)] = "VIREO-BELLII"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("AMPHISPIZA-BILINEATA",df$FILE)] = "AMPHISPIZA-BILINEATA"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("AURIPARUS-FLAVICEPS",df$FILE)] = "AURIPARUS-FLAVICEPS"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("PHAINOPEPLA-NITENS",df$FILE)] = "PHAINOPEPLA-NITENS"
df$SPECIES[!(is.na(df$FILE)) & is.na(df$SPECIES) & grepl("CARDINALIS-SINUATUS",df$FILE)] = "CARDINALIS-SINUATUS"

df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("MELOZONE-FUSCA",df$RUN)] = "MELOZONE-FUSCA"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("POLIOPTILA-MELANURA",df$RUN)] = "POLIOPTILA-MELANURA"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("CAMPYLORHYNCHUS-BRUNNEICAPILLUS",df$RUN)] = "CAMPYLORHYNCHUS-BRUNNEICAPILLUS"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("TOXOSTOMA-CRISSALE",df$RUN)] = "TOXOSTOMA-CRISSALE"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("TOXOSTOMA-CURVIROSTRE",df$RUN)] = "TOXOSTOMA-CURVIROSTRE"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("VIREO-BELLII",df$RUN)] = "VIREO-BELLII"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("AMPHISPIZA-BILINEATA",df$RUN)] = "AMPHISPIZA-BILINEATA"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("AURIPARUS-FLAVICEPS",df$RUN)] = "AURIPARUS-FLAVICEPS"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("PHAINOPEPLA-NITENS",df$RUN)] = "PHAINOPEPLA-NITENS"
df$SPECIES[!(is.na(df$RUN)) & is.na(df$SPECIES) & grepl("CARDINALIS-SINUATUS",df$RUN)] = "CARDINALIS-SINUATUS"

df$FILE[df$FILE==""] = NA

df[is.na(df)] = NA
df[df=="NA"] = NA

df$TIMESTAMP[is.na(df$TIMESTAMP) & grepl("MIGRATE",df$FILE) & !(is.na(df$FILE))]=
  as.character(sapply(df$FILE[is.na(df$TIMESTAMP) & grepl("MIGRATE",df$FILE) & !(is.na(df$FILE))],FUN=function(x){strsplit(x,"-")[[1]][19]}))


## now check cols with same names
##HAP.DIVERSITY.WITHIN, .1, .2, .3, .4
df$DIV_PER_SITE=signif(as.numeric(df$DIV_PER_SITE),6)
df$DIV_PER_SITE.1=signif(as.numeric(df$DIV_PER_SITE.1),6)

df$DIV_PER_SITE[is.na(df$DIV_PER_SITE)] = df$DIV_PER_SITE.1[is.na(df$DIV_PER_SITE)]
df$DIV_PER_SITE.1[is.na(df$DIV_PER_SITE.1)] = df$DIV_PER_SITE[is.na(df$DIV_PER_SITE.1)]
if(sum(colnames(df) %in% "DIV_PER_SITE.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("DIV_PER_SITE.1"))] }## CANNOT RUN TWICE



df$HAPLOTYPE.F_ST=signif(as.numeric(df$HAPLOTYPE.F_ST),6)
df$HAPLOTYPE.F_ST.1=signif(as.numeric(df$HAPLOTYPE.F_ST.1),6)
df$HAPLOTYPE.F_ST.2=signif(as.numeric(df$HAPLOTYPE.F_ST.2),6)
df$HAPLOTYPE.F_ST.3=signif(as.numeric(df$HAPLOTYPE.F_ST.3),6)
df$HAPLOTYPE.F_ST[is.na(df$HAPLOTYPE.F_ST)] = df$HAPLOTYPE.F_ST.1[is.na(df$HAPLOTYPE.F_ST)]
df$HAPLOTYPE.F_ST[is.na(df$HAPLOTYPE.F_ST)] = df$HAPLOTYPE.F_ST.2[is.na(df$HAPLOTYPE.F_ST)]
df$HAPLOTYPE.F_ST[is.na(df$HAPLOTYPE.F_ST)] = df$HAPLOTYPE.F_ST.3[is.na(df$HAPLOTYPE.F_ST)]
df$HAPLOTYPE.F_ST.1[is.na(df$HAPLOTYPE.F_ST.1)] = df$HAPLOTYPE.F_ST[is.na(df$HAPLOTYPE.F_ST.1)]
df$HAPLOTYPE.F_ST.1[is.na(df$HAPLOTYPE.F_ST.1)] = df$HAPLOTYPE.F_ST.2[is.na(df$HAPLOTYPE.F_ST.1)]
df$HAPLOTYPE.F_ST.1[is.na(df$HAPLOTYPE.F_ST.1)] = df$HAPLOTYPE.F_ST.3[is.na(df$HAPLOTYPE.F_ST.1)]
df$HAPLOTYPE.F_ST.2[is.na(df$HAPLOTYPE.F_ST.2)] = df$HAPLOTYPE.F_ST[is.na(df$HAPLOTYPE.F_ST.2)]
df$HAPLOTYPE.F_ST.2[is.na(df$HAPLOTYPE.F_ST.2)] = df$HAPLOTYPE.F_ST.1[is.na(df$HAPLOTYPE.F_ST.2)]
df$HAPLOTYPE.F_ST.2[is.na(df$HAPLOTYPE.F_ST.2)] = df$HAPLOTYPE.F_ST.3[is.na(df$HAPLOTYPE.F_ST.2)]
df$HAPLOTYPE.F_ST.3[is.na(df$HAPLOTYPE.F_ST.3)] = df$HAPLOTYPE.F_ST[is.na(df$HAPLOTYPE.F_ST.3)]
df$HAPLOTYPE.F_ST.3[is.na(df$HAPLOTYPE.F_ST.3)] = df$HAPLOTYPE.F_ST.1[is.na(df$HAPLOTYPE.F_ST.3)]
df$HAPLOTYPE.F_ST.3[is.na(df$HAPLOTYPE.F_ST.3)] = df$HAPLOTYPE.F_ST.2[is.na(df$HAPLOTYPE.F_ST.3)]
if(sum(colnames(df) %in% "HAPLOTYPE.F_ST.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HAPLOTYPE.F_ST.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HAPLOTYPE.F_ST.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HAPLOTYPE.F_ST.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HAPLOTYPE.F_ST.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HAPLOTYPE.F_ST.3"))] }## CANNOT RUN TWICE


df$HUDSON.G_ST=signif(as.numeric(df$HUDSON.G_ST),6)
df$HUDSON.G_ST.1=signif(as.numeric(df$HUDSON.G_ST.1),6)
df$HUDSON.G_ST.2=signif(as.numeric(df$HUDSON.G_ST.2),6)
df$HUDSON.G_ST.3=signif(as.numeric(df$HUDSON.G_ST.3),6)
df$HUDSON.G_ST[is.na(df$HUDSON.G_ST)] = df$HUDSON.G_ST.1[is.na(df$HUDSON.G_ST)]
df$HUDSON.G_ST[is.na(df$HUDSON.G_ST)] = df$HUDSON.G_ST.2[is.na(df$HUDSON.G_ST)]
df$HUDSON.G_ST[is.na(df$HUDSON.G_ST)] = df$HUDSON.G_ST.3[is.na(df$HUDSON.G_ST)]
df$HUDSON.G_ST.1[is.na(df$HUDSON.G_ST.1)] = df$HUDSON.G_ST[is.na(df$HUDSON.G_ST.1)]
df$HUDSON.G_ST.1[is.na(df$HUDSON.G_ST.1)] = df$HUDSON.G_ST.2[is.na(df$HUDSON.G_ST.1)]
df$HUDSON.G_ST.1[is.na(df$HUDSON.G_ST.1)] = df$HUDSON.G_ST.3[is.na(df$HUDSON.G_ST.1)]
df$HUDSON.G_ST.2[is.na(df$HUDSON.G_ST.2)] = df$HUDSON.G_ST[is.na(df$HUDSON.G_ST.2)]
df$HUDSON.G_ST.2[is.na(df$HUDSON.G_ST.2)] = df$HUDSON.G_ST.1[is.na(df$HUDSON.G_ST.2)]
df$HUDSON.G_ST.2[is.na(df$HUDSON.G_ST.2)] = df$HUDSON.G_ST.3[is.na(df$HUDSON.G_ST.2)]
df$HUDSON.G_ST.3[is.na(df$HUDSON.G_ST.3)] = df$HUDSON.G_ST[is.na(df$HUDSON.G_ST.3)]
df$HUDSON.G_ST.3[is.na(df$HUDSON.G_ST.3)] = df$HUDSON.G_ST.1[is.na(df$HUDSON.G_ST.3)]
df$HUDSON.G_ST.3[is.na(df$HUDSON.G_ST.3)] = df$HUDSON.G_ST.2[is.na(df$HUDSON.G_ST.3)]
if(sum(colnames(df) %in% "HUDSON.G_ST.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.G_ST.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HUDSON.G_ST.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.G_ST.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HUDSON.G_ST.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.G_ST.3"))] }## CANNOT RUN TWICE


df$HUDSON.H_ST=signif(as.numeric(df$HUDSON.H_ST),6)
df$HUDSON.H_ST.1=signif(as.numeric(df$HUDSON.H_ST.1),6)
df$HUDSON.H_ST.2=signif(as.numeric(df$HUDSON.H_ST.2),6)
df$HUDSON.H_ST.3=signif(as.numeric(df$HUDSON.H_ST.3),6)
df$HUDSON.H_ST[is.na(df$HUDSON.H_ST)] = df$HUDSON.H_ST.1[is.na(df$HUDSON.H_ST)]
df$HUDSON.H_ST[is.na(df$HUDSON.H_ST)] = df$HUDSON.H_ST.2[is.na(df$HUDSON.H_ST)]
df$HUDSON.H_ST[is.na(df$HUDSON.H_ST)] = df$HUDSON.H_ST.3[is.na(df$HUDSON.H_ST)]
df$HUDSON.H_ST.1[is.na(df$HUDSON.H_ST.1)] = df$HUDSON.H_ST[is.na(df$HUDSON.H_ST.1)]
df$HUDSON.H_ST.1[is.na(df$HUDSON.H_ST.1)] = df$HUDSON.H_ST.2[is.na(df$HUDSON.H_ST.1)]
df$HUDSON.H_ST.1[is.na(df$HUDSON.H_ST.1)] = df$HUDSON.H_ST.3[is.na(df$HUDSON.H_ST.1)]
df$HUDSON.H_ST.2[is.na(df$HUDSON.H_ST.2)] = df$HUDSON.H_ST[is.na(df$HUDSON.H_ST.2)]
df$HUDSON.H_ST.2[is.na(df$HUDSON.H_ST.2)] = df$HUDSON.H_ST.1[is.na(df$HUDSON.H_ST.2)]
df$HUDSON.H_ST.2[is.na(df$HUDSON.H_ST.2)] = df$HUDSON.H_ST.3[is.na(df$HUDSON.H_ST.2)]
df$HUDSON.H_ST.3[is.na(df$HUDSON.H_ST.3)] = df$HUDSON.H_ST[is.na(df$HUDSON.H_ST.3)]
df$HUDSON.H_ST.3[is.na(df$HUDSON.H_ST.3)] = df$HUDSON.H_ST.1[is.na(df$HUDSON.H_ST.3)]
df$HUDSON.H_ST.3[is.na(df$HUDSON.H_ST.3)] = df$HUDSON.H_ST.2[is.na(df$HUDSON.H_ST.3)]
if(sum(colnames(df) %in% "HUDSON.H_ST.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.H_ST.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HUDSON.H_ST.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.H_ST.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HUDSON.H_ST.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.H_ST.3"))] }## CANNOT RUN TWICE


df$KURT_HAPLOTYPE.COUNTS=signif(as.numeric(df$KURT_HAPLOTYPE.COUNTS),6)
df$SKEW_HAPLOTYPE.COUNTS=signif(as.numeric(df$SKEW_HAPLOTYPE.COUNTS),6)



df$KURT_NUC.DIVERSITY.WITHIN=signif(as.numeric(df$KURT_NUC.DIVERSITY.WITHIN),6)
df$KURT_NUC.DIVERSITY.WITHIN.1=signif(as.numeric(df$KURT_NUC.DIVERSITY.WITHIN.1),6)

df$KURT_NUC.DIVERSITY.WITHIN[is.na(df$KURT_NUC.DIVERSITY.WITHIN)] = df$KURT_NUC.DIVERSITY.WITHIN.1[is.na(df$KURT_NUC.DIVERSITY.WITHIN)]
df$KURT_NUC.DIVERSITY.WITHIN.1[is.na(df$KURT_NUC.DIVERSITY.WITHIN.1)] = df$KURT_NUC.DIVERSITY.WITHIN[is.na(df$KURT_NUC.DIVERSITY.WITHIN.1)]
if(sum(colnames(df) %in% "KURT_NUC.DIVERSITY.WITHIN.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("KURT_NUC.DIVERSITY.WITHIN.1"))] }## CANNOT RUN TWICE


df$MEAN_NUC.DIVERSITY.WITHIN=signif(as.numeric(df$MEAN_NUC.DIVERSITY.WITHIN),6)
df$MEAN_NUC.DIVERSITY.WITHIN.1=signif(as.numeric(df$MEAN_NUC.DIVERSITY.WITHIN.1),6)
df$MEAN_NUC.DIVERSITY.WITHIN[is.na(df$MEAN_NUC.DIVERSITY.WITHIN)] = df$MEAN_NUC.DIVERSITY.WITHIN.1[is.na(df$MEAN_NUC.DIVERSITY.WITHIN)]
df$MEAN_NUC.DIVERSITY.WITHIN.1[is.na(df$MEAN_NUC.DIVERSITY.WITHIN.1)] = df$MEAN_NUC.DIVERSITY.WITHIN[is.na(df$MEAN_NUC.DIVERSITY.WITHIN.1)]
if(sum(colnames(df) %in% "MEAN_NUC.DIVERSITY.WITHIN.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("MEAN_NUC.DIVERSITY.WITHIN.1"))] }## CANNOT RUN TWICE


df$MS.HAPLOTYPE.DIVERSITY=signif(as.numeric(df$MS.HAPLOTYPE.DIVERSITY),6)
df$MS.HAPLOTYPE.DIVERSITY2=signif(as.numeric(df$MS.HAPLOTYPE.DIVERSITY2),6)
df$MS.HAPLOTYPE.DIVERSITY[is.na(df$MS.HAPLOTYPE.DIVERSITY)] = df$MS.HAPLOTYPE.DIVERSITY2[is.na(df$MS.HAPLOTYPE.DIVERSITY)]
df$MS.HAPLOTYPE.DIVERSITY2[is.na(df$MS.HAPLOTYPE.DIVERSITY2)] = df$MS.HAPLOTYPE.DIVERSITY[is.na(df$MS.HAPLOTYPE.DIVERSITY2)]
if(sum(colnames(df) %in% "MS.HAPLOTYPE.DIVERSITY2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("MS.HAPLOTYPE.DIVERSITY2"))] }## CANNOT RUN TWICE


df$NEI.G_ST=signif(as.numeric(df$NEI.G_ST),6)
df$NEI.G_ST.1=signif(as.numeric(df$NEI.G_ST.1),6)
df$NEI.G_ST.2=signif(as.numeric(df$NEI.G_ST.2),6)
df$NEI.G_ST.3=signif(as.numeric(df$NEI.G_ST.3),6)
df$NEI.G_ST[is.na(df$NEI.G_ST)] = df$NEI.G_ST.1[is.na(df$NEI.G_ST)]
df$NEI.G_ST[is.na(df$NEI.G_ST)] = df$NEI.G_ST.2[is.na(df$NEI.G_ST)]
df$NEI.G_ST[is.na(df$NEI.G_ST)] = df$NEI.G_ST.3[is.na(df$NEI.G_ST)]
df$NEI.G_ST.1[is.na(df$NEI.G_ST.1)] = df$NEI.G_ST[is.na(df$NEI.G_ST.1)]
df$NEI.G_ST.1[is.na(df$NEI.G_ST.1)] = df$NEI.G_ST.2[is.na(df$NEI.G_ST.1)]
df$NEI.G_ST.1[is.na(df$NEI.G_ST.1)] = df$NEI.G_ST.3[is.na(df$NEI.G_ST.1)]
df$NEI.G_ST.2[is.na(df$NEI.G_ST.2)] = df$NEI.G_ST[is.na(df$NEI.G_ST.2)]
df$NEI.G_ST.2[is.na(df$NEI.G_ST.2)] = df$NEI.G_ST.1[is.na(df$NEI.G_ST.2)]
df$NEI.G_ST.2[is.na(df$NEI.G_ST.2)] = df$NEI.G_ST.3[is.na(df$NEI.G_ST.2)]
df$NEI.G_ST.3[is.na(df$NEI.G_ST.3)] = df$NEI.G_ST[is.na(df$NEI.G_ST.3)]
df$NEI.G_ST.3[is.na(df$NEI.G_ST.3)] = df$NEI.G_ST.1[is.na(df$NEI.G_ST.3)]
df$NEI.G_ST.3[is.na(df$NEI.G_ST.3)] = df$NEI.G_ST.2[is.na(df$NEI.G_ST.3)]
if(sum(colnames(df) %in% "NEI.G_ST.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NEI.G_ST.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "NEI.G_ST.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NEI.G_ST.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "NEI.G_ST.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NEI.G_ST.3"))] }## CANNOT RUN TWICE


df$NUC.F_ST.VS.ALL=signif(as.numeric(df$NUC.F_ST.VS.ALL),6)
df$NUC.F_ST.VS.ALL.1=signif(as.numeric(df$NUC.F_ST.VS.ALL.1),6)
df$NUC.F_ST.VS.ALL.2=signif(as.numeric(df$NUC.F_ST.VS.ALL.2),6)
df$NUC.F_ST.VS.ALL.3=signif(as.numeric(df$NUC.F_ST.VS.ALL.3),6)
df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL)] = df$NUC.F_ST.VS.ALL.1[is.na(df$NUC.F_ST.VS.ALL)]
df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL)] = df$NUC.F_ST.VS.ALL.2[is.na(df$NUC.F_ST.VS.ALL)]
df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL)] = df$NUC.F_ST.VS.ALL.3[is.na(df$NUC.F_ST.VS.ALL)]
df$NUC.F_ST.VS.ALL.1[is.na(df$NUC.F_ST.VS.ALL.1)] = df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL.1)]
df$NUC.F_ST.VS.ALL.1[is.na(df$NUC.F_ST.VS.ALL.1)] = df$NUC.F_ST.VS.ALL.2[is.na(df$NUC.F_ST.VS.ALL.1)]
df$NUC.F_ST.VS.ALL.1[is.na(df$NUC.F_ST.VS.ALL.1)] = df$NUC.F_ST.VS.ALL.3[is.na(df$NUC.F_ST.VS.ALL.1)]
df$NUC.F_ST.VS.ALL.2[is.na(df$NUC.F_ST.VS.ALL.2)] = df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL.2)]
df$NUC.F_ST.VS.ALL.2[is.na(df$NUC.F_ST.VS.ALL.2)] = df$NUC.F_ST.VS.ALL.1[is.na(df$NUC.F_ST.VS.ALL.2)]
df$NUC.F_ST.VS.ALL.2[is.na(df$NUC.F_ST.VS.ALL.2)] = df$NUC.F_ST.VS.ALL.3[is.na(df$NUC.F_ST.VS.ALL.2)]
df$NUC.F_ST.VS.ALL.3[is.na(df$NUC.F_ST.VS.ALL.3)] = df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL.3)]
df$NUC.F_ST.VS.ALL.3[is.na(df$NUC.F_ST.VS.ALL.3)] = df$NUC.F_ST.VS.ALL.1[is.na(df$NUC.F_ST.VS.ALL.3)]
df$NUC.F_ST.VS.ALL.3[is.na(df$NUC.F_ST.VS.ALL.3)] = df$NUC.F_ST.VS.ALL.2[is.na(df$NUC.F_ST.VS.ALL.3)]
if(sum(colnames(df) %in% "NUC.F_ST.VS.ALL.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUC.F_ST.VS.ALL.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "NUC.F_ST.VS.ALL.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUC.F_ST.VS.ALL.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "NUC.F_ST.VS.ALL.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUC.F_ST.VS.ALL.3"))] }## CANNOT RUN TWICE

## SOMEWHERE IN HERE IT FAILS -- though maybe won't fail if run line by line? 
df$NUC.F_ST.VS.ALL.4=signif(as.numeric(df$NUC.F_ST.VS.ALL.4),6)
df$NUC.F_ST.VS.ALL.4[is.na(df$NUC.F_ST.VS.ALL.4)] = df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL.4)]
df$NUC.F_ST.VS.ALL[is.na(df$NUC.F_ST.VS.ALL)] = df$NUC.F_ST.VS.ALL.4[is.na(df$NUC.F_ST.VS.ALL)]
df = df[,-which(names(df) %in% c("NUC.F_ST.VS.ALL.4"))]
if(sum(colnames(df) %in% "NUC.F_ST.VS.ALL.4",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUC.F_ST.VS.ALL.4"))] }## CANNOT RUN TWICE



df$NUCLEOTIDE.F_ST=signif(as.numeric(df$NUCLEOTIDE.F_ST),6)
df$NUCLEOTIDE.F_ST.1=signif(as.numeric(df$NUCLEOTIDE.F_ST.1),6)
df$NUCLEOTIDE.F_ST.2=signif(as.numeric(df$NUCLEOTIDE.F_ST.2),6)
df$NUCLEOTIDE.F_ST.3=signif(as.numeric(df$NUCLEOTIDE.F_ST.3),6)
df$NUCLEOTIDE.F_ST[is.na(df$NUCLEOTIDE.F_ST)] = df$NUCLEOTIDE.F_ST.1[is.na(df$NUCLEOTIDE.F_ST)]
df$NUCLEOTIDE.F_ST[is.na(df$NUCLEOTIDE.F_ST)] = df$NUCLEOTIDE.F_ST.2[is.na(df$NUCLEOTIDE.F_ST)]
df$NUCLEOTIDE.F_ST[is.na(df$NUCLEOTIDE.F_ST)] = df$NUCLEOTIDE.F_ST.3[is.na(df$NUCLEOTIDE.F_ST)]
df$NUCLEOTIDE.F_ST.1[is.na(df$NUCLEOTIDE.F_ST.1)] = df$NUCLEOTIDE.F_ST[is.na(df$NUCLEOTIDE.F_ST.1)]
df$NUCLEOTIDE.F_ST.1[is.na(df$NUCLEOTIDE.F_ST.1)] = df$NUCLEOTIDE.F_ST.2[is.na(df$NUCLEOTIDE.F_ST.1)]
df$NUCLEOTIDE.F_ST.1[is.na(df$NUCLEOTIDE.F_ST.1)] = df$NUCLEOTIDE.F_ST.3[is.na(df$NUCLEOTIDE.F_ST.1)]
df$NUCLEOTIDE.F_ST.2[is.na(df$NUCLEOTIDE.F_ST.2)] = df$NUCLEOTIDE.F_ST[is.na(df$NUCLEOTIDE.F_ST.2)]
df$NUCLEOTIDE.F_ST.2[is.na(df$NUCLEOTIDE.F_ST.2)] = df$NUCLEOTIDE.F_ST.1[is.na(df$NUCLEOTIDE.F_ST.2)]
df$NUCLEOTIDE.F_ST.2[is.na(df$NUCLEOTIDE.F_ST.2)] = df$NUCLEOTIDE.F_ST.3[is.na(df$NUCLEOTIDE.F_ST.2)]
df$NUCLEOTIDE.F_ST.3[is.na(df$NUCLEOTIDE.F_ST.3)] = df$NUCLEOTIDE.F_ST[is.na(df$NUCLEOTIDE.F_ST.3)]
df$NUCLEOTIDE.F_ST.3[is.na(df$NUCLEOTIDE.F_ST.3)] = df$NUCLEOTIDE.F_ST.1[is.na(df$NUCLEOTIDE.F_ST.3)]
df$NUCLEOTIDE.F_ST.3[is.na(df$NUCLEOTIDE.F_ST.3)] = df$NUCLEOTIDE.F_ST.2[is.na(df$NUCLEOTIDE.F_ST.3)]
if(sum(colnames(df) %in% "NUCLEOTIDE.F_ST.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUCLEOTIDE.F_ST.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "NUCLEOTIDE.F_ST.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUCLEOTIDE.F_ST.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "NUCLEOTIDE.F_ST.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("NUCLEOTIDE.F_ST.3"))] }## CANNOT RUN TWICE
## SOMEWHERE IN HERE IT FAILS

if(sum(colnames(df) %in% "PROB_HAPLOTYPE.COUNTS",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PROB_HAPLOTYPE.COUNTS"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "PROB_HAPLOTYPE.COUNTS.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PROB_HAPLOTYPE.COUNTS.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "PROB_MINOR.ALLELE.FREQS",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PROB_MINOR.ALLELE.FREQS"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "PROB_MINOR.ALLELE.FREQS.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PROB_MINOR.ALLELE.FREQS.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "PROB_NUC.DIVERSITY.WITHIN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PROB_NUC.DIVERSITY.WITHIN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "PROB_NUC.DIVERSITY.WITHIN.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PROB_NUC.DIVERSITY.WITHIN.1"))] }## CANNOT RUN TWICE


df$SKEW_NUC.DIVERSITY.WITHIN=signif(as.numeric(df$SKEW_NUC.DIVERSITY.WITHIN),6)
df$SKEW_NUC.DIVERSITY.WITHIN.1=signif(as.numeric(df$SKEW_NUC.DIVERSITY.WITHIN.1),6)
df$SKEW_NUC.DIVERSITY.WITHIN[is.na(df$SKEW_NUC.DIVERSITY.WITHIN)] = df$SKEW_NUC.DIVERSITY.WITHIN.1[is.na(df$SKEW_NUC.DIVERSITY.WITHIN)]
df$SKEW_NUC.DIVERSITY.WITHIN.1[is.na(df$SKEW_NUC.DIVERSITY.WITHIN.1)] = df$SKEW_NUC.DIVERSITY.WITHIN[is.na(df$SKEW_NUC.DIVERSITY.WITHIN.1)]
if(sum(colnames(df) %in% "SKEW_NUC.DIVERSITY.WITHIN.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("SKEW_NUC.DIVERSITY.WITHIN.1"))] }## CANNOT RUN TWICE

df$THETA_ACHAZ.TAJIMA=signif(as.numeric(df$THETA_ACHAZ.TAJIMA),6)
df$THETA_ACHAZ.TAJIMA.1=signif(as.numeric(df$THETA_ACHAZ.TAJIMA.1),6)
df$THETA_ACHAZ.TAJIMA.2=signif(as.numeric(df$THETA_ACHAZ.TAJIMA.2),6)

df$VAR_NUC.DIVERSITY.WITHIN=signif(as.numeric(df$VAR_NUC.DIVERSITY.WITHIN),6)
df$VAR_NUC.DIVERSITY.WITHIN.1=signif(as.numeric(df$VAR_NUC.DIVERSITY.WITHIN.1),6)
df$VAR_NUC.DIVERSITY.WITHIN[is.na(df$VAR_NUC.DIVERSITY.WITHIN)] = df$VAR_NUC.DIVERSITY.WITHIN.1[is.na(df$VAR_NUC.DIVERSITY.WITHIN)]
df$SKEW_NUC.DIVERSITY.WITHIN.1[is.na(df$VAR_NUC.DIVERSITY.WITHIN.1)] = df$VAR_NUC.DIVERSITY.WITHIN[is.na(df$VAR_NUC.DIVERSITY.WITHIN.1)]
if(sum(colnames(df) %in% "VAR_NUC.DIVERSITY.WITHIN.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("VAR_NUC.DIVERSITY.WITHIN.1"))] }## CANNOT RUN TWICE


if(sum(colnames(df) %in% "SNN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("SNN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "SNN.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("SNN.1"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "SNN.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("SNN.2"))] }## CANNOT RUN TWICE


novariance_cols=sapply(1:ncol(df),FUN=function(x){var(df[,x],na.rm=T)})==0
novariance_cols[is.na(novariance_cols)] = FALSE
summary(df[,novariance_cols])

if(sum(colnames(df) %in% "PI.3",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PI.3"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "PI.4",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("PI.4"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "N.UNKNOWNS",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("N.UNKNOWNS"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "P1_NONSYN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("P1_NONSYN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "P1_SYN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("P1_SYN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "P2_NONSYN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("P2_NONSYN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "P2_SYN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("P2_SYN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HUDSON.K_ST.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.K_ST.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "N.GAPS",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("N.GAPS"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "N.POLYALLELIC.SITES",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("N.POLYALLELIC.SITES"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "D_NONSYN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("D_NONSYN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "D_SYN",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("D_SYN"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "FISHER.P.VALUE",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("FISHER.P.VALUE"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HAP.F_ST.VS.ALL.2",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HAP.F_ST.VS.ALL.2"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HAP.F_ST.VS.ALL.4",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HAP.F_ST.VS.ALL.4"))] }## CANNOT RUN TWICE
if(sum(colnames(df) %in% "HUDSON.K_ST.1",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("HUDSON.K_ST.1"))] }## CANNOT RUN TWICE


df$CHROMOSOME[is.na(df$CHROMOSOME) & df$DEMOG=="SCALE1" & !(is.na(df$DEMOG))] = "SIMULATION"
df$MUTATIONPOWER[is.na(df$MUTATIONPOWER) & df$DEMOG=="EMPIRICAL" & !(is.na(df$MUTATIONPOWER))] = "EMPIRICAL"
df$YEAR[is.na(df$YEAR) & df$DEMOG=="EMPIRICAL" & !(is.na(df$YEAR))] = "EMPIRICAL"
df$YEAR[df$YEAR==1 & df$DEMOG=="EMPIRICAL" & !(is.na(df$YEAR))] = "EMPIRICAL"


df$WINDOW[is.na(df$WINDOW) & df$CHROMOSOME=="SIMULATION" & !(is.na(df$CHROMOSOME))] = 0


## metadata if species is not sim 
df$EFFPOPSIZE[is.na(df$EFFPOPSIZE) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$MIGRATION.RATE[is.na(df$MIGRATION.RATE) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$MUTATIONPOWER[is.na(df$MUTATIONPOWER) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$NUM.POP[is.na(df$NUM.POP) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$NUMRUN[is.na(df$NUMRUN) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$YEAR[is.na(df$YEAR) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$TIMESTAMP[is.na(df$TIMESTAMP) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$RUN[is.na(df$RUN) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$SCALE[is.na(df$SCALE) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$SECCON[is.na(df$SECCON) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$RECOMPOWER[is.na(df$RECOMPOWER) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"
df$SIMDEFAULT.[is.na(df$SIMDEFAULT.) & df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES))] = "EMPIRICAL"


df$OVERLAP.SIZE[!(is.na(df$FILE)) & is.na(df$OVERLAP.SIZE) & grepl("_O20000_",df$FILE)] = 20000
df$OVERLAP.SIZE[!(is.na(df$FILE)) & is.na(df$OVERLAP.SIZE) & grepl("_O100000_",df$FILE)] = 100000
df$OVERLAP.SIZE[is.na(df$OVERLAP.SIZE) & df$SPECIES=="SIMULATION" & !(is.na(df$SPECIES))] = 100000
df$WINDOW.SIZE[!(is.na(df$FILE)) & is.na(df$WINDOW.SIZE) & grepl("_W20000_",df$FILE)] = 20000
df$WINDOW.SIZE[!(is.na(df$FILE)) & is.na(df$WINDOW.SIZE) & grepl("_W100000_",df$FILE)] = 100000
df$WINDOW.SIZE[is.na(df$WINDOW.SIZE) & df$SPECIES=="SIMULATION" & !(is.na(df$SPECIES))] = 100000

df$SCALE[is.na(df$SCALE) & df$DEMOG=="SCALE1" & !(is.na(df$DEMOG))] = 1

df[is.na(df$WINDOW),]

df$WINDOW[!(is.na(df$FILE)) & is.na(df$WINDOW) & grepl("GENO.NW",df$FILE)] = 0 

sort(colSums(is.na(df) | df==-100))

## to numeric:
## fu.f_s, fu.li.dd, fu.li.f, genome.haplotype.diversity and .1, hap.diversity.within and .1 and .2, hudson.k_st,
df$FU.F_S=signif(as.numeric(df$FU.F_S),6)
df$FU.LI.D=signif(as.numeric(df$FU.LI.D),6)
df$FU.LI.F=signif(as.numeric(df$FU.LI.F),6)
df$GENOME.HAPLOTYPE.DIVERSITY=signif(as.numeric(df$GENOME.HAPLOTYPE.DIVERSITY),6)
df$GENOME.HAPLOTYPE.DIVERSITY.1=signif(as.numeric(df$GENOME.HAPLOTYPE.DIVERSITY.1),6)
df$HAP.DIVERSITY.WITHIN=signif(as.numeric(df$HAP.DIVERSITY.WITHIN),6)
df$HAP.DIVERSITY.WITHIN.1=signif(as.numeric(df$HAP.DIVERSITY.WITHIN.1),6)
df$HAP.DIVERSITY.WITHIN.2=signif(as.numeric(df$HAP.DIVERSITY.WITHIN.2),6)
df$HUDSON.K_ST=signif(as.numeric(df$HUDSON.K_ST),6)
## hudson.kaplan.rm, kelly.z_ns, kurt_haplotype.counts, kurt_minor.allele.freqs, mean_haplotype.counts, 
df$HUDSON.KAPLAN.RM=signif(as.numeric(df$HUDSON.KAPLAN.RM),6)
df$KELLY.Z_NS=signif(as.numeric(df$KELLY.Z_NS),6)
df$KURT_HAPLOTYPE.COUNTS=signif(as.numeric(df$KURT_HAPLOTYPE.COUNTS),6)
df$KURT_MINOR.ALLELE.FREQS=signif(as.numeric(df$KURT_MINOR.ALLELE.FREQS),6)
df$MEAN_HAPLOTYPE.COUNTS=signif(as.numeric(df$MEAN_HAPLOTYPE.COUNTS),6)
## mean_minor.allele.freqs, n.biallelic.sites, n.segregating.sites, n.sites, nuc.diversity.within and .1 and .2, pi,
df$MEAN_MINOR.ALLELE.FREQS=signif(as.numeric(df$MEAN_MINOR.ALLELE.FREQS),6)
df$N.BIALLELIC.SITES=signif(as.numeric(df$N.BIALLELIC.SITES),6)
df$N.SEGREGATING.SITES=signif(as.numeric(df$N.SEGREGATING.SITES),6)
df$N.SITES=signif(as.numeric(df$N.SITES),6)
df$NUC.DIVERSITY.WITHIN=signif(as.numeric(df$NUC.DIVERSITY.WITHIN),6)
df$NUC.DIVERSITY.WITHIN.1=signif(as.numeric(df$NUC.DIVERSITY.WITHIN.1),6)
df$NUC.DIVERSITY.WITHIN.2=signif(as.numeric(df$NUC.DIVERSITY.WITHIN.2),6)
## pi.1, pi.2, rozas.R_2, rozas.za, rozas.zz, skew_haplotype,counts, skew_minor.allele.freqs, strobek.s, tajima.d, theta_achaz.watterson,
df$PI=signif(as.numeric(df$PI),6)
df$PI.1=signif(as.numeric(df$PI.1),6)
df$PI.2=signif(as.numeric(df$PI.2),6)
df$ROZAS.R_2=signif(as.numeric(df$ROZAS.R_2),6)
df$ROZAS.ZA=signif(as.numeric(df$ROZAS.ZA),6)
df$ROZAS.ZZ=signif(as.numeric(df$ROZAS.ZZ),6)
df$SKEW_HAPLOTYPE.COUNTS=signif(as.numeric(df$SKEW_HAPLOTYPE.COUNTS),6)
df$SKEW_MINOR.ALLELE.FREQS=signif(as.numeric(df$SKEW_MINOR.ALLELE.FREQS),6)
df$STROBECK.S=signif(as.numeric(df$STROBECK.S),6)
df$TAJIMA.D=signif(as.numeric(df$TAJIMA.D),6)
df$THETA_ACHAZ.WATTERSON=signif(as.numeric(df$THETA_ACHAZ.WATTERSON),6)
# theta_fu.li, theta_tajima, theta_waterson, trans.transv.ratio, var_haplotype.counts, var_minor.allele.freqs, wall.b, wall.q
df$THETA_FU.LI=signif(as.numeric(df$THETA_FU.LI),6)
df$THETA_TAJIMA=signif(as.numeric(df$THETA_TAJIMA),6)
df$THETA_WATTERSON=signif(as.numeric(df$THETA_WATTERSON),6)
df$TRANS.TRANSV.RATIO=signif(as.numeric(df$TRANS.TRANSV.RATIO),6)
df$VAR_HAPLOTYPE.COUNTS=signif(as.numeric(df$VAR_HAPLOTYPE.COUNTS),6)
df$VAR_MINOR.ALLELE.FREQS=signif(as.numeric(df$VAR_MINOR.ALLELE.FREQS),6)
df$WALL.B=signif(as.numeric(df$WALL.B),6)
df$WALL.Q=signif(as.numeric(df$WALL.Q),6)


## make sure get rid of -100
## STROBECK.S, ROZAS.R_2, N.VALID.SITES,FU.F_S, FU.F_S.1
df$STROBECK.S[df$STROBECK.S==-100 & !(is.na(df$STROBECK.S))] = NA
df$ROZAS.R_2[df$ROZAS.R_2==-100 & !(is.na(df$ROZAS.R_2))] = NA
df$N.VALID.SITES[df$N.VALID.SITES==-100 & !(is.na(df$N.VALID.SITES))] = NA
df$FU.F_S[df$FU.F_S==-100 & !(is.na(df$FU.F_S))] = NA
df$FU.F_S.1[df$FU.F_S.1==-100 & !(is.na(df$FU.F_S.1))] = NA

if(sum(colnames(df) %in% "N.VALID.SITES",na.rm=T) >= 1){df = df[,-which(names(df) %in% c("N.VALID.SITES"))] }## CANNOT RUN TWICE


df=df[order(df$FILE),]


head(df[grepl("TGUT_",df$FILE) & is.na(df$CHROMOSOME) & !(is.na(df$FILE)),])
df$CHROMOSOME[grepl("TGUT_",df$FILE) & is.na(df$CHROMOSOME) & !(is.na(df$FILE))]=as.character(sapply(df$FILE[grepl("TGUT_",df$FILE) & is.na(df$CHROMOSOME) & !(is.na(df$FILE))],FUN=function(x){strsplit(strsplit(x,"_")[[1]][4],"\\.")[[1]][1]}))

#df$RUN[df$RUN=="SIMULATION" & df$FILE!="SIMULATION" & !(is.na(df$FILE))] = df$FILE[df$RUN=="SIMULATION" & df$FILE!="SIMULATION" & !(is.na(df$FILE))]
df$FILE[df$FILE=="SIMULATION" & df$RUN!="SIMULATION" & !(is.na(df$RUN))] <- df$RUN[df$FILE=="SIMULATION" & df$RUN!="SIMULATION" & !(is.na(df$RUN))]

#df$FILE[df$FILE=="" & df$RUN!="" & !(is.na(df$RUN))] <- df$RUN[df$FILE=="" & df$RUN!="" & !(is.na(df$RUN))]
#df$RUN[df$RUN=="" & df$FILE!="" & !(is.na(df$FILE))] = df$FILE[df$RUN=="" & df$FILE!="" & !(is.na(df$FILE))]


df$WINDOW[df$FILE==1 & df$RUN==1 & is.na(df$WINDOW)]=0
df$YEAR[df$FILE==1 & df$RUN==1 & is.na(df$YEAR)]=6000
df$TIMESTAMP[df$FILE==1 & df$RUN==1 & is.na(df$TIMESTAMP)]=1
df$SIMDEFAULT.[df$FILE==1 & df$RUN==1 & is.na(df$SIMDEFAULT.)]="YES"
df$SECCON[df$FILE==1 & df$RUN==1 & is.na(df$SECCON)]=0
df$RECOMPOWER[df$FILE==1 & df$RUN==1 & is.na(df$RECOMPOWER)]=-8
df$SCALE[df$FILE==1 & df$RUN==1 & is.na(df$SCALE)]=0.02
df$NUM.POP[df$FILE==1 & df$RUN==1 & is.na(df$NUM.POP)]=2
df$MUTATIONPOWER[df$FILE==1 & df$RUN==1 & is.na(df$MUTATIONPOWER)]=-9
df$MISSING[df$FILE==1 & df$RUN==1 & is.na(df$MISSING)]=0
df$EFFPOPSIZE[df$FILE==1 & df$RUN==1 & is.na(df$EFFPOPSIZE)]=4e+05
df$CHROMOSOME[df$FILE==1 & df$RUN==1 & is.na(df$CHROMOSOME)]="SIMULATION"
df$BASES[df$FILE==1 & df$RUN==1 & is.na(df$BASES)]=0
df$OVERLAP.SIZE[df$FILE==1 & df$RUN==1 & is.na(df$OVERLAP.SIZE)]=1e+5
df$WINDOW.SIZE[df$FILE==1 & df$RUN==1 & is.na(df$WINDOW.SIZE)]=1e+5

df$WINDOW[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$WINDOW)]=0
df$WINDOW.SIZE[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$WINDOW.SIZE)]=1e+5
df$TIMESTAMP[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$TIMESTAMP)]="SIMULATION"
df$RECOMPOWER[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$RECOMPOWER)]=-8
df$OVERLAP.SIZE[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$OVERLAP.SIZE)]=1e+5
df$MUTATIONPOWER[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$MUTATIONPOWER)]=-9
df$MISSING[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$MISSING)]=0
df$EFFPOPSIZE[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$EFFPOPSIZE)]=4e+05
df$CHROMOSOME[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$CHROMOSOME)]="SIMULATION"
df$BASES[df$FILE=="SIMULATION" & df$RUN=="SIMULATION" & is.na(df$BASES)]=0


metadata_columns = c("CHROMOSOME","DEMOG","EFFPOPSIZE","FILE","IBD","MIGRATION.RATE","MUTATIONPOWER","NUM.POP","NUMRUN","OVERLAP.SIZE","RECOMPOWER","RUN",
                     "SCALE","SECCON","SIMDEFAULT.","SPECIES","TIMESTAMP","WINDOW","WINDOW.SIZE","YEAR","MISSING","BASES")
dim(unique(df[,metadata_columns]))
dim(unique(df[,-which(names(df) %in% metadata_columns)]))

## get rid of rows with only metadata
df = df[(rowSums(is.na(df[,-which(names(df) %in% metadata_columns)]) | df[,-which(names(df) %in% metadata_columns)]==-100))!=ncol(df[,-which(names(df) %in% metadata_columns)]),]

## try aggregating by file
#v1 <- colnames(df_toagg[1:ncol(df_numericonly)])

dim(unique(df[duplicated(df[,c1]),c1]))

c1 <- c("FILE","RUN")
v1 = colnames(df)[!(colnames(df) %in% c1)]
agg2=aggregate(df[v1],by=df[c1],FUN=function(x){
  vars = unique(x)
  numvars = as.numeric(vars)
  vars=vars[vars!="" & !(is.na(vars))]
  numvars=numvars[numvars!="" & !(is.na(numvars))]
  
  #length(numvars) == length(vars)
  
  if (length(vars)==0){ 
    return(NA)
  } else if (length(vars)>= 2) {
    if(is.numeric(vars) | length(numvars) == length(vars)) {
      vars = as.numeric(vars)
      ## check the dif between max and min. if small enough, just return the mean
      if(max(vars,na.rm=T)-min(vars,na.rm=T) <  1e-8) { 
        return(mean(vars,na.rm=T)) 
      } else {
        if(length(vars)>=10){ 
          return("TOO VARIABLE") 
        } else { 
          return(paste(sort(vars),sep=";",collapse=";")) 
        }
      }
    } else {
      return(paste(sort(vars),sep=";")) 
    }
  } else {return(vars)}
})
head(agg2[3:nrow(agg2),])
#agg=aggregate(df[v1],by=df[c1],FUN=function(x){mean(as.numeric(x),na.rm=T)})

x=rowSums(agg2=="TOO VARIABLE",na.rm=T)>=1

good=agg2[!(x),]
bad=agg2[x,]

confirmgood=sapply(1:nrow(good),FUN=function(i){
  print(i)
  row=(unique(as.character(unlist(good[i,]))))
  return(sum(grepl(",",row),grepl(";",row),na.rm=T)==0)
})

better = good[confirmgood,]
worse = good[!(confirmgood),]

df = df[!(paste(df$FILE,df$RUN) %in% paste(better$FILE,better$RUN)),]
df = rbind(df,better)

worse_num=apply(worse[,3:ncol(worse)],2,FUN=function(x){as.numeric(as.character(x))})
worse_num = cbind(worse[,1:2],worse_num)

for (column_number in 3:ncol(worse_num)){
  print(column_number)
  this_column = worse_num[,column_number]
  ## check if entirely na if so skip
  num_filled=sum(!is.na(this_column),na.rm=T)
  
  if(num_filled != 0){
    for (row_number in 1:nrow(worse_num)){
      this_value = this_column[row_number]
      
      if(is.na(this_value) & !(is.na(worse[row_number,column_number]))){
        combo_value = unlist(worse[row_number,column_number])
        values = strsplit(combo_value,";")[[1]]
        numvalues = as.numeric(values)
        if (sum(is.na(values),na.rm=T)==sum(is.na(numvalues),na.rm=T)) {
          ## this is numeric
          checkdiff = (max(numvalues,na.rm=T)-min(numvalues,na.rm=T))/max(numvalues,na.rm=T)
          
          if(!(is.na(checkdiff))) {
            
            ## if less than 1%
            if (abs(checkdiff)<0.01){
              ## char difs
              if (max(nchar(numvalues),na.rm=T)-min(nchar(numvalues),na.rm=T) == 0) {
                ## average them they give the same information 
                worse[row_number,column_number] = mean(numvalues,na.rm=T)
              } else {
                ## average the ones with the max number of characters
                worse[row_number,column_number] = mean(numvalues[nchar(numvalues) == min(nchar(numvalues),na.rm=T)],na.rm=T)
              }
            }
          }
        }
      }
    }
  }
}

## check worse
confirmworse=sapply(1:nrow(worse),FUN=function(i){
  print(i)
  row=(unique(as.character(unlist(worse[i,]))))
  return(sum(grepl(",",row),grepl(";",row),na.rm=T)==0)
})

okay_again = worse[confirmworse,]
still_bad = worse[!(confirmworse),]

df = df[!(paste(df$FILE,df$RUN) %in% paste(okay_again$FILE,okay_again$RUN)),]
df = rbind(df,okay_again)

still_bad=apply(still_bad,2,as.character)
write.table(still_bad,"temp_worse_df.txt",row.names = F,sep="\t")

still_bad_manualfix = read.table("temp_worse_df_manualfix.txt",header=T,sep="\t")

df = df[!(paste(df$FILE,df$RUN) %in% paste(still_bad_manualfix$FILE,still_bad_manualfix$RUN)),]
df = rbind(df,still_bad_manualfix)

still_bad_manualfix
difficult = still_bad[!(paste(still_bad[,"FILE"],still_bad[,"RUN"]) %in% paste(still_bad_manualfix$FILE,still_bad_manualfix$RUN)),]

## KURT_HAPLOTYPE.COUNTS, SKEW_HAPLOTYPE.COUNTS, 
#y=colSums(bad=="TOO VARIABLE",na.rm=T)
#sort(y)

#worse2=apply(worse,2,as.character)

df=apply(df,2,as.character)
sort(colSums(is.na(df) | df=="-100"))
sort(colSums(is.na(df)))

df[,"N.SITES.CORRECTED"]=as.numeric(df[,"N.SITES.CORRECTED"])
df[,"N.SITES"]=as.numeric(df[,"N.SITES"])



## and now we try to fix the "1" "1" and "SIM" "SIM"
one_one = unique(df[df$RUN=="1" & df$FILE=="1",])
one_one = one_one[order(one_one$DIV_PER_SITE),]
dups = data.frame()
for(row_i in 1:nrow(one_one)){
  j_iterate_over = which(!(row_i:nrow(one_one) %in% as.numeric(unlist(dups))))
  for(row_j in j_iterate_over){
    if (row_j > row_i) {
      if(row_j %% 100 == 0){print(paste(row_i,row_j))}
      rows=one_one[c(row_i,row_j),]
      rows[rows=="NA"]=NA
      completes=rows[,complete.cases(t(rows))]
      if (nrow(unique(completes)) == 1){
        dups = rbind(dups,cbind(row_i,row_j))
      }
    }
  }
}
for (dup_row in sample(nrow(dups):1)){
  print(dup_row)
  row_1 = dups[dup_row,1]
  row_2 = dups[dup_row,2]
  for(i in 1:ncol(one_one)){
    if (is.na(one_one[row_1,i])){
      one_one[row_1,i] = one_one[row_2,i]
    }
    if (is.na(one_one[row_1,i])){
      one_one[row_2,i] = one_one[row_1,i]
    }  
  }
}
one_one = unique(one_one)
dim(one_one)

df = df[!(paste(df$FILE,df$RUN) %in% paste(one_one$FILE,one_one$RUN)),]
df = rbind(df,one_one)


sim_sim = unique(df[df$RUN=="SIMULATION" & df$FILE=="SIMULATION",])
sim_sim = sim_sim[order(sim_sim$DIV_PER_SITE),]

dups = data.frame()
for(row_i in 170:nrow(sim_sim)){
  j_iterate_over = which(!(row_i:nrow(sim_sim) %in% as.numeric(unlist(dups))))
  for(row_j in j_iterate_over){
    if (row_j > row_i) {
      if(row_j %% 100 == 0){print(paste(row_i,row_j))}
      rows=sim_sim[c(row_i,row_j),]
      rows[rows=="NA"]=NA
      completes=rows[,complete.cases(t(rows))]
      if (nrow(unique(completes)) == 1){
        dups = rbind(dups,cbind(row_i,row_j))
      }
    }
  }
}

for (dup_row in sample(nrow(dups):1)){
  print(dup_row)
  row_1 = dups[dup_row,1]
  row_2 = dups[dup_row,2]
  for(i in 1:ncol(sim_sim)){
    if (is.na(sim_sim[row_1,i])){
      sim_sim[row_1,i] = sim_sim[row_2,i]
    }
    if (is.na(sim_sim[row_1,i])){
      sim_sim[row_2,i] = sim_sim[row_1,i]
    }  
  }
}
sim_sim = unique(sim_sim)
dim(sim_sim)

df = df[!(paste(df$FILE,df$RUN) %in% paste(sim_sim$FILE,sim_sim$RUN)),]
df = rbind(df,sim_sim)

df$N.SITES.CORRECTED[is.na(df$N.SITES.CORRECTED) & df$NUMRUN!="EMPIRICAL" & !(is.na(df$N.SITES)) & !(is.na(df$NUMRUN))]=df$N.SITES[is.na(df$N.SITES.CORRECTED) & df$NUMRUN!="EMPIRICAL" & !(is.na(df$N.SITES)) & !(is.na(df$NUMRUN))]

#dim(df[df[,"N.SITES"]>=100000,])

colnames(df)
df=df[,(colSums(is.na(df) | df==-100))!=nrow(df)]
df=df[(rowSums(is.na(df) | df==-100))!=ncol(df),]
df=df[ , order(colnames(df))]

df=as.data.frame(df)

df=df[order(df$BASES,df$CHROMOSOME,df$DEMOG,df$EFFPOPSIZE,df$FILE,df$IBD,df$MIGRATION.RATE,df$MISSING,df$MUTATIONPOWER,df$NUM.POP,df$NUMRUN,
            df$OVERLAP.SIZE,df$RECOMPOWER,df$RUN,df$SCALE,df$SECCON,df$SIMDEFAULT.,df$SPECIES,df$TIMESTAMP,df$WINDOW,df$WINDOW.SIZE,df$YEAR),]

df = unique(df)
dim(df)




write.table(df,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020_CLEANED.txt")

empirical = df[df$SPECIES!="SIMULATION" & !(is.na(df$SPECIES)),]
simulated = df[df$SPECIES=="SIMULATION" & !(is.na(df$SPECIES)),]

#df_empgood = df[,(colSums(is.na(empirical)))!=nrow(empirical)]
#df_simgood = df[,(colSums(is.na(simulated)))!=nrow(simulated)]
df_allgood = df[,(colSums(is.na(empirical)))!=nrow(empirical) & (colSums(is.na(simulated)))!=nrow(simulated)]
df_allgood=df_allgood[,(colSums(is.na(df_allgood)))!=nrow(df_allgood)]
df_allgood=df_allgood[(rowSums(is.na(df_allgood)))!=ncol(df_allgood),]
df_allgood=df_allgood[ , order(colnames(df_allgood))]
df_allgood = unique(df_allgood)
write.table(df_allgood,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020_CLEANED_ALLGOOD.txt")

to_do_pca=c(#"WALL.Q","WALL.B",
  "VAR_NUC.DIVERSITY.WITHIN","VAR_MINOR.ALLELE.FREQS",
  "VAR_HAPLOTYPE.COUNTS",#"THETA_WATTERSON","THETA_TAJIMA","THETA_FU.LI",
  #"THETA_ACHAZ.WATTERSON","THETA_ACHAZ.TAJIMA",
  "TAJIMA.D",#"STROBECK.S",
  "SKEW_NUC.DIVERSITY.WITHIN.1","SKEW_NUC.DIVERSITY.WITHIN","SKEW_MINOR.ALLELE.FREQS",
  "SKEW_HAPLOTYPE.COUNTS",#"ROZAS.R_2","PI","NUC.DIVERSITY.WITHIN",
  "N.SITES","N.SEGREGATING.SITES",
  "N.BIALLELIC.SITES","MEAN_NUC.DIVERSITY.WITHIN","MEAN_MINOR.ALLELE.FREQS","MEAN_HAPLOTYPE.COUNTS",
  #"KURT_NUC.DIVERSITY.WITHIN",
  "KURT_MINOR.ALLELE.FREQS","KURT_HAPLOTYPE.COUNTS",
  "HAP.DIVERSITY.WITHIN","FU.LI.F","FU.LI.D",#"FU.F_S",
  "DIV_PER_SITE")
to_do_pca=c(to_do_pca,"SPECIES","YEAR")

df_pca = df_allgood[,to_do_pca]
df_pca=df_pca[,(colSums(is.na(df_pca)))!=nrow(df_pca)]
df_pca=df_pca[(rowSums(is.na(df_pca)))!=ncol(df_pca),]
#df_pca=df_pca[ , order(colnames(df_pca))]
df_pca = unique(df_pca)
write.table(df_pca,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020_CLEANED_ALLGOOD_PCA.txt")

cc_pca = df_pca[complete.cases(df_pca),]
cc_pca2=apply(cc_pca[,1:20],2,as.numeric)
cc_pca = cbind(cc_pca2,cc_pca[,21:22])

dim(cc_pca)

palette(c("red","blue","green","magenta","goldenrod","orange"))
pca=prcomp(cc_pca[,1:20],center=T,scale.=T)
pca_data = cbind(pca$x,cc_pca)
pca_data$color = 1
pca_data$color[pca_data$YEAR=="2e+06"] = 2
pca_data$color[pca_data$YEAR=="4e+06"] = 2
png("pcatemp.png")
#plot(pca_data$PC1[pca_data$color=="1"],pca_data$PC2[pca_data$color=="1"],
#     col="red",ylim=c(-2,5),xlim=c(-10,20))
plot(pca_data$PC1[pca_data$SPECIES=="SIMULATION"],pca_data$PC2[pca_data$SPECIES=="SIMULATION"],
     col="red",ylim=c(-2,5),xlim=c(-10,20))
#points(pca_data$PC1[pca_data$SPECIES=="SIMULATION" & pca_data$color!="1"],
#       pca_data$PC2[pca_data$SPECIES=="SIMULATION" & pca_data$color!="1"],
#     col=pca_data$color[pca_data$SPECIES=="SIMULATION" & pca_data$color!="1"],ylim=c(-2,5),xlim=c(-10,20))
points(pca_data$PC1[pca_data$SPECIES!="SIMULATION"],pca_data$PC2[pca_data$SPECIES!="SIMULATION"],
       col=rgb(0,0,0,0.2))
dev.off()


pca$rotation

## smaller effpopsizes than 4e-05, other mutpowers, scale does not seem to matter
## recompower standard, missingdata does not seem to matter. having bases matters?

sort(colSums(is.na(df_allgood[df_allgood$SPECIES=="SIMULATION",which(!(colnames(df_allgood) %in% metadata_columns))])))

## variables kept in 33:
var_33 = c("DIV_PER_SITE","FU.LI.D","FU.LI.D.1","FU.LI.F","FU.LI.F.1","HAP.DIVERSITY.WITHIN","KURT_HAPLOTYPE.COUNTS","KURT_MINOR.ALLELE.FREQS",
           "MEAN_HAPLOTYPE.COUNTS","MEAN_MINOR.ALLELE.FREQS","MEAN_NUC.DIVERSITY.WITHIN","N.BIALLELIC.SITES","N.SEGREGATING.SITES","N.SEGREGATING.SITES.1",
           "N.SITES.CORRECTED","PI","SKEW_HAPLOTYPE.COUNTS","SKEW_MINOR.ALLELE.FREQS","SKEW_NUC.DIVERSITY.WITHIN","TAJIMA.D","TAJIMA.D.1","THETA_ACHAZ.TAJIMA",
           "THETA_ACHAZ.TAJIMA.1","THETA_ACHAZ.WATTERSON","THETA_ACHAZ.WATTERSON.1","THETA_FU.LI","THETA_FU.LI.1","THETA_TAJIMA","THETA_WATTERSON",
           "THETA_WATTERSON.1","VAR_HAPLOTYPE.COUNTS","VAR_MINOR.ALLELE.FREQS","VAR_NUC.DIVERSITY.WITHIN")
## kept in 11:
var_11=c("DIV_PER_SITE","FU.LI.D","HAP.DIVERSITY.WITHIN","KURT_HAPLOTYPE.COUNTS","MEAN_HAPLOTYPE.COUNTS","N.SEGREGATING.SITES","N.SITES.CORRECTED",
         "SKEW_MINOR.ALLELE.FREQS","TAJIMA.D.1","THETA_ACHAZ.TAJIMA.1","THETA_FU.LI.1")
## metadata 
metadata_columns = c("CHROMOSOME","DEMOG","EFFPOPSIZE","FILE","IBD","MIGRATION.RATE","MUTATIONPOWER","NUM.POP","NUMRUN","OVERLAP.SIZE","RECOMPOWER","RUN",
                     "SCALE","SECCON","SIMDEFAULT.","SPECIES","TIMESTAMP","WINDOW","WINDOW.SIZE","YEAR","MISSING","BASES")




original=df[df$SIMDEFAULT.=="YES",]
original=original[as.numeric(original$TIMESTAMP)>1 & !(is.na(original$TIMESTAMP)),]
original$TIMESTAMP=as.numeric(original$TIMESTAMP)
original = original[original$TIMESTAMP<=1564172831,]
original = original[,c(var_33,metadata_columns)]
original=unique(original)
## 1551823140 = Tuesday, March 5, 2019 9:59:00 PM 
## 1558968497 = Monday, May 27, 2019 2:48:17 PM 
## 1560514797= Friday, June 14, 2019 12:19:57 PM 
## 1564172831 = , oldest 1mil runs
## 1582335852 = Saturday, February 22, 2020 1:44:12 AM
## want timestanmps that are less than 1580000000
head(original[original$YEAR=="120000" & original$NUMRUN==1,])

table(original[,c("TIMESTAMP","YEAR")])

#df = smartbind(df,dfxy)

original=original[,c(var_33,"NUMRUN","YEAR")]
original=unique(original)
table(original[,c("NUMRUN","YEAR")])

org_33 = original[complete.cases(original),]
table(org_33[,c("NUMRUN","YEAR")])

dfxy = read.table("/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/stats/MERGED_snakes_combo.txt",header=T,
                  na.strings=c("NA","NAN","-100",""),fill=T)
dfxy=dfxy[!(grepl("MISSINGFALSE",dfxy$FILE)),]
dfxy=dfxy[!(grepl("MISSINGTRUE.MONOMORPHICTRUE",dfxy$FILE)),]

dfxy_2 = gtools::smartbind(dfxy,org_33)
dfxy_2$N.SITES.CORRECTED=dfxy_2$N.SITES
dfxy = dfxy_2[1:nrow(dfxy),c(var_33)]
dfxy_2 = dfxy_2[1:nrow(dfxy),]

sum(complete.cases(dfxy))

rotation = pca$rotation

snake_pca=predict(object=pca,newdata=dfxy)
snake_pca = cbind(snake_pca,dfxy_2)



org_33_pca = org_33[,1:33]

palette(c("black","grey","red","pink","blue","lightblue","goldenrod","goldenrod1"))
pca=prcomp(org_33_pca,center=T,scale.=T)
pca_data = cbind(pca$x,org_33)
pca_data$color = 1
pca_data$color[pca_data$YEAR=="2e+06"] = 2
pca_data$color[pca_data$YEAR=="4e+06"] = 2
pdf("pcatemp_org.pdf",height=3,width=8)
par(mfrow=c(1,3))
plot(snake_pca$PC1,snake_pca$PC2,col="lightgreen",ylab="PC2",xlab="PC1")
for (i in 1:8) {
  d = pca_data[pca_data$NUMRUN==i,]
  
  car::ellipse(center = colMeans( pca_data[pca_data$NUMRUN==i,c("PC2","PC1")],na.rm = T), 
               shape = cov( pca_data[pca_data$NUMRUN==i,c("PC2","PC1")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=1,lwd=c(1,1,2,2,3,3,4,4)[i],
               radius = sqrt(qchisq(.95, df=2)),col = "black",
               lty=c(1,2,1,2,1,2,1,2)[i])
}
points(snake_pca$PC1[!(grepl("LOCUS",snake_pca$FILE))],
       snake_pca$PC2[!(grepl("LOCUS",snake_pca$FILE))],col="darkgreen")

plot(snake_pca$PC1[!(grepl("LOCUS",snake_pca$FILE))],
       snake_pca$PC2[!(grepl("LOCUS",snake_pca$FILE))],col="darkgreen",type="n",
     ylab="PC2",xlab="PC1")
points(pca_data$PC1,pca_data$PC2,col="green",cex=0.3)
for (i in 1:8) {
  d = pca_data[pca_data$NUMRUN==i,]
  
  car::ellipse(center = colMeans( pca_data[pca_data$NUMRUN==i,c("PC2","PC1")],na.rm = T), 
               shape = cov( pca_data[pca_data$NUMRUN==i,c("PC2","PC1")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=1,lwd=c(1,1,2,2,3,3,4,4)[i],
               radius = sqrt(qchisq(.95, df=2)),col = "black",
               lty=c(1,2,1,2,1,2,1,2)[i])
}
points(snake_pca$PC1[grepl("LAMP",snake_pca$FILE)],
       snake_pca$PC2[grepl("LAMP",snake_pca$FILE)],col="red",pch=0)
points(snake_pca$PC1[grepl("ATROX",snake_pca$FILE)],
       snake_pca$PC2[grepl("ATROX",snake_pca$FILE)],col="blue",pch=0)
points(snake_pca$PC1[grepl("CATE",snake_pca$FILE)],
       snake_pca$PC2[grepl("CATE",snake_pca$FILE)],col="goldenrod",pch=0)
points(snake_pca$PC1[grepl("SCUT",snake_pca$FILE)],
       snake_pca$PC2[grepl("SCUT",snake_pca$FILE)],col="magenta",pch=0)

plot(pca_data$PC1,pca_data$PC2,col="green",cex=0.3,
     ylim=c(-30,30),xlim=c(-30,30),ylab="PC2",xlab="PC1")
for (i in 1:8) {
  d = pca_data[pca_data$NUMRUN==i,]
  
  car::ellipse(center = colMeans( pca_data[pca_data$NUMRUN==i,c("PC2","PC1")],na.rm = T), 
               shape = cov( pca_data[pca_data$NUMRUN==i,c("PC2","PC1")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=1,lwd=c(1,1,2,2,3,3,4,4)[i],
               radius = sqrt(qchisq(.95, df=2)),col = "black",
               lty=c(1,2,1,2,1,2,1,2)[i])
}
points(snake_pca$PC1,snake_pca$PC2,col="lightgreen")
points(snake_pca$PC1[!(grepl("LOCUS",snake_pca$FILE))],
       snake_pca$PC2[!(grepl("LOCUS",snake_pca$FILE))],col="darkgreen")


dev.off()


View(snake_pca[!(grepl("LOCUS",snake_pca$FILE)),])
summary(org_33)

pdf("sim_snake_boxplots_df.pdf")
par(mfrow=c(2,2))
for(i in 1:ncol(df)){
  
  if(!(is.na(unique(as.numeric(df[,i])))) & colnames(df)[i] %in% colnames(dfxy) ){
    if( length(as.numeric(dfxy[,colnames(df)[i]]))>0 ){
    boxplot(list(as.numeric(df[,i]),
                 as.numeric(dfxy[,colnames(df)[i]])
                 #as.numeric(dfxy[!(grepl("LOCUS",dfxy$FILE)),colnames(df)[i]]),
                 #as.numeric(dfxy[(grepl("LOCUS",dfxy$FILE)),colnames(df)[i]])
                 ),
            ylab=colnames(df)[i],
            names=c("Sims","Snakes"))
    
   
    #  print(i)
    #  print(colnames(df)[i])
    #  test=t.test(as.numeric(df[,i]),as.numeric(dfxy[,colnames(df)[i]]))
    #  print(test$p.value)
    }
    
    
  }
  

}
dev.off()

pdf("sim_snake_hist_df.pdf",width=6,height=3)
par(mfrow=c(2,5),mar=c(4,2,2,0))
for(i in 1:ncol(df)){
  
  if(!(is.na(unique(as.numeric(df[,i])))) & colnames(df)[i] %in% colnames(dfxy) ){
    if( length(as.numeric(dfxy[,colnames(df)[i]]))>0 ){
      
      
      for(j in c(1:8,"SIMULATION","EMPIRICAL")){
        
        
        mini=min(c(as.numeric(df[,i]),as.numeric(dfxy[,colnames(df)[i]])),na.rm=T)
        maxi=max(c(as.numeric(df[,i]),as.numeric(dfxy[,colnames(df)[i]])),na.rm=T)
        
        hist((as.numeric(df[df$NUMRUN==as.character(j),i])),
             xlab=colnames(df)[i],col=as.numeric(j),main=j,ylab="",
             #xlim=c(mini,maxi),breaks=20,
             add=F)
        abline(v=mean(as.numeric(dfxy[,colnames(df)[i]]),na.rm=T),
               col="red",lty=3,lwd=3)
        
        abline(v=mean(as.numeric(dfxy[,colnames(df)[i]]),na.rm=T)-sd(as.numeric(dfxy[,colnames(df)[i]]),na.rm=T),
               col="blue",lty=3,lwd=3)
        
        abline(v=mean(as.numeric(dfxy[,colnames(df)[i]]),na.rm=T)+sd(as.numeric(dfxy[,colnames(df)[i]]),na.rm=T),
               col="blue",lty=3,lwd=3)
        
      }
      

      
      
      #  print(i)
      #  print(colnames(df)[i])
      #  test=t.test(as.numeric(df[,i]),as.numeric(dfxy[,colnames(df)[i]]))
      #  print(test$p.value)
    }
    
    
  }
  
  
}
dev.off()



## keeping snakes
keep_snake_var = c("FU.LI.D","FU.LI.F","MEAN_MINOR.ALLELE.FREQS","MEAN_NUC.DIVERSITY.WITHIN",
                   "SKEW_MINOR.ALLELE.FREQS","SKEW_NUC.DIVERSITY.WITHIN","TAJIMA.D",
                   "VAR_MINOR.ALLELE.FREQS","VAR_NUC.DIVERSITY.WITHIN")
maybe_snake_var = c("DIV_PER_SITE","HAP.DIVERSITY.WITHIN","KURT_MINOR.ALLELE.FREQS","N.BIALLELIC.SITES",
                    "N.SEGREGATING.SITES","THETA_FU.LI","THETA_TAJIMA","THETA_WATTERSON")
no_snake_var = c("KURT_HAPLOTYPE.COUNTS","MEAN_HAPLOTYPE.COUNTS","N.SITES.CORRECTED","PI",
                 "SKEW_HAPLOTYPE.COUNTS","THETA_ACHAZ.TAJIMA","THETA_ACHAZ.WATTERSON","VAR_HAPLOTYPE.COUNTS")


df_pca3 = df[,c(keep_snake_var,"NUMRUN")]
df_pca3=df_pca3[complete.cases(df_pca3),]
pca3 = prcomp(df_pca3[,-(ncol(df_pca3))],center = T,scale. = T)
pca3_data = pca3$x
pdf("numrun.snakes.pdf")
plot(pca3_data[,1:2],col="grey",cex=0.1)
points(pca3_data[df_pca3$NUMRUN %in% c(1:8),1:2],cex=0.3,
     col=as.numeric(as.factor(df_pca3$NUMRUN[df_pca3$NUMRUN %in% c(1:8)])))
for (i in 1:8) {
  demog=c(1:8)[i]
  car::ellipse(center = colMeans( pca3_data[df_pca3$NUMRUN==demog,c("PC1","PC2")],na.rm = T), 
               shape = cov( pca3_data[df_pca3$NUMRUN==demog,c("PC1","PC2")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=1,
               radius = sqrt(qchisq(.95, df=2)),col = "black")
}
dev.off()



df_pca2=df[,c(keep_snake_var,"DEMOG")]
df_pca2=df_pca2[complete.cases(df_pca2),]
pca2 = prcomp(df_pca2[,-(ncol(df_pca2))],center = T,scale. = T)
snake2 = predict(pca2,dfxy[,keep_snake_var])
snake2=cbind(snake2,dfxy)
#snake2=snake2[complete.cases(snake2),]
pca2_data = pca2$x
pdf("snake_smaller_pca.pdf")
plot(pca2_data[,1:2],cex=0.3)
points(pca2_data[df_pca2$DEMOG=="EMPIRICAL",1:2],cex=0.3,col="blue")
points(snake2[(grepl("LOCUS",snake2$FILE)),1:2],col="darkred",cex=0.6)
points(snake2[!(grepl("LOCUS",snake2$FILE)),1:2],col="red",pch=16)
legend("top",bty="n",legend=c("Sim","Birds","Snake Loci","Snake Genome"),
       col=c("black","blue","darkred","red"),
       fill=c("black","blue","darkred","red"))
dev.off()

pdf("snake_smaller_pca_data.pdf")
plot(pca2_data[df_pca2$DEMOG %in% c("GENEFLOW","ISOLATION","PANMIXIA","SECCON"),1:2],cex=0.3,
     col=as.numeric(as.factor(df_pca2$DEMOG[df_pca2$DEMOG %in% c("GENEFLOW","ISOLATION","PANMIXIA","SECCON")])))
for (i in 1:4) {
  demog=c("GENEFLOW","ISOLATION","PANMIXIA","SECCON")[i]
  car::ellipse(center = colMeans( pca2_data[df_pca2$DEMOG==demog,c("PC1","PC2")],na.rm = T), 
               shape = cov( pca2_data[df_pca2$DEMOG==demog,c("PC1","PC2")],use="pairwise.complete.obs"),
               center.pch=i,center.cex=1,
               radius = sqrt(qchisq(.95, df=2)),col = "black")
}
dev.off()


small_variables = c("DIV_PER_SITE","HAP.DIVERSITY.WITHIN","HAP.DIVERSITY.WITHIN.1",
                    "HUDSON.KAPLAN.RM","N.BIALLELIC.SITES","N.SEGREGATING.SITES",
                    "N.SITES","THETA_WATTERSON","FU.LI.D")
small_pca = df[,c(small_variables,"NUMRUN","YEAR")]
small_pca = small_pca[complete.cases(small_pca),]
small_snake = dfxy[!(grepl("LOCUS",dfxy$FILE)),c(small_variables,"FILE")]
small_snake = small_snake[complete.cases(small_snake),]

dim(small_pca)

pca_small = prcomp(small_pca[,1:9],center = T,scale. = T)
pca_small_data = pca_small$x
pca_snake = predict(pca_small,small_snake[,1:9])
plot(pca_small_data[,1:2],xlim=c(-30,65),ylim=c(-10,150)
     )
points(pca_snake[,1:2],col="red")


df_ag_33 = df_allgood[,c(metadata_columns,var_33)]
df_ag_11 = df_allgood[,c(metadata_columns,var_11)]

df_ag_33=df_ag_33[,(colSums(is.na(df_ag_33)))!=nrow(df_ag_33)]
df_ag_33=df_ag_33[(rowSums(is.na(df_ag_33)))!=ncol(df_ag_33),]
df_ag_33=df_ag_33[ , order(colnames(df_ag_33))]
df_ag_33 = unique(df_ag_33)
write.table(df_ag_33,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020_CLEANED_ALLGOOD_33.txt")

df_ag_11=df_ag_11[,(colSums(is.na(df_ag_11)))!=nrow(df_ag_11)]
df_ag_11=df_ag_11[(rowSums(is.na(df_ag_11)))!=ncol(df_ag_11),]
df_ag_11=df_ag_11[ , order(colnames(df_ag_11))]
df_ag_11 = unique(df_ag_11)
write.table(df_ag_11,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/SPECIES/STATS/MERGED_ALL_SHEETS_TOGETHER_9SEPT2020_CLEANED_ALLGOOD_11.txt")




df3=apply(df2,2,as.numeric)
summary(df3[,var_11])


cols_numeric=sapply(1:ncol(df),FUN=function(x){class(df[,x])!="character"})

cor_matrix = cor(df[,cols_numeric],use="pairwise.complete.obs")
diag(cor_matrix) = NA

cor_matrix2=cor_matrix[(colSums(is.na(cor_matrix)))!=nrow(cor_matrix),(colSums(is.na(cor_matrix)))!=nrow(cor_matrix)]
diag(cor_matrix2) = NA

for (i in 1:nrow(cor_matrix)){
  for(j in 1:ncol(cor_matrix)){
    if (i < j){
      value = cor_matrix[i,j]
      if( !(is.na(value)) & abs(value) == 1) {
        print(paste(rownames(cor_matrix)[i],colnames(cor_matrix)[j],sep=","))
      }
    }
  }
}

pdf("temp.pdf")
corrplot::corrplot(cor_matrix2,method="number",tl.cex=0.2,na.label="square",na.label.col="grey",
                   order="alphabet",cl.cex=0.2,number.cex=0.2)
dev.off()




pdf("simulations_data_differences.pdf")
par(mfrow=c(2,2),mar=c(4,4,0,0))
for(i in rev(c(1:5,7:53,55:113,115:193,195:212))){
  print("##########"); print(colnames(df)[i]); print(i)
  temp = df[,c(i,6,54,213,114)]
  temp[,1] = as.numeric(as.character(temp[,1]))
  temp = temp[temp$DEMOG!="EMPIRICAL",]
  temp = temp[complete.cases(temp),]
  
  

  
  if(length(unique(as.numeric(temp[,1]))) >1 & length(unique(as.character(temp[,2])))>1 
     & length(unique(as.character(temp[,3])))>1 & length(unique(as.character(temp[,4])))>1){
    #mod=aov(as.numeric(temp[,1])~as.character(temp[,2])+as.character(temp[,3])+as.character(temp[,4]))
    #print(summary(mod))
    #TukeyHSD(mod)
    
    boxplot(temp[,1]~temp$DEMOG,xlab="DEMOG",ylab=colnames(temp)[1])
    boxplot(temp[,1]~temp$IBD,xlab="IBD",ylab=colnames(temp)[1])
    boxplot(temp[,1]~temp$YEAR,xlab="YEAR",ylab=colnames(temp)[1])
    boxplot(temp[,1]~temp$NUMRUN,xlab="NUMRUN",ylab=colnames(temp)[1])
    
  }
}
dev.off()





