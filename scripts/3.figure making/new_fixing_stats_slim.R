x=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/MERGED_basics_bunch_of_files_together_15OCT2020.txt",
             header=T,sep=" ")
card=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/cardinalis vcf/MERGED_cardcard_12final_14oct2020.txt",
                header=T,sep="\t")
snek=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/MERGED_snakes_14oct2020.txt",
                header=T,sep="\t")

merged=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/MERGED_card_snake_sims_15OCT2020.txt",header=T,sep="\t")

lists = list(x,card,snek)
merged = plyr::rbind.fill(lists)
merged = unique(merged)
dim(merged)

## remove hudson gst .3, hst.3, kst.3
## remove pi.2, NUC.DIVERSITY.WITHIN.2, NUCLEOTIDE.F_ST.3
## remove NUC.F_ST.VS.ALL.1, NUC.F_ST.VS.ALL.3, NUC.F_ST.VS.ALL.4

##cl and cl.window are identical

cor(merged[,c("NUC.F_ST.VS.ALL","NUC.F_ST.VS.ALL.1",
              "NUC.F_ST.VS.ALL.3","NUC.F_ST.VS.ALL.4")],use="pairwise.complete.obs")

plot(merged[,c("TAJIMA.D","TAJIMA.D.1")])



merged$NUC.F_ST.VS.ALL.3[is.na(merged$NUC.F_ST.VS.ALL.3)] = merged$NUC.F_ST.VS.ALL.4[is.na(merged$NUC.F_ST.VS.ALL.3)]
merged$NUC.F_ST.VS.ALL.4[is.na(merged$NUC.F_ST.VS.ALL.4)] = merged$NUC.F_ST.VS.ALL.3[is.na(merged$NUC.F_ST.VS.ALL.4)]

merged$NUC.F_ST.VS.ALL.3[is.na(merged$NUC.F_ST.VS.ALL.3)] = merged$NUC.F_ST.VS.ALL[is.na(merged$NUC.F_ST.VS.ALL.3)]
merged$NUC.F_ST.VS.ALL[is.na(merged$NUC.F_ST.VS.ALL)] = merged$NUC.F_ST.VS.ALL.3[is.na(merged$NUC.F_ST.VS.ALL)]

merged$NUC.F_ST.VS.ALL.4[is.na(merged$NUC.F_ST.VS.ALL.4)] = merged$NUC.F_ST.VS.ALL[is.na(merged$NUC.F_ST.VS.ALL.4)]
merged$NUC.F_ST.VS.ALL[is.na(merged$NUC.F_ST.VS.ALL)] = merged$NUC.F_ST.VS.ALL.4[is.na(merged$NUC.F_ST.VS.ALL)]

merged$NUC.F_ST.VS.ALL.1[is.na(merged$NUC.F_ST.VS.ALL.1)] = merged$NUC.F_ST.VS.ALL[is.na(merged$NUC.F_ST.VS.ALL.1)]
merged$NUC.F_ST.VS.ALL[is.na(merged$NUC.F_ST.VS.ALL)] = merged$NUC.F_ST.VS.ALL.1[is.na(merged$NUC.F_ST.VS.ALL)]


merged$NUCLEOTIDE.F_ST[is.na(merged$NUCLEOTIDE.F_ST)] = merged$NUCLEOTIDE.F_ST.3[is.na(merged$NUCLEOTIDE.F_ST)]
merged$NUCLEOTIDE.F_ST.3[is.na(merged$NUCLEOTIDE.F_ST.3)] = merged$NUCLEOTIDE.F_ST[is.na(merged$NUCLEOTIDE.F_ST.3)]
merged$NEI.G_ST[is.na(merged$NEI.G_ST)] = merged$NEI.G_ST.3[is.na(merged$NEI.G_ST)]
merged$NEI.G_ST.3[is.na(merged$NEI.G_ST.3)] = merged$NEI.G_ST[is.na(merged$NEI.G_ST.3)]


merged$NUC.DIVERSITY.WITHIN[is.na(merged$NUC.DIVERSITY.WITHIN)] = merged$NUC.DIVERSITY.WITHIN.2[is.na(merged$NUC.DIVERSITY.WITHIN)]
merged$NUC.DIVERSITY.WITHIN.2[is.na(merged$NUC.DIVERSITY.WITHIN.2)] = merged$NUC.DIVERSITY.WITHIN[is.na(merged$NUC.DIVERSITY.WITHIN.2)]

merged$PI.1[is.na(merged$PI.1)] = merged$PI.2[is.na(merged$PI.1)]
merged$PI.2[is.na(merged$PI.2)] = merged$PI.1[is.na(merged$PI.2)]


merged$HUDSON.K_ST[is.na(merged$HUDSON.K_ST)] = merged$HUDSON.K_ST.3[is.na(merged$HUDSON.K_ST)]
merged$HUDSON.K_ST.3[is.na(merged$HUDSON.K_ST.3)] = merged$HUDSON.K_ST[is.na(merged$HUDSON.K_ST.3)]
merged$HUDSON.G_ST[is.na(merged$HUDSON.G_ST)] = merged$HUDSON.G_ST.3[is.na(merged$HUDSON.G_ST)]
merged$HUDSON.G_ST.3[is.na(merged$HUDSON.G_ST.3)] = merged$HUDSON.G_ST[is.na(merged$HUDSON.G_ST.3)]
merged$HUDSON.H_ST[is.na(merged$HUDSON.H_ST)] = merged$HUDSON.H_ST.3[is.na(merged$HUDSON.H_ST)]
merged$HUDSON.H_ST.3[is.na(merged$HUDSON.H_ST.3)] = merged$HUDSON.H_ST[is.na(merged$HUDSON.H_ST.3)]


merged$YEAR[is.na(merged$YEAR) & grepl("CARDCARD",merged$FILE)] = "CARDINALIS"
merged$YEAR[is.na(merged$YEAR) & grepl("LOCUS",merged$FILE)] = "SCUTULATUS"
merged$YEAR[is.na(merged$YEAR) & grepl("GETULA",merged$FILE)] = "GETULA"
merged$YEAR[is.na(merged$YEAR) & grepl("CATENIFER",merged$FILE)] = "CATENIFER"
merged$YEAR[is.na(merged$YEAR) & grepl("SCUTULATUS",merged$FILE)] = "SCUTULATUS"
merged$YEAR[is.na(merged$YEAR) & grepl("ATROX",merged$FILE)] = "ATROX"
merged$YEAR[is.na(merged$YEAR) & grepl("SCUTUALTUS",merged$FILE)] = "SCUTULATUS"
merged$YEAR[is.na(merged$YEAR) & grepl("CARDTEST",merged$FILE)] = "CARDINALIS"

merged$DEMOG[is.na(merged$DEMOG) & grepl("CARDCARD",merged$FILE)] = "CARDINALIS"
merged$DEMOG[is.na(merged$DEMOG) & grepl("LOCUS",merged$FILE)] = "SCUTULATUS"
merged$DEMOG[is.na(merged$DEMOG) & grepl("GETULA",merged$FILE)] = "GETULA"
merged$DEMOG[is.na(merged$DEMOG) & grepl("CATENIFER",merged$FILE)] = "CATENIFER"
merged$DEMOG[is.na(merged$DEMOG) & grepl("SCUTULATUS",merged$FILE)] = "SCUTULATUS"
merged$DEMOG[is.na(merged$DEMOG) & grepl("ATROX",merged$FILE)] = "ATROX"
merged$DEMOG[is.na(merged$DEMOG) & grepl("SCUTUALTUS",merged$FILE)] = "SCUTULATUS"
merged$DEMOG[is.na(merged$DEMOG) & grepl("CARDTEST",merged$FILE)] = "CARDINALIS"

merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL1",merged$FILE)] = "PANMIXIA"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL2",merged$FILE)] = "PANMIXIA"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL3",merged$FILE)] = "ISOLATION"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL4",merged$FILE)] = "ISOLATION"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL5",merged$FILE)] = "GENEFLOW"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL6",merged$FILE)] = "GENEFLOW"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL7",merged$FILE)] = "SECONDARY"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MODEL8",merged$FILE)] = "SECONDARY"

merged$DEMOG[is.na(merged$DEMOG) & grepl("POP-1",merged$FILE)] = "PANMIXIA"
merged$DEMOG[is.na(merged$DEMOG) & grepl("SECCON-1",merged$FILE)] = "SECONDARY"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MIGRATE-0-",merged$FILE)] = "ISOLATION"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MIGRATE-0.0-",merged$FILE)] = "ISOLATION"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MIGRATE-2",merged$FILE)] = "GENEFLOW"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MIGRATE-0.1",merged$FILE)] = "GENEFLOW"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MIGRATE-0.01",merged$FILE)] = "GENEFLOW"
merged$DEMOG[is.na(merged$DEMOG) & grepl("MIGRATE-0.001",merged$FILE)] = "GENEFLOW"

merged$IBD[is.na(merged$IBD) & grepl("IBD-0-",merged$FILE)] = 0
merged$IBD[is.na(merged$IBD) & grepl("IBD-1-",merged$FILE)] = 1

merged$IBD[is.na(merged$IBD) & grepl("MODEL1",merged$FILE)] = 0
merged$IBD[is.na(merged$IBD) & grepl("MODEL3",merged$FILE)] = 0
merged$IBD[is.na(merged$IBD) & grepl("MODEL5",merged$FILE)] = 0
merged$IBD[is.na(merged$IBD) & grepl("MODEL7",merged$FILE)] = 0

merged$IBD[is.na(merged$IBD) & grepl("MODEL2",merged$FILE)] = 1
merged$IBD[is.na(merged$IBD) & grepl("MODEL4",merged$FILE)] = 1
merged$IBD[is.na(merged$IBD) & grepl("MODEL6",merged$FILE)] = 1
merged$IBD[is.na(merged$IBD) & grepl("MODEL8",merged$FILE)] = 1

merged$IBD[is.na(merged$IBD) & grepl("CARDCARD",merged$FILE)] = "CARDINALIS"
merged$IBD[is.na(merged$IBD) & grepl("LOCUS",merged$FILE)] = "SCUTULATUS"
merged$IBD[is.na(merged$IBD) & grepl("GETULA",merged$FILE)] = "GETULA"
merged$IBD[is.na(merged$IBD) & grepl("CATENIFER",merged$FILE)] = "CATENIFER"
merged$IBD[is.na(merged$IBD) & grepl("SCUTULATUS",merged$FILE)] = "SCUTULATUS"
merged$IBD[is.na(merged$IBD) & grepl("ATROX",merged$FILE)] = "ATROX"
merged$IBD[is.na(merged$IBD) & grepl("SCUTUALTUS",merged$FILE)] = "SCUTULATUS"
merged$IBD[is.na(merged$IBD) & grepl("CARDTEST",merged$FILE)] = "CARDINALIS"


merged$PROPORTION.BIALLELIC.SITES = merged$N.BIALLELIC.SITES / merged$N.SITES
merged$PROPORTION.GAPS = merged$N.GAPS / merged$N.SITES
merged$PROPORTION.POLYALLELIC.SITES = merged$N.POLYALLELIC.SITES / merged$N.SITES
merged$PROPORTION.SEGREGATING.SITES = merged$N.SEGREGATING.SITES / merged$N.SITES
merged$PROPORTION.SEGREGATING.SITES.1 = merged$N.SEGREGATING.SITES.1 / merged$N.SITES
merged$PROPORTION.SEGREGATING.SITES.2 = merged$N.SEGREGATING.SITES.2 / merged$N.SITES
merged$PROPORTION.UNKNOWNS = merged$N.UNKNOWNS / merged$N.SITES
merged$PROPORTION.VALID.SITES = merged$N.VALID.SITES / merged$N.SITES




merged$CL.POP.1[is.na(merged$CL.POP.1)] = merged$CL.POP1[is.na(merged$CL.POP.1)]
merged$CL.POP1[is.na(merged$CL.POP1)] = merged$CL.POP.1[is.na(merged$CL.POP1)]
merged$CL.POP.2[is.na(merged$CL.POP.2)] = merged$CL.POP2[is.na(merged$CL.POP.2)]
merged$CL.POP2[is.na(merged$CL.POP2)] = merged$CL.POP.2[is.na(merged$CL.POP2)]

merged$CLR.POP.1[is.na(merged$CLR.POP.1)] = merged$CLR.POP1[is.na(merged$CLR.POP.1)]
merged$CLR.POP1[is.na(merged$CLR.POP1)] = merged$CLR.POP.1[is.na(merged$CLR.POP1)]
merged$CLR.POP.2[is.na(merged$CLR.POP.2)] = merged$CLR.POP2[is.na(merged$CLR.POP.2)]
merged$CLR.POP2[is.na(merged$CLR.POP2)] = merged$CLR.POP.2[is.na(merged$CLR.POP2)]


merged$DIV_PER_SITE[is.na(merged$DIV_PER_SITE)] = merged$DIV_PER_SITE.1[is.na(merged$DIV_PER_SITE)]
merged$DIV_PER_SITE.1[is.na(merged$DIV_PER_SITE.1)] = merged$DIV_PER_SITE[is.na(merged$DIV_PER_SITE.1)]

merged$FAY.WU.H[is.na(merged$FAY.WU.H)] = merged$FAY.WU.H.1[is.na(merged$FAY.WU.H)]
merged$FAY.WU.H[is.na(merged$FAY.WU.H)] = merged$FAY.WU.H.2[is.na(merged$FAY.WU.H)]

merged$FAY.WU.H.1[is.na(merged$FAY.WU.H.1)] = merged$FAY.WU.H[is.na(merged$FAY.WU.H.1)]
merged$FAY.WU.H.1[is.na(merged$FAY.WU.H.1)] = merged$FAY.WU.H.2[is.na(merged$FAY.WU.H.1)]

merged$FAY.WU.H.2[is.na(merged$FAY.WU.H.2)] = merged$FAY.WU.H[is.na(merged$FAY.WU.H.2)]
merged$FAY.WU.H.2[is.na(merged$FAY.WU.H.2)] = merged$FAY.WU.H.1[is.na(merged$FAY.WU.H.2)]



merged$KURT_CL.POP.1[is.na(merged$KURT_CL.POP.1)] = merged$KURT_CL.POP1[is.na(merged$KURT_CL.POP.1)]
merged$KURT_CL.POP1[is.na(merged$KURT_CL.POP1)] = merged$KURT_CL.POP.1[is.na(merged$KURT_CL.POP1)]
merged$KURT_CL.POP.2[is.na(merged$KURT_CL.POP.2)] = merged$KURT_CL.POP2[is.na(merged$KURT_CL.POP.2)]
merged$KURT_CL.POP2[is.na(merged$KURT_CL.POP2)] = merged$KURT_CL.POP.2[is.na(merged$KURT_CL.POP2)]

merged$MEAN_CL.POP.1[is.na(merged$MEAN_CL.POP.1)] = merged$MEAN_CL.POP1[is.na(merged$MEAN_CL.POP.1)]
merged$MEAN_CL.POP1[is.na(merged$MEAN_CL.POP1)] = merged$MEAN_CL.POP.1[is.na(merged$MEAN_CL.POP1)]
merged$MEAN_CL.POP.2[is.na(merged$MEAN_CL.POP.2)] = merged$MEAN_CL.POP2[is.na(merged$MEAN_CL.POP.2)]
merged$MEAN_CL.POP2[is.na(merged$MEAN_CL.POP2)] = merged$MEAN_CL.POP.2[is.na(merged$MEAN_CL.POP2)]

merged$SKEW_CL.POP.1[is.na(merged$SKEW_CL.POP.1)] = merged$SKEW_CL.POP1[is.na(merged$SKEW_CL.POP.1)]
merged$SKEW_CL.POP1[is.na(merged$SKEW_CL.POP1)] = merged$SKEW_CL.POP.1[is.na(merged$SKEW_CL.POP1)]
merged$SKEW_CL.POP.2[is.na(merged$SKEW_CL.POP.2)] = merged$SKEW_CL.POP2[is.na(merged$SKEW_CL.POP.2)]
merged$SKEW_CL.POP2[is.na(merged$SKEW_CL.POP2)] = merged$SKEW_CL.POP.2[is.na(merged$SKEW_CL.POP2)]

merged$VAR_CL.POP.1[is.na(merged$VAR_CL.POP.1)] = merged$VAR_CL.POP1[is.na(merged$VAR_CL.POP.1)]
merged$VAR_CL.POP1[is.na(merged$VAR_CL.POP1)] = merged$VAR_CL.POP.1[is.na(merged$VAR_CL.POP1)]
merged$VAR_CL.POP.2[is.na(merged$VAR_CL.POP.2)] = merged$VAR_CL.POP2[is.na(merged$VAR_CL.POP.2)]
merged$VAR_CL.POP2[is.na(merged$VAR_CL.POP2)] = merged$VAR_CL.POP.2[is.na(merged$VAR_CL.POP2)]


write.table(merged,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/MERGED_card_snake_sims_15OCT2020.txt",
            row.names = F,sep="\t")



merged=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/MERGED_card_snake_sims_15OCT2020_trimmed.txt",header=T,sep="\t")


write.table(cor(merged[,c(1:2,4:5,7:31,33:176,178:382)],use="pairwise.complete.obs"),"test.corrplot.stats.card.snake.sim.txt",sep="\t")
write.table(merged,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/MERGED_card_snake_sims_15OCT2020_trimmed.txt",
            row.names = F,sep="\t")

pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/boxplots2.0.pdf")
for(i in 1:ncol(merged)){
  boxplot_data = as.numeric(merged[,i])
  if(length(unique(boxplot_data))>1){
    boxplot(boxplot_data~merged$DEMOG,las=2,main=colnames(merged)[i])
  }
}
dev.off()


## subset only emp and sim both

cols_to_keep = c("DEMOG",
                 "DIV_PER_SITE",
                 "FILE",
                 "FU.LI.D.1",
                 "FU.LI.D.2",
                 "FU.LI.D",
                 "FU.LI.F.1",
                 "FU.LI.F.2",
                 "FU.LI.F",
                 "HAP.DIVERSITY.WITHIN.1",
                 "HAP.DIVERSITY.WITHIN.2",
                 "HAP.DIVERSITY.WITHIN.3",
                 "HAP.DIVERSITY.WITHIN",
                 "HAP.F_ST.VS.ALL.1",
                 "HAP.F_ST.VS.ALL.2",
                 "HAP.F_ST.VS.ALL",
                 "HUDSON.G_ST",
                 "HUDSON.H_ST",
                 "HUDSON.K_ST",
                 "HUDSON.KAPLAN.RM.1",
                 "HUDSON.KAPLAN.RM.2",
                 "HUDSON.KAPLAN.RM",
                 "IBD",
                 "KELLY.Z_NS.1",
                 "KELLY.Z_NS",
                 "KURT_FU.LI.D.1.WINDOW",
                 "KURT_FU.LI.D.WINDOW",
                 "KURT_FU.LI.F.1.WINDOW",
                 "KURT_FU.LI.F.WINDOW",
                 "KURT_HUDSON.G_ST.WINDOW",
                 "KURT_HUDSON.H_ST.WINDOW",
                 "KURT_HUDSON.K_ST.WINDOW",
                 "KURT_HUDSON.KAPLAN.RM.1.WINDOW",
                 "KURT_HUDSON.KAPLAN.RM.WINDOW",
                 "KURT_KELLY.Z_NS.1.WINDOW",
                 "KURT_KELLY.Z_NS.WINDOW",
                 "KURT_N.SEGREGATING.SITES.1.WINDOW",
                 "KURT_N.SEGREGATING.SITES.WINDOW",
                 "KURT_NEI.G_ST.WINDOW",
                 "KURT_NUC.DIVERSITY.WITHIN.1",
                 "KURT_NUC.DIVERSITY.WITHIN",
                 "KURT_ROZAS.ZA.1.WINDOW",
                 "KURT_ROZAS.ZA.WINDOW",
                 "KURT_ROZAS.ZZ.1.WINDOW",
                 "KURT_ROZAS.ZZ.WINDOW",
                 "KURT_TAJIMA.D.1.WINDOW",
                 "KURT_TAJIMA.D.WINDOW",
                 "KURT_THETA_ACHAZ.TAJIMA.1.WINDOW",
                 "KURT_THETA_ACHAZ.TAJIMA.WINDOW",
                 "KURT_THETA_ACHAZ.WATTERSON.1.WINDOW",
                 "KURT_THETA_ACHAZ.WATTERSON.WINDOW",
                 "KURT_THETA_FU.LI.1.WINDOW",
                 "KURT_THETA_FU.LI.WINDOW",
                 "KURT_THETA_TAJIMA.1.WINDOW",
                 "KURT_THETA_TAJIMA.WINDOW",
                 "KURT_THETA_WATTERSON.1.WINDOW",
                 "KURT_THETA_WATTERSON.WINDOW",
                 "KURT_WALL.B.1.WINDOW",
                 "KURT_WALL.B.WINDOW",
                 "KURT_WALL.Q.1.WINDOW",
                 "KURT_WALL.Q.WINDOW",
                 "MEAN_FU.LI.D.1.WINDOW",
                 "MEAN_FU.LI.D.WINDOW",
                 "MEAN_FU.LI.F.1.WINDOW",
                 "MEAN_FU.LI.F.WINDOW",
                 "MEAN_HUDSON.G_ST.WINDOW",
                 "MEAN_HUDSON.H_ST.WINDOW",
                 "MEAN_HUDSON.K_ST.WINDOW",
                 "MEAN_HUDSON.KAPLAN.RM.1.WINDOW",
                 "MEAN_HUDSON.KAPLAN.RM.WINDOW",
                 "MEAN_KELLY.Z_NS.1.WINDOW",
                 "MEAN_KELLY.Z_NS.WINDOW",
                 "MEAN_N.SEGREGATING.SITES.1.WINDOW",
                 "MEAN_N.SEGREGATING.SITES.WINDOW",
                 "MEAN_NEI.G_ST.WINDOW",
                 "MEAN_NUC.DIVERSITY.WITHIN.1",
                 "MEAN_NUC.DIVERSITY.WITHIN",
                 "MEAN_ROZAS.ZA.1.WINDOW",
                 "MEAN_ROZAS.ZA.WINDOW",
                 "MEAN_ROZAS.ZZ.1.WINDOW",
                 "MEAN_ROZAS.ZZ.WINDOW",
                 "MEAN_TAJIMA.D.1.WINDOW",
                 "MEAN_TAJIMA.D.WINDOW",
                 "MEAN_THETA_ACHAZ.TAJIMA.1.WINDOW",
                 "MEAN_THETA_ACHAZ.TAJIMA.WINDOW",
                 "MEAN_THETA_ACHAZ.WATTERSON.1.WINDOW",
                 "MEAN_THETA_ACHAZ.WATTERSON.WINDOW",
                 "MEAN_THETA_FU.LI.1.WINDOW",
                 "MEAN_THETA_FU.LI.WINDOW",
                 "MEAN_THETA_TAJIMA.1.WINDOW",
                 "MEAN_THETA_TAJIMA.WINDOW",
                 "MEAN_THETA_WATTERSON.1.WINDOW",
                 "MEAN_THETA_WATTERSON.WINDOW",
                 "MEAN_WALL.B.1.WINDOW",
                 "MEAN_WALL.B.WINDOW",
                 "MEAN_WALL.Q.1.WINDOW",
                 "MEAN_WALL.Q.WINDOW",
                 "N.BIALLELIC.SITES",
                 "N.POLYALLELIC.SITES",
                 "N.SEGREGATING.SITES.1",
                 "N.SEGREGATING.SITES.2",
                 "N.SEGREGATING.SITES",
                 "N.SITES",
                 "NEI.G_ST.1",
                 "NEI.G_ST",
                 "NUC.DIVERSITY.WITHIN.1",
                 "NUC.DIVERSITY.WITHIN.2",
                 "NUC.DIVERSITY.WITHIN.3",
                 "NUC.DIVERSITY.WITHIN",
                 "NUC.F_ST.VS.ALL",
                 "P1_NONSYN",
                 "PI.1",
                 "PI.2",
                 "PI.3",
                 "PI.4",
                 "PI",
                 "PROB_NUC.DIVERSITY.WITHIN.1",
                 "PROPORTION.BIALLELIC.SITES",
                 "PROPORTION.POLYALLELIC.SITES",
                 "PROPORTION.SEGREGATING.SITES.1",
                 "PROPORTION.SEGREGATING.SITES.2",
                 "PROPORTION.SEGREGATING.SITES",
                 "ROZAS.R_2.1",
                 "ROZAS.R_2.2",
                 "ROZAS.R_2",
                 "ROZAS.ZA.1",
                 "ROZAS.ZA",
                 "ROZAS.ZZ.1",
                 "ROZAS.ZZ",
                 "SKEW_FU.LI.D.1.WINDOW",
                 "SKEW_FU.LI.D.WINDOW",
                 "SKEW_FU.LI.F.1.WINDOW",
                 "SKEW_FU.LI.F.WINDOW",
                 "SKEW_HUDSON.G_ST.WINDOW",
                 "SKEW_HUDSON.H_ST.WINDOW",
                 "SKEW_HUDSON.K_ST.WINDOW",
                 "SKEW_HUDSON.KAPLAN.RM.1.WINDOW",
                 "SKEW_HUDSON.KAPLAN.RM.WINDOW",
                 "SKEW_KELLY.Z_NS.1.WINDOW",
                 "SKEW_KELLY.Z_NS.WINDOW",
                 "SKEW_N.SEGREGATING.SITES.1.WINDOW",
                 "SKEW_N.SEGREGATING.SITES.WINDOW",
                 "SKEW_NEI.G_ST.WINDOW",
                 "SKEW_NUC.DIVERSITY.WITHIN.1",
                 "SKEW_NUC.DIVERSITY.WITHIN",
                 "SKEW_ROZAS.ZA.1.WINDOW",
                 "SKEW_ROZAS.ZA.WINDOW",
                 "SKEW_ROZAS.ZZ.1.WINDOW",
                 "SKEW_ROZAS.ZZ.WINDOW",
                 "SKEW_TAJIMA.D.1.WINDOW",
                 "SKEW_TAJIMA.D.WINDOW",
                 "SKEW_THETA_ACHAZ.TAJIMA.1.WINDOW",
                 "SKEW_THETA_ACHAZ.TAJIMA.WINDOW",
                 "SKEW_THETA_ACHAZ.WATTERSON.1.WINDOW",
                 "SKEW_THETA_ACHAZ.WATTERSON.WINDOW",
                 "SKEW_THETA_FU.LI.1.WINDOW",
                 "SKEW_THETA_FU.LI.WINDOW",
                 "SKEW_THETA_TAJIMA.1.WINDOW",
                 "SKEW_THETA_TAJIMA.WINDOW",
                 "SKEW_THETA_WATTERSON.1.WINDOW",
                 "SKEW_THETA_WATTERSON.WINDOW",
                 "SKEW_WALL.B.1.WINDOW",
                 "SKEW_WALL.B.WINDOW",
                 "SKEW_WALL.Q.1.WINDOW",
                 "SKEW_WALL.Q.WINDOW",
                 "TAJIMA.D.1",
                 "TAJIMA.D.2",
                 "TAJIMA.D",
                 "THETA_ACHAZ.TAJIMA.1",
                 "THETA_ACHAZ.TAJIMA.2",
                 "THETA_ACHAZ.TAJIMA",
                 "THETA_ACHAZ.WATTERSON.1",
                 "THETA_ACHAZ.WATTERSON.2",
                 "THETA_ACHAZ.WATTERSON",
                 "THETA_FAY.WU.2",
                 "THETA_FU.LI.1",
                 "THETA_FU.LI.2",
                 "THETA_FU.LI",
                 "THETA_TAJIMA.1",
                 "THETA_TAJIMA.2",
                 "THETA_TAJIMA",
                 "THETA_WATTERSON.1",
                 "THETA_WATTERSON.2",
                 "THETA_WATTERSON",
                 "TRANS.TRANSV.RATIO",
                 "VAR_FU.LI.D.1.WINDOW",
                 "VAR_FU.LI.D.WINDOW",
                 "VAR_FU.LI.F.1.WINDOW",
                 "VAR_FU.LI.F.WINDOW",
                 "VAR_HUDSON.G_ST.WINDOW",
                 "VAR_HUDSON.H_ST.WINDOW",
                 "VAR_HUDSON.K_ST.WINDOW",
                 "VAR_HUDSON.KAPLAN.RM.1.WINDOW",
                 "VAR_HUDSON.KAPLAN.RM.WINDOW",
                 "VAR_KELLY.Z_NS.1.WINDOW",
                 "VAR_KELLY.Z_NS.WINDOW",
                 "VAR_N.SEGREGATING.SITES.1.WINDOW",
                 "VAR_N.SEGREGATING.SITES.WINDOW",
                 "VAR_NEI.G_ST.WINDOW",
                 "VAR_NUC.DIVERSITY.WITHIN.1",
                 "VAR_NUC.DIVERSITY.WITHIN",
                 "VAR_ROZAS.ZA.1.WINDOW",
                 "VAR_ROZAS.ZA.WINDOW",
                 "VAR_ROZAS.ZZ.1.WINDOW",
                 "VAR_ROZAS.ZZ.WINDOW",
                 "VAR_TAJIMA.D.1.WINDOW",
                 "VAR_TAJIMA.D.WINDOW",
                 "VAR_THETA_ACHAZ.TAJIMA.1.WINDOW",
                 "VAR_THETA_ACHAZ.TAJIMA.WINDOW",
                 "VAR_THETA_ACHAZ.WATTERSON.1.WINDOW",
                 "VAR_THETA_ACHAZ.WATTERSON.WINDOW",
                 "VAR_THETA_FU.LI.1.WINDOW",
                 "VAR_THETA_FU.LI.WINDOW",
                 "VAR_THETA_TAJIMA.1.WINDOW",
                 "VAR_THETA_TAJIMA.WINDOW",
                 "VAR_THETA_WATTERSON.1.WINDOW",
                 "VAR_THETA_WATTERSON.WINDOW",
                 "VAR_WALL.B.1.WINDOW",
                 "VAR_WALL.B.WINDOW",
                 "VAR_WALL.Q.1.WINDOW",
                 "VAR_WALL.Q.WINDOW",
                 "WALL.B.1",
                 "WALL.B",
                 "WALL.Q.1",
                 "WALL.Q",
                 "YEAR")

merged_subset = merged[,cols_to_keep]

write.table(cor(merged_subset[,c(2,4:22,24:225)],use="pairwise.complete.obs"),"test.corrplot.stats.card.snake.sim.txt",sep="\t")

write.table(merged_subset,"/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/MERGED_card_snake_sims_15OCT2020_trimmed_subset.txt",
            row.names = F,sep="\t")


merged_subset=read.table("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/MERGED_card_snake_sims_15OCT2020_trimmed_subset.txt",header=T,sep="\t")

sort(colSums(is.na(merged_subset)))
sort(colSums(is.na(merged_subset[merged_subset$DEMOG %in% c("ISOLATION","SECONDARY","PANMIXIA","GENEFLOW"),])))
sort(colSums(is.na(merged_subset[merged_subset$DEMOG %in% c("CARDINALIS"),])))
sort(colSums(is.na(merged_subset[merged_subset$DEMOG %in% c("SCUTULATUS","CATENIFER","ATROX","GETULA"),])))


## more trimmed 
cols_to_keep_pca = c("DEMOG","DIV_PER_SITE","FILE","FU.LI.D","FU.LI.D.1","FU.LI.F","FU.LI.F.1","HAP.DIVERSITY.WITHIN","HAP.F_ST.VS.ALL","HUDSON.G_ST","HUDSON.H_ST","HUDSON.K_ST","HUDSON.KAPLAN.RM","HUDSON.KAPLAN.RM.1","IBD","KELLY.Z_NS","KELLY.Z_NS.1","KURT_NUC.DIVERSITY.WITHIN","MEAN_NUC.DIVERSITY.WITHIN","N.BIALLELIC.SITES","N.POLYALLELIC.SITES","N.SEGREGATING.SITES","N.SEGREGATING.SITES.1","N.SITES","NEI.G_ST","NUC.DIVERSITY.WITHIN","NUC.F_ST.VS.ALL","P1_NONSYN","PI","PROPORTION.BIALLELIC.SITES","PROPORTION.POLYALLELIC.SITES","PROPORTION.SEGREGATING.SITES","PROPORTION.SEGREGATING.SITES.1","ROZAS.ZA","ROZAS.ZA.1","ROZAS.ZZ","ROZAS.ZZ.1","SKEW_NUC.DIVERSITY.WITHIN","TAJIMA.D","TAJIMA.D.1","THETA_ACHAZ.TAJIMA","THETA_ACHAZ.TAJIMA.1","THETA_ACHAZ.WATTERSON","THETA_ACHAZ.WATTERSON.1","THETA_FU.LI","THETA_FU.LI.1","THETA_TAJIMA","THETA_TAJIMA.1","THETA_WATTERSON","THETA_WATTERSON.1","TRANS.TRANSV.RATIO","VAR_NUC.DIVERSITY.WITHIN","WALL.B.1","YEAR")
merged_pca = merged_subset[,cols_to_keep_pca]
write.table(cor(merged_pca[,c(2,4:14,16:53)],use="pairwise.complete.obs"),"test.corrplot.stats.card.snake.sim.txt",sep="\t")



merged_pca_complete = merged_pca[complete.cases(merged_pca),]
pca=prcomp(merged_pca_complete[,c(2,4,6:14,16:27,29:50,52,53)],center=T,scale. = T)


colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorBlindBlack8_transparent  <- c("#00000011", "#E69F0011", "#56B4E911", "#009E7311", 
                       "#F0E44211", "#0072B211", "#D55E0011", "#CC79A711")



pdf("PCA_new_sumstats.pdf")
plot(pca$x[,1:2],
     col=colorBlindBlack8[as.numeric(as.factor(merged_pca_complete$DEMOG))],
     pch=as.numeric(as.factor(merged_pca_complete$DEMOG)),cex=0.75)
legend(x="bottomleft",
       legend=levels(as.factor(merged_pca_complete$DEMOG)),
       cex=0.75,col=colorBlindBlack8,pch=1:8)
plot(pca$x[,1:2],
     col="grey",xlim=c(-30,10),ylim=c(-20,10),
     pch=as.numeric(as.factor(merged_pca_complete$DEMOG)),cex=0.75)
legend(x="bottomleft",
       legend=levels(as.factor(merged_pca_complete$DEMOG)),
       cex=0.75,col=colorBlindBlack8,pch=1:8)
for (i in 1:8) {
  d = pca$x[,1:2]
  
  included=(as.numeric(as.factor(merged_pca_complete$DEMOG))==i)
  
  if(sum(included)>1){
    
    car::ellipse(center = colMeans( d[included,],na.rm = T), 
                 shape = cov( d[included,],use="pairwise.complete.obs"),
                 center.pch=i,center.cex=2,
                 radius = sqrt(qchisq(.95, df=2)),col = colorBlindBlack8[i])
  }

}
plot(pca$x[,1:2],
     col="grey",xlim=c(-30,10),ylim=c(-20,10),
     pch=as.numeric(as.factor(merged_pca_complete$DEMOG)),cex=0.75)
legend(x="bottomleft",
       legend=levels(as.factor(merged_pca_complete$DEMOG))[c(1:3,7)],
       cex=0.75,col=colorBlindBlack8[c(1:3,7)],pch=c(1:3,7))
for (i in c(1:3,7)) {
  d = pca$x[,1:2]
  
  included=(as.numeric(as.factor(merged_pca_complete$DEMOG))==i)
  
  if(sum(included)>1){
    
    car::ellipse(center = colMeans( d[included,],na.rm = T), 
                 shape = cov( d[included,],use="pairwise.complete.obs"),
                 center.pch=i,center.cex=2,
                 radius = sqrt(qchisq(.95, df=2)),col = colorBlindBlack8[i])
  }
  
}
plot(pca$x[,1:2],
     col="grey",xlim=c(-30,10),ylim=c(-20,10),
     pch=as.numeric(as.factor(merged_pca_complete$DEMOG)),cex=0.75)
legend(x="bottomleft",
       legend=levels(as.factor(merged_pca_complete$DEMOG))[c(4:6,8)],
       cex=0.75,col=colorBlindBlack8[c(4:6,8)],pch=c(4:6,8))
for (i in c(4:6,8)) {
  d = pca$x[,1:2]
  
  included=(as.numeric(as.factor(merged_pca_complete$DEMOG))==i)
  
  if(sum(included)>1){
    
    car::ellipse(center = colMeans( d[included,],na.rm = T), 
                 shape = cov( d[included,],use="pairwise.complete.obs"),
                 center.pch=i,center.cex=2,
                 radius = sqrt(qchisq(.95, df=2)),col = colorBlindBlack8[i])
  }
  
}

par(mfrow=c(2,2))
for(i in 1:8){
  d = pca$x[,1:2]
  included=(as.numeric(as.factor(merged_pca_complete$DEMOG))==i)
  
  plot(d,xlim=c(-30,10),ylim=c(-20,10),cex=0.75,col="lightgrey",main=levels(as.factor(merged_pca_complete$DEMOG))[i])
  points(d[included,],col=colorBlindBlack8[i],pch=i)
  
}

dev.off()

