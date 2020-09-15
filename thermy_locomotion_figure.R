filename="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/review/cfb_ch1_multidata2.csv"
y = read.table(filename,header=T,sep=",",fill=T,stringsAsFactors = F)
smally = y[,c("LOCOMOTION.DISPERSAL","ENDO.ECTO","STRUCTURE")]

colors=c("#d7191c","#fdae61","#ffffbf","#abdda4","#2b83ba")
ramp=colorRampPalette(colors)

pdf("black_background_ecto_loco.pdf",
    h=4,w=6)
par(mfrow=c(1,2),mar=c(2,0,1,0),
    bg="black",col.axis="white",
    col.lab="white",
    col.sub="white",
    col.main="white",col="white",
    fg="white")
e=table(smally[,c("ENDO.ECTO","STRUCTURE")])
esum=colSums(e)
e[1,] = e[1,]/esum
e[2,] = e[2,]/esum
barplot(e,col=c("white","grey"),width=esum,yaxt="n")
d=table(smally[,c("LOCOMOTION.DISPERSAL","STRUCTURE")])
dsum=colSums(d)
d[1,] = d[1,]/dsum
d[2,] = d[2,]/dsum
d[3,] = d[3,]/dsum
d[4,] = d[4,]/dsum
d[5,] = d[5,]/dsum
d[6,] = d[6,]/dsum
barplot(d,col=ramp(6),width=dsum,yaxt="n")
dev.off()
