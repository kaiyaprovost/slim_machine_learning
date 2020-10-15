data=read.table("/Users/kprovost/Dropbox (AMNH)/temp_ch1.txt",sep="\t",header=T)

data$age.width = data$oldest.age-data$youngest.age

data = data[order(data$oldest.age,data$youngest.age,data$mean.age),]

hist(data$age.width)

par(mar=c(4,0,0,0))
palette(c("black","red","blue"))
plot(data$oldest.age,1:length(data$oldest.age),
     xlim=c(0,25),ylim=c(0,50),
     col=(5-as.numeric(as.factor(data$Explicitly.Dated.))),
     yaxt="n",ylab="",xlab="Age",cex=0.5)
points(data$youngest.age,1:length(data$oldest.age),cex=0.5,
       col=(5-as.numeric(as.factor(data$Explicitly.Dated.))))
points(data$mean.age,1:length(data$oldest.age),
       col=(5-as.numeric(as.factor(data$Explicitly.Dated.))))
segments(x0=data$oldest.age,x1=data$youngest.age,
         y0=1:length(data$oldest.age),
         col=(5-as.numeric(as.factor(data$Explicitly.Dated.))))
abline(v=c(1,2.58,5.333,23.03),col="grey",lty=2)
legend("bottomright",col=c("black","red"),lty=1,
       pch=1,legend=c("Explicitly Estimated","Not Explicitly Estimated"))
