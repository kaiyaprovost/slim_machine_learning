lines = c()

#file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/logtile_X_train.**NN_expand_1604531680.temp_1.out.txt"
#lines= c(lines,readLines(file))

#file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/NN_Results_Empirical/logtile_X_train.**NN_expand_1604531680.temp_2.out.txt"
#lines = c(lines,readLines(file))

file="/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/logfile_X_train.**NN_expand_1604531680.temp_1.out.txt"
lines = c(lines,readLines(file))

lines=gsub(",", " ", lines, ignore.case = FALSE, perl = FALSE,
           fixed = FALSE, useBytes = FALSE)
lines=lines[!(grepl("ConvergenceWarning",lines))]
lines=lines[(grepl("Iteration",lines))]

df = read.table(text=lines,sep=" ")
loss = as.numeric(df[,6])
iteration = as.numeric(df[,2])
color=iteration

for(i in 1:length(iteration)){
  if(i == 1){
    color[i]=1
  } else {
    if(iteration[i-1] < iteration[i]){
      color[i] = color[i-1]
    } else {
      color[i] = color[i-1]+1
    }
  }
}

library(viridis)
palette(viridis::viridis(max(color)))
#palette(rainbow(max(color)))
pdf("/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/loss_vs_epoch.pdf")
plot(iteration,loss,type="n",xlab="Epoch",ylab="Loss")
for(i in 1:max(color)){
  lines(iteration[color==i],loss[color==i],col=i,lwd=2,lty= (i %% 3)+1)
}
#plot(iteration,loss,col=color,cex=0.25,pch=color-1)
dev.off()
