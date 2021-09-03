#!/usr/local/bin/Rscript

args<-commandArgs(TRUE)
cat("Now generating statistical plots!","\n")

Stat <- read.delim(args[2], header=TRUE)
foo <- read.delim(args[3], header=FALSE)


#setwd(args[1])

pdf('Graphical-Results.pdf')

myhist <- hist(Stat[,2],breaks=50,plot=FALSE)
multiplier <- myhist$counts / myhist$density
mydensity <- density(Stat[,2])
mydensity$y <- mydensity$y * multiplier[1]
plot(myhist,ylim=c(0,max(myhist$counts)*1.1),main="Distribution of read lengths", xlab="Length of reads",col="lightblue")
lines(mydensity,col="purple")
abline(v = mean(Stat[,2]), col = "red")
abline(v = median(Stat[,2]), col = "red",lty="dashed")
abline(v = mean(Stat[,2])-sd(Stat[,2]),col="darkgreen",lty="dotdash")
abline(v = mean(Stat[,2])+sd(Stat[,2]),col="darkgreen",lty="dotdash")
legend("topright", c("Mean", "Median","Stdev","Density"), col = c("red","red","darkgreen","purple"),lty = c(1,2,4,1),cex=0.8)


myhist <- hist(Stat[,3],breaks=50,plot=FALSE)
multiplier <- myhist$counts / myhist$density
mydensity <- density(Stat[,3])
mydensity$y <- mydensity$y * multiplier[1]
plot(myhist,ylim=c(0,max(myhist$counts)*1.1),main="Distribution of GC content", xlab="GC content [%]",col="lightblue")
lines(mydensity,col="purple")
abline(v = mean(Stat[,3]), col = "red")
abline(v = median(Stat[,3]), col = "red",lty="dashed")
abline(v = mean(Stat[,3])-sd(Stat[,3]),col="darkgreen",lty="dotdash")
abline(v = mean(Stat[,3])+sd(Stat[,3]),col="darkgreen",lty="dotdash")
legend("topright", c("Mean", "Median","Stdev","Density"), col = c("red","red","darkgreen","purple"),lty = c(1,2,4,1),cex=0.8)


myhist <- hist(Stat[,4],breaks=50,plot=FALSE)
multiplier <- myhist$counts / myhist$density
mydensity <- density(Stat[,4])
mydensity$y <- mydensity$y * multiplier[1]
plot(myhist,ylim=c(0,max(myhist$counts)*1.1),main="Distribution of qualities", xlab="Mean quality value",col="lightblue")
lines(mydensity,col="purple")
abline(v = mean(Stat[,4]), col = "red")
abline(v = median(Stat[,4]), col = "red",lty="dashed")
abline(v = mean(Stat[,4])-sd(Stat[,4]),col="darkgreen",lty="dotdash")
abline(v = mean(Stat[,4])+sd(Stat[,4]),col="darkgreen",lty="dotdash")
legend("topright", c("Mean", "Median","Stdev","Density"), col = c("red","red","darkgreen","purple"),lty = c(1,2,4,1),cex=0.8)



nf <- layout(mat = matrix(c(1,2),2,1),  height = c(1,3))
par(mar=c(0.1, 5.1, 2.1, 2.1),mgp=c(3,1,0))

Count<-hist(Stat[,2],plot=FALSE)
Count$counts<-cumsum(Count$counts)
Count$counts<-Count$counts/max(Count$counts)*100
plot(rev(Count$counts),type="l",ylab="Length [%]",frame=TRUE,xaxt="n",col="red",lwd=3,main="Quality depending on position in read")
legend("topright",paste("n =",length(Stat[,2])),cex=0.8,bty="n")

par(mar=c(7.1, 5.1, 0.1, 2.1),mgp=c(3,1,0))

names<-c("1","2","3","4","5-9","10-14","15-19","20-24","25-30","31-39","40-49","50-74","75-99","100-149","150-199","200-249","250-299","300-399","400-499","500-599","600-699","700-799","800-899","900-999",">1000")
uniqueData <- unique(foo[,1])
data.list <- vector("list", length(uniqueData))
for(i in 1:length(uniqueData)){
  index <- which(foo[,1]==uniqueData[i])
  data.list[[i]] <- unlist(apply(foo[index,2:3], 1, function(x)(rep(x[1], x[2]))))
}
h <- boxplot(data.list, plot=TRUE,xaxt="n",outline=FALSE,ylab="Quality value")
par(las=3,mgp=c(5,1,0))
title(xlab="Position in reads")
axis(1,at=1:length(h$stats[3,]),labels=names[1:length(h$stats[3,])],xlab="Position in reads")
lines(h$stats[3,],col="red")

invisible(dev.off())

cat("Program ended successfully!\n")