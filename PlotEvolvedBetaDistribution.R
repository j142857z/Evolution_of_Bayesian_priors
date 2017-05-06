
args<-commandArgs(TRUE)

outputFile <- args[1]

realAlpha <- as.double(args[2])
realBeta <- as.double(args[3])

alphaFile1 <- args[4]
betaFile1 <- args[5]

alphaFile2 <- args[6]
betaFile2 <- args[7]

alphaFile3 <- args[8]
betaFile3 <- args[9]

alphaFile4 <- args[10]
betaFile4 <- args[11]

yAxisLim <- as.double(args[12])



alphas1 = as.matrix(read.csv(alphaFile1,sep=",",head=FALSE))
betas1 = as.matrix(read.csv(betaFile1,sep=",",head=FALSE))

alphas2 = as.matrix(read.csv(alphaFile2,sep=",",head=FALSE))
betas2 = as.matrix(read.csv(betaFile2,sep=",",head=FALSE))

alphas3 = as.matrix(read.csv(alphaFile3,sep=",",head=FALSE))
betas3 = as.matrix(read.csv(betaFile3,sep=",",head=FALSE))

alphas4 = as.matrix(read.csv(alphaFile4,sep=",",head=FALSE))
betas4 = as.matrix(read.csv(betaFile4,sep=",",head=FALSE))

evolvedAlpha1 <- mean(alphas1[nrow(alphas1),])
evolvedBeta1 <- mean(betas1[nrow(betas1),])

evolvedAlpha2 <- mean(alphas2[nrow(alphas2),])
evolvedBeta2 <- mean(betas2[nrow(betas2),])

evolvedAlpha3 <- mean(alphas3[nrow(alphas3),])
evolvedBeta3 <- mean(betas3[nrow(betas3),])

evolvedAlpha4 <- mean(alphas4[nrow(alphas4),])
evolvedBeta4 <- mean(betas4[nrow(betas4),])


pdf(outputFile)

fontSize <- 2

par(mar=c(4,5,2,2))

x=seq(0,1,length=100)
y1=dbeta(x,realAlpha,realBeta)
y2=dbeta(x,evolvedAlpha1,evolvedBeta1)
y3=dbeta(x,evolvedAlpha2,evolvedBeta2)
y4=dbeta(x,evolvedAlpha3,evolvedBeta3)
y5=dbeta(x,evolvedAlpha4,evolvedBeta4)

# plot (c(0,1),c(0,2.9),type="n", xlab="",ylab="PDF", cex.lab=fontSize, cex.axis=fontSize, ylim=c(min(y1,y2,y3,y4,y5),max(y1,y2,y3,y4,y5)))
plot (c(0,1),c(0,2.9),type="n", xlab="",ylab="PDF", cex.lab=fontSize, cex.axis=fontSize, ylim=c(0,yAxisLim))

lines(x,y1,col="black",lwd=2.5, lty=2)
lines(x,y2,col="blue",lwd=2.5)
lines(x,y3,col="green",lwd=2.5)
lines(x,y4,col="red",lwd=2.5)
lines(x,y5,col="yellow",lwd=2.5)

dev.off()