
args<-commandArgs(TRUE)

alphaFile <- args[1]
betaFile <- args[2]
outputFile <- args[3]


# pdf(outputFile, width = 12, height = 12)
png(outputFile, width = 800, height = 800)

fontSize <- 2.5
pointSize <- 0.1

alphas = as.matrix(read.csv(alphaFile,sep=",",head=FALSE))

betas = as.matrix(read.csv(betaFile,sep=",",head=FALSE))

generations <- matrix(, nrow = nrow(alphas), ncol = ncol(alphas))

for(i in 1:nrow(generations)){
  generations[i,]<-i
}

par(mar=c(5,6,2,2))

# plot(generations,alphas, ylim=c(0, max(c(max(alphas), max(betas)))), type='p', xlab='Generations', ylab=expression(paste("Population (", alpha, ", ", beta, ")")), pch=16, col='red', cex=pointSize, cex.lab=fontSize, cex.axis=fontSize)

plot(generations,alphas, ylim=c(0, 6), type='p', xlab='Generations', ylab=expression(paste("Population (", alpha, ", ", beta, ")")), pch=16, col='red', cex=pointSize, cex.lab=fontSize, cex.axis=fontSize)

points(generations, y = betas, type = "p", cex=pointSize, pch=16, col='blue')

dev.off()