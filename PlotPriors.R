
args<-commandArgs(TRUE)

environmentStatesFile <- args[1]
priorsFile <- args[2]
outputFile <- args[3]


# pdf(outputFile, width = 12, height = 12)
png(outputFile, width = 800, height = 800)

fontSize <- 2.5
pointSize <- 0.1

environmentStates = as.matrix(read.csv(environmentStatesFile,sep=",",head=FALSE))

priors = as.matrix(read.csv(priorsFile,sep=",",head=FALSE))

generations <- matrix(, nrow = nrow(priors), ncol = ncol(priors))

for(i in 1:nrow(generations)){
  generations[i,]<-i
}

par(mar=c(5,6,2,2))

plot(generations,priors, type='p', xlab='Generations', ylab='Priors', pch=16, col='red', cex=pointSize, cex.lab=fontSize, cex.axis=fontSize, ylim=c(0,1))

# points(1:nrow(environmentStates), y = environmentStates[,3], type = "p", cex=pointSize)
points(x = environmentStates[,1], y = environmentStates[,3], type = "p", cex=pointSize)

dev.off()