args<-commandArgs(TRUE)

outputFile <- args[1]

byValue <- as.double(args[2])

minAlpha <- as.double(args[3])
maxAlpha <- as.double(args[4])

minBeta <- as.double(args[5])
maxBeta <- as.double(args[6])

alphas <- seq(minAlpha,maxAlpha,by=byValue)
betas <- seq(minBeta,maxBeta,by=byValue)


calculateBetaEntropy <- function(alphaArg, betaArg){
  entropy <- log(beta(alphaArg,betaArg)) - (alphaArg - 1)*digamma(alphaArg)-(betaArg - 1)*digamma(betaArg)+(alphaArg+betaArg-2)*digamma(alphaArg+betaArg)
  return(entropy)
}


dataGrid <- NULL

xBy <- 0.01
xMin <- xBy
xMax <- 0.99

for (alpha in alphas) {
  for (beta in betas) {
    entropy <- calculateBetaEntropy(alpha,beta)
    dataGrid <- rbind(dataGrid,c(alpha,beta,entropy))
  }
}

write.table(dataGrid, file = outputFile, sep = ",", row.names=FALSE, col.names=FALSE)