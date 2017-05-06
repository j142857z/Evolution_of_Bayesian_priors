args = commandArgs(trailingOnly=TRUE)

inputFilename = args[1];
outputFilename = args[2];

require(scatterplot3d)

data <- read.csv(file=inputFilename, head=FALSE, sep=",")

colnames(data) = c("n", "alpha", "beta", "entropy", "population");

x=data$n;
y=data$entropy;
z=data$population;

pdf(outputFilename);

scatterplot3d(x=x, y=y, z=z, xlab="Learning length", ylab="Entropy", zlab="Bayesian population", pch=16, highlight.3d=FALSE, color="blue", cex.symbols=0.6);

dev.off();
