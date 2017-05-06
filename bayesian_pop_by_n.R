require(ggplot2)

data <- read.csv (file="BayesianPopulationByEntropyAndAlphaBetaBelow1.csv", head=FALSE, sep=",")

colnames(data) = c("n", "alpha", "beta", "entropy", "population");

data = data[data[, "alpha"] < data[, "beta"],];

means = data[, "alpha"] / (data[, "alpha"] + data[, "beta"]);

x=data$n;
y=means;
z=data$population;

pdf("Population.pdf");

ggplot(data, aes(x=x, y=y, z=z)) + geom_density2d() + ggtitle("") +  labs(x="n",y=expression(mu))

dev.off();