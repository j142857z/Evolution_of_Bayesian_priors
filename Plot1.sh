#!/bin/bash

./t1 100 100 2 10000 15 95 1 10

Rscript PlotPriors.R EnvironmentStatesByGeneration.csv PriorsByGeneration.csv Priors1.png
Rscript PlotGenomes.R PopulationAlphaByGeneration.csv PopulationBetaByGeneration.csv Genomes1.png

mv PopulationAlphaByGeneration.csv alphas1.csv
mv PopulationBetaByGeneration.csv betas1.csv

./t1 100 100 2 10000 15 95 1 10

Rscript PlotPriors.R EnvironmentStatesByGeneration.csv PriorsByGeneration.csv Priors2.png
Rscript PlotGenomes.R PopulationAlphaByGeneration.csv PopulationBetaByGeneration.csv Genomes2.png

mv PopulationAlphaByGeneration.csv alphas2.csv
mv PopulationBetaByGeneration.csv betas2.csv

./t1 100 100 2 10000 15 95 1 10

Rscript PlotPriors.R EnvironmentStatesByGeneration.csv PriorsByGeneration.csv Priors3.png
Rscript PlotGenomes.R PopulationAlphaByGeneration.csv PopulationBetaByGeneration.csv Genomes3.png

mv PopulationAlphaByGeneration.csv alphas3.csv
mv PopulationBetaByGeneration.csv betas3.csv

./t1 100 100 2 10000 15 95 1 10

Rscript PlotPriors.R EnvironmentStatesByGeneration.csv PriorsByGeneration.csv Priors4.png
Rscript PlotGenomes.R PopulationAlphaByGeneration.csv PopulationBetaByGeneration.csv Genomes4.png

mv PopulationAlphaByGeneration.csv alphas4.csv
mv PopulationBetaByGeneration.csv betas4.csv


Rscript PlotEvolvedBetaDistribution.R EvolvedBetaDistribution.pdf 15 95 alphas1.csv betas1.csv alphas2.csv betas2.csv alphas3.csv betas3.csv alphas4.csv betas4.csv 12.5

# rm EnvironmentStatesByGeneration.csv PopulationAlphaByGeneration.csv PopulationBetaByGeneration.csv PriorsByGeneration.csv ObjectiveAlpha.csv ObjectiveBeta.csv
