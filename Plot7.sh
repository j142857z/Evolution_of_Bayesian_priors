#!/bin/bash

Rscript PlotEntropyOfBetaDistributions.R Entropy2.csv 1 10 100 10 100

matlab -nodisplay -nosplash -r "OutputSurfacePlot('Entropy2.csv', 'Entropy2.pdf', 10, 100, 10, 100); exit"

rm Entropy2.csv