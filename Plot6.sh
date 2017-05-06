#!/bin/bash

Rscript PlotEntropyOfBetaDistributions.R Entropy1.csv 0.01 0.01 1 0.01 1

matlab -nodisplay -nosplash -r "OutputSurfacePlot('Entropy1.csv', 'Entropy1.pdf', 0, 1, 0, 1); exit"

rm Entropy1.csv