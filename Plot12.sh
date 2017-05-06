#!/bin/bash

./t4 100 100 2 10000 2 2 1 10 50000

matlab -nodisplay -nosplash -r "OutputHistogramPlot('BayesianPopulationPercentages.csv', 'HistPlot.pdf'); exit"