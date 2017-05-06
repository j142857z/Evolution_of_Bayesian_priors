#!/bin/bash

# ./t3 100 100 2 1000 100 1 10

./EvolutionaryCompetitionByEnvironment 100 100 2 1000 100 1 10 10 100 10 100 37

matlab -nodisplay -nosplash -r "OutputSurfacePlot('fileOfAvgBayesianPopulationByHyperparameterPair.csv', 'SurfPlot.pdf', 10, 100, 10, 100); exit"
