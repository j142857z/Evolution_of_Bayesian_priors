#!/bin/bash

# ./t2 100 100 2 1000 100 1 10

./EvolutionaryCompetitionByEnvironment 100 100 2 1000 100 1 10 0.025 1 0.025 1 40

matlab -nodisplay -nosplash -r "OutputSurfacePlot('fileOfAvgBayesianPopulationByHyperparameterPair.csv', 'SurfPlot.pdf', 0, 1, 0, 1); exit"
