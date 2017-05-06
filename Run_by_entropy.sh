#!/bin/bash

populationSize='100'
maxAlleleValue='100'
minNumberOfObservations='2'
maxNumberOfObservations='11'
numberOfGenerations='1000'
numberOfRepetitions='100'
numberOfGenerationsPerChange='1'
numberOfProbabilityValuesPerGeneration='10'

minAlpha='10'
maxAlpha='100'
minBeta='10'
maxBeta='100'
numberOfHyperparameterPoints='10'

./EvolutionaryCompetitionByEnvironmentAndLearningLength $populationSize $maxAlleleValue $minNumberOfObservations $maxNumberOfObservations $numberOfGenerations $numberOfRepetitions $numberOfGenerationsPerChange $numberOfProbabilityValuesPerGeneration $minAlpha $maxAlpha $minBeta $maxBeta $numberOfHyperparameterPoints

Rscript bayesian_pop_scat.R "BayesianPopulationByEntropyAndAlphaBetaAbove1.csv"  "Population_by_n_2.pdf"
