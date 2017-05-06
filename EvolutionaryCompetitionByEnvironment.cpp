#include "common.h"

/**
 * This program produces an ordered list of (alpha,beta) parameters and runs the
 * simulation with each pair a certain number of times. Each repetition of the experiment (with
 * each set of hyperparameters) the average population of bayesians is calculated (average in the
 * course of all generations) then the average from each repetition is calculated for each pair
 * of hyperparameters. Then the program outputs a text file with these averages per hyperparameter
 * pairs.
 *
 * The format of the output file is as follows.
 *
 * alpha_1,beta_1,avg(alpha_1,beta_1)
 * alpha_1,beta_2,avg(alpha_1,beta_2)
 * ...
 * alpha_m,beta_n,avg(alpha_1,beta_n)
 *
 * The list of hyperparameters can be changed only from the code in this file. There is no way from
 * receiving this list as an argument to the program.
 *
 *
 * Also, for each of hyperparameters the program outputs a file
 * fileOfBayesianAverageFrequencywith(alpha,beta).csv in the format.
 *
 * (Avg bayesian population in repetition 1);(Avg bayesian population in repetition 2);(Avg bayesian population in repetition 3);...;(Avg bayesian population in repetition N)
 *
 * where N is the number of repetitions given as argument to the program.
 * */

/**
 * This experiment simulates competition between frequentist and bayesian learners. All individuals
 * observe the same environment, i.e., observations are not independent and the probability values
 * they observe are always for all.
 * */

 /**
  * The arguments must be provided in the following order.
  * 
  * 1. The size of the population. As a positive integer.
  * 
  * 2. The maximum value of alpha and beta that will be used as allele for all individuals. The
  * alpha and beta values of each individual will be in the interval (0,maxAlleleValue]. As a
  * float value.
  * 
  * 3. The number of learning trials each individual carries out when estimating a single
  * environment state. As a positive integer.
  * 
  * 4. 
  */
int main(int argc, char **argv) {
	if (argc < 8) {
		cout
				<< "<Population-Size> <Max-Allele-Value> <Maximum-Learning-Trials> "
				<< "<Num-Of-Generations> <NumberOfRepetitions> <GenerationsPerChange> <ProbabilitiesByGeneration>"
				<< endl;
		return -1;
	}
	srand((unsigned) time(0));

	int populationSize = atoi(argv[1]);
	int maxAlleleValue = atof(argv[2]);
	learningLength = atoi(argv[3]);
	int numberOfGenerations = atoi(argv[4]);
	int numberOfRepetitions = atoi(argv[5]);
	numberOfGenerationsPerChange = atoi(argv[6]);
	unsigned int numberOfProbabilityValuesPerGeneration = atoi(argv[7]);
	
	float minAlpha = atof(argv[8]);
	float maxAlpha = atof(argv[9]);
	float minBeta = atof(argv[10]);
	float maxBeta = atof(argv[11]);
	float numberOfHyperparameterPoints = atoi(argv[12]);
	
	float spaceBetweenAlphaPoints = (maxAlpha - minAlpha) / (numberOfHyperparameterPoints - 1);
	float spaceBetweenBetaPoints = (maxBeta - minBeta) / (numberOfHyperparameterPoints - 1);

	vector<hyperparameterType> listOfEnvironmentHyperparameterPairs = getListOfHyperparameterPairs(minAlpha, maxAlpha, spaceBetweenAlphaPoints, minBeta, maxBeta, spaceBetweenBetaPoints);

	map<hyperparameterType, float> bayesianAveragePopulationByEnvironmentHyperparameterPair;
	map<hyperparameterType, map<float, float> > frequencyOfBayesianAveragesByEnvironmentHyperparameterPair;

	ofstream fileOfAvgBayesianPopulationByHyperparameterPair;
	fileOfAvgBayesianPopulationByHyperparameterPair.open("fileOfAvgBayesianPopulationByHyperparameterPair.csv");

	for (int i = 0; i < listOfEnvironmentHyperparameterPairs.size(); i++) {

		float environmentAlpha = listOfEnvironmentHyperparameterPairs[i].first;
		float environmentBeta = listOfEnvironmentHyperparameterPairs[i].second;

		std::vector<long> avgBayesianPopulationByRepetition;
		std::vector<long> avgFrequentistPopulationsByRepetition;
		for (int repetition = 1; repetition <= numberOfRepetitions;
				repetition++) {
			runEvolutionarySimulation(environmentAlpha, environmentBeta, populationSize, maxAlleleValue, numberOfGenerations, numberOfGenerationsPerChange, numberOfProbabilityValuesPerGeneration, avgBayesianPopulationByRepetition, avgFrequentistPopulationsByRepetition, frequencyOfBayesianAveragesByEnvironmentHyperparameterPair);
		}

		float avgFrequentistPopulationEvolvedInAllRepetitions =
				avgFrequentistPopulationsByRepetition.size() != 0 ?
						(((float) std::accumulate(
								avgFrequentistPopulationsByRepetition.begin(),
								avgFrequentistPopulationsByRepetition.end(), 0))
								/ ((float) avgFrequentistPopulationsByRepetition.size())) :
						0;
		float avgBayesianPopulationEvolvedInAllRepetitions =
				avgBayesianPopulationByRepetition.size() != 0 ?
						(((float) std::accumulate(
								avgBayesianPopulationByRepetition.begin(),
								avgBayesianPopulationByRepetition.end(), 0))
								/ ((float) avgBayesianPopulationByRepetition.size())) :
						0;

		bayesianAveragePopulationByEnvironmentHyperparameterPair.insert(
				std::make_pair(std::make_pair(environmentAlpha, environmentBeta),
						avgBayesianPopulationEvolvedInAllRepetitions));
		
		fileOfAvgBayesianPopulationByHyperparameterPair << environmentAlpha << "," << environmentBeta << ","
				<< avgBayesianPopulationEvolvedInAllRepetitions << endl;

	}
	fileOfAvgBayesianPopulationByHyperparameterPair.close();

	
	for (map<hyperparameterType, map<float, float> >::iterator frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter =
			frequencyOfBayesianAveragesByEnvironmentHyperparameterPair.begin();
			frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter
					!= frequencyOfBayesianAveragesByEnvironmentHyperparameterPair.end();
			frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter++) {
		float thisAlpha =
				(frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter->first).first;
		float thisBeta =
				(frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter->first).second;

		stringstream fileNameStream;
		fileNameStream << "fileOfBayesianAverageFrequencywith(" << thisAlpha << ","
				<< thisBeta << ").csv";
		ofstream fileOfFrequenciesOfBayesianAveragesByEnvHyperparameterPair;
		fileOfFrequenciesOfBayesianAveragesByEnvHyperparameterPair.open(
				fileNameStream.str().c_str());

		for (map<float, float>::iterator diter =
				(frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter->second).begin();
				diter
						!= (frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter->second).end();
				diter++) {
			float pop = diter->first;
			float freq = diter->second;

			for (int count = 1; count <= freq; count++) {
				fileOfFrequenciesOfBayesianAveragesByEnvHyperparameterPair << pop;
				fileOfFrequenciesOfBayesianAveragesByEnvHyperparameterPair << ",";
			}

			diter++;
			if (diter
					!= (frequencyOfBayesianAveragesByEnvironmentHyperparameterPairIter->second).end()) {
			}
			diter--;
		}

		fileOfFrequenciesOfBayesianAveragesByEnvHyperparameterPair.close();
	}
	return 0;
}