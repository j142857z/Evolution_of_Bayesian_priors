#include "common.h"

int main(int argc, char **argv) {
  
	srand((unsigned) time(0));

	int populationSize = atoi(argv[1]);
	int maxAlleleValue = atof(argv[2]);
	int minNumberOfObservations = atof(argv[3]);
	int maxNumberOfObservations = atof(argv[4]);
	int numberOfGenerations = atoi(argv[5]);
	int numberOfRepetitions = atoi(argv[6]);
	numberOfGenerationsPerChange = atoi(argv[7]);
	unsigned int numberOfProbabilityValuesPerGeneration = atoi(argv[8]);
	
	float minAlpha = atof(argv[9]);
	float maxAlpha = atof(argv[10]);
	float minBeta = atof(argv[11]);
	float maxBeta = atof(argv[12]);
	float numberOfHyperparameterPoints = atoi(argv[13]);
	
	float spaceBetweenAlphaPoints = (maxAlpha - minAlpha) / (numberOfHyperparameterPoints - 1);
	float spaceBetweenBetaPoints = (maxBeta - minBeta) / (numberOfHyperparameterPoints - 1);

	vector<hyperparameterType> listOfEnvironmentHyperparameterPairs = getListOfHyperparameterPairs(minAlpha, maxAlpha, spaceBetweenAlphaPoints, minBeta, maxBeta, spaceBetweenBetaPoints);

	map<hyperparameterType, float> bayesianAveragePopulationByEnvironmentHyperparameterPair;
	map<hyperparameterType, map<float, float> > frequencyOfBayesianAveragesByEnvironmentHyperparameterPair;

	ofstream fileOfAvgBayesianPopulationByHyperparameterPair;
	fileOfAvgBayesianPopulationByHyperparameterPair.open("fileOfAvgBayesianPopulationByHyperparameterPair.csv");
	
	for(learningLength=minNumberOfObservations; learningLength<=maxNumberOfObservations; learningLength++){
	  for (int i = 0; i < listOfEnvironmentHyperparameterPairs.size(); i++) {
	    float environmentAlpha = listOfEnvironmentHyperparameterPairs[i].first;
	    float environmentBeta = listOfEnvironmentHyperparameterPairs[i].second;

	    std::vector<long> avgBayesianPopulationByRepetition;
	    std::vector<long> avgFrequentistPopulationsByRepetition;
	    for (int repetition = 1; repetition <= numberOfRepetitions; repetition++) {
	      runEvolutionarySimulation(environmentAlpha, environmentBeta, populationSize, maxAlleleValue, numberOfGenerations, numberOfGenerationsPerChange, numberOfProbabilityValuesPerGeneration, avgBayesianPopulationByRepetition, avgFrequentistPopulationsByRepetition, frequencyOfBayesianAveragesByEnvironmentHyperparameterPair);
	    }

	    float avgFrequentistPopulationEvolvedInAllRepetitions = avgFrequentistPopulationsByRepetition.size() != 0 ? (((float) std::accumulate(avgFrequentistPopulationsByRepetition.begin(), avgFrequentistPopulationsByRepetition.end(), 0)) / ((float) avgFrequentistPopulationsByRepetition.size())) : 0;

	    float avgBayesianPopulationEvolvedInAllRepetitions = avgBayesianPopulationByRepetition.size() != 0 ? (((float) std::accumulate(avgBayesianPopulationByRepetition.begin(), avgBayesianPopulationByRepetition.end(), 0)) / ((float) avgBayesianPopulationByRepetition.size())) : 0;

	    bayesianAveragePopulationByEnvironmentHyperparameterPair.insert(std::make_pair(std::make_pair(environmentAlpha, environmentBeta), avgBayesianPopulationEvolvedInAllRepetitions));
	    
	    double environmentEntropy = calculateEntropy(environmentAlpha, environmentBeta);
	    
	    fileOfAvgBayesianPopulationByHyperparameterPair << learningLength << ","  << environmentAlpha << "," << environmentBeta << "," << environmentEntropy << "," << avgBayesianPopulationEvolvedInAllRepetitions << endl;
	  }
	}
	fileOfAvgBayesianPopulationByHyperparameterPair.close();

	return 0;
}