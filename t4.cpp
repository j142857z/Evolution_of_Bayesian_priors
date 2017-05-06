#include "common.h"

/**
 * Example call to the program.
 *
 * ./t8 100 100 2 1000 4 4 1 10 10
 *
 * */

float calculateFitnessOfBayesianOrFrequentistGenomeWithCostForBayesianAbility(
		GAGenome& g);

int main(int argc, char **argv) {
	if (argc < 9) {
		cout
				<< "<Population-Size> <Max-Allele-Value> <Maximum-Learning-Trials> "
				<< "<Num-Of-Generations> <Bayesian-cost> <Alpha> <Beta> <GenerationsPerChange> <ProbabilitiesByGeneration>"
				<< endl;
		return -1;
	}
	srand((unsigned) time(0));

	int populationSize = atoi(argv[1]);
	int maxAlleleValue = atoi(argv[2]);
	learningLength = atoi(argv[3]);
	int numberOfGenerations = atoi(argv[4]);
	alpha = atof(argv[5]);
	betaHyperparameter = atof(argv[6]);
	numberOfGenerationsPerChange = atoi(argv[7]);
	unsigned int numberOfProbabilityValuesPerGeneration = atoi(argv[8]);
	unsigned int numberOfReps = atoi(argv[9]);
	
	ofstream populationByLearningTypeFile;
	stringstream ss;
	ss << "BayesianPopulationPercentages.csv";
	populationByLearningTypeFile.open(ss.str().c_str());

	for(int rep=1; rep<=numberOfReps; rep++){
	  environment = new BetaEnvironment(alpha, betaHyperparameter,
			  numberOfGenerationsPerChange,
			  numberOfProbabilityValuesPerGeneration);
	  environment->change();

	  GAAlleleSet<float> learnerTypeAlleleSet;
	  learnerTypeAlleleSet.add(FREQUENTIST);
	  learnerTypeAlleleSet.add(BAYESIAN);

	  GAAlleleSetArray<float> alleleSetArray;
	  alleleSetArray.add(FLT_MIN, maxAlleleValue); //Allele encoding the alpha estimates
	  alleleSetArray.add(FLT_MIN, maxAlleleValue); //Allele encoding the beta estimates
	  alleleSetArray.add(learnerTypeAlleleSet); //Allele encoding the type of learning

	  GARealGenome genomePrototype(alleleSetArray,
			  calculateFitnessOfBayesianOrFrequentistGenomeWithCostForBayesianAbility);
	  genomePrototype.crossover(GARealUniformCrossover);
	  genomePrototype.mutator(GARealGaussianMutatorCustom);
	  genomePrototype.initializer(initialiseGenomeWithLearningType);

	  GASimpleGA ga(genomePrototype);
	  ga.populationSize(populationSize);
	  ga.initialize();
	  ga.nGenerations(numberOfGenerations);
	  ga.pMutation(MUTATIONRATE);
	  ga.pCrossover(CROSSOVERRATE);

	  std::vector<long> avgBayesians;
	  std::vector<long> avgFrequentists;

	  for (int i = 1; i <= numberOfGenerations; i++) {
		  ga.step();
		  environment->change();

		  int numberOfBayesiansObservedInThisGeneration = 0;
		  int numberOfFrequentistsObservedInThisGeneration = 0;
		  for (int k = 0; k < ga.population().size(); k++) {
			  genomePrototype = ga.population().individual(k);
			  if (genomePrototype.gene(2) == BAYESIAN) {
				  numberOfBayesiansObservedInThisGeneration++;
			  } else if (genomePrototype.gene(2) == FREQUENTIST) {
				  numberOfFrequentistsObservedInThisGeneration++;
			  } else {
				  cout << "This situation should never occur" << endl;
				  return -2;
			  }
		  }
		  avgBayesians.push_back(numberOfBayesiansObservedInThisGeneration);
		  avgFrequentists.push_back(numberOfFrequentistsObservedInThisGeneration);
	  }
	  populationByLearningTypeFile
			  << (((float) std::accumulate(avgBayesians.begin(),
					  avgBayesians.end(), 0)) / ((float) avgBayesians.size()))
			  << ","<<"\n";
	}
	populationByLearningTypeFile.close();
	return 0;
}

float calculateBayesianEstimateAndReturnSquaredErrorWithCost(GAGenome& g,
		float probability, int &successfulTrials, int &failedTrials) {
	GARealGenome & genome = (GARealGenome &) g;
	for (float numberOfLearningAttempts = 1.0;
			numberOfLearningAttempts <= learningLength;
			numberOfLearningAttempts = numberOfLearningAttempts + 1.0) {
		if (numberOfLearningAttempts > 0.0) {
			if (util->runBernoulliTrial(probability)) {
				successfulTrials++;
			} else {
				failedTrials++;
			}
		}
	}
	return pow(
			(calculateBayesianEstimate(genome, successfulTrials, failedTrials)
					- probability), 2);
}

float calculateBayesianFitnessWithSumOfSquaredErrorsWithCost(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	float score = calculateFitnessOfBayesianGenomeWithSumOfSquaredErrors(g);
	return score >= 0 ? score : 0;
}

float calculateFitnessOfBayesianOrFrequentistGenomeWithCostForBayesianAbility(
		GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	if (genome[2] == BAYESIAN) {
		return calculateBayesianFitnessWithSumOfSquaredErrorsWithCost(g);
	} else if (genome[2] == FREQUENTIST) {
		return calculateFitnessOfFrequentistGenomeWithSumOfSquaredErrors(g);
	} else {
		return 0;
	}
}

/**
 * GAPopulation already has one member function to write the whole population to a text file. But
 * genomes are written one per line.
 * */
void writePopulationToFile(GAPopulation population, STD_OSTREAM& ofstream) {

}
