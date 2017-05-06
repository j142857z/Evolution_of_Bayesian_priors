/**
 * Example run where 100 individuals, with genomes encoding positive real values up to 100 and who
 * are allowed to take a maximum of 2 observations, evolve for 10,000 generations in an
 * environment where the environment is of size 1 and changes every 1,000 generations with values
 * sampled from a beta distribution with parameters 0.001 and 0.001.
 *
 * ./t2 100 100 2 10000 0.001 0.001 1000 1
 * */

#include "common.h"

int main(int argc, char **argv) {
	if (argc < 9) {
		cout
				<< "<Population-Size> <Max-Allele-Value> <Maximum-Learning-Trials> "
				<< "<Num-Of-Generations> <Alpha> <Beta> <GenerationsPerChange> <ProbabilitiesByGeneration>..."
				<< endl;
		return -1;
	}
// 	srand((unsigned) time(0));
	srand (time(NULL));

	int populationSize = atoi(argv[1]);
	int maxAlleleValue = atoi(argv[2]);
	learningLength = atoi(argv[3]);
	int numberOfGenerations = atoi(argv[4]);
	alpha = atof(argv[5]);
	betaHyperparameter = atof(argv[6]);
	numberOfGenerationsPerChange = atoi(argv[7]);
	unsigned int numberOfProbabilityValuesPerGeneration = atoi(argv[8]);

	/**
	 * The environment must be initialised and moved forward to the first generation because it is
	 * required by the fitness function in the tests done by the genetic algorithm object upon
	 * creation. These tests take place before the algorithm starts iterating.
	 * */
	environment = new BetaEnvironment(alpha, betaHyperparameter,
			numberOfGenerationsPerChange,
			numberOfProbabilityValuesPerGeneration);
	environment->change();

	GAAlleleSetArray<float> alleleSetArray;
	alleleSetArray.add(FLT_MIN, maxAlleleValue);
	alleleSetArray.add(FLT_MIN, maxAlleleValue);

	GARealGenome genomePrototype(alleleSetArray,
			calculateFitnessOfBayesianGenomeWithSumOfSquaredErrors);
	genomePrototype.crossover(GARealUniformCrossover);
	genomePrototype.mutator(GARealGaussianMutatorCustom);
	genomePrototype.initializer(initialiseGenomeWithNoLearningType);

	GASimpleGA ga(genomePrototype);
	ga.populationSize(populationSize);
	ga.initialize();
	ga.nGenerations(numberOfGenerations);
	ga.pMutation(MUTATIONRATE);
	ga.pCrossover(CROSSOVERRATE);

	/**
	 * The output file where the prior-to-learning estimates are written.
	 */
	ofstream estimatesOfSuccessProbabilityByGenerationFile;
	estimatesOfSuccessProbabilityByGenerationFile.open(
			estimatesOfSuccessProbabilityByGenerationFileName);

	/**
	 * The output file where the real values of the environment probability are written.
	 */
	ofstream objectiveProbabilitiesFile;
	objectiveProbabilitiesFile.open(objectiveProbabilitiesFileName);

	ofstream objectiveAlphaFile;
	objectiveAlphaFile.open(objectiveAlphaFileName);
	ofstream objectiveBetaFile;
	objectiveBetaFile.open(objectiveBetaFileName);
	ofstream estimatesOfAlphaByGenerationFile;
	estimatesOfAlphaByGenerationFile.open(estimatesOfAlphaByGenerationFileName);
	ofstream estimatesOfBetaByGenerationFile;
	estimatesOfBetaByGenerationFile.open(estimatesOfBetaByGenerationFileName);

	objectiveAlphaFile << 1 << "," << numberOfGenerations << "," << alpha;
	objectiveBetaFile << 1 << "," << numberOfGenerations << ","
			<< betaHyperparameter;

	for (int i = 1; i <= numberOfGenerations; i++) {
		ga.step();
		environment->change();

		vector<float> currentP = environment->getCurrentListOfProbabilities();
		for (vector<float>::iterator currentPIter = currentP.begin();
				currentPIter != currentP.end(); currentPIter++) {
			objectiveProbabilitiesFile << i << "," << i << "," << *currentPIter
					<< endl;
		}

		float alphaEstimates[ga.population().size()];
		float betaEstimates[ga.population().size()];

		/**
		 * Scan all genomes...
		 */
		for (int k = 0; k < ga.population().size(); k++) {
			genomePrototype = ga.population().individual(k);

			alphaEstimates[k] = genomePrototype.gene(0);
			betaEstimates[k] = genomePrototype.gene(1);

			/**
			 * ...calculate the prior-to-learning estimate of each one and add it to the file.
			 */
			estimatesOfSuccessProbabilityByGenerationFile
					<< calculateBayesianEstimate(genomePrototype, 0, 0);

			estimatesOfAlphaByGenerationFile << genomePrototype.gene(0);
			estimatesOfBetaByGenerationFile << genomePrototype.gene(1);

			/**
			 * If k is not the last genome, add a pipe to the current last line of the file...
			 */
			if (k != (ga.population().size() - 1)) {
				estimatesOfSuccessProbabilityByGenerationFile << ",";
				estimatesOfAlphaByGenerationFile << ",";
				estimatesOfBetaByGenerationFile << ",";
			}

			/**
			 * ...otherwise jump to the next line but only if this is not the last generation.
			 */
			else if (k == (ga.population().size() - 1)
					&& i != numberOfGenerations) {
				estimatesOfSuccessProbabilityByGenerationFile << "" << endl;
				estimatesOfAlphaByGenerationFile << "" << endl;
				estimatesOfBetaByGenerationFile << "" << endl;
			}
		}
	}
	estimatesOfSuccessProbabilityByGenerationFile.close();
	objectiveProbabilitiesFile.close();
	objectiveAlphaFile.close();
	objectiveBetaFile.close();
	estimatesOfAlphaByGenerationFile.close();
	estimatesOfBetaByGenerationFile.close();

	return 0;
}
