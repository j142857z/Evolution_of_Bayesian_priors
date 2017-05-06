#ifndef COMMON_H_
#define COMMON_H_

#include <ga/ga.h>
#include <ga/GARealGenome.h>
#include <boost/random.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <float.h>
#include <map>
#include <math.h>
#include "betaenvironment.h"



/**
 * Allele value indicating a genome is a bayesian learner.
 * */
const float BAYESIAN = 0;
const float FREQUENTIST = 1;

/**
 * A mutation rate equal to 0.01 yields more stable, separate populations. But the 0.05 rate
 * is the one that produces the best results in the other types of experiments.
 */
const float MUTATIONRATE = 0.05;
const float CROSSOVERRATE = 0.95;

//vector<float> * LOLDenvironment = new vector<float>;

/**
 * The environment must be initialised appropriately by each main function but it must remain
 * visible to the methods in this file.
 * */
Environment * environment;

/**
 * The name of the output file where the prior-to-learning estimates are written.
 */
// char estimatesOfSuccessProbabilityByGenerationFileName[] =
// 		"EstimatesOfSuccessProbabilityByGeneration.dat";
char estimatesOfSuccessProbabilityByGenerationFileName[] =
		"PriorsByGeneration.csv";


/**
 * The name of the output file where the real values of the LOLDenvironment probability are written.
 */
// char objectiveProbabilitiesFileName[] = "ObjectiveProbabilities.dat";
char objectiveProbabilitiesFileName[] = "EnvironmentStatesByGeneration.csv";

/**
 * The name of the output file where the real values of the LOLDenvironment probability are written.
 */
char objectiveAlphaFileName[] = "ObjectiveAlpha.csv";

/**
 * The name of the output file where the real values of the LOLDenvironment probability are written.
 */
char objectiveBetaFileName[] = "ObjectiveBeta.csv";

/**
 * The name of the output file where the prior-to-learning estimates are written.
 */
char estimatesOfAlphaByGenerationFileName[] = "PopulationAlphaByGeneration.csv";

/**
 * The name of the output file where the prior-to-learning estimates are written.
 */
char estimatesOfBetaByGenerationFileName[] = "PopulationBetaByGeneration.csv";

/**
 * The name of the output file where the prior-to-learning estimates are written.
 */
char estimatesOfAlphaAVGByGenerationFileName[] =
		"EstimatesOfAlphaAVGByGeneration.dat";

/**
 * The name of the output file where the prior-to-learning estimates are written.
 */
char estimatesOfBetaAVGByGenerationFileName[] =
		"EstimatesOfBetaAVGByGeneration.dat";

/**
 * The number of times the environment probability changes during the learning process and fitness evaluation.
 * If zero then all individuals carry out a single learning process with the current LOLDenvironment probability.
 */
float numberOfGenerationsPerChange = 0;

/**
 * The probability of ocurrence of the event the population is expected to learn. It is initialised as 0.5
 * because it is used in the fitness function implementation which in turn is used by the GA library during
 * initialisation, prior to actual evolution. This situation does not affect the simulation since the variable
 * is changed properly when evolution is going to begin.
 */
float environmentProbability = 0.5;

/**
 * The alpha hyperparameter of the Beta distribution that generates the LOLDenvironment probability.
 */
float alpha = 1;

/**
 * The beta hyperparameter of the Beta distribution that generates the LOLDenvironment probability.
 */
float betaHyperparameter = 1;

/**
 * The maximum difference admitted between the environment probability and an individual's own estimate for both
 * values to be considered sufficiently equal and the individual to be deemed to have been able to learn.
 */
float maximumProbabilityEstimateError = 0.025;

/**
 * The maximum number of learning trials each individual is entitled to attempt in order to estimate the
 * value of the LOLDenvironment probability.
 */
float learningLength = 1000;
/**
 * The random number generator used for generating initial allele values.
 */
boost::mt19937 rng(std::time(0));

/**
 * The Normal distribution used for generating initial allele values.
 */
boost::normal_distribution<> normalDistribution(0.0, 0.05);

/**
 * The function used for generating Normal-distributed initial allele values.
 */
boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > generateNormalPseudorandomNumber(
		rng, normalDistribution);

/**
 * In order to run the complete beta function, parameter x must be 1.
 */
double beta(double a, double b, double x);

double digamma(double z);

float calculateBayesianEstimate(GARealGenome& genome, int successfulTrials,
		int failedTrials);

float evaluateBayesianGenome(GAGenome& g, float probability,
		int &successfulTrials, int &failedTrials);

float evaluateBayesianGenome(GAGenome& g, float probability);

float evaluateFrequentistGenome(GAGenome& g, float probability,
		int &successfulTrials, int &failedTrials);

float evaluateFrequentistGenome(GAGenome& g, float probability);

Util * util = new Util;

/*
 * Fitness function to be used as part of a GARealGenome prototype. It evaluates a bayesian genome
 * by doing the following.
 *
 * 1. The list of probabilities is taken from the environment. This list is assumed to be the same
 * for all agents in the same generations.
 *
 * 2. For each probability the agent makes a maximum of [learningLength] observations
 * of the event E or until the agent's bayesian estimate is assumed to be sufficiently close. This
 * is measured with the variable [maximumProbabilityEstimateError].
 *
 * 3. The agent get's a partial fitness for evaluating each probability. The final score is
 * is calculated as the average of the partial measures.
 * */
float calculateFitnessOfBayesianGenome(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	int successes = 0;
	int failures = 0;
	float summedFitness = 0.0;
	vector<float> values = environment->getCurrentListOfProbabilities();
	for (int k = 0; k < values.size(); k++) {
		summedFitness += evaluateBayesianGenome(g, values.at(k), successes,
				failures);
	}
	return summedFitness / (float) values.size();
}

/*
 * It evaluates a frequentist genome by doing the following.
 *
 * 1. The list of probabilities is taken from the environment. This list is assumed to be the same
 * for all agents in the same generations.
 *
 * 2. For each probability the agent makes a maximum of [learningLength] observations
 * of the event E or until the agent's frequentist estimate is assumed to be sufficiently close.
 * This is measured with the variable [maximumProbabilityEstimateError].
 *
 * 3. The agent get's a partial fitness for evaluating each probability. The final score is
 * is calculated as the average of the partial measures.
 *
 * This function is not expected to be used as part of a GARealGenome prototype. But it could be
 * used as such. The reason for not doing so is that no programs are expected to run with a
 * population of only frequentists.
 * */
float calculateFitnessOfFrequentistGenome(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	int successes = 0;
	int failures = 0;
	float summedFitness = 0.0;
	vector<float> values = environment->getCurrentListOfProbabilities();
	for (int k = 0; k < values.size(); k++) {
		summedFitness += evaluateFrequentistGenome(g, values.at(k), successes,
				failures);
	}
	return summedFitness / (float) values.size();

}

/**
 * This function evaluates the (partial) fitness of a bayesian agent in terms of how it manages to
 * estimate one probability by repeating observations of the corresponding E experiment.
 *
 * fitness = 1 - t/tmax
 *
 * Where tmax = learningLength.
 *
 * This is an auxiliary function. Not suitable for a GARealGenome prototype.
 *
 * This function is expected to be used to evaluate a partial fitness of a single genome. Because
 * of this, this function received two references corresponding to the number of observations
 * (success and failures) made by that genome up to that point.
 */
float evaluateBayesianGenome(GAGenome& g, float probability,
		int &successfulTrials, int &failedTrials) {
	GARealGenome & genome = (GARealGenome &) g;
	for (float numberOfLearningAttempts = 0.0;
			numberOfLearningAttempts <= learningLength;
			numberOfLearningAttempts = numberOfLearningAttempts + 1.0) {
		if (numberOfLearningAttempts > 0.0) {
			if (util->runBernoulliTrial(probability)) {
				successfulTrials++;
			} else {
				failedTrials++;
			}
		}
		if (fabs(
				calculateBayesianEstimate(genome, successfulTrials,
						failedTrials) - probability)
				<= maximumProbabilityEstimateError) {
			return (1.0
					- (numberOfLearningAttempts / learningLength));
		}
	}
	return 0.0;
}

/**
 * This function evaluates the (partial) fitness of a frequentist agent in terms of how it manages
 * to estimate one probability by repeating observations of the corresponding E experiment.
 *
 * fitness = 1 - t/tmax
 *
 * Where tmax = learningLength.
 *
 * This is an auxiliary function. Not suitable for a GARealGenome prototype.
 *
 * This function is expected to be used to evaluate a partial fitness of a single genome. Because
 * of this, this function received two references corresponding to the number of observations
 * (success and failures) made by that genome up to that point.
 */
float evaluateFrequentistGenome(GAGenome& g, float probability,
		int &successfulTrials, int &failedTrials) {
	GARealGenome & genome = (GARealGenome &) g;
	for (float numberOfLearningAttempts = 1.0;
			numberOfLearningAttempts <= learningLength;
			numberOfLearningAttempts = numberOfLearningAttempts + 1.0) {
		if (util->runBernoulliTrial(probability)) {
			successfulTrials++;
		} else {
			failedTrials++;
		}
		float estimate = ((float) successfulTrials
				/ ((float) successfulTrials + (float) failedTrials));
		if (fabs(estimate - probability) <= maximumProbabilityEstimateError) {
			return (1.0
					- (numberOfLearningAttempts / learningLength));
		}
	}
	return 0.0;
}

/**
 * Prototype fitness function to evaluate both bayesians and frequentists.
 *
 * This function works only with genomes whose third gene codes the type of learning: bayesian or
 * frequentist.
 *
 * Neither bayesians nor frequentists pay any type of cost.
 *
 * Both types of agents are evaluated with a Hinton-like mathematical function, in terms of the
 * number of observations an agent required to estimate individually each of the probabilities
 * of the environment in the current generation. Total fitness is the average of the partial
 * scores after estimating each probability.
 *
 * We can't pass the the list of environment values as a parameter because we need the function to
 * comply with the requirements of the GALib library. Then the environment variable must be set and
 * updated accordingly before calling this function or the step function of the genetic algorithm
 * object.
 * */
float calculateFitnessOfBayesianOrFrequentistGenome(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	if (genome[2] == BAYESIAN) {
		return calculateFitnessOfBayesianGenome(g);
	} else if (genome[2] == FREQUENTIST) {
		return calculateFitnessOfFrequentistGenome(g);
	} else {
		return 0;
	}
}

void initialiseGenomeWithNoLearningType(GAGenome &genome) {
	float initialAlpha = 1 + generateNormalPseudorandomNumber();
	float initialBeta = 1 + generateNormalPseudorandomNumber();
	if (initialAlpha < 0) {
		initialAlpha = FLT_MIN;
	}
	if (initialBeta < 0) {
		initialBeta = FLT_MIN;
	}
	((GARealGenome &) genome)[0] = initialAlpha;
	((GARealGenome &) genome)[1] = initialBeta;
}

void initialiseGenomeWithLearningType(GAGenome &genome) {
	initialiseGenomeWithNoLearningType(genome);
	/**
	 * There are only two possible types of learners and initialisation can be made with the flip
	 * of a coin.
	 * */
	((GARealGenome &) genome)[2] = FREQUENTIST;
	if (util->runBernoulliTrial(0.5)) {
		((GARealGenome &) genome)[2] = BAYESIAN;
	}
}

/**
 * Function used to mutate a genome using Gaussian values. Each gene is mutated with the
 * probability passed as argument.
 *
 * Gaussian values are generated with the global function generateNormalPseudorandomNumber.
 *
 * This function is an almost complete copy of the implementation provided in the GARealGenome
 * class source file (GARealGenome.C). The only difference is that this copy uses the
 * generateNormalPseudorandomNumber function because the original implementation uses a Normal
 * distribution with fixed parameters (0,1). With this local implementation we can adjust the
 * distribution with different, custom parameters.
 */
int GARealGaussianMutatorCustom(GAGenome& g, float mutationProbability) {
	GA1DArrayAlleleGenome<float> &child =
			DYN_CAST(GA1DArrayAlleleGenome<float> &, g);
			register int n, i;
			if (mutationProbability <= 0.0)
			return (0);

			float nMut = mutationProbability * (float) (child.length());
			int length = child.length() - 1;
			if (nMut < 1.0) { // we have to do a flip test on each element
				nMut = 0;
				for (i = length; i >= 0; i--) {
					float value = child.gene(i);
					if (GAFlipCoin(mutationProbability)) {
						if (child.alleleset(i).type() == GAAllele::ENUMERATED
								|| child.alleleset(i).type() == GAAllele::DISCRETIZED)
						value = child.alleleset(i).allele();
						else if (child.alleleset(i).type() == GAAllele::BOUNDED) {
							/**
							 * The following line is one of two changes made to the original function.
							 * The other one is below.
							 * */
							value += generateNormalPseudorandomNumber();
							value = GAMax(child.alleleset(i).lower(), value);
							value = GAMin(child.alleleset(i).upper(), value);
						}
						child.gene(i, value);
						nMut++;
					}
				}
			} else { // only mutate the ones we need to
				for (n = 0; n < nMut; n++) {
					int idx = GARandomInt(0, length);
					float value = child.gene(idx);
					if (child.alleleset(idx).type() == GAAllele::ENUMERATED
							|| child.alleleset(idx).type() == GAAllele::DISCRETIZED)
					value = child.alleleset(idx).allele();
					else if (child.alleleset(idx).type() == GAAllele::BOUNDED) {
						/**
						 * The following line is one of two changes made to the original function.
						 * The other one is above.
						 * */
						value += generateNormalPseudorandomNumber();
						value = GAMax(child.alleleset(idx).lower(), value);
						value = GAMin(child.alleleset(idx).upper(), value);
					}
					child.gene(idx, value);
				}
			}
			return ((int) nMut);
		}

	/**
	 * This function calculates the bayesian estimate of a probability given a number of observed
	 * occurrences and non-occurrences of event A.
	 *
	 * This function assumes the genome passed as argument is a bayesian. No check is done. Use of this
	 * function on a frequentist genome will lead to incorrect results.
	 * */
float calculateBayesianEstimate(GARealGenome& genome, int successfulTrials,
		int failedTrials) {
	return (successfulTrials + genome.gene(0))
			/ (successfulTrials + failedTrials + genome.gene(0) + genome.gene(1));
}

/**
 * This function receives a (bayesian) genome with the number of (Bernoulli) successes and failures
 * it has observed in the current generation along with a probability value to estimate.
 *
 * The function returns the squared error between the estimate after
 * [learningLength] and the actual value of p.
 * */
float calculateBayesianEstimateAndReturnSquaredError(GAGenome& g,
		float probability, int &successfulTrials, int &failedTrials) {
	GARealGenome & genome = (GARealGenome &) g;
	for (float numberOfLearningAttempts = 1.0;
			numberOfLearningAttempts <= learningLength;
			numberOfLearningAttempts = numberOfLearningAttempts + 1.0) {
		if (util->runBernoulliTrial(probability)) {
			successfulTrials++;
		} else {
			failedTrials++;
		}
	}
	return pow(
			(calculateBayesianEstimate(genome, successfulTrials, failedTrials)
					- probability), 2);
}

float calculateFrequentistEstimateAndReturnSquaredError(GAGenome& g,
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
	float frequentistEstimate = ((float) ((float) successfulTrials)
			/ ((float) successfulTrials + (float) failedTrials));
	return pow((frequentistEstimate - probability), 2);
}

/*
 * Fitness function to be used as part of a GARealGenome prototype. It evaluates a bayesian genome
 * by doing the following.
 *
 * 1. The list of probabilities is taken from the environment. This list is assumed to be the same
 * for all agents in the same generations.
 *
 * 2. For each probability the agent makes a total of [learningLength] observations
 * of the event E.
 *
 * 3. The sum of squared errors (phi) between each value of p and the corresponding estimate is
 * calculated. The final fitness score is
 *
 * fitness = 1 / (1 + phi)
 * */
float calculateFitnessOfBayesianGenomeWithSumOfSquaredErrors(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	int successes = 0;
	int failures = 0;
	float sumOfSquaredErrors = 0.0;
	vector<float> values = environment->getCurrentListOfProbabilities();
	for (int k = 0; k < values.size(); k++) {
		sumOfSquaredErrors += calculateBayesianEstimateAndReturnSquaredError(g,
				values.at(k), successes, failures);
	}
	return 1 / (1 + sumOfSquaredErrors);
}

/*
 * Fitness function that can be used as part of a GARealGenome prototype. It evaluates a
 * frequentist genome by doing the following.
 *
 * 1. The list of probabilities is taken from the environment. This list is assumed to be the same
 * for all agents in the same generations.
 *
 * 2. For each probability the agent makes a total of [learningLength] observations
 * of the event E.
 *
 * 3. The sum of squared errors (phi) between each value of p and the corresponding estimate is
 * calculated. The final fitness score is
 *
 * fitness = 1 / (1 + phi)
 * */
float calculateFitnessOfFrequentistGenomeWithSumOfSquaredErrors(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	int successes = 0;
	int failures = 0;
	float sumOfSquaredErrors = 0.0;
	vector<float> values = environment->getCurrentListOfProbabilities();
	for (int k = 0; k < values.size(); k++) {
		sumOfSquaredErrors += calculateFrequentistEstimateAndReturnSquaredError(
				g, values.at(k), successes, failures);
	}
	return 1 / (1 + sumOfSquaredErrors);
}

float calculateBayesianFitnessAsSumOfSquaredErrors(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	float score = calculateFitnessOfBayesianGenomeWithSumOfSquaredErrors(g);
	return score >= 0 ? score : 0;
}

typedef pair<float, float> hyperparameterType;

vector<hyperparameterType> getListOfHyperparameterPairs(float minAlpha, float maxAlpha, float spaceBetweenAlphaPoints, float minBeta, float maxBeta, float spaceBetweenBetaPoints) {
  vector<hyperparameterType> listOfEnvironmentHyperparameterPairs;
  float xAlpha = minAlpha;
  while (xAlpha <= maxAlpha) {
    float xBeta = minBeta;
    while (xBeta <= maxBeta) {
      listOfEnvironmentHyperparameterPairs.push_back(make_pair(xAlpha, xBeta));
      xBeta += spaceBetweenBetaPoints;
    }
    xAlpha += spaceBetweenAlphaPoints;
  }  
  return listOfEnvironmentHyperparameterPairs;
}

float calculateFitnessOfBayesianOrFrequentistGenomeAsSumOfSquaredErrors(GAGenome& g) {
	GARealGenome & genome = (GARealGenome &) g;
	if (genome[2] == BAYESIAN) {
		return calculateBayesianFitnessAsSumOfSquaredErrors(g);
	} else if (genome[2] == FREQUENTIST) {
		return calculateFitnessOfFrequentistGenomeWithSumOfSquaredErrors(g);
	} else {
		return 0;
	}
}

void runEvolutionarySimulation(float environmentAlpha, float environmentBeta, int populationSize, float maxAlleleValue, int numberOfGenerations, int numberOfGenerationsPerChange, int numberOfProbabilityValuesPerGeneration, std::vector<long> &avgBayesianPopulationByRepetition, std::vector<long> &avgFrequentistPopulationsByRepetition, map<hyperparameterType, map<float, float> > &frequencyOfBayesianAveragesByEnvironmentHyperparameterPair) {
/**
			 * The environment must be initialised and moved forward to the first generation because it is
			 * required by the fitness function in the tests done by the genetic algorithm object upon
			 * creation. These tests take place before the algorithm starts iterating.
			 * */
			environment = new BetaEnvironment(environmentAlpha, environmentBeta,
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
					calculateFitnessOfBayesianOrFrequentistGenomeAsSumOfSquaredErrors);
			genomePrototype.crossover(GARealUniformCrossover);
			genomePrototype.mutator(GARealGaussianMutatorCustom);
			genomePrototype.initializer(initialiseGenomeWithLearningType);

			GASimpleGA ga(genomePrototype);
			ga.populationSize(populationSize);
			ga.initialize();
			ga.nGenerations(numberOfGenerations);
			ga.pMutation(MUTATIONRATE);
			ga.pCrossover(CROSSOVERRATE);

			std::vector<long> bayesianPopulationByGeneration;
			std::vector<long> frequentistPopulationByGeneration;

			for (int i = 1; i <= numberOfGenerations; i++) {
				ga.step();
				environment->change();

				int numberOfBayesiansObservedInThisGeneration = 0;
				int numberOfFrequentistsObservedInThisGeneration = 0;
				for (int k = 0; k < ga.population().size(); k++) {
					genomePrototype = ga.population().individual(k);
					if (genomePrototype.gene(2) == BAYESIAN) {
						numberOfBayesiansObservedInThisGeneration++;
					} else { // genomePrototype.gene(2) == FREQUENTIST
						numberOfFrequentistsObservedInThisGeneration++;
					}
				}
				bayesianPopulationByGeneration.push_back(
						numberOfBayesiansObservedInThisGeneration);
				frequentistPopulationByGeneration.push_back(
						numberOfFrequentistsObservedInThisGeneration);
			}
			float avgBayesiansInThisRepetition = (((float) std::accumulate(
					bayesianPopulationByGeneration.begin(),
					bayesianPopulationByGeneration.end(), 0))
					/ ((float) bayesianPopulationByGeneration.size()));
			
			avgBayesianPopulationByRepetition.push_back(
					avgBayesiansInThisRepetition);
			
			float avgFrequentistsInThisRepetition = (((float) std::accumulate(
					frequentistPopulationByGeneration.begin(),
					frequentistPopulationByGeneration.end(), 0))
					/ ((float) frequentistPopulationByGeneration.size()));
			
			avgFrequentistPopulationsByRepetition.push_back(
					avgFrequentistsInThisRepetition);

			frequencyOfBayesianAveragesByEnvironmentHyperparameterPair[std::make_pair(
					environmentAlpha, environmentBeta)][floor(
					avgBayesiansInThisRepetition)]++;
}

double calculateEntropy(float environmentAlpha, float environmentBeta){
  double entropy = log(beta(environmentAlpha, environmentBeta, 1.0)) - (environmentAlpha - 1)*digamma(environmentAlpha) -(environmentBeta - 1)*digamma(environmentBeta) + (environmentAlpha + environmentBeta - 2)*digamma(environmentAlpha + environmentBeta);  
  return entropy;
}

#endif /* COMMON_H_ */
