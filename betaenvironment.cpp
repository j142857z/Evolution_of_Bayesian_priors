/*
 * Environment.cpp
 *
 *  Created on: 20 Sep 2012
 *      Author: j142857z
 */

#include "betaenvironment.h"

/**
 * The environment changes every [frequencyOfChanges] generations, starting with the first and each
 * generation [numberOfProbabilityChangesPerGeneration] values are available.
 * */
BetaEnvironment::BetaEnvironment(float alpha, float beta,
		unsigned int generationsPerChange,
		unsigned int numberOfProbabilityValuesPerGeneration) {
	this->alpha = alpha;
	this->beta = beta;
	currentProbability = 0;
	this->numberOfProbabilityChangesPerGeneration =
			numberOfProbabilityValuesPerGeneration;
	current = new vector<float>;
	generations = 0;
	this->generationsPerChange = generationsPerChange;
}

BetaEnvironment::~BetaEnvironment() {
	delete current;
}

void BetaEnvironment::change() {
	generations++;
	if ((generations % generationsPerChange) == 1 || generations == 1 || generationsPerChange == 1) {
		Util util;
		current->clear();
		for (int i = 1; i <= numberOfProbabilityChangesPerGeneration; i++) {
			double newValue = util.generateBeta(alpha, beta);
			current->push_back(newValue);
		}
	}
}

vector<float> BetaEnvironment::getCurrentListOfProbabilities() {
	return *current;
}

