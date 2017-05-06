/*
 * Environment.h
 *
 *  Created on: 20 Sep 2012
 *      Author: j142857z
 */

#ifndef BETAENVIRONMENT_H_
#define BETAENVIRONMENT_H_

#include <vector>
#include "environment.h"
#include "util.h"

using namespace std;

class BetaEnvironment : public Environment {
public:
	BetaEnvironment(float alpha, float beta, unsigned int generationsPerChange,
			unsigned int numberOfProbabilityChangesPerGeneration);
	virtual ~BetaEnvironment();
	void change();
	vector<float> getCurrentListOfProbabilities();
private:
	float alpha;
	float beta;
	float numberOfProbabilityChangesPerGeneration;
	int generationsPerChange;
	int currentProbability;
	vector<float> * current;
	int generations;
};

#endif /* BETAENVIRONMENT_H_ */
