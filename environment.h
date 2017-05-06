/*
 * environment.h
 *
 *  Created on: 20 Sep 2012
 *      Author: j142857z
 */

#ifndef ENVIRONMENT_H_
#define ENVIRONMENT_H_

#include <vector>

using namespace std;

class Environment {
public:
	virtual void change()=0;
	virtual vector<float> getCurrentListOfProbabilities()=0;
protected:
	unsigned int currentGeneration;
};

#endif /* ENVIRONMENT_H_ */
