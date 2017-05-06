#include "util.h"

Util::Util() {
}

Util::~Util() {
}

float Util::generateRandomNumber() {
	return ((float) rand()) / (float) RAND_MAX;
}

float Util::generateRandomNumber(float lowerLimit, float upperLimit) {
	return lowerLimit + (this->generateRandomNumber()) * (upperLimit
			- lowerLimit);
}

bool Util::runBernoulliTrial(float p) {
	if (generateRandomNumber() <= p) {
		return true;
	}
	return false;
}

/**
 * This function generates Beta-distributed pseudorandom numbers with the given hyperparameters.
 * The seed must be updated with a call to the srand function.
 *
 * Solution taken from:
 * http://goo.gl/JaAxR
 * */
double Util::generateBeta(double alphaArg, double betaArg) {
	double randFromUnif = generateRandomNumber(0, 1);
	beta_distribution<> betaDistribution(alphaArg, betaArg);
	return quantile(betaDistribution, randFromUnif);
}

/**
 * An auxiliary function to calculate the mean of a series of values stored in an array.
 * */
double Util::getMean(double values[], int length) {
	double sum = 0.0;
	for (int i = 0; i < length; i++) {
		sum += values[i];
	}
	return sum / length;
}

/**
 * An auxiliary function to calculate the standard deviation of a series of values stored in an
 * array.
 * */
double Util::getVariance(double values[], int length, double avg) {
	double sumOfSquareDifferences = 0.0;
	for (int i = 0; i < length; i++) {
		sumOfSquareDifferences += pow((values[i] - avg), 2.0);
	}
	return sqrt(sumOfSquareDifferences / length);
}

std::vector<std::string> Util::split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while(std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> Util::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    return split(s, delim, elems);
}
