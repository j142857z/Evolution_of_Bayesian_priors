#ifndef UTIL_H_
#define UTIL_H_

#include <boost/math/distributions/beta.hpp>

using namespace boost::math;

class Util {
public:
	Util();
	virtual ~Util();
	float generateRandomNumber();
	float generateRandomNumber(float lowerLimit, float upperLimit);
	bool runBernoulliTrial(float p);
	double generateBeta(double alphaArg, double betaArg);
	double getMean(double values[], int length);
	double getVariance(double values[], int length, double avg);
	std::vector<std::string> split(const std::string &s, char delim,
			std::vector<std::string> &elems);
	std::vector<std::string> split(const std::string &s, char delim);

	int round(float d) {
		return (int) floor(d + 0.5);
	}
};

#endif /* UTIL_H_ */
