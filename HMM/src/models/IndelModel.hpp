/*
 * IndelModel.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef INDELMODEL_HPP_
#define INDELMODEL_HPP_

#include <vector>

using namespace std;

namespace EBC
{

class IndelModel
{
protected:
	double gapExtensionProbability;

	double gapOpeningProbability;

	//divergence time
	double time;

	unsigned int paramsNumber;

	vector<double> parameterHiBounds;
	vector<double> parameterLoBounds;

	//bool logMode;

public:
	IndelModel(unsigned int);

	virtual double calculateGapOpening(double time) = 0;
	virtual double calculateGapExtension(double time) = 0;


	//set parameters - time + the rest of parameters
	virtual void setParameters(double*) = 0;

	virtual void setParameters(vector<double>)=0;

	void setTime(double t)
	{
		this->time = t;
	}

	virtual void summarize()=0;

	virtual void calculate()=0;

	double getGapExtensionProbability() const
	{
		return gapExtensionProbability;
	}

	double getGapOpeningProbability() const
	{
		return gapOpeningProbability;
	}

	virtual double* getParameters()=0;

	unsigned int getParamsNumber() const
	{
		return paramsNumber;
	}

	inline double getHiBound(unsigned int pos)
	{
		return parameterHiBounds[pos];
	}

	inline double getLoBound(unsigned int pos)
	{
		return parameterLoBounds[pos];
	}
};

} /* namespace EBC */
#endif /* INDELMODEL_HPP_ */
