/*
 * IndelModel.hpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef INDELMODEL_HPP_
#define INDELMODEL_HPP_

namespace EBC
{

class IndelModel
{
protected:
	double gapExtensionProbability;

	double gapOpeningProbability;

	unsigned int paramsNumber;

	//bool logMode;

public:
	IndelModel();

	//set parameters - time + the rest of parameters
	virtual void setParameters(double*) = 0;

	double getGapExtensionProbability() const
	{
		return gapExtensionProbability;
	}

	double getGapOpeningProbability() const
	{
		return gapOpeningProbability;
	}

	unsigned int getParamsNumber() const
	{
		return paramsNumber;
	}
};

} /* namespace EBC */
#endif /* INDELMODEL_HPP_ */
