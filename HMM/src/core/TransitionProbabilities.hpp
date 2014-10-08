/*
 * TransitionProbabilities.hpp
 *
 *  Created on: Oct 7, 2014
 *      Author: root
 */

#ifndef TRANSITIONPROBABILITIES_HPP_
#define TRANSITIONPROBABILITIES_HPP_

#include "models/IndelModel.hpp"

namespace EBC
{

class TransitionProbabilities
{
protected:

	double gapOpening;
	double gapExtension;

	double time;

	IndelModel* indelModel;

public:
	TransitionProbabilities(IndelModel* im);

	virtual ~TransitionProbabilities();

	void calculate();

	double getGapExtension() const
	{
		return gapExtension;
	}

	double getGapOpening() const
	{
		return gapOpening;
	}

	void setTime(double time)
	{
		this->time = time;
	}
};

} /* namespace EBC */

#endif /* TRANSITIONPROBABILITIES_HPP_ */
