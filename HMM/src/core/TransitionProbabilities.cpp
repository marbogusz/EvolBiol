/*
 * TransitionProbabilities.cpp
 *
 *  Created on: Oct 7, 2014
 *      Author: root
 */

#include <core/TransitionProbabilities.hpp>

namespace EBC
{

TransitionProbabilities::TransitionProbabilities(IndelModel* im) : indelModel(im)
{
}

TransitionProbabilities::~TransitionProbabilities()
{
}

void TransitionProbabilities::calculate()
{
	this->gapExtension = indelModel->calculateGapExtension(this->time);
	this->gapOpening = indelModel->calculateGapOpening(this->time);
}

} /* namespace EBC */
