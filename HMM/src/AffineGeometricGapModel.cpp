/*
 * AffineGeometricGapModel.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "AffineGeometricGapModel.hpp"
#include <cmath>

namespace EBC
{

AffineGeometricGapModel::AffineGeometricGapModel(double lambda, double t, double gapExtenstionProbs)
{
	calculateGeometricProbability(lambda, t);
	this->gapExtensionProbability = gapExtenstionProbs;
	//time, lambda and gap extension
	this->paramsNumber = 3;
}

AffineGeometricGapModel::AffineGeometricGapModel()
{
	this->paramsNumber = 3;
}

AffineGeometricGapModel::~AffineGeometricGapModel()
{
	// TODO Auto-generated destructor stub
}

void AffineGeometricGapModel::calculateGeometricProbability(double lambda, double t)
{
	//FIXME - this might not be making sense!
	this->gapOpeningProbability = log(1-exp(lambda*t));
}

void AffineGeometricGapModel::setParameters(double* params)
{
	//FIXME - maybe vector instead ?
	//time goes first
	this->calculateGeometricProbability(params[1], params[0]);
	this->gapExtensionProbability = params[2];

}

} /* namespace EBC */
