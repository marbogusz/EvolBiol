/*
 * AffineGeometricGapModel.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "AffineGeometricGapModel.hpp"
#include <cmath>
#include <iostream>
#include "Maths.hpp"

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
	//logMode = true;
	this->paramsNumber = 3;
}

AffineGeometricGapModel::~AffineGeometricGapModel()
{
	// TODO Auto-generated destructor stub
}

void AffineGeometricGapModel::calculateGeometricProbability(double lambda, double t)
{
	//FIXME - this might not be making sense!
	double exponent = 1-exp(-1*lambda*t);
	this->gapOpeningProbability = exponent;
}

void AffineGeometricGapModel::setParameters(double* params)
{
	//FIXME - maybe vector instead ?
	//time goes first
	this->calculateGeometricProbability(Maths::logistic(params[1]),
			Maths::logistic(params[0]));
	this->gapExtensionProbability = Maths::logistic(params[2]);

	//std::cerr << " Opening " << this->gapOpeningProbability << "\t" << "Extension " << this->gapExtensionProbability << std::endl;

}

} /* namespace EBC */
