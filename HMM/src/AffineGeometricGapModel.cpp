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

using namespace std;

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

	this->lambda = params[1];
	this->time = params[0];

	this->calculateGeometricProbability(params[1],params[0]);
	this->gapExtensionProbability = params[2];

	//std::cerr << " Opening " << this->gapOpeningProbability << "\t" << "Extension " << this->gapExtensionProbability << std::endl;

}

void AffineGeometricGapModel::summarize()
{
	cout << endl << " Affine geometric gap model parameters : " << endl;
	cout << "lambda\t\ttime\text. prob\topen. prob" << endl;
	cout << lambda <<"\t" << time << "\t" <<  gapExtensionProbability << "\t" << gapOpeningProbability << endl;
}

} /* namespace EBC */
