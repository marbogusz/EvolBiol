/*
 * AffineGeometricGapModel.cpp
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#include "models/NegativeBinomialGapModel.hpp"
#include <cmath>
#include <iostream>
#include "core/Maths.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{

NegativeBinomialGapModel::NegativeBinomialGapModel() : IndelModel(Definitions::NBIndelParamCount)
{
	//logMode = true;
	this->parameterLoBounds[0] = this->parameterLoBounds[1] = 0.000001;
	this->parameterHiBounds[0] = 0.15;  //lambda limitation due to the model!
	this->parameterHiBounds[1] = 0.9;
}

NegativeBinomialGapModel::~NegativeBinomialGapModel()
{
	// TODO Auto-generated destructor stub
}

double NegativeBinomialGapModel::calculateGapOpening(double time)
{
	return 1-exp(-1*lambda*time);
}

double NegativeBinomialGapModel::calculateGapExtension(double time)
{
	return this->gapExtensionProbability;
}

void NegativeBinomialGapModel::setParameters(vector<double> vc)
{
	this->lambda = vc[0];
	this->gapExtensionProbability = vc[1];
}

void NegativeBinomialGapModel::calculate()
{
	calculateGeometricProbability(this->lambda, this->time);
}


void NegativeBinomialGapModel::calculateGeometricProbability(double lambda, double t)
{
	//FIXME - this might not be making sense!
	double exponent = 1-exp(-1*lambda*t);
	this->gapOpeningProbability = exponent;
}

void NegativeBinomialGapModel::setParameters(double* params)
{

	this->lambda = params[0];

	this->gapExtensionProbability = params[1];

	//std::cerr << " Opening " << this->gapOpeningProbability << "\t" << "Extension " << this->gapExtensionProbability << std::endl;

}

void NegativeBinomialGapModel::summarize()
{
	cout << endl << " Affine geometric gap model parameters : " << endl;
	cout << "lambda\text. prob" << endl;
	cout << lambda << "\t" <<  gapExtensionProbability << endl;
}

} /* namespace EBC */
