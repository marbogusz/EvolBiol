/*
 * AminoacidSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "AminoacidSubstitutionModel.hpp"

namespace EBC
{

AminoacidSubstitutionModel::AminoacidSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha, Definitions::aaModelDefinition modelDef) :
	SubstitutionModelBase(dict,alg,alpha), eigenDecomposed(false)
{
	this->paramsNumber = 1;
	this->parameters = new double[1];

	this->maxRate = modelDef.maxRate;
	this->piFreqs = new double[matrixSize];

	std::copy(modelDef.aaRates,modelDef.aaRates + this->matrixFullSize, this->qMatrix);
	std::copy(modelDef.aaFreqs,modelDef.aaFreqs + this->matrixSize, this->piFreqs);

	this->setDiagonalMeans();
	this->doEigenDecomposition();
	eigenDecomposed=true;
}


void AminoacidSubstitutionModel::calculatePt()
{
	//One parameter model
	this->time = parameters[0];
	if(this->eigenDecomposed == false)
	{
		this->setDiagonalMeans();
		this->doEigenDecomposition();
		eigenDecomposed=true;
	}
	calculateGammaPtMatrices();
}

AminoacidSubstitutionModel::~AminoacidSubstitutionModel()
{
}

void AminoacidSubstitutionModel::summarize()
{
	std::cout << "AA Substitution model summary " << std::endl << "Divergence time: " << this->time << std::endl;
}

} /* namespace EBC */


