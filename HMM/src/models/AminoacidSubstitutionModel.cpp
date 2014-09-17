/*
 * AminoacidSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "models/AminoacidSubstitutionModel.hpp"

namespace EBC
{

AminoacidSubstitutionModel::AminoacidSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha, Definitions::aaModelDefinition modelDef) :
	SubstitutionModelBase(dict,alg,alpha, Definitions::AAParamCount), eigenDecomposed(false)
{
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
	if(this->eigenDecomposed == false)
	{
		this->setDiagonalMeans();
		this->doEigenDecomposition();
		eigenDecomposed=true;
	}
	//calculateGammaPtMatrices();
}

AminoacidSubstitutionModel::~AminoacidSubstitutionModel()
{
}

void AminoacidSubstitutionModel::summarize()
{
	std::cout << "AA Substitution model summary " << std::endl << "Divergence time: " << this->time << std::endl;
}

} /* namespace EBC */


