/*
 * NucleotideSubstitutionMode.cpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#include "NucleotideSubstitutionModel.hpp"

namespace EBC
{

NucleotideSubstitutionModel::NucleotideSubstitutionModel(Dictionary* dict, Maths* alg, unsigned int alpha) :
	SubstitutionModelBase(dict,alg,alpha)
{
}


void NucleotideSubstitutionModel::calculatePt()
{
	this->buildSmatrix();
	this->setDiagonalMeans();
	this->doEigenDecomposition();
	//qMatrix, roots, u and v matrices

	calculateGammaPtMatrices();
}

NucleotideSubstitutionModel::~NucleotideSubstitutionModel()
{
}

} /* namespace EBC */


