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

void NucleotideSubstitutionModel::doEigenDecomposition()
{
	std::fill(roots, roots+matrixSize, 0);
	std::fill(uMatrix, uMatrix+matrixFullSize, 0);
	std::fill(vMatrix, vMatrix+matrixFullSize, 0);
	std::fill(squareRoots, squareRoots+matrixFullSize, 0);
	this->algebra->eigenQREV(qMatrix, piFreqs, matrixSize, roots, uMatrix, vMatrix, squareRoots);
}

void NucleotideSubstitutionModel::calculatePt()
{
	this->buildSmatrix();
	this->setDiagonalMeans();
	this->doEigenDecomposition();
	algebra->expLambdaT(roots, this->time, matrixSize);
	algebra->matrixByDiagonalMultiply(uMatrix, roots, matrixSize);
	this->pMatrix = this->algebra->matrixMultiply(uMatrix, vMatrix, matrixSize);
}

NucleotideSubstitutionModel::~NucleotideSubstitutionModel()
{
}

} /* namespace EBC */


