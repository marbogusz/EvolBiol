/*
 * PMatrix.cpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#include <core/PMatrix.hpp>

namespace EBC
{

PMatrix::PMatrix(SubstitutionModelBase* m) : model(m),  matrixSize(m->getMatrixSize()) ,time(0),
		rateCategories(model->getRateCategories()), ptMatrices(rateCategories, nullptr)
{
	this->matrixFullSize = matrixSize*matrixSize;
}

PMatrix::~PMatrix()
{
	for (int i = 0; i < ptMatrices.size(); i++)
	{
		delete [] ptMatrices[i];
	}
}

void PMatrix::setTime(double t)
{
		this->time = t;
}

void PMatrix::summarize()
{
	cout << "P(t) matrix summary :" << endl;
	cout << "Divergence time : " << time << endl;

	cout << "Pairwise Site patterns " << endl;
}

} /* namespace EBC */


