/*
 * PMatrix.hpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#ifndef PMATRIX_HPP_
#define PMATRIX_HPP_

#include "models/SubstitutionModelBase.hpp"
#include "core/ProgramException.hpp"
#include <vector>
#include <array>

using namespace std;

namespace EBC
{

class PMatrix
{
protected:

	SubstitutionModelBase* model;

	unsigned int matrixSize;

	double time;

	unsigned int rateCategories;

	vector<double*> ptMatrices;

	unsigned int matrixFullSize;

public:
	PMatrix(SubstitutionModelBase* m);
	virtual ~PMatrix();

	void setTime(double t);

	virtual void calculate()=0;

	inline double getEquilibriumFreq(unsigned int xi)
	{
		return model->getEquilibriumFrequencies(xi);
	}

	inline double getLogEquilibriumFreq(unsigned int xi)
	{
		return model->getLogEquilibriumFrequencies(xi);
	}

	void summarize();

};

} /* namespace EBC */

#endif /* PMATRIX_HPP_ */
