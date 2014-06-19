/*
 * SubstitutionModelBase.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author:mbogusz
 */

#ifndef S_MODEL_BASE_H_
#define S_MODEL_BASE_H_

#include "Dictionary.hpp"
#include "Definitions.hpp"
#include "Maths.hpp"
#include "HmmException.hpp"
#include <cmath>
#include <vector>

namespace EBC
{

class SubstitutionModelBase
{

protected:

	vector<double> parameterHiBounds;
	vector<double> parameterLoBounds;

	unsigned int matrixSize;

	unsigned int matrixFullSize;

	//number of discrete gamma rate categories
	unsigned int rateCategories;

	//gamma distribution alpha parameter
	double alpha;

	double* gammaFrequencies;

	double* gammaRates;

	//equilibrium frequencies
	double* piFreqs;

	//q-rates
	double* qMatrix;

	double* vMatrix;

	double* uMatrix;

	//FIXME - remove and treat separately!
	double* pMatrix;

	//square roots matrix for eigen decomposition
	double* squareRoots;

	//roots vector for eigen decomposition
	double* roots;

	//ML parameters - model parameters + time!
	double* parameters;

	Dictionary* dictionary;

	Maths* maths;

	unsigned int paramsNumber;

	double meanRate;

	//current divergence time
	double time;

	//Allocate the memory;
	void allocateMatrices();

	void destroyMatrices();

	void doEigenDecomposition();

	void setDiagonalMeans();

	void calculateGamma();

	void calculateGammaPtMatrices();

public:

	SubstitutionModelBase(Dictionary*, Maths*i, unsigned int, unsigned int);

	virtual ~SubstitutionModelBase();

	virtual void summarize()=0;

	virtual void calculatePt()=0;

	void setObservedFrequencies(double* observedFrequencies);

	double getPXiYi(unsigned int xi, unsigned int yi);

	double getQXi(unsigned int xi);

	virtual void setParameters(vector<double>&)=0;

	void setTime(double t)
	{
		this->time = t;
	}
//getters and setters

	inline double getHiBound(unsigned int pos)
	{
		return parameterHiBounds[pos];
	}

	inline double getLoBound(unsigned int pos)
	{
		return parameterLoBounds[pos];
	}

	void setAlpha(double a)
	{
		this->alpha = a;
		calculateGamma();
	}

	double getAlpha()
	{
		return alpha;
	}


	void setParameters(const double* params)
	{
		for (unsigned int i=0;i<paramsNumber;i++)
			this->parameters[i] = params[i];
			//this->parameters[i] = fabs(params[i]);
	}


	unsigned int getParamsNumber() const 
	{
		return paramsNumber;
	}

	double* getParameters() 
	{
		return parameters;
	}

};

} /* namespace EBC */
#endif /* MODEL_H_ */