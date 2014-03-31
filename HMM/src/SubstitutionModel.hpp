/*
 * SubstitutionModel.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef MODEL_H_
#define MODEL_H_

#include "Dictionary.hpp"
#include "Definitions.hpp"
#include "Maths.hpp"
#include <cmath>

namespace EBC
{

class SubstitutionModel
{

protected:

	int matrixSize;

	int matrixFullSize;

	//equilibrium frequencies
	double* piFreqs;
	//FIXME -  currently observed freqs are equal to calculated ones
	double* observedFrequencies;

	//q-rates
	double* qMatrix;

	double* vMatrix;

	double* uMatrix;

	double* pMatrix;

	//square roots matrix for eigen decomposition
	double* squareRoots;

	//roots vector for eigen decomposition
	double* roots;

	//ML parameters
	double* parameters;

	Dictionary* dictionary;

	Maths* algebra;

	unsigned int paramsNumber;

	double meanRate;

	double time;

	void outputPtPi();

private:
	//Allocate the memory;
	void allocateMatrices();

	void destroyMatrices();

public:

	//logarithmic mode - log frequencies, distances
	//and parameters!
	bool logMode;

	SubstitutionModel(Dictionary*, Maths*);

	~SubstitutionModel();

	virtual void summarize()=0;

	void setObservedFrequencies(double* observedFrequencies);

	void setDiagMeans();

	void doEigenDecomposition();

	//returns a P(t) matrix
	double* calculatePt(double time);

	void calculatePt();

	virtual void setParametersInMatrix()=0;

	double getPXiYi(unsigned int xi, unsigned int yi);

	double getQXi(unsigned int xi);

//getters and setters

	void setParameters(const double* params)
	{
		for (int i=0;i<paramsNumber;i++)
			this->parameters[i] = logMode?exp(params[i]):params[i];
			//this->parameters[i] = fabs(logMode?log(params[i]):params[i]);
	}


	unsigned int getParamsNumber() const {
		return paramsNumber;
	}

	double* getParameters() {
		return parameters;
	}

};

} /* namespace EBC */
#endif /* MODEL_H_ */
