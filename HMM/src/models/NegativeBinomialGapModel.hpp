/*
 * NegativeBinomialGapModel.h
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef NBGAPMODEL_H_
#define NBGAPMODEL_H_

#include "models/IndelModel.hpp"


namespace EBC
{

class NegativeBinomialGapModel: public EBC::IndelModel
{
protected:
	double lambda;

public:

	NegativeBinomialGapModel();

	double calculateGapOpening(double time);
	double calculateGapExtension(double time);

	virtual ~NegativeBinomialGapModel();

	void calculateGeometricProbability(double lambda, double t);

	void setParameters(double* params);

	void setParameters(vector<double>);

	void calculate();

	void summarize();
};

} /* namespace EBC */
#endif /* AFFINEGEOMETRICGAPMODEL_H_ */
