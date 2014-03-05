/*
 * AffineGeometricGapModel.h
 *
 *  Created on: Feb 25, 2014
 *      Author: root
 */

#ifndef AFFINEGEOMETRICGAPMODEL_H_
#define AFFINEGEOMETRICGAPMODEL_H_

#include "IndelModel.hpp"

namespace EBC
{

class AffineGeometricGapModel: public EBC::IndelModel
{
protected:
	double lambda;
	double time;

public:

	AffineGeometricGapModel(double lambda, double t, double gapExtenstionProbs);

	AffineGeometricGapModel();

	virtual ~AffineGeometricGapModel();

	void calculateGeometricProbability(double lambda, double t);

	void setParameters(double* params);

	void summarize();
};

} /* namespace EBC */
#endif /* AFFINEGEOMETRICGAPMODEL_H_ */
