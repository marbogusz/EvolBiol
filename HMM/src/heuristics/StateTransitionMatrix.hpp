/*
 * StateTransitionMatrix.hpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef STATETRANSITIONMATRIX_HPP_
#define STATETRANSITIONMATRIX_HPP_

#include "models/IndelModel.hpp"

namespace EBC
{


//This is a 3x3 state transition matrix
class StateTransitionMatrix
{
protected:

	IndelModel* indelModel;

	unsigned int matrixSize;

	//transition probabilities
	double d;
	double e;

	double matrixData[9];
	double pis[3];

public:
	virtual ~StateTransitionMatrix();

	//get equilibrium frequencies
	void calculatePIs();

	void calculateParameters();

	StateTransitionMatrix();

	double getTransition(unsigned int i, unsigned int j);

	double getPi(unsigned int);
};

} /* namespace EBC */

#endif /* STATETRANSITIONMATRIX_HPP_ */
