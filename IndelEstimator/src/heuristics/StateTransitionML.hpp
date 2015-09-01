/*
 * StateTransitionMatrix.hpp
 *
 *  Created on: Oct 6, 2014
 *      Author: root
 */

#ifndef STATETRANSITIONMATRIX_HPP_
#define STATETRANSITIONMATRIX_HPP_

#include "models/IndelModel.hpp"
#include "core/Definitions.hpp"
#include "core/TransitionProbabilities.hpp"
#include "core/SequenceElement.hpp"

namespace EBC
{


//This is a 3x3 state transition matrix
class StateTransitionML
{
protected:

	TransitionProbabilities* tpb;
	double time;

	unsigned int matrixSize;

	//transition probabilities
	double g;
	double e;

	double md[Definitions::stateCount][Definitions::stateCount];
	//match,insert, delete
	unsigned int counts[Definitions::stateCount][Definitions::stateCount];
	double pis[Definitions::stateCount];

	unsigned char isGap;

	bool useStateEq;

	Definitions::StateId firstState;

public:
	virtual ~StateTransitionML();

	//get equilibrium frequencies
	void calculatePIs();

	void calculateParameters();

	StateTransitionML(IndelModel* im, double time, unsigned char, bool);

	void addSample(vector<unsigned char>*, vector<unsigned char>* s2);
	void addSample(vector<SequenceElement*>*, vector<SequenceElement*>* s2);

	double getLnL();
};

} /* namespace EBC */

#endif /* STATETRANSITIONMATRIX_HPP_ */
