/*
 * PairHmmDeletionState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRHMMDELETIONSTATE_HPP_
#define PAIRHMMDELETIONSTATE_HPP_

#include "PairHmmState.hpp"

namespace EBC
{

class PairHmmDeletionState: public EBC::PairHmmState
{
protected:

	void initializeData();

	double gapOpen;
	double gapExtension;

public:
	PairHmmDeletionState(unsigned int x, unsigned int y, double gapOpening, double gapExtension);


	virtual ~PairHmmDeletionState();
};

} /* namespace EBC */
#endif /* PAIRHMMDELETIONSTATE_HPP_ */
