/*
 * PairHmmInsertionState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRHMMINSERTIONSTATE_HPP_
#define PAIRHMMINSERTIONSTATE_HPP_

#include "PairHmmState.hpp"

namespace EBC
{

class PairHmmInsertionState: public EBC::PairHmmState
{
protected:

	void initializeData();

	double gapOpen;
	double gapExtension;

public:
	PairHmmInsertionState(unsigned int x, unsigned int y, double gapOpening, double gapExtension);


	virtual ~PairHmmInsertionState();
};

} /* namespace EBC */
#endif /* PAIRHMMINSERTIONSTATE_HPP_ */
