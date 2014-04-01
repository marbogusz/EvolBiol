/*
 * ReducedPairHmmInsertState.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#ifndef REDUCEDPAIRHMMINSERTSTATE_HPP_
#define REDUCEDPAIRHMMINSERTSTATE_HPP_

#include "DpReducedMatrix.hpp"
#include "PairHmmStateBase.hpp"

namespace EBC
{

class ReducedPairHmmInsertState: public DpReducedMatrix, public PairHmmStateBase
{
public:
	ReducedPairHmmInsertState(unsigned int xSize, unsigned int ySize);

	virtual ~ReducedPairHmmInsertState();

	void initializeData();

	inline void nextRow()
	{
		double* tmp = previousRow;
		previousRow = currentRow;
		currentRow = tmp;
	}

};

} /* namespace EBC */

#endif /* REDUCEDPAIRHMMINSERTSTATE_HPP_ */
