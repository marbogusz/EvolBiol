/*
 * ReducedPairHmmMatchState.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#ifndef REDUCEDPAIRHMMMATCHSTATE_HPP_
#define REDUCEDPAIRHMMMATCHSTATE_HPP_

#include "DpReducedMatrix.hpp"
#include "PairHmmStateBase.hpp"

namespace EBC
{

class ReducedPairHmmMatchState: public DpReducedMatrix, public PairHmmStateBase
{
public:
	ReducedPairHmmMatchState(unsigned int xSize, unsigned int ySize);

	virtual ~ReducedPairHmmMatchState();

	void initializeData();

	inline void nextRow()
	{
		//DEBUG("MC");
		//DEBUGV(currentRow, ySize);
		double* tmp = previousRow;
		previousRow = currentRow;
		currentRow = tmp;
		currentRow[0] = minVal;
	}

};

} /* namespace EBC */

#endif /* REDUCEDPAIRHMMMATCHSTATE_HPP_ */
