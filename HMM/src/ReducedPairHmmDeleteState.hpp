/*
 * ReducedPairHmmDeleteState.hpp
 *
 *  Created on: Mar 31, 2014
 *      Author: root
 */

#ifndef REDUCEDPAIRHMMDELETESTATE_HPP_
#define REDUCEDPAIRHMMDELETESTATE_HPP_

#include "DpReducedMatrix.hpp"
#include "PairHmmStateBase.hpp"

namespace EBC
{

class ReducedPairHmmDeleteState: public DpReducedMatrix, public PairHmmStateBase
{
public:
	ReducedPairHmmDeleteState(unsigned int xSize, unsigned int ySize);

	virtual ~ReducedPairHmmDeleteState();

	void initializeData();

	inline void nextRow()
	{
		//DEBUG("DC");
		//DEBUGV(currentRow, ySize);
		double* tmp = previousRow;
		previousRow = currentRow;
		currentRow = tmp;
	}

};

} /* namespace EBC */

#endif /* REDUCEDPAIRHMMMATCHSTATE_HPP_ */
