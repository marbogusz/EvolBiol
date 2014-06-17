/*
 * PairwiseHmmInsertionState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRWISEHMMINSERTIONSTATE_HPP_
#define PAIRWISEHMMINSERTIONSTATE_HPP_

#include "PairwiseHmmStateBase.hpp"

namespace EBC
{

class PairwiseHmmInsertState: public EBC::PairwiseHmmStateBase
{

public:

	void initializeData();

	PairwiseHmmInsertState(unsigned int x, unsigned int y);

	PairwiseHmmInsertState(DpMatrixBase*);

	void setDirection(unsigned int i, unsigned int j);


	virtual ~PairwiseHmmInsertState();
};

} /* namespace EBC */
#endif /* PAIRHMMINSERTIONSTATE_HPP_ */
