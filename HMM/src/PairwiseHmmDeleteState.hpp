/*
 * PairwiseHmmDeleteState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRWISEHMMDELETIONSTATE_HPP_
#define PAIRWISEHMMDELETIONSTATE_HPP_

#include "PairwiseHmmStateBase.hpp"

namespace EBC
{

class PairwiseHmmDeleteState: public EBC::PairwiseHmmStateBase
{
protected:

	void initializeData();

public:
	PairwiseHmmDeleteState(unsigned int x, unsigned int y);

	PairwiseHmmDeleteState(DpMatrixBase*);

	void setDirection(unsigned int i, unsigned int j);

	virtual ~PairwiseHmmDeleteState();
};

} /* namespace EBC */
#endif /* PAIRHMMDELETIONSTATE_HPP_ */
