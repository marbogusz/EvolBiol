/*
 * PairwiseHmmMatchState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRWISEHMMMATCHSTATE_HPP_
#define PAIRWISEHMMMATCHSTATE_HPP_

#include "PairwiseHmmStateBase.hpp"
#include "DpMatrixBase.hpp"

namespace EBC
{

class PairwiseHmmMatchState: public EBC::PairwiseHmmStateBase
{
protected:

	void initializeData();

public:
	PairwiseHmmMatchState(unsigned int x, unsigned int y);

	PairwiseHmmMatchState(DpMatrixBase*);

	void setDirection(unsigned int i, unsigned int j);

	virtual ~PairwiseHmmMatchState();
};

} /* namespace EBC */
#endif /* PAIRWISEHMMMATCHSTATE_HPP_ */
