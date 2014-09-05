/*
 * PairwiseHmmMatchState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRWISEHMMMATCHSTATE_HPP_
#define PAIRWISEHMMMATCHSTATE_HPP_

#include "hmm/PairwiseHmmStateBase.hpp"
#include "hmm/DpMatrixBase.hpp"

namespace EBC
{

class PairwiseHmmMatchState: public EBC::PairwiseHmmStateBase
{

public:
	void initializeData();

	PairwiseHmmMatchState(unsigned int x, unsigned int y);

	PairwiseHmmMatchState(DpMatrixBase*);

	void setDirection(unsigned int i, unsigned int j);

	virtual ~PairwiseHmmMatchState();
};

} /* namespace EBC */
#endif /* PAIRWISEHMMMATCHSTATE_HPP_ */
