//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

/*
 * PairwiseHmmInsertionState.hpp
 *
 *  Created on: Feb 26, 2014
 *      Author: root
 */

#ifndef PAIRWISEHMMINSERTIONSTATE_HPP_
#define PAIRWISEHMMINSERTIONSTATE_HPP_

#include "hmm/PairwiseHmmStateBase.hpp"

namespace EBC
{

class PairwiseHmmInsertState: public EBC::PairwiseHmmStateBase
{

public:

	void initializeData(double lnl, bool backwards=false);

	PairwiseHmmInsertState(unsigned int x, unsigned int y);

	PairwiseHmmInsertState(DpMatrixBase*);

	void setDirection(unsigned int i, unsigned int j);


	virtual ~PairwiseHmmInsertState();
};

} /* namespace EBC */
#endif /* PAIRHMMINSERTIONSTATE_HPP_ */
