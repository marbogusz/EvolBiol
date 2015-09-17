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
 * HMMPathSample.h
 *
 *  Created on: Jan 26, 2015
 *      Author: marcin
 */

#ifndef SAMPLING_HMMPATHSAMPLE_HPP_
#define SAMPLING_HMMPATHSAMPLE_HPP_

#include "core/Definitions.hpp"

#include "hmm/PairwiseHmmStateBase.hpp"

namespace EBC {

class HMMPathSample {

protected:
	//state transition count matrix
	unsigned int stateTrans[Definitions::stateCount][Definitions::stateCount];
	//perfrmance friendly constant size array for both AA and nucleotides
	unsigned int sitePatterns[Definitions::aminoacidCount+1][Definitions::aminoacidCount+1];
	//ised to calculate lnl of the sequence
	PairwiseHmmStateBase* firstState;

public:
	HMMPathSample();

	inline void addTransition(unsigned char state1, unsigned char state2)
	{
		stateTrans[state1][state2] += 1;
	}

	inline void addSitePattern(unsigned char state1, unsigned char state2)
	{
		sitePatterns[state1][state2] += 1;
	}

	inline unsigned int getSitePattern(unsigned char state1, unsigned char state2)
	{
		return sitePatterns[state1][state2];
	}

	inline unsigned int getTransition(unsigned char state1, unsigned char state2)
	{
		return stateTrans[state1][state2];
	}

	inline PairwiseHmmStateBase* getFirstState()
	{
		return firstState;
	}

	void setLastState(PairwiseHmmStateBase* state)
	{
		this->firstState=state;
	}

	bool operator < (const HMMPathSample& pt) const
	{
		//does not matter
		return true;
	}

};

} /* namespace EBC */

#endif /* SAMPLING_HMMPATHSAMPLE_HPP_ */
