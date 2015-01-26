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
	unsigned int sitePatterns[Definitions::aminoacidCount][Definitions::aminoacidCount];
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
