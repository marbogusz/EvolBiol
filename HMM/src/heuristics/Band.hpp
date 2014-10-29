/*
 * Band.hpp
 *
 *  Created on: Oct 29, 2014
 *      Author: root
 */

#ifndef HEURISTICS_BAND_HPP_
#define HEURISTICS_BAND_HPP_

#include <vector>
#include "hmm/EvolutionaryPairHMM.hpp"

using namespace std;

namespace EBC
{

class Band
{
protected:

	vector<pair<unsigned int, unsigned int>> matchBand;
	vector<pair<unsigned int, unsigned int>> insertBand;
	vector<pair<unsigned int, unsigned int>> deleteBand;

	double posteriorLikelihoodLimit;
	double posteriorLikelihoodDelta;

	//Match state
	PairwiseHmmStateBase* M;
	//Insert state
	PairwiseHmmStateBase* X;
	//Delete state
	PairwiseHmmStateBase* Y;


public:
	Band(unsigned int size);
	virtual ~Band();

	inline void setMatchRangeAt(unsigned int pos, unsigned int start, unsigned int end)
	{
		matchBand[pos] = std::make_pair(start,end);
	}
	inline void setInsertRangeAt(unsigned int pos, unsigned int start, unsigned int end)
	{
		insertBand[pos] = std::make_pair(start,end);
	}
	inline void setDeleteRangeAt(unsigned int pos, unsigned int start, unsigned int end)
	{
		deleteBand[pos] = std::make_pair(start,end);
	}

	inline std::pair<unsigned int, unsigned int> getMatchRangeAt(unsigned int pos)
	{
		return matchBand[pos];
	}
	inline std::pair<unsigned int, unsigned int> getInsertRangeAt(unsigned int pos)
	{
		return insertBand[pos];
	}
	inline std::pair<unsigned int, unsigned int> getDeleteRangeAt(unsigned int pos)
	{
		return deleteBand[pos];
	}

	void processPosteriorProbabilities(EvolutionaryPairHMM* hmm);
};

} /* namespace EBC */

#endif /* HEURISTICS_BAND_HPP_ */
