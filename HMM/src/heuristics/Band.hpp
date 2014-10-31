/*
 * Band.hpp
 *
 *  Created on: Oct 29, 2014
 *      Author: root
 */

#ifndef HEURISTICS_BAND_HPP_
#define HEURISTICS_BAND_HPP_

#include <vector>

using namespace std;

namespace EBC
{

class Band
{
protected:

	vector<pair<int, int>> matchBand;
	vector<pair<int, int>> insertBand;
	vector<pair<int, int>> deleteBand;


public:
	Band(unsigned int size);
	virtual ~Band();

	inline void setMatchRangeAt(unsigned int pos, int start, int end)
	{
		matchBand[pos] = std::make_pair(start,end);
	}
	inline void setInsertRangeAt(unsigned int pos, int start, int end)
	{
		insertBand[pos] = std::make_pair(start,end);
	}
	inline void setDeleteRangeAt(unsigned int pos, int start, int end)
	{
		deleteBand[pos] = std::make_pair(start,end);
	}

	inline std::pair<int, int> getMatchRangeAt(unsigned int pos)
	{
		return matchBand[pos];
	}
	inline std::pair<int, int> getInsertRangeAt(unsigned int pos)
	{
		return insertBand[pos];
	}
	inline std::pair<int, int> getDeleteRangeAt(unsigned int pos)
	{
		return deleteBand[pos];
	}
};

} /* namespace EBC */

#endif /* HEURISTICS_BAND_HPP_ */
