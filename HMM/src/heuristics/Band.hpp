/*
 * Band.hpp
 *
 *  Created on: Oct 29, 2014
 *      Author: root
 */

#ifndef HEURISTICS_BAND_HPP_
#define HEURISTICS_BAND_HPP_

#include <vector>
#include "core/FileLogger.hpp"

using namespace std;

namespace EBC
{

class Band
{
protected:

	vector<pair<int, int> > matchBand;
	vector<pair<int, int> > insertBand;
	vector<pair<int, int> > deleteBand;


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

	inline void output()
	{
		for(int i =0; i< matchBand.size(); i++)
		{
			FileLogger::DebugLogger() << "M/X/Y bands " << i << "\t" << getMatchRangeAt(i).first <<"\t" << getMatchRangeAt(i).second
							<< "\t" << getInsertRangeAt(i).first <<"\t" << getInsertRangeAt(i).second
							<< "\t" << getDeleteRangeAt(i).first <<"\t" << getDeleteRangeAt(i).second << "\n";

		}
	}

	const vector<pair<int, int> >& getDeleteBand() const {
		return deleteBand;
	}

	const vector<pair<int, int> >& getInsertBand() const {
		return insertBand;
	}

	const vector<pair<int, int> >& getMatchBand() const {
		return matchBand;
	}
};

} /* namespace EBC */

#endif /* HEURISTICS_BAND_HPP_ */
