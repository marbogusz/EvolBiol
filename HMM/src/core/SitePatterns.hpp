/*
 * SitePatterns.hpp
 *
 *  Created on: Sep 17, 2014
 *      Author: root
 */

#ifndef SITEPATTERNS_HPP_
#define SITEPATTERNS_HPP_

#include <map>
#include <array>

using namespace std;

namespace EBC
{

template <std::size_t N>
class SitePatterns
{
protected:
	map<array<unsigned int,N>, unsigned int> patterns;

public:
	inline void addPattern(array<unsigned int,N> pat)
	{
			patterns[pat]++;
	}

	inline double getTotalLikelihood()
	{
		double result = 0;
		for (auto pt : patterns)
		{
			result += pt.second.first * pt.second.second;
		}

		return result;
	}

	SitePatterns();
	virtual ~SitePatterns();
};


} /* namespace EBC */

#endif /* SITEPATTERNS_HPP_ */
