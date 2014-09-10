/*
 * DistanceMatrix.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: root
 */

#include "core/DistanceMatrix.hpp"
#include "core/Definitions.hpp"
#include <random>

namespace EBC
{


EBC::DistanceMatrix::DistanceMatrix(int size) :taxas(size)
{

}

double EBC::DistanceMatrix::getDistance(unsigned int i, unsigned int j)
{
	double dst = 0;

	dst = this->distances[make_pair(i,j)];

	//DEBUG("Distance matrix getting distance");

	return dst;
}

void DistanceMatrix::buildMap()
{
}

void DistanceMatrix::addDistance(unsigned int s1, unsigned int s2,
		double distance)
{
	distances[make_pair(s1,s2)] = distance;
	distances[make_pair(s2,s1)] = distance;

	revdistances.insert(make_pair(distance,(make_pair(s1,s2))));
}

pair<unsigned int, unsigned int>& EBC::DistanceMatrix::getPairWithinDistance(
		double lo, double hi)
{
	default_random_engine generator;
	uniform_int_distribution<int> distribution(0,revdistances.size()-1);
	//iterate over keys
	auto itlow= revdistances.lower_bound(lo);
	auto ithi= revdistances.upper_bound(hi);
	auto end = revdistances.end();
	auto begin = revdistances.begin();

	if (itlow != end)
	{
		//we're in the bracket
		unsigned int dist = std::distance(itlow,ithi);
		uniform_int_distribution<int> dist2(0,dist);
		std::advance(itlow, dist2(generator));
		return (*itlow).second;

	}
	//return the lowest anyway ?
	else
	{
		std::advance(begin, distribution(generator) );
		return (*begin).second;
	}

}

unsigned int DistanceMatrix::getThirdLeafWithinDistance(double targetLen, unsigned int l1, unsigned int l2)
{
	unsigned int leaves = 0;
	double tempSum;
	map<double, unsigned int> mappings;
	for (leaves =0; leaves < this->taxas; leaves++)
	{
		tempSum = getDistance(l1,leaves) + getDistance(l2,leaves);
		mappings[tempSum] = leaves;
	}
	//range iterator over mappings!
	//find the closest to target len
	auto itlow= mappings.lower_bound(0.7*targetLen);
	auto ithi= mappings.upper_bound(1.3*targetLen);
	auto end = mappings.end();

	if (itlow == end)
		return (*(end--)).second;
	else if (ithi == end)
		return (*(mappings.begin())).second;
	else return (*itlow).second;
}

} /* namespace EBC */
