/*
 * DistanceMatrix.cpp
 *
 *  Created on: Sep 2, 2014
 *      Author: root
 */

#include "DistanceMatrix.hpp"

namespace EBC
{


EBC::DistanceMatrix::DistanceMatrix(int size) :taxas(size), revdistances((size*(size-1))/2), distances(size*(size-1))
{

}

double EBC::DistanceMatrix::getDistance(unsigned int i, unsigned int j)
{
	double dst = 0;

	dst = this->distances[make_pair(i,j)];

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

	revdistances[distance] = make_pair(s1,s2);
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

	if (itlow != end)
	{
		//we're in the bracket
		unsigned int dist = std::distance(itlow,ithi);
		uniform_int_distribution<int> dist2(0,dist);
		return *(itlow+dist2(generator));

	}
	//return the lowest anyway ?
	else return *(end-distribution(generator));

}

unsigned int DistanceMatrix::getThirdLeafWithinDistance(unsigned int l1, unsigned int l2, double targetLen)
{
	unsigned int leaves = 0;
	double tempSum;
	map<double, unsigned int> mappings;
	for (leaves =0; leaves < this->taxas leaves++)
	{
		tempSum = getDistance(a,leaves) + getDistance(b,leaves);
		mappings[tempSum]
	}
	//range itetator!

}

} /* namespace EBC */
