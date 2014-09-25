/*
 * DistanceMatrix.hpp
 *
 *  Created on: Sep 2, 2014
 *      Author: root
 */

#ifndef DISTANCEMATRIX_HPP_
#define DISTANCEMATRIX_HPP_

#include <map>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

namespace EBC
{

class DistanceMatrix
{
private:

	//map the length to the pair of sequences
	multimap <double, pair<unsigned int, unsigned int> > revdistances;

	//pairwise distances ordered, double the size of dictionary!
	map<pair<unsigned int, unsigned int>,double> distances;

	unsigned int taxas;

	void buildMap();


public:
	DistanceMatrix(int size);

	void addDistance(unsigned int s1, unsigned int s2, double distance);

	double getDistance(unsigned int s1, unsigned int s2);

	pair<unsigned int, unsigned int>& getPairWithinDistance(double lo, double hi);

	unsigned int getThirdLeafWithinDistance(double targetLen, unsigned int l1, unsigned int l2);


};

} /* namespace EBC */

#endif /* DISTANCEMATRIX_HPP_ */
