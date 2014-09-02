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

namespace EBC
{

class DistanceMatrix
{
private:

	map<string, unsigned int> dictionary;

	double** data;

public:
	DistanceMatrix(int size);
};

} /* namespace EBC */

#endif /* DISTANCEMATRIX_HPP_ */
