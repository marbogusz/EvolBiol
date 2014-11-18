/*
 * Band.cpp
 *
 *  Created on: Oct 29, 2014
 *      Author: root
 */

#include <heuristics/Band.hpp>

namespace EBC
{

Band::Band(unsigned int size) : matchBand(size), insertBand(size), deleteBand(size)
{
}

Band::Band(unsigned int len1, unsigned int len2, double coverage) :  matchBand(len1+1), insertBand(len1+1), deleteBand(len1+1)
{
	double ratio = static_cast<double>(len1+1)/(len2+1);
	double idealWidth = coverage * (len1+1);
	unsigned int calcHalfSize = static_cast<unsigned int>(idealWidth/2);
	unsigned int halfSize = calcHalfSize < Definitions::minBandDelta ? Definitions::minBandDelta : calcHalfSize;
	unsigned int estDiagonal;
	int min, max;

	//deal with the first column separately
	this->setMatchRangeAt(0,-1,-1);
	this->setInsertRangeAt(0,0,halfSize);
	this->setDeleteRangeAt(0,-1,-1);

	for (int col = 1; col <= len2; col++)
	{
		estDiagonal = static_cast<unsigned int>(col*ratio);

		min = (estDiagonal-halfSize) < 0 ? 0 : estDiagonal-halfSize;
		max = (estDiagonal+halfSize) > len1?  len1 : estDiagonal+halfSize;

		this->setMatchRangeAt(col,min,max);
		this->setInsertRangeAt(col,min,max);
		this->setDeleteRangeAt(col,min,max);
	}
}

Band::~Band()
{
}

} /* namespace EBC */


