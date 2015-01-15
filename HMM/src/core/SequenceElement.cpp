/*
 * SequenceElement.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: mbogusz
 */

#include "core/SequenceElement.hpp"

namespace EBC
{



SequenceElement::SequenceElement(bool isGap, short index, short* alternativeIndexes, string sy) : symbol(sy)
{
	this->isGap = isGap;
	//matrix index of -1 means gap or multiple options!
	this->matrixIndex = index;
	//FIXME - implememnt the following
	//vector<short> alternativeIndexes
}

SequenceElement::SequenceElement() : SequenceElement(false,0, NULL, ""){}

bool SequenceElement::operator== (SequenceElement &cP1, SequenceElement &cP2)
{
	return cP1.matrixIndex == cP2.matrixIndex;
}

} /* namespace EBC */
