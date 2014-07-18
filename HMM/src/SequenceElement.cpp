/*
 * SequenceElement.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: mbogusz
 */

#include "SequenceElement.hpp"

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

} /* namespace EBC */
