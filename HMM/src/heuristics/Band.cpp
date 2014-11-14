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

Band::~Band()
{
}

} /* namespace EBC */
