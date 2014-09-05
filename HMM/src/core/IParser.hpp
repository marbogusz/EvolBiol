/*
 * IParser.hpp
 *
 *  Created on: Oct 7, 2013
 *      Author: mbogusz
 */

#ifndef IPARSER_H_
#define IPARSER_H_

#include "core/Definitions.hpp"
#include <string>

using namespace std;

namespace EBC
{

class IParser
{
public:
	virtual string getNextSequence() = 0;
	virtual unsigned int getSequenceCount() = 0;
	virtual string getSequenceAt(unsigned int) = 0;
};

} /* namespace EBC */
#endif /* IPARSER_H_ */
