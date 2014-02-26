/*
 * CommandReader.hpp
 *
 *  Created on: Oct 7, 2013
 *      Author: mbogusz
 */

#ifndef COMMANDREADER_H_
#define COMMANDREADER_H_

#include "IParser.hpp"
#include "ParseException.hpp"

namespace EBC
{

class CommandReader
{
private:

	int count;			//parameter count
	char** args;		//command line arguments

public:
	CommandReader(int argc, char** argv);
	IParser* getParser() throw (ProgramException&);
};

} /* namespace EBC */
#endif /* COMMANDREADER_H_ */
