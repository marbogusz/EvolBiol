/*
 * CommandReader.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: mbogusz
 */

#include "CommandReader.hpp"
#include "FileParser.hpp"
#include "TextInputParser.hpp"
#include <sstream>
#include <cstring>

using namespace std;

namespace EBC
{

CommandReader::CommandReader(int argc, char** argv) : count(argc), args(argv)
{

}

IParser* CommandReader::getParser() throw (ProgramException&)
{
	stringstream exSs;
	exSs << "Usage : " << this->args[0] << " -[f|c] [file name| sequence strings]";

	if (strcmp(this->args[1], "-f") == 0)
	{
		if (this->count != 3)
		{

			throw ProgramException(exSs.str());
		}
		else
		{
			return new FileParser(this->args[2]);
		}
	}
	/*else if(strcmp(this->args[1], "-c") == 0)
	{
		if (this->count != 4)
		{
			throw ParseException(exSs.str());
		}
		else
		{
			return new TextInputParser(this->args[2], this->args[3], true);
		}
	}
	*/
	else
	{
		throw ProgramException(exSs.str());
	}
}

} /* namespace EBC */
