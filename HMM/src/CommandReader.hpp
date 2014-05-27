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
#include <dlib/cmd_line_parser.h>

namespace EBC
{

class CommandReader
{
private:

	int count;			//parameter count
	char** args;		//command line arguments

public:

	dlib::command_line_parser parser;

	CommandReader(int argc, char** argv);
	IParser* getParser() throw (ProgramException&);

	inline bool isViterbi()
	{
		return parser.option("V");
	}

	inline bool isForward()
	{
		return parser.option("F");
	}

	vector<double> getIndelParams();

	vector<double> getSubstParams();

	bool isOutputViterbiAlignment()
	{
		return parser.option("ov");
	}

	Definitions::ModelType getModelType()
	{
		if (parser.option("rev"))
		{
			return Definitions::ModelType::GTR;
		}
		if (parser.option("hky"))
		{
			return Definitions::ModelType::HKY85;
		}
	}

	string getInputFileName()
	{
		return parser.option("in").argument();
	}

	Definitions::OptimizationType getOptimizationType()
	{
			unsigned int opt = get_option(parser,"o",1);
			return (Definitions::OptimizationType) opt;
	}

	bool getBanding()
	{
		unsigned int opt = get_option(parser,"b",0);
		return (opt == 1);
	}

};

} /* namespace EBC */
#endif /* COMMANDREADER_H_ */
