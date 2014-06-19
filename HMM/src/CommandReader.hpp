/*
 * CommandReader.hpp
 *
 *  Created on: Oct 7, 2013
 *      Author: mbogusz
 */

#ifndef COMMANDREADER_H_
#define COMMANDREADER_H_

#include "IParser.hpp"
#include "HmmException.hpp"
#include <dlib/cmd_line_parser.h>

namespace EBC
{

class CommandReader
{
public:

	dlib::command_line_parser parser;

	CommandReader(int argc, char** argv);
	IParser* getParser() throw (HmmException&);

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
		if (parser.option("lg"))
		{
			return Definitions::ModelType::LG;
		}
	}

	string getInputFileName()
	{
		return parser.option("in").argument();
	}

	Definitions::OptimizationType getOptimizationType()
	{
			unsigned int opt = get_option(parser,"o",0);
			return (Definitions::OptimizationType) opt;
	}

	bool getBanding()
	{
		unsigned int opt = get_option(parser,"b",0);
		return (opt == 1);
	}

	unsigned int getBandFactor()
	{
		unsigned int opt = get_option(parser,"bf",5);
		return opt;
	}

	double getDistance()
	{
		double opt = get_option(parser,"d",-1.0);
		return opt;
	}

	double getAlpha()
	{
		return get_option(parser,"initAlpha",0.5);
	}

	double getCategories()
	{
		return get_option(parser,"rateCat",0);
	}

	bool estimateAlpha()
	{
		int res = get_option(parser,"estimateAlpha",0);
		return res == 1;
	}
	Definitions::SequenceType getSequenceType()
	{
		if (parser.option("rev") || parser.option("hky"))
			return Definitions::SequenceType::Nucleotide;
		else return Definitions::SequenceType::Aminoacid;
		//FIXME codons
	}
};

} /* namespace EBC */
#endif /* COMMANDREADER_H_ */
