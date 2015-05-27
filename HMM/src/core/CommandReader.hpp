/*
 * CommandReader.hpp
 *
 *  Created on: Oct 7, 2013
 *      Author: mbogusz
 */

#ifndef COMMANDREADER_H_
#define COMMANDREADER_H_

#include "core/IParser.hpp"
#include "core/HmmException.hpp"
#include "core/Definitions.hpp"
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

	inline bool isMLE()
	{
		return parser.option("M");
	}

	//2 below are optional and new
	inline bool isFdist()
	{
		return parser.option("X");
	}
	inline bool isVdist()
	{
		return parser.option("Y");
	}

	inline bool isFixedAlignment()
	{
		return parser.option("fa");
	}

	vector<double> getIndelParams();

	vector<double> getSubstParams();

	bool isOutputViterbiAlignment()
	{
		return parser.option("ov");
	}

	Definitions::AlgorithmType getAlgorithmType()
		{
			if (parser.option("V"))
			{
				return Definitions::AlgorithmType::Viterbi;
			}
			if (parser.option("F"))
			{
				return Definitions::AlgorithmType::Forward;
			}
			//default;
			return Definitions::AlgorithmType::MLE;
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
		if (parser.option("m0"))
		{
			return Definitions::ModelType::M0;
		}
		//default;
		return Definitions::ModelType::GTR;
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


	FileLogger::logType getLoggingLevel()
	{
		if (parser.option("lDD"))
			return FileLogger::L_DMP;
		if (parser.option("lE"))
			return FileLogger::L_ERR;
		if (parser.option("lW"))
			return FileLogger::L_WARN;
		if (parser.option("lI"))
			return FileLogger::L_INF;
		if (parser.option("lD"))
			return FileLogger::L_DBG;
		//info by default!
		return
			FileLogger::L_INF;
	}

	bool getBanding()
	{
		if (parser.option("b"))
			return true;
		else
			return false;
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
		return get_option(parser,"rateCat",1);
	}

	bool estimateAlpha()
	{
		int res = get_option(parser,"estimateAlpha",1);
		return res == 1;
	}
	Definitions::SequenceType getSequenceType()
	{
		if (parser.option("rev") || parser.option("hky"))
			return Definitions::SequenceType::Nucleotide;
		else if (parser.option("lg"))
			return Definitions::SequenceType::Aminoacid;
		else if (parser.option("m0"))
			return Definitions::SequenceType::Codon;
	}
};

} /* namespace EBC */
#endif /* COMMANDREADER_H_ */
