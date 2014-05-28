/*
 * CommandReader.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: Marcin Bogusz
 */

#include "CommandReader.hpp"
#include "FileParser.hpp"
#include "TextInputParser.hpp"
#include <sstream>
#include <cstring>
#include <cstdlib>

using namespace std;

namespace EBC
{

CommandReader::CommandReader(int argc, char** argv) : count(argc), args(argv)
{
	try
	{
		parser.add_option("V", "Run Viterbi algorithm using user parameters");
		parser.add_option("F", "Run Forward algorithm");
		parser.add_option("in","This option takes one argument which specifies the name of the file we want to analyze",1);
		parser.add_option("rev", "REV Substitution Model");
		parser.add_option("hky", "HKY85 Substitution Model");
		parser.add_option("i","indel parameters (NB probability and rate)",2);
		parser.add_option("d","evolutionary distance",1);
		parser.set_group_name("Miscellaneous Options");
		parser.add_option("h","Display this help message.");
		parser.add_option("b","Set banding 0/1, default is 0",1);
		parser.add_option("bf","Set band factor, default is 5%",1);
		parser.add_option("o","Set optimizer, 0- BFGS, 1- BOBYQA default is 1",1);
		parser.add_option("param_rev","GTR model parameters",5);
		parser.add_option("param_hky","HKY85 model parameters",1);
		parser.add_option("ov","Output viterbi alignment for estimated parameters");

		parser.parse(argc,argv);

		const char* one_time_opts[] = {"V", "F", "in", "i","d" ,"h","b","o", "ov"};
		parser.check_one_time_options(one_time_opts);

		parser.check_incompatible_options("V", "F");
		parser.check_incompatible_options("rev", "hky");
		//parser.check_incompatible_options("d", "F");

		const char* f_sub_opts[] = {"b","o","ov"};
		const char* rev_sub_opts[] = {"param_rev"};
		const char* hky_sub_opts[] = {"param_hky"};
		const char* band_sub_opts[] = {"bf"};
		parser.check_sub_options("F", f_sub_opts);
		parser.check_sub_options("rev", rev_sub_opts);
		parser.check_sub_options("hky", hky_sub_opts);
		parser.check_sub_options("b", band_sub_opts);

		parser.check_option_arg_range("param_hky", 0.0000001, 20.0);
		parser.check_option_arg_range("param_rev", 0.0000001, 10.0);
		parser.check_option_arg_range("i", 0.0000001, 1.0);
		parser.check_option_arg_range("d", 0.0000001, 3.0);
		parser.check_option_arg_range("bf", 1, 100);

		if (!parser.option("V") && !parser.option("F"))
		{
		    cout << "Usage: HMM (-F|-V) --in input_file (rev|hky) [param_rev .... | param_hky ...] [i indel parameters] [d distance] [b] [o=0|1] [ov]\n";
		    parser.print_options();
			throw ProgramException("Specify which algorithm you want to run!\n");
		}

		parser.check_option_arg_range("o", 0, 1);
		parser.check_option_arg_range("b", 0, 1);
		parser.check_option_arg_range("bf", 0, 100);
		if (parser.option("h"))
		{
			// display all the command line options
		    cout << "Usage: HMM (-F|-V) --in input_file (rev|hky) [param_rev .... | param_hky ...] [i indel parameters] [d distance] [b] [o=0|1] [ov]\n";
		    parser.print_options();
		}
	}
	catch (exception& e)
	{
	        throw ProgramException(e.what());
	}
}

vector<double> CommandReader::getSubstParams()
{
	int i;
	vector<double> vec;
	if (parser.option("hky"))
	{
		if(parser.option("param_hky"))
		{
			for (i=0; i< 1; i++)
			{
				DEBUG("hky parameter " << i <<  ": " << parser.option("param_hky").argument(i));
				vec.push_back(atof(parser.option("param_hky").argument(i).c_str()));
			}
		}
	}
	else if (parser.option("rev"))
	{
		if (parser.option("param_rev"))
		{
			for (i=0; i< 5; i++)
			{
				DEBUG("Rev parameter " << i <<  ": " << parser.option("param_rev").argument(i));
				vec.push_back(atof(parser.option("param_rev").argument(i).c_str()));
			}
		}
	}
	else throw ProgramException("Model not specified");

	return vec;
}

vector<double> CommandReader::getIndelParams()
{
	int i;
	vector<double> vec;
	if (parser.option("i"))
	{
		for (i=0; i< 2; i++)
		{
			DEBUG("indel parameter " << i << ": " << parser.option("i").argument(i));
			vec.push_back(atof(parser.option("i").argument(i).c_str()));
		}
	}
	return vec;
}

IParser* CommandReader::getParser() throw (ProgramException&)
{
	if (parser.option("in"))
	{
		return new FileParser((string(parser.option("in").argument())).c_str());
	}
	else
		throw ProgramException("input file not specified");
}

} /* namespace EBC */
