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
		parser.add_option("o","Set optimizer, default is bfgs",1);
		parser.add_option("param_rev","GTR model parameters",5);
		parser.add_option("param_hky","HKY85 model parameters",1);
		parser.add_option("ov","Output viterbi alignment for estimated parameters");

		parser.parse(argc,argv);

		const char* one_time_opts[] = {"V", "F", "in", "m", "s", "i","d","h","b","o"};
		parser.check_one_time_options(one_time_opts);

		parser.check_incompatible_options("V", "F");
		parser.check_incompatible_options("rev", "hky");
		parser.check_incompatible_options("d", "F");

		const char* f_sub_opts[] = {"b","o","ov"};
		const char* rev_sub_opts[] = {"param_rev"};
		const char* hky_sub_opts[] = {"param_hky"};
		parser.check_sub_options("F", f_sub_opts);
		parser.check_sub_options("rev", rev_sub_opts);
		parser.check_sub_options("hky", hky_sub_opts);

		parser.check_option_arg_range("param_hky", 0.0000001, 20.0);
		parser.check_option_arg_range("param_rev", 0.0000001, 10.0);
		parser.check_option_arg_range("i", 0.0000001, 1.0);
		parser.check_option_arg_range("d", 0.0000001, 2.0);


		//if (parser.option("h"))
		//{
			// display all the command line options
		//    cout << "Usage: HMM (-F|-V) --in input_file model poarameters misc options\n";
		 //   parser.print_options();
		//}
	}
	catch (exception& e)
	{
	        throw ProgramException(e.what());
	}
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
