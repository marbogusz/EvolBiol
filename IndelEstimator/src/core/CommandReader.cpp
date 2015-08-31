/*
 * CommandReader.cpp
 *
 *  Created on: Oct 7, 2013
 *      Author: Marcin Bogusz
 */

#include "core/CommandReader.hpp"
#include "core/FileParser.hpp"
#include "core/TextInputParser.hpp"
#include <sstream>
#include <cstring>
#include <cstdlib>

using namespace std;

namespace EBC
{

CommandReader::CommandReader(int argc, char** argv)
{
	try
	{
		parser.add_option("A", "Run indel analysis using aligned data");
		parser.add_option("R", "Run indel analysis using raw sequence data");
		parser.add_option("in","This option takes one argument which specifies the name of the file we want to analyze",1);
		parser.add_option("tree","This option takes one argument which specifies the name of the tree file (newick)",1);
		parser.add_option("rev", "REV Substitution Model");
		parser.add_option("hky", "HKY85 Substitution Model");
		parser.add_option("lg", "Le & Gasquel AA Substitution Model");
		parser.set_group_name("Miscellaneous Options");
		parser.add_option("h","Display this help message.");

		parser.add_option("rateCat", "Specify gamma rate categories, default is 5",1);
		parser.add_option("initAlpha", "Specify initial alpha parameter, default is 0.5",1 );
		parser.add_option("estimateAlpha", "Specify to estimate alpha 0|1, default is 1",1 );

		parser.add_option("lE", "log error");
		parser.add_option("lW", "log warning");
		parser.add_option("lI", "log info");
		parser.add_option("lD", "log debug");
		parser.add_option("lDD", "log dump");

		parser.parse(argc,argv);

		const char* one_time_opts[] = {"A", "R", "in", "tree"};
		parser.check_one_time_options(one_time_opts);

		parser.check_incompatible_options("A", "R");
		parser.check_incompatible_options("rev", "hky");

		const char* r_sub_opts[] = {"rev","hky","lg", "rateCat", "initAlpha", "estimateAlpha"};
		parser.check_sub_options("R", r_sub_opts);

		parser.check_option_arg_range("initAlpha", 0.0000001, 1000.0);

		if (!(parser.option("A") || parser.option("R")) && !parser.option("in") && !parser.option("tree") )

		{
		    cout << "Usage: indestimate -{A|R} --in sequence_file --tree tree_file\n";

		    parser.print_options();
			throw ProgramException("Specify which either A or R option along with sequence and tree files\n");
		}

		if (parser.option("h"))
		{
			// display all the command line options
		    cout << "Usage: indestimate -{A|R} --in sequence_file --tree tree_file\n";
		    parser.print_options();
		}
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
