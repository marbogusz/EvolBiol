/*
 * FileParser.hpp
 *
 *  Created on: Oct 1, 2013
 *      Author: mbogusz
 */

#ifndef FILEPARSER_H_
#define FILEPARSER_H_

#include <fstream>
#include "ParseException.hpp"
#include "IParser.hpp"
#include <vector>

using namespace std;

namespace EBC
{

class FileParser : public IParser
{
private:

	string filename;
	ifstream infile;
	bool endReached;
	vector<string> sequences;
	vector<string>::iterator it;

public:

	FileParser(char* filename);
	string getNextSequence();

	unsigned int getSequenceCount();

	string getSequenceAt(unsigned int);


	~FileParser();
};

} /* namespace EBC */
#endif /* FILEPARSER_H_ */
