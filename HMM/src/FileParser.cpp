/*
 * FileParser.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: mbogusz
 */

#include "FileParser.hpp"
#include <iostream>


namespace EBC
{

FileParser::FileParser(const char* filename)
{
	this->filename = filename;
	//this->infile.exceptions(ifstream::failbit | ifstream::badbit);
	this->infile.open(this->filename.c_str(), fstream::in);

	if(!this->infile)
	{
		throw ProgramException("Can't open the file");
	}
	string tmp;

	while (std::getline(infile, tmp,'\n'))
	{
		//DEBUG("FILE SEQUENCE " << tmp );
		sequences.push_back(tmp);
	}
	infile.close();
	it=sequences.begin();
}

string FileParser::getNextSequence()
{
	//DEBUG("Get next sequence it: " << *it);
	return *it++;
}

unsigned int FileParser::getSequenceCount()
{
	return sequences.size();
}

string FileParser::getSequenceAt(unsigned int position)
{
	return sequences.at(position);
}

FileParser::~FileParser()
{
	if(this->infile.is_open())
	{
		this->infile.close();
	}
}

} /* namespace EBC */
