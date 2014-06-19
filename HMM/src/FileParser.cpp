/*
 * FileParser.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: mbogusz
 */

#include "FileParser.hpp"
#include <iostream>
#include <algorithm>


namespace EBC
{

FileParser::FileParser(const char* filename)
{
	this->filename = filename;
	//this->infile.exceptions(ifstream::failbit | ifstream::badbit);
	this->infile.open(this->filename.c_str(), fstream::in);

	if(!this->infile)
	{
		throw HmmException("Can't open the file");
	}
	string tmp;
	string seq = "";
	bool justStarted = true;

	while (std::getline(infile, tmp,'\n'))
	{
		//wait for the first sequence description
		if(justStarted)
		{
			if (!isDefinitionLine(tmp))
				continue;
			else
			{
				justStarted = false;
			}
		}
		else
		{
			if (!isDefinitionLine(tmp))
			{
				seq += tmp;
			}
			else
			{
				trimWsChars(seq);
				sequences.push_back(seq);
				seq = "";
			}
		}
	}

	if(seq != "")
	{
		trimWsChars(seq);
		sequences.push_back(seq);
		seq = "";
	}
	infile.close();
	it=sequences.begin();
}

bool FileParser::isDefinitionLine(string& s)
{
	std::size_t found = s.find(">");
	return (found !=std::string::npos);
}

void FileParser::trimWsChars(string& s)
{
		s.erase(std::remove_if( s.begin(), s.end(), [](char c){ return (c =='\r' || c =='\t' || c == ' ' || c == '\n');}), s.end() );
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


