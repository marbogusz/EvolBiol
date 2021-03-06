/*
 * FileParser.cpp
 *
 *  Created on: Oct 1, 2013
 *      Author: mbogusz
 */

#include "core/FileParser.hpp"
#include <iostream>
#include <algorithm>


namespace EBC
{

FileParser::FileParser(const char* filename)
{
	//FIXME - check the file structure to avoid errors from phylip, nexus and other non-fasta files

	sequences = new vector<string>();
	names = new vector<string>();

	this->filename = filename;
	//this->infile.exceptions(ifstream::failbit | ifstream::badbit);
	this->infile.open(this->filename.c_str(), fstream::in);

	if(!this->infile)
	{
		throw ProgramException(string("Can't open the file : ") + string(filename) + string("\n"));
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
				names->push_back(getSequenceName(tmp));
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
				names->push_back(getSequenceName(tmp));
				trimWsChars(seq);
				sequences->push_back(seq);
				seq = "";
			}
		}
	}

	if(seq != "")
	{
		trimWsChars(seq);
		sequences->push_back(seq);
		seq = "";
	}
	infile.close();
	it=sequences->begin();
	itN=names->begin();
}

bool FileParser::isDefinitionLine(string& s)
{
	std::size_t found = s.find(">");
	return (found !=std::string::npos);
}

string FileParser::getSequenceName(string& s)
{
	//string name(s);
	s.erase(std::remove_if( s.begin(), s.end(), [](char c){ return (c =='>' || c =='\r' || c =='\t' || c == ' ' || c == '\n');}), s.end() );
	DEBUG("Found sequence named " << s);
	return s;//name;
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

string FileParser::getNextName()
{
	//DEBUG("Get next sequence it: " << *it);
	return *itN++;
}

unsigned int FileParser::getSequenceCount()
{
	return sequences->size();
}

string FileParser::getSequenceAt(unsigned int position)
{
	return sequences->at(position);
}

string FileParser::getSequenceNameAt(unsigned int position)
{
	return names->at(position);
}

FileParser::~FileParser()
{
	if(this->infile.is_open())
	{
		this->infile.close();
	}
}

} /* namespace EBC */


