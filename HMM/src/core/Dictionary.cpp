/*
 * Dictionary.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#include "core/Dictionary.hpp"
#include "core/Definitions.hpp"
#include <algorithm>

namespace EBC
{

void Dictionary::setAlphabet(char dict[], unsigned short size)
{
	unsigned short i;
	for(i=0; i<size; i++)
	{
		this->alphabet.push_back(string(1,dict[i]));
		this->translator.insert(std::make_pair(string(1,dict[i]),i));
	}
	//add gap
	this->translator.insert(std::make_pair(string(1,'-'),size));
	this->alphabetSize = size;
}

void Dictionary::outputAlphabet()
{
	cout << "Model dictionary: " << endl;
	vector<string>::iterator it = alphabet.begin();
	while (it != alphabet.end())
	{
		cout << *it << "\t";
		it++;
	}
	cout << endl;
}

string Dictionary::getSymbolAt(short index)
{
	return alphabet[index];
}


short Dictionary::getSymbolIndex(string& symbol)
{
	return translator[symbol];
}

short Dictionary::getSymbolIndex(char symbol)
{
	string tmpSearchstring(1,symbol);
	return translator[tmpSearchstring];
}

vector<SequenceElement> Dictionary::translate(string& sequence, bool disregardIndels)
{
	//FIXME - deal with indels!
	vector<SequenceElement> translatedVector;
	translatedVector.reserve((sequence.size()));
	short currentEl;

	DEBUG("Transled: ");
	for(string::iterator it = sequence.begin(); it < sequence.end(); it++)
	{
		currentEl = getSymbolIndex(*it);
		if (currentEl == alphabetSize && disregardIndels)
			continue;
		DEBUGN(currentEl);
		//translatedVector.push_back(SequenceElement((currentEl == alphabetSize), currentEl,NULL, getSymbolAt(currentEl)));
		translatedVector.push_back(SequenceElement((currentEl == alphabetSize), currentEl,NULL, getSymbolAt(currentEl)));
	}
	DEBUGN(std::endl);
	return translatedVector;

}

const char Dictionary::nucleotides[] = {'T', 'C', 'A', 'G'};
//const char Dictionary::nucleotides[] = {'A', 'C', 'G', 'T'};
const char Dictionary::aminoacids[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}


NucleotideDictionary::NucleotideDictionary()
{
	this->setAlphabet((char*)Dictionary::nucleotides,4);
}


AminoacidDictionary::AminoacidDictionary()
{
	this->setAlphabet((char*)Dictionary::aminoacids,20);
}



}

