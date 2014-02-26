/*
 * Dictionary.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#include "Dictionary.hpp"
#include "Definitions.hpp"
#include <algorithm>

namespace EBC
{

void Dictionary::setAlphabet(char dict[], unsigned short size)
{
	unsigned short i;
	for(i=0; i<size; i++)
	{
		this->alphabet.push_back(string(1,dict[i]));
	}
	this->alphabetSize = size;
}

string Dictionary::getSymbolAt(short index)
{
	return alphabet[index];
}

const char Dictionary::nucleotides[] = {'T', 'C', 'A', 'G'};
const char Dictionary::aminoacids[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}


NucleotideDictionary::NucleotideDictionary()
{
	this->setAlphabet((char*)Dictionary::nucleotides,4);
}

NucleotideDictionary::~NucleotideDictionary()
{


}


//Fixme - implement something faster based on map collection
short NucleotideDictionary::getSymbolIndex(string symbol)
{
	return find(alphabet.begin(), alphabet.end(), symbol) - alphabet.begin();
}

short NucleotideDictionary::getSymbolIndex(char symbol)
{
	string tmpSearchstring(1,symbol);
	return find(alphabet.begin(), alphabet.end(), tmpSearchstring) - alphabet.begin();
}

vector<SequenceElement> NucleotideDictionary::translate(string& sequence, bool disregardIndels)
{
	//FIXME - deal with indels!
	vector<SequenceElement> translatedVector;
	short currentEl;

	DEBUG("Transled: ");
	for(string::iterator it = sequence.begin(); it < sequence.end(); it++)
	{
		currentEl = getSymbolIndex(*it);
		DEBUGN(currentEl);
		translatedVector.push_back(SequenceElement(false, currentEl,NULL));
	}
	DEBUGN(std::endl);
	return translatedVector;

}

}

