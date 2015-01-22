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
	//this->alphabet.reserve(size+1);

	//includes gap
	alphabet.append(dict,size+1);

	for(unsigned short i=0; i<=size; i++)
	{
		//this->alphabet.push_back(string(1,dict[i]));
		this->translator.insert(std::make_pair(alphabet[i],new SequenceElement(i==gapId, i, NULL, alphabet[i])));
	}

	//alphabet size does not include gap e.g. size is 4 for nucleotides
	this->alphabetSize = size;
}

void Dictionary::outputAlphabet()
{
	cout << "Model dictionary: " << endl;

	cout << alphabet << endl;
}

char Dictionary::getSymbolAt(unsigned short index)
{
	return alphabet[index];
}


unsigned char Dictionary::getSymbolIndex(char symbol)
{
	return (translator[symbol])->getMatrixIndex();
}

SequenceElement* Dictionary::getSequenceElement(char symbol)
{
	return translator[symbol];
}

vector<SequenceElement*>* Dictionary::translate(string& sequence, bool disregardIndels)
{

	vector<SequenceElement*> *translatedVector = new vector<SequenceElement*>(sequence.size());
	unsigned short currentEl;

//	DEBUG("Transled: ");
	unsigned int pos = 0;
	for(string::iterator it = sequence.begin(); it < sequence.end(); it++)
	{
		(*translatedVector)[pos] = getSequenceElement(*it);
		pos++;
		//currentEl = getSymbolIndex(*it);
		//if (currentEl == alphabetSize && disregardIndels)
		//	continue;
//		DEBUGN(currentEl);
		//translatedVector.push_back(SequenceElement((currentEl == alphabetSize), currentEl,NULL, getSymbolAt(currentEl)));
		//translatedVector.push_back(SequenceElement((currentEl == alphabetSize), currentEl,NULL, getSymbolAt(currentEl)));
	}
//	DEBUGN(std::endl);
	return translatedVector;

}

const char Dictionary::nucleotides[] = {'T', 'C', 'A', 'G', '-'};
//const char Dictionary::nucleotides[] = {'A', 'C', 'G', 'T'};
const char Dictionary::aminoacids[] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-'};

const char Dictionary::gapChar = '-';

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}


NucleotideDictionary::NucleotideDictionary()
{
	gapId = 4;
	this->setAlphabet((char*)Dictionary::nucleotides,4);

}


AminoacidDictionary::AminoacidDictionary()
{
	gapId = 20;
	this->setAlphabet((char*)Dictionary::aminoacids,20);

}



}

