//==============================================================================
// Pair-HMM phylogenetic tree estimator
// 
// Copyright (c) 2015 Marcin Bogusz.
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
//     This program is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses>.
//==============================================================================

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
		this->translator.insert(std::make_pair(tolower(alphabet[i]),new SequenceElement(i==gapId, i, NULL, alphabet[i])));
	}

	//alphabet size does not include gap e.g. size is 4 for nucleotides
	this->alphabetSize = size;
}

void Dictionary::outputAlphabet()
{
	cout << "Model dictionary: " << endl;

	cout << alphabet << endl;
}

char Dictionary::getSymbolAt(unsigned char index)
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

vector<SequenceElement*>* Dictionary::translate(string& sequence, bool removeGaps)
{

	vector<SequenceElement*> *translatedVector = new vector<SequenceElement*>();
	translatedVector->reserve(sequence.size());
	unsigned short currentEl;

//	DEBUG("Transled: ");
	unsigned int pos = 0;
	for(string::iterator it = sequence.begin(); it < sequence.end(); it++)
	{
		auto el = getSequenceElement(*it);

		if(el->isIsGap() && removeGaps)
			continue;
		translatedVector->push_back(el);
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

