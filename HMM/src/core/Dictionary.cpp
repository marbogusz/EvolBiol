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

void Dictionary::setAlphabet(const string dict[], unsigned short size)
{
	unsigned short i;

	for (i=0; i<=size;i++){
		alphabet.push_back(dict[i]);
	}

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
	for(auto sym : alphabet){
		cout << sym << ",";
	}
}

string& Dictionary::getSymbolAt(unsigned char i)
{
	return (alphabet[i]);
}

SequenceElement* Dictionary::getSequenceElement(string& symbol)
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
		string cstrg = string(1, *it);
		(*translatedVector)[pos] = getSequenceElement(cstrg);
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

const string Dictionary::nucleotides[] = {"T", "C", "A", "G", "-"};
//const char Dictionary::nucleotides[] = {'A', 'C', 'G', 'T'};
const string Dictionary::aminoacids[] = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","-"};

const string Dictionary::gapChar = "-";

unsigned short Dictionary::getAlphabetSize()
{
	return this->alphabetSize;
}


NucleotideDictionary::NucleotideDictionary()
{
	gapId = 4;
	this->setAlphabet(Dictionary::nucleotides,4);

}


AminoacidDictionary::AminoacidDictionary()
{
	gapId = 20;
	this->setAlphabet(Dictionary::aminoacids,20);

}



}

