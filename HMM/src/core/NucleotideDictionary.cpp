/*
 * Dictionary.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#include "core/NucleotideDictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{

NucleotideDictionary::NucleotideDictionary(Definitions::FrequencyScheme fs) : Dictionary(fs)
{
	gapId = 4;
	this->setAlphabet(Dictionary::nucleotides,4);
	this->translator.insert(std::make_pair("U",new SequenceElement(false, 0, NULL, "U")));

	elementFrequencies = new double[alphabetSize];

	for (unsigned int i=0; i<alphabetSize; i++)
				this->elementFrequencies[i] = 0;

	if (fs == Definitions::FrequencyScheme::Equal){
		setEqualFrequencies();
	}



}

}//Namespace definition
