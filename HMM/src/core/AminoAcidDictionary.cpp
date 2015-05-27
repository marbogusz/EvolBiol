/*
 * Dictionary.cpp
 *
 *  Created on: Jan 14, 2014
 *      Author: root
 */

#include "core/AminoAcidDictionary.hpp"
#include "core/Definitions.hpp"
#include "core/HmmException.hpp"
#include <algorithm>

namespace EBC
{


AminoacidDictionary::AminoacidDictionary(Definitions::FrequencyScheme fs) : Dictionary(fs)
{
	gapId = 20;
	this->setAlphabet(Dictionary::aminoacids,20);

	elementFrequencies = new double[alphabetSize];

	for (unsigned int i=0; i<alphabetSize; i++)
			this->elementFrequencies[i] = 0;

	if (fs == Definitions::FrequencyScheme::Equal){
		setEqualFrequencies();
	}

}


}//Namespace definition
