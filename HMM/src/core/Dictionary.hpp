/*
 * Dictionary.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"

using namespace std;

namespace EBC
{
	class Dictionary
	{
	protected:
		unsigned short alphabetSize;
		
		unsigned char gapId;

		Definitions::FrequencyScheme fScheme;

		//equilibruim frequencies based on various strategies
		double* elementFrequencies;	

		vector<string> alphabet;
		
		map<string,SequenceElement*> translator;

		void setEqualFrequencies();

		unsigned int simpleCount;

	public:

		static const string nucleotides[5];
		static const string aminoacids[21];
		static const string codons[65];
		static const int geneticCode[64];
		static const string gapChar;

		virtual vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);

		virtual unsigned short getAlphabetSize();

		Dictionary(Definitions::FrequencyScheme fs);

		SequenceElement* getSequenceElement(string& symbol);

		virtual string& getSymbolAt(unsigned char i);

		virtual void outputAlphabet();

		virtual double* getElementFrequencies();

		inline unsigned char getGapID()
		{
			return gapId;
		}

	protected:
		virtual void setAlphabet(const string alphabet[], unsigned short size);

	};
}


#endif /* DICTIONARY_H_ */
