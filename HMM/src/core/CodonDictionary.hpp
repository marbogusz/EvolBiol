/*
 * Dictionary.hpp
 *
 *  Created on: Jan 13, 2014
 *      Author: root
 */

#ifndef CODONDICTIONARY_H_
#define CODONDICTIONARY_H_

#include <vector>
#include <string>
#include <map>
#include "core/SequenceElement.hpp"
#include "core/Definitions.hpp"
#include "core/Dictionary.hpp"

using namespace std;

namespace EBC
{

	class CodonDictionary : public Dictionary
	{
	protected:
		int reducedGenCode[64];

		bool isPurine(char base);
		bool isPyramidine(char base);

		void calculateFreqeuencies();

	public:
		CodonDictionary(Definitions::FrequencyScheme fs);

		vector<SequenceElement*>* translate(string &sequence, bool disregardIndels = false);

		//void setGeneticCode(const int[]);

		int getAminoacidId(unsigned int codonId);

		unsigned int getNumberOfDifferentPositions(unsigned int codon1, unsigned int codon2,
				bool& synonymous, bool& transition);
	};
}


#endif /* DICTIONARY_H_ */
