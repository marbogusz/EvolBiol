/*
 * Definitions.hpp
 *
 *  Created on: Jan 14, 2014
 *      Author: mbogusz
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <iostream>

#define DEBUG_BUILD 1

#ifdef DEBUG_BUILD
#  define DEBUG(x) do { std::cerr << "DBG! " << x << std::endl; } while (0)
#  define DEBUGN(x) do { std::cerr << x ; } while (0)
#  define DEBUGV(x,n) do { for (int i=0;i<n;i++) std::cerr << x[i] << "\t"; std::cerr<<std::endl;} while (0)
#else
#  define DEBUG(x) do {} while (0)
#  define DEBUGN(x) do {} while (0)
#  define DEBUGV(x,n) do {} while (0)
#endif

namespace EBC
{

class Definitions
{
	enum SequenceType {Nucleotide, Aminoacid, Codon};
	//TODO - CODON lookup table


};

} /* namespace EBC */
#endif /* DEFINITIONS_H_ */
