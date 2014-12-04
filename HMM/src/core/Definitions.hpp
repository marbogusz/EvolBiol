/*
 * Definitions.hpp
 *
 *  Created on: Jan 14, 2014
 *      Author: mbogusz
 */

#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

#include <iostream>
#include "core/FileLogger.hpp"

//#define DEBUG_BUILD 1

#define DUMP(x) do { FileLogger::DumpLogger() << "   [DUMP]\t" << x << "\n"; } while (0)
//#define DUMP(x) do {} while (0)
#define DEBUG(x) do { FileLogger::DebugLogger() << "  [DEBUG]\t" << x << "\n"; } while (0)
#define INFO(x) do { FileLogger::InfoLogger() << " [INFO]\t"  << x<< "\n"; } while (0)
#define WARN(x) do { FileLogger::WarningLogger() << "! [WARN]\t"  << x << "\n"; } while (0)
#define ERROR(x) do { FileLogger::ErrorLogger() << "!!! [ERROR]\t"  << x << "\n"; } while (0)

#  define DEBUGN(x) do {} while (0)
#  define DEBUGV(x,n) do {} while (0)
#  define DEBUGM(x,n) do {} while (0)

/*
#ifdef DEBUG_BUILD
#  define DEBUG(x) do { std::cerr<< "DBG! " << x << std::endl; } while (0)
#  define DEBUGN(x) do { std::cerr << x ; } while (0)
#  define DEBUGV(x,n) do { for (int i=0;i<n;i++) std::cerr << x[i] << "\t"; std::cerr << std::endl;} while (0)
#  define DEBUGM(x,h,w) do { for (int i=0;i<h;i++){std::cerr << endl; for(int j=0;j<w;j++) std::cerr << x[i][j] << "\t";} std::cerr <<std::endl;} while (0)
#else
#  define DEBUG(x) do {} while (0)
#  define DEBUGN(x) do {} while (0)
#  define DEBUGV(x,n) do {} while (0)
#  define DEBUGM(x,n) do {} while (0)
#endif
*/
namespace EBC
{



class Definitions
{
public:

	constexpr static const double defaultGapPenalty = 0.5;

	constexpr static const unsigned int stateCount = 3;

	constexpr static const unsigned int modelEstimatorPathSamples = 5;

	constexpr static const unsigned int heuristicsTreeSize = 3;

	constexpr static const double divergenceBound = 3.5;
	constexpr static const double initialLambda = 0.02;
	constexpr static const double initialEpsilon = 0.5;

	constexpr static const unsigned int minBandDelta = 7; //This makes the min band width of 15 characters

	constexpr static const double minMatrixLikelihood = -1000000.0;



	struct aaModelDefinition
	{
		double maxRate;
		double aaRates[400];
		double aaFreqs[20];
	};

	enum SequenceType {Nucleotide, Aminoacid, Codon};
	//TODO - CODON lookup table
	enum ModelType {GTR, HKY85, LG};

	enum OptimizationType {BFGS, BOBYQA};

	enum AlgorithmType {Forward, Viterbi, MLE};

	enum DpMatrixType {Full, Limited};

	enum StateId {Match, Insert , Delete};

	static const unsigned int HKY85ParamCount;

	static const unsigned int GTRParamCount;

	static const unsigned int AAParamCount;

	static const unsigned int NBIndelParamCount;

	static const unsigned int nucleotideCount;

	static const unsigned int aminoacidCount;

	static const double blosum62[20][20];

	static const double blosum62gapOpening;

	static const double blosum62gapExtension;

	static aaModelDefinition aaLgModel;

};

} /* namespace EBC */
#endif /* DEFINITIONS_H_ */
