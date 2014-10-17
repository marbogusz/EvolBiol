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
#  define DEBUG(x) do { std::cout << "DBG! " << x << std::endl; } while (0)
#  define DEBUGN(x) do { std::cout << x ; } while (0)
#  define DEBUGV(x,n) do { for (int i=0;i<n;i++) std::cout << x[i] << "\t"; std::cout << std::endl;} while (0)
#  define DEBUGM(x,h,w) do { for (int i=0;i<h;i++){std::cout << endl; for(int j=0;j<w;j++) std::cout << x[i][j] << "\t";} std::cout <<std::endl;} while (0)
#else
#  define DEBUG(x) do {} while (0)
#  define DEBUGN(x) do {} while (0)
#  define DEBUGV(x,n) do {} while (0)
#  define DEBUGM(x,n) do {} while (0)
#endif

namespace EBC
{



class Definitions
{
public:

	constexpr static const double defaultGapPenalty = 0.5;

	constexpr static const unsigned int stateCount = 3;
	constexpr static const unsigned int heuristicsTreeSize = 3;
	constexpr static const double divergenceBound = 3.5;

	constexpr static const double initialLambda = 0.02;
	constexpr static const double initialEpsilon = 0.5;

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
	/*constexpr static const double LGaaModelPis[] = {0.079066,  0.055941,  0.041977,  0.053052,0.012937, 0.040767,  0.071586,  0.057337, 0.022355,
			0.062157,  0.099081, 0.064600, 0.022951, 0.042302, 0.044040, 0.061197, 0.053287, 0.012066, 0.034155, 0.069147};

	constexpr static const double LGaaMaxRate = 10.649107;


	        aaWAG[ 0][ 0] = 0.0000000; aaWAG[ 1][ 0] = 0.5515710; aaWAG[ 2][ 0] = 0.5098480; aaWAG[ 3][ 0] = 0.7389980; aaWAG[ 4][ 0] = 1.0270400;
	        aaWAG[ 5][ 0] = 0.9085980; aaWAG[ 6][ 0] = 1.5828500; aaWAG[ 7][ 0] = 1.4167200; aaWAG[ 8][ 0] = 0.3169540; aaWAG[ 9][ 0] = 0.1933350;
	        aaWAG[10][ 0] = 0.3979150; aaWAG[11][ 0] = 0.9062650; aaWAG[12][ 0] = 0.8934960; aaWAG[13][ 0] = 0.2104940; aaWAG[14][ 0] = 1.4385500;
	        aaWAG[15][ 0] = 3.3707900; aaWAG[16][ 0] = 2.1211100; aaWAG[17][ 0] = 0.1131330; aaWAG[18][ 0] = 0.2407350; aaWAG[19][ 0] = 2.0060100;
	        aaWAG[ 0][ 1] = 0.5515710; aaWAG[ 1][ 1] = 0.0000000; aaWAG[ 2][ 1] = 0.6353460; aaWAG[ 3][ 1] = 0.1473040; aaWAG[ 4][ 1] = 0.5281910;
	        aaWAG[ 5][ 1] = 3.0355000; aaWAG[ 6][ 1] = 0.4391570; aaWAG[ 7][ 1] = 0.5846650; aaWAG[ 8][ 1] = 2.1371500; aaWAG[ 9][ 1] = 0.1869790;
	        aaWAG[10][ 1] = 0.4976710; aaWAG[11][ 1] = 5.3514200; aaWAG[12][ 1] = 0.6831620; aaWAG[13][ 1] = 0.1027110; aaWAG[14][ 1] = 0.6794890;
	        aaWAG[15][ 1] = 1.2241900; aaWAG[16][ 1] = 0.5544130; aaWAG[17][ 1] = 1.1639200; aaWAG[18][ 1] = 0.3815330; aaWAG[19][ 1] = 0.2518490;
	        aaWAG[ 0][ 2] = 0.5098480; aaWAG[ 1][ 2] = 0.6353460; aaWAG[ 2][ 2] = 0.0000000; aaWAG[ 3][ 2] = 5.4294200; aaWAG[ 4][ 2] = 0.2652560;
	        aaWAG[ 5][ 2] = 1.5436400; aaWAG[ 6][ 2] = 0.9471980; aaWAG[ 7][ 2] = 1.1255600; aaWAG[ 8][ 2] = 3.9562900; aaWAG[ 9][ 2] = 0.5542360;
	        aaWAG[10][ 2] = 0.1315280; aaWAG[11][ 2] = 3.0120100; aaWAG[12][ 2] = 0.1982210; aaWAG[13][ 2] = 0.0961621; aaWAG[14][ 2] = 0.1950810;
	        aaWAG[15][ 2] = 3.9742300; aaWAG[16][ 2] = 2.0300600; aaWAG[17][ 2] = 0.0719167; aaWAG[18][ 2] = 1.0860000; aaWAG[19][ 2] = 0.1962460;
	        aaWAG[ 0][ 3] = 0.7389980; aaWAG[ 1][ 3] = 0.1473040; aaWAG[ 2][ 3] = 5.4294200; aaWAG[ 3][ 3] = 0.0000000; aaWAG[ 4][ 3] = 0.0302949;
	        aaWAG[ 5][ 3] = 0.6167830; aaWAG[ 6][ 3] = 6.1741600; aaWAG[ 7][ 3] = 0.8655840; aaWAG[ 8][ 3] = 0.9306760; aaWAG[ 9][ 3] = 0.0394370;
	        aaWAG[10][ 3] = 0.0848047; aaWAG[11][ 3] = 0.4798550; aaWAG[12][ 3] = 0.1037540; aaWAG[13][ 3] = 0.0467304; aaWAG[14][ 3] = 0.4239840;
	        aaWAG[15][ 3] = 1.0717600; aaWAG[16][ 3] = 0.3748660; aaWAG[17][ 3] = 0.1297670; aaWAG[18][ 3] = 0.3257110; aaWAG[19][ 3] = 0.1523350;
	        aaWAG[ 0][ 4] = 1.0270400; aaWAG[ 1][ 4] = 0.5281910; aaWAG[ 2][ 4] = 0.2652560; aaWAG[ 3][ 4] = 0.0302949; aaWAG[ 4][ 4] = 0.0000000;
	        aaWAG[ 5][ 4] = 0.0988179; aaWAG[ 6][ 4] = 0.0213520; aaWAG[ 7][ 4] = 0.3066740; aaWAG[ 8][ 4] = 0.2489720; aaWAG[ 9][ 4] = 0.1701350;
	        aaWAG[10][ 4] = 0.3842870; aaWAG[11][ 4] = 0.0740339; aaWAG[12][ 4] = 0.3904820; aaWAG[13][ 4] = 0.3980200; aaWAG[14][ 4] = 0.1094040;
	        aaWAG[15][ 4] = 1.4076600; aaWAG[16][ 4] = 0.5129840; aaWAG[17][ 4] = 0.7170700; aaWAG[18][ 4] = 0.5438330; aaWAG[19][ 4] = 1.0021400;
	        aaWAG[ 0][ 5] = 0.9085980; aaWAG[ 1][ 5] = 3.0355000; aaWAG[ 2][ 5] = 1.5436400; aaWAG[ 3][ 5] = 0.6167830; aaWAG[ 4][ 5] = 0.0988179;
	        aaWAG[ 5][ 5] = 0.0000000; aaWAG[ 6][ 5] = 5.4694700; aaWAG[ 7][ 5] = 0.3300520; aaWAG[ 8][ 5] = 4.2941100; aaWAG[ 9][ 5] = 0.1139170;
	        aaWAG[10][ 5] = 0.8694890; aaWAG[11][ 5] = 3.8949000; aaWAG[12][ 5] = 1.5452600; aaWAG[13][ 5] = 0.0999208; aaWAG[14][ 5] = 0.9333720;
	        aaWAG[15][ 5] = 1.0288700; aaWAG[16][ 5] = 0.8579280; aaWAG[17][ 5] = 0.2157370; aaWAG[18][ 5] = 0.2277100; aaWAG[19][ 5] = 0.3012810;
	        aaWAG[ 0][ 6] = 1.5828500; aaWAG[ 1][ 6] = 0.4391570; aaWAG[ 2][ 6] = 0.9471980; aaWAG[ 3][ 6] = 6.1741600; aaWAG[ 4][ 6] = 0.0213520;
	        aaWAG[ 5][ 6] = 5.4694700; aaWAG[ 6][ 6] = 0.0000000; aaWAG[ 7][ 6] = 0.5677170; aaWAG[ 8][ 6] = 0.5700250; aaWAG[ 9][ 6] = 0.1273950;
	        aaWAG[10][ 6] = 0.1542630; aaWAG[11][ 6] = 2.5844300; aaWAG[12][ 6] = 0.3151240; aaWAG[13][ 6] = 0.0811339; aaWAG[14][ 6] = 0.6823550;
	        aaWAG[15][ 6] = 0.7049390; aaWAG[16][ 6] = 0.8227650; aaWAG[17][ 6] = 0.1565570; aaWAG[18][ 6] = 0.1963030; aaWAG[19][ 6] = 0.5887310;
	        aaWAG[ 0][ 7] = 1.4167200; aaWAG[ 1][ 7] = 0.5846650; aaWAG[ 2][ 7] = 1.1255600; aaWAG[ 3][ 7] = 0.8655840; aaWAG[ 4][ 7] = 0.3066740;
	        aaWAG[ 5][ 7] = 0.3300520; aaWAG[ 6][ 7] = 0.5677170; aaWAG[ 7][ 7] = 0.0000000; aaWAG[ 8][ 7] = 0.2494100; aaWAG[ 9][ 7] = 0.0304501;
	        aaWAG[10][ 7] = 0.0613037; aaWAG[11][ 7] = 0.3735580; aaWAG[12][ 7] = 0.1741000; aaWAG[13][ 7] = 0.0499310; aaWAG[14][ 7] = 0.2435700;
	        aaWAG[15][ 7] = 1.3418200; aaWAG[16][ 7] = 0.2258330; aaWAG[17][ 7] = 0.3369830; aaWAG[18][ 7] = 0.1036040; aaWAG[19][ 7] = 0.1872470;
	        aaWAG[ 0][ 8] = 0.3169540; aaWAG[ 1][ 8] = 2.1371500; aaWAG[ 2][ 8] = 3.9562900; aaWAG[ 3][ 8] = 0.9306760; aaWAG[ 4][ 8] = 0.2489720;
	        aaWAG[ 5][ 8] = 4.2941100; aaWAG[ 6][ 8] = 0.5700250; aaWAG[ 7][ 8] = 0.2494100; aaWAG[ 8][ 8] = 0.0000000; aaWAG[ 9][ 8] = 0.1381900;
	        aaWAG[10][ 8] = 0.4994620; aaWAG[11][ 8] = 0.8904320; aaWAG[12][ 8] = 0.4041410; aaWAG[13][ 8] = 0.6793710; aaWAG[14][ 8] = 0.6961980;
	        aaWAG[15][ 8] = 0.7401690; aaWAG[16][ 8] = 0.4733070; aaWAG[17][ 8] = 0.2625690; aaWAG[18][ 8] = 3.8734400; aaWAG[19][ 8] = 0.1183580;
	        aaWAG[ 0][ 9] = 0.1933350; aaWAG[ 1][ 9] = 0.1869790; aaWAG[ 2][ 9] = 0.5542360; aaWAG[ 3][ 9] = 0.0394370; aaWAG[ 4][ 9] = 0.1701350;
	        aaWAG[ 5][ 9] = 0.1139170; aaWAG[ 6][ 9] = 0.1273950; aaWAG[ 7][ 9] = 0.0304501; aaWAG[ 8][ 9] = 0.1381900; aaWAG[ 9][ 9] = 0.0000000;
	        aaWAG[10][ 9] = 3.1709700; aaWAG[11][ 9] = 0.3238320; aaWAG[12][ 9] = 4.2574600; aaWAG[13][ 9] = 1.0594700; aaWAG[14][ 9] = 0.0999288;
	        aaWAG[15][ 9] = 0.3194400; aaWAG[16][ 9] = 1.4581600; aaWAG[17][ 9] = 0.2124830; aaWAG[18][ 9] = 0.4201700; aaWAG[19][ 9] = 7.8213000;
	        aaWAG[ 0][10] = 0.3979150; aaWAG[ 1][10] = 0.4976710; aaWAG[ 2][10] = 0.1315280; aaWAG[ 3][10] = 0.0848047; aaWAG[ 4][10] = 0.3842870;
	        aaWAG[ 5][10] = 0.8694890; aaWAG[ 6][10] = 0.1542630; aaWAG[ 7][10] = 0.0613037; aaWAG[ 8][10] = 0.4994620; aaWAG[ 9][10] = 3.1709700;
	        aaWAG[10][10] = 0.0000000; aaWAG[11][10] = 0.2575550; aaWAG[12][10] = 4.8540200; aaWAG[13][10] = 2.1151700; aaWAG[14][10] = 0.4158440;
	        aaWAG[15][10] = 0.3447390; aaWAG[16][10] = 0.3266220; aaWAG[17][10] = 0.6653090; aaWAG[18][10] = 0.3986180; aaWAG[19][10] = 1.8003400;
	        aaWAG[ 0][11] = 0.9062650; aaWAG[ 1][11] = 5.3514200; aaWAG[ 2][11] = 3.0120100; aaWAG[ 3][11] = 0.4798550; aaWAG[ 4][11] = 0.0740339;
	        aaWAG[ 5][11] = 3.8949000; aaWAG[ 6][11] = 2.5844300; aaWAG[ 7][11] = 0.3735580; aaWAG[ 8][11] = 0.8904320; aaWAG[ 9][11] = 0.3238320;
	        aaWAG[10][11] = 0.2575550; aaWAG[11][11] = 0.0000000; aaWAG[12][11] = 0.9342760; aaWAG[13][11] = 0.0888360; aaWAG[14][11] = 0.5568960;
	        aaWAG[15][11] = 0.9671300; aaWAG[16][11] = 1.3869800; aaWAG[17][11] = 0.1375050; aaWAG[18][11] = 0.1332640; aaWAG[19][11] = 0.3054340;
	        aaWAG[ 0][12] = 0.8934960; aaWAG[ 1][12] = 0.6831620; aaWAG[ 2][12] = 0.1982210; aaWAG[ 3][12] = 0.1037540; aaWAG[ 4][12] = 0.3904820;
	        aaWAG[ 5][12] = 1.5452600; aaWAG[ 6][12] = 0.3151240; aaWAG[ 7][12] = 0.1741000; aaWAG[ 8][12] = 0.4041410; aaWAG[ 9][12] = 4.2574600;
	        aaWAG[10][12] = 4.8540200; aaWAG[11][12] = 0.9342760; aaWAG[12][12] = 0.0000000; aaWAG[13][12] = 1.1906300; aaWAG[14][12] = 0.1713290;
	        aaWAG[15][12] = 0.4939050; aaWAG[16][12] = 1.5161200; aaWAG[17][12] = 0.5157060; aaWAG[18][12] = 0.4284370; aaWAG[19][12] = 2.0584500;
	        aaWAG[ 0][13] = 0.2104940; aaWAG[ 1][13] = 0.1027110; aaWAG[ 2][13] = 0.0961621; aaWAG[ 3][13] = 0.0467304; aaWAG[ 4][13] = 0.3980200;
	        aaWAG[ 5][13] = 0.0999208; aaWAG[ 6][13] = 0.0811339; aaWAG[ 7][13] = 0.0499310; aaWAG[ 8][13] = 0.6793710; aaWAG[ 9][13] = 1.0594700;
	        aaWAG[10][13] = 2.1151700; aaWAG[11][13] = 0.0888360; aaWAG[12][13] = 1.1906300; aaWAG[13][13] = 0.0000000; aaWAG[14][13] = 0.1614440;
	        aaWAG[15][13] = 0.5459310; aaWAG[16][13] = 0.1719030; aaWAG[17][13] = 1.5296400; aaWAG[18][13] = 6.4542800; aaWAG[19][13] = 0.6498920;
	        aaWAG[ 0][14] = 1.4385500; aaWAG[ 1][14] = 0.6794890; aaWAG[ 2][14] = 0.1950810; aaWAG[ 3][14] = 0.4239840; aaWAG[ 4][14] = 0.1094040;
	        aaWAG[ 5][14] = 0.9333720; aaWAG[ 6][14] = 0.6823550; aaWAG[ 7][14] = 0.2435700; aaWAG[ 8][14] = 0.6961980; aaWAG[ 9][14] = 0.0999288;
	        aaWAG[10][14] = 0.4158440; aaWAG[11][14] = 0.5568960; aaWAG[12][14] = 0.1713290; aaWAG[13][14] = 0.1614440; aaWAG[14][14] = 0.0000000;
	        aaWAG[15][14] = 1.6132800; aaWAG[16][14] = 0.7953840; aaWAG[17][14] = 0.1394050; aaWAG[18][14] = 0.2160460; aaWAG[19][14] = 0.3148870;
	        aaWAG[ 0][15] = 3.3707900; aaWAG[ 1][15] = 1.2241900; aaWAG[ 2][15] = 3.9742300; aaWAG[ 3][15] = 1.0717600; aaWAG[ 4][15] = 1.4076600;
	        aaWAG[ 5][15] = 1.0288700; aaWAG[ 6][15] = 0.7049390; aaWAG[ 7][15] = 1.3418200; aaWAG[ 8][15] = 0.7401690; aaWAG[ 9][15] = 0.3194400;
	        aaWAG[10][15] = 0.3447390; aaWAG[11][15] = 0.9671300; aaWAG[12][15] = 0.4939050; aaWAG[13][15] = 0.5459310; aaWAG[14][15] = 1.6132800;
	        aaWAG[15][15] = 0.0000000; aaWAG[16][15] = 4.3780200; aaWAG[17][15] = 0.5237420; aaWAG[18][15] = 0.7869930; aaWAG[19][15] = 0.2327390;
	        aaWAG[ 0][16] = 2.1211100; aaWAG[ 1][16] = 0.5544130; aaWAG[ 2][16] = 2.0300600; aaWAG[ 3][16] = 0.3748660; aaWAG[ 4][16] = 0.5129840;
	        aaWAG[ 5][16] = 0.8579280; aaWAG[ 6][16] = 0.8227650; aaWAG[ 7][16] = 0.2258330; aaWAG[ 8][16] = 0.4733070; aaWAG[ 9][16] = 1.4581600;
	        aaWAG[10][16] = 0.3266220; aaWAG[11][16] = 1.3869800; aaWAG[12][16] = 1.5161200; aaWAG[13][16] = 0.1719030; aaWAG[14][16] = 0.7953840;
	        aaWAG[15][16] = 4.3780200; aaWAG[16][16] = 0.0000000; aaWAG[17][16] = 0.1108640; aaWAG[18][16] = 0.2911480; aaWAG[19][16] = 1.3882300;
	        aaWAG[ 0][17] = 0.1131330; aaWAG[ 1][17] = 1.1639200; aaWAG[ 2][17] = 0.0719167; aaWAG[ 3][17] = 0.1297670; aaWAG[ 4][17] = 0.7170700;
	        aaWAG[ 5][17] = 0.2157370; aaWAG[ 6][17] = 0.1565570; aaWAG[ 7][17] = 0.3369830; aaWAG[ 8][17] = 0.2625690; aaWAG[ 9][17] = 0.2124830;
	        aaWAG[10][17] = 0.6653090; aaWAG[11][17] = 0.1375050; aaWAG[12][17] = 0.5157060; aaWAG[13][17] = 1.5296400; aaWAG[14][17] = 0.1394050;
	        aaWAG[15][17] = 0.5237420; aaWAG[16][17] = 0.1108640; aaWAG[17][17] = 0.0000000; aaWAG[18][17] = 2.4853900; aaWAG[19][17] = 0.3653690;
	        aaWAG[ 0][18] = 0.2407350; aaWAG[ 1][18] = 0.3815330; aaWAG[ 2][18] = 1.0860000; aaWAG[ 3][18] = 0.3257110; aaWAG[ 4][18] = 0.5438330;
	        aaWAG[ 5][18] = 0.2277100; aaWAG[ 6][18] = 0.1963030; aaWAG[ 7][18] = 0.1036040; aaWAG[ 8][18] = 3.8734400; aaWAG[ 9][18] = 0.4201700;
	        aaWAG[10][18] = 0.3986180; aaWAG[11][18] = 0.1332640; aaWAG[12][18] = 0.4284370; aaWAG[13][18] = 6.4542800; aaWAG[14][18] = 0.2160460;
	        aaWAG[15][18] = 0.7869930; aaWAG[16][18] = 0.2911480; aaWAG[17][18] = 2.4853900; aaWAG[18][18] = 0.0000000; aaWAG[19][18] = 0.3147300;
	        aaWAG[ 0][19] = 2.0060100; aaWAG[ 1][19] = 0.2518490; aaWAG[ 2][19] = 0.1962460; aaWAG[ 3][19] = 0.1523350; aaWAG[ 4][19] = 1.0021400;
	        aaWAG[ 5][19] = 0.3012810; aaWAG[ 6][19] = 0.5887310; aaWAG[ 7][19] = 0.1872470; aaWAG[ 8][19] = 0.1183580; aaWAG[ 9][19] = 7.8213000;
	        aaWAG[10][19] = 1.8003400; aaWAG[11][19] = 0.3054340; aaWAG[12][19] = 2.0584500; aaWAG[13][19] = 0.6498920; aaWAG[14][19] = 0.3148870;
	        aaWAG[15][19] = 0.2327390; aaWAG[16][19] = 1.3882300; aaWAG[17][19] = 0.3653690; aaWAG[18][19] = 0.3147300; aaWAG[19][19] = 0.0000000;

	        wagPi[ 0] = 0.08662790;
	            wagPi[ 1] = 0.04397200;
	            wagPi[ 2] = 0.03908940;
	            wagPi[ 3] = 0.05704510;
	            wagPi[ 4] = 0.01930780;
	            wagPi[ 5] = 0.03672810;
	            wagPi[ 6] = 0.05805890;
	            wagPi[ 7] = 0.08325180;
	            wagPi[ 8] = 0.02443130;
	            wagPi[ 9] = 0.04846600;
	            wagPi[10] = 0.08620970;
	            wagPi[11] = 0.06202860;
	            wagPi[12] = 0.01950273;
	            wagPi[13] = 0.03843190;
	            wagPi[14] = 0.04576310;
	            wagPi[15] = 0.06951790;
	            wagPi[16] = 0.06101270;
	            wagPi[17] = 0.01438590;
	            wagPi[18] = 0.03527420;
	            wagPi[19] = 0.07089560;

	     */

};

} /* namespace EBC */
#endif /* DEFINITIONS_H_ */
