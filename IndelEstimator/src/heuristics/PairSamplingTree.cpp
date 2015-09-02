/*
 * TripletSamplingTree.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include <heuristics/PairSamplingTree.hpp>
#include <stack>
#include <random>
#include <iostream>
#include <cmath>
#include "core/Definitions.hpp"

namespace EBC
{

PairSamplingTree::PairSamplingTree(DistanceMatrix* dm) : distMat(dm)
{
	DEBUG("Creating PairSamplingTree");
}

PairSamplingTree::~PairSamplingTree()
{
	for (auto nod : nodes)
		delete nod;
}

Node* PairSamplingTree::mostRecentAncestor(Node* n1, Node* n2)
{
	Node* tmpNode1 = n1;
	Node* tmpNode2 = n2;

	while(tmpNode1 != NULL)
	{
		tmpNode2 = n2;
		while(tmpNode2 != NULL)
		{
			if (*tmpNode1 == *tmpNode2)
				return tmpNode1;
			tmpNode2 = tmpNode2->parent;
		}
	tmpNode1 = tmpNode1->parent;
	}
	return NULL;
}

double PairSamplingTree::distanceToParent(Node* n1, Node* par)
{
	double distance = 0;
	Node* tmp = n1;

	if (*n1 == *par)
		return 0;

	while (tmp != par)
	{
		distance += tmp->distanceToParent;
		tmp = tmp->parent;
	}
	return distance;
}

vector<array<unsigned int, 2> > PairSamplingTree::sampleFromTree()
{

	vector<array<unsigned int, 2> > result;

	//FIXME - remove this magic number
	pair<double,double> idealRange = make_pair(0.35,0.75);
	pair<double,double> secondaryRange = make_pair(0.2,0.85);

	Node *firstNd, *secondNd;
	unsigned int pairsNo = 0;
	//copy

	distMat->getSize();

	for (unsigned int i = 0;  i < distMat->getSize(); i++)
	{
		leafNodes[i] = nullptr;
	}

	availableNodes = leafNodes;

	bool found;
	unsigned int s1,s2;

	double tmpd1, tmpd2;

	//get all the pairs within a good range
	found = false;

	vector<pair<unsigned int, unsigned int> > vecPairs = distMat->getPairsWithinDistance(idealRange.first, idealRange.second);

	if (vecPairs.size() == 0){
		vecPairs.push_back(distMat->getPairWithinDistance(idealRange.first, idealRange.second));
		DEBUG("Pair sampling tree : no pairs within desired distance found");
	}

	for(auto pr : vecPairs)
		result.push_back({pr.first,pr.second});

	return result;

	/*
	//FIXME - magagic numbers
	while(pairsNo < 20)
	{
		vecPairs.push_back(distMat->getPairWithinDistance(idealRange.first, idealRange.second));
	}


	for(auto pr : vecPairs){
		//FIXME - magagic numbers
		if (pairNo > 20)
			break;
		s1 = pr.first;
		s2 = pr.second;

		if (availableNodes.find(s1) == availableNodes.end()
				|| availableNodes.find(s2) == availableNodes.end())
			continue;

		//break vec pairs
		if (treeNo >= 3)
			break;

		for (auto nd : availableNodes)
		{

			if (nd.first == s1 || nd.first == s2)
				continue;
			tmpd1 =  distMat->getDistance(nd.first, s1);
			tmpd2 =  distMat->getDistance(nd.first, s2);

			if (isWithinRange(tmpd1, idealRange)) {
				result.push_back({{s2,s1,nd.first}});
				DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else if(isWithinRange(tmpd2, idealRange)){
				result.push_back({{s1,s2,nd.first}});
				DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else{
				continue;
			}
		}
	}
	if(!found){
		auto pr = vecPairs[0]; //get the best one
		s1 = pr.first;
		s2 = pr.second;

		for (auto nd : availableNodes)
		{

			if (nd.first == s1 || nd.first == s2)
				continue;
			tmpd1 =  distMat->getDistance(nd.first, s1);
			tmpd2 =  distMat->getDistance(nd.first, s2);


			if (isWithinRange(tmpd1, secondaryRange)) {
				result.push_back({{s2,s1,nd.first}});
				DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else if(isWithinRange(tmpd2, secondaryRange)){
				result.push_back({{s1,s2,nd.first}});
				DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << nd.first);
				availableNodes.erase(s1);
				availableNodes.erase(s2);
				availableNodes.erase(nd.first);
				distMat->invalidate(pr);
				found = true;
				treeNo++;
				break;
			}
			else{
				continue;
			}
		}
		if(!found){
			//nothing was found. this means that either all the distances are very small or big
			//small = get the biggest, big - get the smallest
			unsigned int id = 0;
			double tdist = 0;
			//bool largeDistances = false;

			//if (distMat->getDistance(s1, s2) > secondaryRange.second){
			//	largeDistances = true;
			//}
			double currentDist = 100;

			availableNodes.erase(s1);
			availableNodes.erase(s2);

			for(auto nd : availableNodes)
			{
				//deviation from the sweetspot
				//FIXME - magicnumber for the sweetspot
				tdist = abs(0.5 - distMat->getDistance(nd.first, s1)) + abs(0.5 - distMat->getDistance(nd.first, s2));

					if (tdist < currentDist){
						currentDist = tdist;
						id = nd.first;
					}

			}

			tmpd1 = distMat->getDistance(s1, id);
			tmpd2 = distMat->getDistance(s2, id);

			if (tmpd1 > distMat->getDistance(s1,s2)){
				if(tmpd1 < tmpd2){
					result.push_back({{s2,s1,id}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << id);
				}
				else{
					result.push_back({{s1,s2,id}});
					DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << id);
				}
			}
			else{
				if(tmpd1 < tmpd2){
					result.push_back({{s1,s2,id}});
					DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << id);
				}
				else{
					result.push_back({{s2,s1,id}});
					DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << id);
				}
			}

//			result.push_back({{s1,s2,id}});
//			DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << id);
		}

	}



	return result;
	*/
}

} /* namespace EBC */
