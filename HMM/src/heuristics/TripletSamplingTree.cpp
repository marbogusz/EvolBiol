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
 * TripletSamplingTree.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include "heuristics/TripletSamplingTree.hpp"
//#include <regex>
#include <stack>
#include <random>
#include <iostream>
#include <cmath>
#include "core/Definitions.hpp"

namespace EBC
{

TripletSamplingTree::TripletSamplingTree(GuideTree& gt) : distMat(gt.getDistanceMatrix())
{
	this->averageLeafbranch = 0;
	this->leafBranchSd = 0;
	DEBUG("Creating TripletSamplingTree");

	idealTreeSize = 1.0;
	this->fromNewick(gt.getNewickTree());

}

TripletSamplingTree::~TripletSamplingTree()
{
	for (auto nod : nodes)
		delete nod;
}

Node* TripletSamplingTree::mostRecentAncestor(Node* n1, Node* n2)
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

double TripletSamplingTree::distanceToParent(Node* n1, Node* par)
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


//TODO - this is no longer in use really.
void TripletSamplingTree::fromNewick(const string& newick)
{
	//parse Newick!
	//start with a bracket, end with a bracket and a semicolon
	//TODO - check if the newick string is well formed
/*
	bool endReached = false;
	unsigned int i = 0;
	unsigned int ids = 0;
	unsigned int sequenceNo;
	regex regDescDist("(([0-9]+):([0-9\\.]+)).*");
	regex regInterDist(":([0-9\\.]+).*");
	Node *tmpNode, *tmpParent, *tmpCurrent;
	stack<Node*> workNodes;

	DEBUG("K-mer newick tree : " << newick);
	//cout << newick << endl;

	string nodeName;
	double currentDistance;

	while (!endReached && i < newick.size())
	{
		if (std::isspace(newick[i]))
		{
			i++;
		}
		else if(newick[i] == '(')
		{
			DUMP("newick (");
			tmpNode = new Node(++ids);
			nodes.push_back(tmpNode);
			workNodes.push(tmpNode);
			i++;
		}
		else if(newick[i] == ',')
		{
			DUMP("newick ,");
			workNodes.pop();
			if (workNodes.empty())
			{
				//push new father
				tmpNode = new Node(++ids);
				nodes.push_back(tmpNode);
				workNodes.push(tmpNode);
			}

			tmpParent = workNodes.top();
			tmpParent->setChild(tmpCurrent);
			tmpCurrent->setParent(tmpParent);

			tmpNode = new Node(++ids);
			nodes.push_back(tmpNode);
			workNodes.push(tmpNode);

			i++;
		}
		else if(newick[i] == ')')
		{

			DUMP("newick )");
			tmpCurrent = workNodes.top();
			workNodes.pop();
			if (workNodes.size() > 0)
			{
				tmpParent = workNodes.top();
				tmpParent->setChild(tmpCurrent);
				tmpCurrent->setParent(tmpParent);
			}


			i++;
		}
		else if(newick[i] == ';')
		{
			DUMP("newick ;");
			endReached = true;
			i++;
		}
		else
		{
			smatch basematch;
			string sstr = newick.substr(i);

			//FIXME - Change ASAP this piece of code is really dangerous - get the sequence number from the seq dictionary!!!
			if (regex_match(sstr, basematch, regDescDist))
			{
				sequenceNo = std::stoi(basematch[2].str());
				currentDistance =  stof(basematch[3].str());
				tmpCurrent = workNodes.top();
				tmpCurrent->setSequenceNumber(sequenceNo);
				tmpCurrent->setDistance(currentDistance);
				tmpCurrent->setLeaf();
				this->leafNodes[sequenceNo] = tmpCurrent;
				this->averageLeafbranch += tmpCurrent->distanceToParent;
				i += basematch[1].str().size();
			}
			else if (regex_match(sstr, basematch, regInterDist))
			{
				currentDistance =   stof(basematch[1].str());
				tmpCurrent = workNodes.top();
				tmpCurrent->setDistance(currentDistance);
				i += basematch[1].str().size() + 1;
			}
		}
	}
	this->averageLeafbranch /= leafNodes.size();

	DEBUG("Newick tree parsed");

	for (auto node : leafNodes)
	{
		this->leafBranchSd += (std::pow(node.second->distanceToParent - averageLeafbranch,2));
	}
	this->leafBranchSd /= leafNodes.size();
	this->leafBranchSd = std::sqrt(leafBranchSd);
	DEBUG("Newick tree parsed");
*/
}

vector<array<unsigned int, 3> > TripletSamplingTree::sampleFromDM()
{
	vector<array<unsigned int, 3> > result;

	unsigned int s1, s2, s3;
	double remainingDistance = 2.0 * this->idealTreeSize;
	auto pair = distMat->getPairWithinDistance(0.5*this->idealTreeSize, 0.8*this->idealTreeSize);
	s1 = pair.first;
	s2 = pair.second;

	remainingDistance -= distMat->getDistance(s1, s2);

	s3 = distMat->getThirdLeafWithinDistance(remainingDistance, s1, s2);

	result.push_back({{s1,s2,s3}});
	DEBUG("Triplet tree DM sampled values : " << s1 << ", " << s2 << ", " << s3);
	cout << "Triplet tree DM sampled values : " << s1 << ", " << s2 << ", " << s3 << endl;
	return result;
}

vector<array<unsigned int, 3> > TripletSamplingTree::sampleFromTree()
{
	vector<array<unsigned int, 3> > result;

	//FIXME - remove this magic number
	pair<double,double> idealRange = make_pair(0.35,0.75);
	pair<double,double> secondaryRange = make_pair(0.15,0.85);

	unsigned int treeNo = 0;
	unsigned int maxTrees = 1;


	for (unsigned int i = 0;  i < distMat->getSize(); i++)
	{
		leafNodes[i] = nullptr;
	}
	availableNodes = leafNodes;

	bool found = false;
	unsigned int s1,s2;

	double tmpd1, tmpd2;

	auto vecPairsDesired = distMat->getPairsWithinDistance(idealRange.first, idealRange.second);
	auto vecPairsSecondary = distMat->getPairsWithinDistance(secondaryRange.first, secondaryRange.second);

	if (vecPairsDesired.size() == 0 && vecPairsSecondary.size() == 0){
		vecPairsSecondary.push_back(distMat->getPairWithinDistance(idealRange.first, idealRange.second));
		DUMP("Triplet sampling tree : no pairs within desired distance found, geting the closest one");
	}
	//Check ideal range
	for(auto pr : vecPairsDesired){
		if (treeNo >maxTrees)
			break;
		s1 = pr.first;
		s2 = pr.second;

		if (availableNodes.find(s1) == availableNodes.end()
				|| availableNodes.find(s2) == availableNodes.end())
			continue;

		for (auto nd : availableNodes)
		{

			if (nd.first == s1 || nd.first == s2)
				continue;
			tmpd1 =  distMat->getDistance(nd.first, s1);
			tmpd2 =  distMat->getDistance(nd.first, s2);

			if (isWithinRange(tmpd1, idealRange) && isWithinRange(tmpd1, idealRange)) {
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
			else{
				continue;
			}
		}
	}
	//not found
	if(!found){
		for(auto pr : vecPairsSecondary){
			if (treeNo >maxTrees)
				break;
			s1 = pr.first;
			s2 = pr.second;

			if (availableNodes.find(s1) == availableNodes.end()
					|| availableNodes.find(s2) == availableNodes.end())
				continue;

			for (auto nd : availableNodes)
			{

				if (nd.first == s1 || nd.first == s2)
					continue;
				tmpd1 =  distMat->getDistance(nd.first, s1);
				tmpd2 =  distMat->getDistance(nd.first, s2);

				if (isWithinRange(tmpd1, secondaryRange) && isWithinRange(tmpd1, secondaryRange)) {
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
				else{
					continue;
				}
			}
		}
	}
	//Still not found?
	if(!found){
		for(auto pr : vecPairsDesired){
			if (treeNo >maxTrees)
				break;
			s1 = pr.first;
			s2 = pr.second;

			if (availableNodes.find(s1) == availableNodes.end()
					|| availableNodes.find(s2) == availableNodes.end())
				continue;

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
	}
	//still nothing
	if(!found){
		for(auto pr : vecPairsSecondary){
			if (treeNo >maxTrees)
				break;
			s1 = pr.first;
			s2 = pr.second;

			if (availableNodes.find(s1) == availableNodes.end()
					|| availableNodes.find(s2) == availableNodes.end())
				continue;


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
		}
	}
	//Nothing within those ranges - just get something
	if(!found){
		//check secondary range
		auto pr = vecPairsSecondary[0]; //get the best one
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
			treeNo++;
		}

	}
	INFO(treeNo << " triplets sampled");
	return result;
}

/*
vector<array<unsigned int, 3> > TripletSamplingTree::sampleFromTree()
{
	vector<array<unsigned int, 3> > result;

	//FIXME - remove this magic number
	pair<double,double> idealRange = make_pair(0.35,0.75);
	pair<double,double> secondaryRange = make_pair(0.2,0.85);

	Node *firstNd, *secondNd;
	unsigned int treeNo = 0;
	//copy

	distMat->getSize();

	for (unsigned int i = 0;  i < distMat->getSize(); i++)
	{
		leafNodes[i] = nullptr;
	}

	availableNodes = leafNodes;


	//randomly select 1 leaf
	double treeSize;
	this->idealTreeSize = 3 * this->averageLeafbranch;
	bool found;
	unsigned int s1,s2;

	double tmpd1, tmpd2;

	//get all the pairs within a good range
	found = false;

	auto vecPairs = distMat->getPairsWithinDistance(idealRange.first, idealRange.second);

	if (vecPairs.size() == 0){
		vecPairs.push_back(distMat->getPairWithinDistance(idealRange.first, idealRange.second));
		DUMP("Triplet sampling tree : no pairs within desired distance found");
	}
	for(auto pr : vecPairs){
		if (treeNo >2)
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
}
*/
} /* namespace EBC */
