/*
 * TripletSamplingTree.cpp
 *
 *  Created on: Sep 1, 2014
 *      Author: root
 */

#include "heuristics/TripletSamplingTree.hpp"
#include <regex>
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



void TripletSamplingTree::fromNewick(const string& newick)
{
	//parse Newick!
	//start with a bracket, end with a bracket and a semicolon
	//TODO - check if the newick string is well formed
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
	pair<double,double> secondaryRange = make_pair(0.2,0.85);

	Node *firstNd, *secondNd;
	unsigned int treeNo = 0;
	//copy
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

		firstNd = leafNodes[pr.first];
		secondNd  = leafNodes[pr.second];
		s1 = pr.first;
		s2 = pr.second;



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
		firstNd = leafNodes[pr.first];
		secondNd  = leafNodes[pr.second];
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
			unsigned int id = 0;

			while(id == s1 || id == s2)
				id ++;
			result.push_back({{s1,s2,id}});
		}

	}
	/*
	while(availableNodes.size() >=3 && treeNo < 3)
	{
		auto pairN = distMat->getPairWithinDistance(idealRange.first, idealRange.second);
		firstNd = leafNodes[pairN.first];
		s1 = pairN.first;
		s2 = pairN.second;
		secondNd  = leafNodes[pairN.second];



		found = false;
		for (auto nd : availableNodes)
		{

			if (nd.second->nodeId == firstNd->nodeId || nd.second->nodeId == secondNd->nodeId)
				continue;
			tmpd1 =  distMat->getDistance(nd.second->nodeId, firstNd->nodeId);
			tmpd2 =  distMat->getDistance(nd.second->nodeId, secondNd->nodeId);

			if (isWithinRange(tmpd1, idealRange)) {
				result.push_back({{s2,s1,nd.first}});
				DEBUG("Sampled triplet : " << s2 << "\t\t" << s1 << "\t\t" << nd.first);
				availableNodes.erase(pairN.first);
				availableNodes.erase(pairN.second);
				availableNodes.erase(nd.first);
				found = true;
				treeNo++;
				break;
			}
			else if(isWithinRange(tmpd2, idealRange)){
				result.push_back({{s1,s2,nd.first}});
				DEBUG("Sampled triplet : " << s1 << "\t\t" << s2 << "\t\t" << nd.first);
				availableNodes.erase(pairN.first);
				availableNodes.erase(pairN.second);
				availableNodes.erase(nd.first);
				found = true;
				treeNo++;
				break;
			}
			else{
				continue;
			}
		}
		//min 1 triplet, max 3 triplets
		if (!found && treeNo < 1)
		{
			// grab any unused leaf node
			result.push_back({{s1,s2,availableNodes.begin()->first}});
			availableNodes.erase(availableNodes.begin());
		}
		treeNo++;

	}
	*/
	//distance to root

	//get another one with distance between 0.3 and 0.7
	//what about root node - root node has no parent set - go through list and check ?????
	//find their common ancestor
	//based on the ancestor, calculate the distance to it and target the next taxon based on the distance to root + common to root
	//invalidate the leaves laying on the common path.


	return result;
}

} /* namespace EBC */
