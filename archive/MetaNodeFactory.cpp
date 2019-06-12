#include "MetaNodeFactory.h"


std::map<int, Roots::MetaNode*> Roots::MetaNodeFactory::vertNodeMap = std::map<int, Roots::MetaNode*>();
std::vector<Roots::MetaNode> Roots::MetaNodeFactory::nodes = std::vector<Roots::MetaNode>();


int Roots::MetaNodeFactory::maxId = -1;

namespace Roots{
	std::istream& MetaNodeFactory::LoadNode(std::istream& in, Skeleton *srcSkel)
	{
		std::string line;
		std::getline(in, line);

		std::vector<std::string> words;
		boost::split(words, line, boost::is_any_of(" "));
		int idx = boost::lexical_cast<int>(words[0]);
		int srcVert = boost::lexical_cast<int>(words[1]);
		std::vector<int> neighbors;
		for (int i = 2; i < words.size(); ++i)
		{
			neighbors.push_back(boost::lexical_cast<int>(words[i]));
		}

		findOrCreateNode(srcSkel, idx, srcVert, neighbors);

		return in;
	}


	MetaNode* MetaNodeFactory::findOrCreateNode(Skeleton* srcSkel, int srcVert)
	{
		//if the meta node isn't mapped, create it
		if (vertNodeMap.count(srcVert) == 0)
		{
			int id;
			if (MetaNodeFactory::maxId < 0)
			{
				id = 0;
				MetaNodeFactory::maxId = 0;
			}
			else
			{
				++MetaNodeFactory::maxId;
				id = MetaNodeFactory::maxId;
			}

			MetaNode node = MetaNode(srcSkel, srcVert, id);
			nodes.push_back(node);
			vertNodeMap[srcVert] = &nodes[id];
		}

		return vertNodeMap[srcVert];
	}

	MetaNode* MetaNodeFactory::findOrCreateNode(Skeleton* srcSkeleton, int idx, int srcVert, std::vector<int> neighbors)
	{
		if (maxId < idx)
		{
			maxId = idx;
			nodes.resize(maxId + 1, MetaNode());

			nodes[idx] = MetaNode(srcSkeleton, idx, srcVert, neighbors);
			vertNodeMap[srcVert] = &nodes[idx];
			nodes.resize(idx);
		}

		return vertNodeMap[srcVert];
	}

	void MetaNodeFactory::RemoveNode(MetaNode toRemove)
	{
		vertNodeMap.erase(toRemove.mSrcVert);
		for (int laterId = toRemove.mIdx; laterId < nodes.size(); ++laterId)
		{
			nodes[laterId].mIdx -= 1;
		}
		for (int i = 0; i < nodes.size(); ++i)
		{
			for (int k = 0; k < nodes[i].mNeighbors.size(); ++k)
			{
				if (nodes[i].mNeighbors[k] > toRemove.mIdx)
				{
					nodes[i].mNeighbors[k] -= 1;
				}
				else if (nodes[i].mNeighbors[k] == toRemove.mIdx)
				{
					nodes[i].mNeighbors.erase(nodes[i].mNeighbors.begin() + k);
				}
			}
		}
		nodes.erase(nodes.begin() + toRemove.mIdx);
		--maxId;
	}
}
