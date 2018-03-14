#include "MetaGraph.h"

#include <fstream>

namespace Roots
{
	MetaGraph::MetaGraph() 
	{ 
		mSkeleton = Skeleton();
		mTopLevelEdges = {};
		BuildFlatEdges();
	}

	MetaGraph::MetaGraph(Skeleton skel, int rootNode, int nextNode)
	{
		mSkeleton = skel;
		mTopLevelEdges = {};
		int node1 = rootNode;
		int node2 = nextNode;
		bool hasNeighbor = false;
		if (node2 == node1)
		{
			//if the input nodes are uninitialized/the same, then check if 
			//the first input node has any neighbors and reassign the second
			//node to be that neighbor.
			if (mSkeleton.mNeighbors[node1].size() > 0)
			{
				node2 = mSkeleton.mNeighbors[node1][0];
				hasNeighbor = true;
			}
		}
		else
		{
			for(int i = 0; i < mSkeleton.mNeighbors[node1].size(); ++i)
			{
				if(mSkeleton.mNeighbors[node1][i] == node2)
				{
					hasNeighbor = true;
				}
			}
			if(!hasNeighbor)
			{
				if (mSkeleton.mNeighbors[node1].size() > 0)
				{
					node2 = mSkeleton.mNeighbors[node1][0];
					hasNeighbor = true;
				}
			}
		}

		MetaEdge *root;

		if (hasNeighbor)
		{
			root = MetaEdgeFactory::CreateEdge(&mSkeleton, node1, node2);
			mTopLevelEdges.push_back(root);
		}
		else
		{
			MetaNodeFactory::findOrCreateNode(&mSkeleton, node1);
		}

		for (int i = 0; i < MetaEdge::vertsInMetaEdge.size(); ++i)
		{
			if (!MetaEdge::vertsInMetaEdge[i])
			{
				int node1 = i;
				hasNeighbor = false;
				if (mSkeleton.mNeighbors[node1].size() > 0)
				{
					node2 = mSkeleton.mNeighbors[node1][0];
					hasNeighbor = true;
				}
				
				if (hasNeighbor)
				{
					root = MetaEdgeFactory::CreateEdge(&mSkeleton, node1, node2);
					mTopLevelEdges.push_back(root);
				}
				else
				{
					MetaNodeFactory::findOrCreateNode(&mSkeleton, node1);
				}
			}
		}

		BuildFlatEdges();
	}

	void MetaGraph::BuildFlatEdges()
	{
		mFlatEdges = std::vector<MetaEdge*>(MetaEdgeFactory::maxId + 1);
		std::vector<MetaEdge*> toDealWith = mTopLevelEdges;
		while (toDealWith.size() > 0)
		{
			std::vector<MetaEdge*> dealWithChildren = {};
			for each(MetaEdge* edge in toDealWith)
			{
				for each(MetaEdge* edgeChild in edge->GetChildren())
				{
					dealWithChildren.push_back(edgeChild);
				}
				mFlatEdges[edge->GetId()] = edge;
			}
			toDealWith = dealWithChildren;
		}
	}

	void MetaGraph::SaveToFile(std::string filename)
	{
		std::ofstream fileStream;
		fileStream.open(filename);

		fileStream << mSkeleton;

		fileStream << MetaNodeFactory::nodes.size() << std::endl;

		for (int i = 0; i < MetaNodeFactory::nodes.size(); ++i)
		{
			fileStream << MetaNodeFactory::nodes[i];
		}

		for each(MetaEdge *topEdge in mTopLevelEdges)
		{
			fileStream << topEdge->GetId() << " ";
		}
		fileStream << std::endl;

		BuildFlatEdges();
		for (int i = 0; i < mFlatEdges.size(); ++i)
		{
			fileStream << mFlatEdges[i];
		}
	}

	MetaGraph::MetaGraph(std::string filename)
	{
		mSkeleton = Skeleton();

		std::ifstream fileStream;
		fileStream.open(filename);

		fileStream >> mSkeleton;

		std::string line;
		std::getline(fileStream, line);

		int numNodes = boost::lexical_cast<int>(line);
		for (int node = 0; node < numNodes; ++node)
		{
			//MetaNode::LoadNode(fileStream, &mSkeleton);
		}

		std::vector<int> topLevelEdgeIndices;
		std::getline(fileStream, line);
		std::vector<std::string> edgeIndiceStrings;
		boost::split(edgeIndiceStrings, line, boost::is_any_of(" "));
		for each(std::string edgeIndex in edgeIndiceStrings)
		{
			topLevelEdgeIndices.push_back(boost::lexical_cast<int>(edgeIndex));
		}

		mTopLevelEdges = {};
		std::vector<std::string> edgeLines = {};
		while (std::getline(fileStream, line))
		{
			edgeLines.push_back(line);
		}

		for each(int topLevelEdgeIndex in topLevelEdgeIndices)
		{
			mTopLevelEdges.push_back(MetaEdgeFactory::recursiveLoadEdges(edgeLines, topLevelEdgeIndex, &mSkeleton, false));
		}

		BuildFlatEdges();
	}

	Skeleton* MetaGraph::GetSkeleton() { return &mSkeleton; }
	std::vector<MetaEdge*> MetaGraph::GetTopLevelEdges() { return mTopLevelEdges; }
	std::vector<MetaEdge*> MetaGraph::GetFlatEdges() { return mFlatEdges; }

	//Json::Value MetaGraph::ToJson()
	//{
	//	Json::Value metaJson;
	//	metaJson["MetaGraph"] = Json::Value();
	//	Json::Value graphJson = metaJson["MetaGraph"];
	//	
	//	graphJson["Skeleton"] = mSkeleton.ToJson();
	//	
	//	graphJson["MetaNodes"] = Json::Value(Json::ValueType::arrayValue);
	//	for each(auto vertNodePair in MetaEdge::vertNodeMap)
	//	{
	//		vertNodePair.second.addToJsonArray(graphJson["MetaNodes"]);
	//	}
	//	for each(MetaNode node in mNodes)
	//	{
	//		node.addToJsonArray(graphJson["MetaNodes"]);
	//	}

	//	graphJson["MetaEdges"] = Json::Value(Json::ValueType::arrayValue);
	//	for each(MetaEdge *topEdge in mTopLevelEdges)
	//	{
	//		topEdge->recursiveAddToJson(graphJson["MetaEdges"]);
	//	}

	//	return metaJson;
	//}

	//void MetaGraph::loadFromJson(Json::Value toLoad)
	//{
	//	//todo check for existence of a metagraph in json.
	//	//elsewise, create default metagraph from skeleton json
	//	Json::Value graphJson = toLoad["MetaGrah"];

	//	
	//	mSkeleton.LoadFromJson(graphJson["Skeleton"]);
	//}
}