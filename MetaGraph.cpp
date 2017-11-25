#include "MetaGraph.h"


namespace Roots
{
	MetaGraph::MetaGraph() 
	{ 
		mSkeleton = Skeleton();
		mNodes = {};
		mTopLevelEdges = {};
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
			root = new MetaEdge(&mSkeleton, node1, node2);
			mTopLevelEdges.push_back(root);
		}
		else
		{
			mNodes.push_back(MetaNode(&mSkeleton, node1));
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
					root = new MetaEdge(&mSkeleton, node1, node2);
					mTopLevelEdges.push_back(root);
				}
				else
				{
					mNodes.push_back(MetaNode(&mSkeleton, node1));
				}
			}
		}
	}

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