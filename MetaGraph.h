#pragma once

#include "MetaEdge.h"
#include "Skeleton.h"
#include "SkeletonEdge.h"
//#include "Serializable.h"

namespace Roots
{

	class MetaGraph
	{
	public:
		//Serializable Implementation >> //
		/*MetaGraph(Json::Value json);
		Json::Value ToJson();*/
		// << Serializable Implementation//

		MetaGraph();
		MetaGraph(Skeleton skel, int rootNode = 0, int nextNode = 0);
		
		//void loadFromJson(Json::Value toLoad);


	private:
		Skeleton mSkeleton;
		std::vector<MetaNode> mNodes;
		std::vector<MetaEdge*> mTopLevelEdges;

	};

}