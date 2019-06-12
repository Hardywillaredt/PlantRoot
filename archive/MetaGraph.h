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
		MetaGraph(std::string filename);

		void BuildFlatEdges();

		void SaveToFile(std::string filename);
		
		//void loadFromJson(Json::Value toLoad);

		Skeleton* GetSkeleton();
		std::vector<MetaEdge*> GetTopLevelEdges();
		std::vector<MetaEdge*> GetFlatEdges();

	private:
		Skeleton mSkeleton;
		std::vector<MetaEdge*> mTopLevelEdges;
		std::vector<MetaEdge*> mFlatEdges;

	};

}
