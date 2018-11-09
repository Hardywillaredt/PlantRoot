#pragma once

#include "MetaNode.h"
#include <map>

namespace Roots
{

	struct MetaNodeFactory
	{
		static int maxId;

		static std::map<int, MetaNode*> vertNodeMap;
		static std::vector<MetaNode> nodes;
		static MetaNode* findOrCreateNode(Skeleton *srcSkeleton, int srcVert);
		static MetaNode* findOrCreateNode(Skeleton *srcSkeleton, int idx, int srcVert, std::vector<int> neighbors);
		static void RemoveNode(MetaNode toRemove);

		std::istream& LoadNode(std::istream& in, Skeleton *srcSkel);
	};
}
