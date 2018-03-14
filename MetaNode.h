#pragma once

#include "Skeleton.h"

#include "boost/python.hpp"
#include <map>


namespace Roots
{
	struct MetaNode
	{

		MetaNode();
		MetaNode(Skeleton *srcSkeleton, int srcVert, int idx);
		MetaNode(Skeleton *srcSkeleton, int idx, int srcVert, std::vector<int> neighbors);


		friend std::ostream& operator<<(std::ostream& out, const MetaNode& node);
		//static std::istream& LoadNode(std::istream& in, Skeleton *srcSkel);

		Skeleton *mSrcSkeleton;
		int mSrcVert;
		//node id (first node created should be )
		int mIdx;
		//vector of neighboring node id's
		std::vector<int> mNeighbors;
	};
}