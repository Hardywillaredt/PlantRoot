#include "MetaNode.h"
#include <algorithm>

#include <iostream>
#include <string>


namespace Roots
{

	MetaNode::MetaNode()
	{
		mSrcSkeleton = nullptr;
		mIdx = 0;
		mSrcVert = 0;
		mNeighbors = {};
	}

	MetaNode::MetaNode(Skeleton *srcSkeleton, int srcVert, int idx)
	{
		mSrcSkeleton = srcSkeleton;
		mSrcVert = srcVert;
		mIdx = idx;
		mNeighbors = {};
	}

	MetaNode::MetaNode(Skeleton* srcSkeleton, int idx, int srcVert, std::vector<int> neighbors)
	{
		mSrcSkeleton = srcSkeleton;
		mSrcVert = srcVert;
		mIdx = idx;
		mNeighbors = neighbors;
	}



	std::ostream& operator<<(std::ostream& out, const MetaNode& node)
	{
		out << node.mIdx << " " << node.mSrcVert << " ";
		for each(int neighbor in node.mNeighbors)
		{
			out << neighbor << " ";
		}
		out << std::endl;
		return out;
	}
}
