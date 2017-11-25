//#include "HierarchicalRoot.h"
//#include <typeinfo>
//
//
//
//namespace Roots
//{
//	HierarchicalRoot::HierarchicalRoot(edgePtrList edges, HierarchicalRoot *parent, int order, HierarchicalRootType type)
//	{
//		if (edges.size() > 0)
//		{
//			mOrderedEdges = tidyUpRoot(edges, mIsValid, endVert1, endVert2);
//		}
//		for each(SkeletonEdge * edge in mOrderedEdges)
//		{
//			edge->RegisterListener(this);
//		}
//		
//		if (parent != nullptr)
//		{
//			mOrder = parent->getOrder() + 1;
//		}
//		else if (order != -1)
//		{
//			mOrder = order;
//		}
//		else
//		{
//			mOrder = 0;
//		}
//
//		mType = type;
//	}
//
//	int HierarchicalRoot::getOrder()
//	{
//		return mOrder;
//	}
//
//	void HierarchicalRoot::SpeakerUpdated(void *updated, std::string className)
//	{
//		if (className == SkeletonEdge::ClassName())
//		{
//			SkeletonEdge* updatedEdge = (SkeletonEdge*)updated;
//			//do nothing?...
//		}
//		else
//		{
//			//no other speaker types
//		}
//	}
//
//	void HierarchicalRoot::SpeakerDestroyed(void *destroyed, std::string className)
//	{
//		if (mOrderedEdges.size() <= 0)
//		{
//			return;
//		}
//		if (className == SkeletonEdge::ClassName() && destroyed != nullptr)
//		{
//			SkeletonEdge* destroyedEdge = (SkeletonEdge*)destroyed;
//			
//			if (mOrderedEdges[0][0] != destroyedEdge[0] && mOrderedEdges[mOrderedEdges.size() - 1][0] != destroyedEdge[0])
//			{
//				mIsValid = false;
//			}
//
//			for (int i = 0; i < mOrderedEdges.size(); ++i)
//			{
//				if (mOrderedEdges[i][0] == destroyedEdge[0])
//				{
//					mOrderedEdges.erase(mOrderedEdges.begin() + i, mOrderedEdges.begin() + 1);
//				}
//			}
//			if (mOrderedEdges.size() < 1)
//			{
//				mIsValid = false;
//			}
//		}
//	}
//}