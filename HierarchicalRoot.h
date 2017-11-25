//#pragma once
//#include "Skeleton.h"
//#include "Listener.h"
////////// NO LONGER A THING//////////////////
//
//namespace Roots
//{
//	/*
//	Annotated root type, as determined by algorithmic or user labelling.
//	*/
//	enum HierarchicalRootType
//	{
//		NoType=0,
//		Type2,
//		Type3,
//
//
//		NumTypes
//	};
//
//	class HierarchicalRoot : public Listener
//	{
//	public:
//		HierarchicalRoot(edgePtrList edges = edgePtrList(), HierarchicalRoot *parent = nullptr, int order = -1, HierarchicalRootType = NoType);
//		
//		virtual void SpeakerUpdated(void *updated, std::string className);
//		virtual void SpeakerDestroyed(void *deleted, std::string className);
//
//		int getOrder();
//
//	private:
//		HierarchicalRoot *mParent;
//		std::vector<HierarchicalRoot*> mChildren;
//		edgePtrList mOrderedEdges;
//		int mOrder;
//		HierarchicalRootType mType;
//		bool mIsValid;
//		int endVert1, endVert2;
//	};
//}