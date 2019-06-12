#pragma once

#include "Geometry.h"
#include "SkeletonEdge.h"
#include "Skeleton.h"
#include "MetaNodeFactory.h"

#include "boost/python.hpp"


#include <map>

namespace Roots
{

	bool stringToBool(std::string in);

	std::string boolToString(bool isTrue);

	class MetaEdge
	{
	public:
		//map of skeleton vertex ids to their corresponding metanodes.  If there is not an entry for a given index, then there is
		//not a metanode defined for that skeleton vertex.
		
		//set of indicators whether the vertex at each index in the skeleton vertex list is contained in a MetaEdge
		static std::vector<bool> vertsInMetaEdge;

		static std::string metaEdgeString;


		MetaEdge(int idx = -1);
		MetaEdge(Skeleton *srcSkel, int idx, int startVertex, int secondVertex, MetaEdge *parent=nullptr);
		MetaEdge(const MetaEdge &toCopy);
		MetaEdge(Skeleton *srcSkel, int aIdx, int order, bool aIsPartOfLoop, int parentSkelVertex, int childSkelVertex, std::vector<int> skelVertices, std::vector<MetaEdge*> aChildren);



		bool validateVertexPair(int vert1, int vert2, Skeleton *srcSkel, int &alternateVert2);
		

	public:
		//locally defined items (these are loaded and saved)
		std::vector<int> mVertices;
		int mOrder;
		std::vector<int> mChildrenIndices;
		int mParentIndex;
		int mParentNodeIndex, mChildNodeIndex;
		bool mIsPartOfLoop;
		int mIdx;

		//items defined by reference to an object (not loaded/saved)
		//the meta edge only owns the children pointers, do not delete any of the others on deconstruction
		Skeleton *mSrcSkel;
		MetaNode *mParentNode, *mChildNode;
		MetaEdge *mParent;
		std::vector<MetaEdge*> mChildren;
		std::vector<SkeletonEdge*> mEdges;
		

	public:
		std::vector<int> GetVertices();
		std::vector<SkeletonEdge*> GetEdges();
		Skeleton* GetSkeleton();
		MetaNode* GetParentNode();
		MetaNode* GetChildNode();
		std::vector<MetaEdge*> GetChildren();
		MetaEdge* GetParent();
		int GetId();
		int GetOrder();
		bool isInLoop();
		
		void SetParent(MetaEdge *parent);

		friend std::ostream& operator<<(std::ostream& out, const MetaEdge& edge);
		

		//given a set of lines corresponding to all the metaedges in the metagraph, recursively load each meta edge and its children
		//from those lines.  Build only the tree to which this meta edge is connected, other trees must built recursively 
		//static MetaEdge* recursiveLoad(std::vector<std::string> lines, int lineToLoad, Skeleton* srcSkeleton, bool isTopDown = false);

	};

	
	struct MetaEdgeFactory
	{
		static int maxId;
		static std::vector<MetaEdge> allEdges;


		static MetaEdge* CreateEdge();
		static MetaEdge* CreateEdge(Skeleton *srcSkel, int startVertex, int secondVertex, MetaEdge *parent = nullptr);
		static MetaEdge* CreateEdge(Skeleton *srcSkel, int aIdx, int order, bool aIsPartOfLoop, int parentSkelVertex, int childSkelVertex, std::vector<int> skelVertices, std::vector<MetaEdge*> aChildren);

		static MetaEdge* recursiveLoadEdges(std::vector<std::string> lines, int lineToLoad, Skeleton* srcSkeleton, bool isTopDown = false);

	private:

		/*
		Builds all descendants of this edge at the child node, and if it
		is the absolute parent edge, then also builds out a hierarchy from
		the parent node (as this is actually a child node to this edge)
		*/
		static void buildDescendants(MetaEdge *parentEdge);

		/*
		Build all descendants of this edge at the given vertex.  Detects
		any loops in the MetaEdge structure by determining if any of the
		neighbors of that vertex are already part of a MetaEdge
		(that isn't this meta edge)
		*/
		static void buildDescendantsAtNode(MetaEdge *parentEdge, int vertIndex);
	};

}
