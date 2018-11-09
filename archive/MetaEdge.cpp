#include "MetaEdge.h"
#include <algorithm>


std::vector<bool> Roots::MetaEdge::vertsInMetaEdge = std::vector<bool>();

std::vector<Roots::MetaEdge> Roots::MetaEdgeFactory::allEdges = std::vector<Roots::MetaEdge>();
int Roots::MetaEdgeFactory::maxId = -1;
std::string Roots::MetaEdge::metaEdgeString = "metaEdge";


namespace
{
	//meta edge strings for json mapping
	const std::string VerticesString = "Vertices";
	const std::string IndexString = "Index";
	const std::string OrderString = "Order";
	const std::string ChildrenIndicesString = "ChildrenIndices";
	const std::string ParentIndexString = "ParentIndex";
	const std::string ParentNodeIndexString = "ParentNodeIndex";
	const std::string ChildNodeIndexString = "ChildNodeIndex";
	const std::string IsPartOfLoopString = "IsPartOfLoop";
}


namespace Roots
{
	

	bool stringToBool(std::string in)
	{
		int isTrue = boost::lexical_cast<int>(in);
		return isTrue == 1 ? false : true;
	}

	std::string boolToString(bool isTrue)
	{
		return isTrue ? "1" : "0";
	}



	MetaEdge::MetaEdge(int idx)
	{
		mVertices = {};
		mIdx = idx;
		mOrder = -1;
		mIsPartOfLoop = false;
		mSrcSkel = nullptr;
		mParentNode = nullptr;
		mParentNodeIndex = -1;
		mChildNode = nullptr;
		mChildNodeIndex = -1;
		mParent = nullptr;
		mParentIndex = -1;
		mChildren = {};
		mChildrenIndices = {};
		
	}

	MetaEdge::MetaEdge(Skeleton *srcSkel, int idx, int startVertex, int secondVertex, MetaEdge *parent)
	{
		mVertices = {};
		mIdx = idx;
		
		mParent = parent;
		if(mParent != nullptr)
		{
			mParentIndex = mParent->GetId();
			mOrder = mParent->GetOrder() + 1;
		}
		else
		{
			mParentIndex = -1;
			mOrder = 0;
		}
		mSrcSkel = srcSkel;
		mParentNode = MetaNodeFactory::findOrCreateNode(mSrcSkel, startVertex);
		mParentNodeIndex = mParentNode->mIdx;
		mParent = parent;
		mChildren = {};
		mChildrenIndices = {};
		mIsPartOfLoop = false;
		
		
		if (mIdx == 0)
		{
			vertsInMetaEdge = std::vector<bool>(mSrcSkel->mNumVertices, false);
		}
		
		
		int vPrev = startVertex;
		int vNext = secondVertex;
		int alternateVert2 = -1;
		bool secondVertValid = validateVertexPair(startVertex, secondVertex, mSrcSkel, alternateVert2);
		if(secondVertValid)
		{
			vNext = secondVertex;
		}
		else if(alternateVert2 >= 0)
		{
			vNext = alternateVert2;
		}
		else
		{
			//std::cout << "The supplied vertex pair cannot be built into a MetaEdge, and there are no alternative vertices ";
			//std::cout << "neighboring the starting vertex which are not already part of a MetaEdge" << std::endl;
			return;
		}
		
		mVertices.push_back(vPrev);

		while (mSrcSkel->mNeighbors[vNext].size() == 2)
		{
			mVertices.push_back(vNext);

			for (int i = 0; i < 2; ++i)
			{
				if (mSrcSkel->mNeighbors[vNext][i] != vPrev)
				{
					vPrev = vNext;
					vNext = mSrcSkel->mNeighbors[vNext][i];
				}
			}
		}

		mVertices.push_back(vNext);

		for each(int vertId in mVertices)
		{
			vertsInMetaEdge[vertId] = true;
		}

		for(int i = 0; i < mVertices.size(); ++i)
		{
			SkeletonEdge* memberEdge = mSrcSkel->GetEdge(mVertices[i], mVertices[i+1]);
			if(memberEdge != nullptr)
			{
				mEdges.push_back(memberEdge);
			}
		}

		mChildNode = MetaNodeFactory::findOrCreateNode(mSrcSkel, vNext);
		mChildNodeIndex = mChildNode->mIdx;

		mChildNode->mNeighbors.push_back(mParentNode->mIdx);
		mParentNode->mNeighbors.push_back(mChildNode->mIdx);
	}

	MetaEdge::MetaEdge(const MetaEdge &toCopy)
	{
		mOrder = toCopy.mOrder;
		mIdx = toCopy.mIdx;
		mVertices = toCopy.mVertices;
		mSrcSkel = toCopy.mSrcSkel;
	}

	MetaEdge::MetaEdge(Skeleton *srcSkel, int aIdx, int order, bool aIsPartOfLoop, int parentSkelVertex, int childSkelVertex, std::vector<int> skelVertices, std::vector<MetaEdge*> aChildren)
	{
		mSrcSkel = srcSkel;
		mIdx = aIdx;
		mOrder = order;
		mIsPartOfLoop = aIsPartOfLoop;
		mParentNode = MetaNodeFactory::findOrCreateNode(mSrcSkel, parentSkelVertex);
		mChildNode = MetaNodeFactory::findOrCreateNode(mSrcSkel, childSkelVertex);

		mVertices = skelVertices;
		if (vertsInMetaEdge.size() < mSrcSkel->mNumVertices)
		{
			vertsInMetaEdge.resize(mSrcSkel->mNumVertices);
		}
		mEdges = {};
		for (int i = 0; i < mVertices.size() - 1; ++i)
		{
			int vert1 = mVertices[i];
			int vert2 = mVertices[i + 1];
			mEdges.push_back(mSrcSkel->GetEdge(vert1, vert2));
		}
		for each(int vert in mVertices)
		{
			vertsInMetaEdge[vert] = true;
		}

		mChildren = aChildren;
		mChildrenIndices = {};
		for each(MetaEdge *child in mChildren)
		{
			mChildrenIndices.push_back(child->GetId());
			child->SetParent(this);
		}

	}

	bool MetaEdge::validateVertexPair(int vert1, int vert2, Skeleton *srcSkel, int &alternateVert2)
	{
		bool isValid = false;
		if(vert1 != vert2)
		{
			bool isNeighbor = false;
			for each(int neighbor in srcSkel->mNeighbors[vert1])
			{
				if(neighbor == vert2)
				{
					isNeighbor = true;
				}
			}
			if(!isNeighbor)
			{
				isValid = false;
			}
		}
		else
		{
			isValid = false;
		}

		alternateVert2 = -1;
		if(!isValid)
		{
			for each(int neighbor in srcSkel->mNeighbors[vert1])
			{
			
				if(!vertsInMetaEdge[neighbor])
				{
					alternateVert2 = neighbor;
				}
			}
		}
		return isValid;
	}




	//MetaEdge* MetaEdge::recursiveLoad(std::vector<std::string> lines, int lineToLoad, Skeleton* srcSkeleton, bool isTopDown)
	//{
	//	std::string myLine = lines[lineToLoad];
	//	std::vector<std::string> data = {};
	//	boost::split(data, myLine, boost::is_any_of(" "));

	//	int dataI = 0;

	//	int idx = boost::lexical_cast<int>(data[dataI]); ++dataI;
	//	
	//	int order = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	bool isPartOfLoop = stringToBool(data[dataI]); ++dataI;

	//	int parentEdgeId = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	if (parentEdgeId >= 0)
	//	{
	//		return recursiveLoad(lines, parentEdgeId, srcSkeleton);
	//	}
	//	
	//	int parentNodeId = boost::lexical_cast<int>(data[dataI]); ++dataI;
	//	
	//	int parentNodeSrcVert = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	int childNodeId = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	int childNodeSrcVert = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	int numSkelVerts = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	std::vector<int> skelVerts(numSkelVerts);

	//	for (int skelI = 0; skelI < numSkelVerts; ++skelI, ++dataI)
	//	{
	//		skelVerts[skelI] = boost::lexical_cast<int>(data[dataI]);
	//	}

	//	int numChildren = boost::lexical_cast<int>(data[dataI]); ++dataI;

	//	std::vector<MetaEdge*> children = {};

	//	for (int childI = 0; childI < numChildren; ++childI, ++dataI)
	//	{
	//		int childLine = boost::lexical_cast<int>(data[dataI]);
	//		children.push_back(MetaEdge::recursiveLoad(lines, childLine, srcSkeleton, true));
	//	}

	//	return new MetaEdge(srcSkeleton, idx, order, isPartOfLoop, parentNodeSrcVert, childNodeSrcVert, skelVerts, children);
	//}

	std::ostream& operator<<(std::ostream& out, const MetaEdge& edge)
	{
		out << edge.mIdx << " " << edge.mOrder << " " << boolToString(edge.mIsPartOfLoop) << " ";

		//we need to know if this is a top level edge or not when loading, 
		//so provide the index of the parent if it is not a top level edge
		//or -1 if it is a top level edge
		if (edge.mOrder == 0)
		{
			out << -1 << " ";
		}
		else
		{
			out << edge.mParent->GetId() << " ";
		}
		out << edge.mParentNode->mIdx << " " << edge.mParentNode->mSrcVert << " ";
		out << edge.mChildNode->mIdx << " " << edge.mChildNode->mSrcVert << " ";
		out << edge.mVertices.size() << " ";
		for (int i = 0; i < edge.mVertices.size(); ++i)
		{
			out << edge.mVertices[i] << " ";
		}
		out << edge.mChildren.size() << " ";
		for (int i = 0; i < edge.mChildren.size(); ++i)
		{
			out << edge.mChildren[i]->GetId() << " ";
		}
		out << std::endl;
		return out;
	}


	std::vector<int> MetaEdge::GetVertices() { return mVertices; }
	std::vector<SkeletonEdge*> MetaEdge::GetEdges() { return mEdges; }
	Skeleton* MetaEdge::GetSkeleton() { return mSrcSkel; }
	MetaNode* MetaEdge::GetParentNode() { return mParentNode; }
	MetaNode* MetaEdge::GetChildNode() { return mChildNode; }
	std::vector<MetaEdge*> MetaEdge::GetChildren() { return mChildren; }
	MetaEdge* MetaEdge::GetParent() { return mParent; }
	int MetaEdge::GetId() { return mIdx; }
	int MetaEdge::GetOrder() { return mOrder; }
	bool MetaEdge::isInLoop() { return mIsPartOfLoop; }

	void MetaEdge::SetParent(MetaEdge *parent) { mParent = parent; return; }


	MetaEdge* MetaEdgeFactory::CreateEdge()
	{
		if (maxId < 0)
		{
			maxId = 0;
		}
		else
		{
			++maxId;
		}
		
		allEdges.push_back(MetaEdge(maxId));
		return &allEdges[maxId];
	}
	MetaEdge* MetaEdgeFactory::CreateEdge(Skeleton *srcSkel, int startVertex, int secondVertex, MetaEdge *parent)
	{
		if (maxId < 0)
		{
			maxId = 0;
		}
		else
		{
			++maxId;
		}
		allEdges.push_back(MetaEdge(srcSkel, maxId, startVertex, secondVertex, parent));
		buildDescendants(&allEdges[maxId]);
		return &allEdges[maxId];
	}

	MetaEdge* MetaEdgeFactory::CreateEdge(Skeleton *srcSkel, int aIdx, int order, bool aIsPartOfLoop, int parentSkelVertex, int childSkelVertex, std::vector<int> skelVertices, std::vector<MetaEdge*> aChildren)
	{
		if (aIdx > maxId)
		{
			maxId = aIdx;
			allEdges.resize(maxId + 1, MetaEdge());
		}
		allEdges[aIdx] = MetaEdge(srcSkel, aIdx, order, aIsPartOfLoop, parentSkelVertex, childSkelVertex, skelVertices, aChildren);
		return &allEdges[aIdx];
	}

	MetaEdge* MetaEdgeFactory::recursiveLoadEdges(std::vector<std::string> lines, int lineToLoad, Skeleton* srcSkeleton, bool isTopDown)
	{
		std::string myLine = lines[lineToLoad];
		std::vector<std::string> data = {};
		boost::split(data, myLine, boost::is_any_of(" "));

		int dataI = 0;
		int idx = boost::lexical_cast<int>(data[dataI]); ++dataI;
		int order = boost::lexical_cast<int>(data[dataI]); ++dataI;
		bool isPartOfLoop = stringToBool(data[dataI]); ++dataI;
		int parentEdgeId = boost::lexical_cast<int>(data[dataI]); ++dataI;

		if (parentEdgeId >= 0)
		{
			return recursiveLoadEdges(lines, parentEdgeId, srcSkeleton);
		}

		int parentNodeId = boost::lexical_cast<int>(data[dataI]); ++dataI;
		int parentNodeSrcVert = boost::lexical_cast<int>(data[dataI]); ++dataI;
		int childNodeId = boost::lexical_cast<int>(data[dataI]); ++dataI;
		int childNodeSrcVert = boost::lexical_cast<int>(data[dataI]); ++dataI;
		int numSkelVerts = boost::lexical_cast<int>(data[dataI]); ++dataI;

		std::vector<int> skelVerts(numSkelVerts);
		for (int skelI = 0; skelI < numSkelVerts; ++skelI, ++dataI)
		{
			skelVerts[skelI] = boost::lexical_cast<int>(data[dataI]);
		}

		int numChildren = boost::lexical_cast<int>(data[dataI]); ++dataI;
		std::vector<MetaEdge*> children = {};
		for (int childI = 0; childI < numChildren; ++childI, ++dataI)
		{
			int childLine = boost::lexical_cast<int>(data[dataI]);
			children.push_back(MetaEdgeFactory::recursiveLoadEdges(lines, childLine, srcSkeleton, true));
		}

		MetaEdge *result = CreateEdge(srcSkeleton, idx, order, isPartOfLoop, parentNodeSrcVert, childNodeSrcVert, skelVerts, children);
		for each(MetaEdge * child in children)
		{
			child->SetParent(result);
		}

		return result;
	}



	void MetaEdgeFactory::buildDescendants(MetaEdge *parentEdge)
	{
		buildDescendantsAtNode(parentEdge, parentEdge->mChildNode->mSrcVert);
		//if this is the primary root, then we also need to build out children from its 'parent' node.
		if (parentEdge->mOrder == 0)
		{
			buildDescendantsAtNode(parentEdge, parentEdge->mParentNode->mSrcVert);
		}
	}

	void MetaEdgeFactory::buildDescendantsAtNode(MetaEdge *parentEdge, int vertIndex)
	{
		for each(int neighborVert in parentEdge->mSrcSkel->mNeighbors[vertIndex])
		{
			if (!MetaEdge::vertsInMetaEdge[neighborVert])
			{
				MetaEdge *child = CreateEdge(parentEdge->GetSkeleton(), parentEdge->mParentNode->mSrcVert, neighborVert, parentEdge);
				buildDescendants(child);
				parentEdge->mChildren.push_back(child);
				child->SetParent(parentEdge);
			}
			else
			{
				//if the neighbor is in a metaEdge already, and that neighbor is neither the second or second last
				//vertex on this edge, then there must be a loop.  Set the is loop indicator.
				if (neighborVert != parentEdge->mVertices[1] && neighborVert != parentEdge->mVertices[parentEdge->mVertices.size() - 2])
				{
					//std::cout << "Loop Found" << std::endl;
					parentEdge->mIsPartOfLoop = true;
				}
			}
		}
	}

}
