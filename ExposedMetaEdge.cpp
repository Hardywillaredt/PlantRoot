//#include "ExposedMetaEdge.h"
//
//#include <fstream>
//
//
//list TestClass::GetTestVecD()
//{
//	list result = list();
//	for (int i = 0; i < myDVec.size(); ++i)
//	{
//		result.append<float>(myDVec[i]);
//	}
//	return result;
//}
//void TestClass::SetTestVecD(boost::python::list& inList)
//{
//	myDVec = std::vector<float>(0);
//	for (int i = 0; i < len(inList); ++i)
//	{
//		myDVec.push_back(extract<float>(inList[i]));
//	}
//
//}
//
//
//PyRootAttributes::PyRootAttributes(Roots::RootAttributes *srcAttributes)
//{
//	data = list();
//	mSrcAttributes = srcAttributes;
//	reload();
//}
//
//PyRootAttributes::PyRootAttributes()
//{
//	data = list();
//	mSrcAttributes = nullptr;
//}
//
//
//void PyRootAttributes::reload()
//{
//	if (mSrcAttributes != nullptr)
//	{
//		data = list();
//		for (int i = 0; i < Roots::RootAttributes::NumAttributes; ++i)
//		{
//			data.append<float>(mSrcAttributes->operator[](i));
//		}
//	}
//}
//
//
//PyPoint3d::PyPoint3d(Point3d p)
//{
//	x = p.x;
//	y = p.y;
//	z = p.z;
//}
//
//PyPoint3d::PyPoint3d()
//{
//	x = 0;
//	y = 0;
//	z = 0;
//}
//
//PySkeleton::PySkeleton()
//{
//	mVertices = list();
//	mEdges = list();
//	mAttributes = list();
//}
//
//PySkeleton::PySkeleton(Roots::Skeleton *srcSkeleton)
//{
//	mSrcSkeleton = *srcSkeleton;
//	reload();
//}
//
//
//PySkeleton::PySkeleton(std::string skelFile)
//{
//	Log::WriteLine("loading pyskeleton from file");
//
//	mSrcSkeleton = Roots::Skeleton(skelFile);
//	reload();
//}
//
//void PySkeleton::reload()
//{
//
//	Roots::vertList verts = mSrcSkeleton.getVertices();
//	mVertices = list();
//	mVertArray = list();
//	for (int i = 0; i < verts.size(); ++i)
//	{
//		mVertices.append<PyPoint3d>(PyPoint3d(verts[i]));
//		mVertArray.append<float>(verts[i].x);
//		mVertArray.append<float>(verts[i].y);
//		mVertArray.append<float>(verts[i].z);
//
//
//	}
//
//	std::vector<Roots::edgeList> edges = mSrcSkeleton.getEdges();
//	mEdges = list();
//	mAttributes = list();
//	for (int i = 0; i < edges.size(); ++i)
//	{
//		for (int k = 0; k < edges[i].size(); ++k)
//		{
//			mEdges.append<int>(edges[i][k].v0);
//			mEdges.append<int>(edges[i][k].v1);
//			mAttributes.append<PyRootAttributes>(PyRootAttributes(&edges[i][k].attributes));
//		}
//	}
//
//	float cx, cy, cz, r;
//	mSrcSkeleton.GetBoundingSphere(cx, cy, cz, r);
//	mSphereCenter = PyPoint3d(Point3d(cx, cy, cz));
//	mSphereR = r;
//}
//
//void PySkeleton::MoveCenterTo(float cx, float cy, float cz)
//{
//	Point3d newCenter = Point3d(cx, cy, cz);
//	mSrcSkeleton.MoveCenterTo(newCenter);
//	reload();
//}
//
//void PySkeleton::ResetCenter()
//{
//	mSrcSkeleton.ResetCenter();
//	reload();
//}
//
//
//PyMetaNode::PyMetaNode(Roots::MetaNode srcNode)
//{
//	mIdx = srcNode.mIdx;
//	mVertex = srcNode.mSrcVert;
//	mNeighbors = list();
//	for (int i = 0; i < srcNode.mNeighbors.size(); ++i)
//	{
//		mNeighbors.append<int>(srcNode.mNeighbors[i]);
//	}
//}
//
//PyMetaNode::PyMetaNode()
//{
//	mIdx = -1;
//	mVertex = -1;
//	mNeighbors = list();
//}
//
//PyMetaEdge::PyMetaEdge(Roots::MetaEdge *srcMetaEdge)
//{
//	mSrcEdge = srcMetaEdge;
//	reload();
//}
//
//PyMetaEdge::PyMetaEdge()
//{
//	mSrcEdge = nullptr;
//	mOrder = -1;
//	mIndex = -1;
//	mParentIndex = -1;
//	mParentNodeIndex = -1;
//	mChildNodeIndex = -1;
//	mChildIndices = {};
//	mVertIndices = {};
//	mEdges = {};
//	mAttributes = {};
//}
//
//
//void PyMetaEdge::reload()
//{
//	if (mSrcEdge != nullptr)
//	{
//		mOrder = mSrcEdge->GetOrder();
//		mIndex = mSrcEdge->GetId();
//		if (mOrder > 0)
//		{
//			mParentIndex = mSrcEdge->GetParent()->GetId();
//		}
//		else
//		{
//			mParentIndex = -1;
//		}
//
//		Roots::MetaNode *parentNode = mSrcEdge->GetParentNode();
//		if (parentNode != nullptr)
//		{
//			mParentNodeIndex = parentNode->mIdx;
//		}
//		else
//		{
//			std::cerr << "Parent node not defined for metaedge " << std::endl;
//		}
//		
//
//		Roots::MetaNode *childNode = mSrcEdge->GetChildNode();
//		if (childNode != nullptr)
//		{
//			mChildNodeIndex = childNode->mIdx;
//		}
//		else
//		{
//			std::cerr << "Child node not defined for metaedge " << std::endl;
//		}
//		
//		mChildIndices = list();
//		for each(Roots::MetaEdge *child in mSrcEdge->GetChildren())
//		{
//			if (child != nullptr)
//			{
//				mChildIndices.append<int>(child->GetId());
//			}
//			else
//			{
//				std::cerr << "Child edge is not defined for metaedge " << std::endl;
//			}
//		}
//
//		mVertIndices = list();
//		std::vector<int> verts = mSrcEdge->GetVertices();
//		for (int i = 0; i < verts.size(); ++i)
//		{
//			mVertIndices.append<int>(verts[i]);
//		}
//
//		std::vector<Roots::SkeletonEdge*> edges = mSrcEdge->GetEdges();
//		mEdges = list();
//		mAttributes = list();
//		for (int i = 0; i < edges.size(); ++i)
//		{
//			mEdges.append<int>(edges[i]->v0);
//			mEdges.append<int>(edges[i]->v1);
//			mAttributes.append<PyRootAttributes>(PyRootAttributes(&edges[i]->attributes));
//		}
//	}
//}
//
//
//PyMetaGraph::PyMetaGraph(std::string srcFile)
//{
//	mSrcGraph = Roots::MetaGraph(srcFile);
//	reload();
//}
//
//PyMetaGraph::PyMetaGraph(Roots::MetaGraph *srcGraph)
//{
//	mSrcGraph = *srcGraph;
//	reload();
//}
//
//void PyMetaGraph::reload()
//{
//
//	mPySkel = PySkeleton(mSrcGraph.GetSkeleton());
//	std::vector<Roots::MetaEdge*> flatEdges = mSrcGraph.GetFlatEdges();
//	mPyEdges = list();
//	for (int i = 0; i < flatEdges.size(); ++i)
//	{
//		mPyEdges.append<PyMetaEdge>(PyMetaEdge(flatEdges[i]));
//	}
//
//	mNodes = list();
//
//		
//	for (int i = 0; i < Roots::MetaNodeFactory::nodes.size(); ++i)
//	{
//		mNodes.append<PyMetaNode>(PyMetaNode(Roots::MetaNodeFactory::nodes[i]));
//	}
//	
//}
//
