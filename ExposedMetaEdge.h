#include "MetaGraph.h"
#include "boost/python.hpp"
#include "boost\python\suite\indexing\map_indexing_suite.hpp"
#include "boost\python\suite\indexing\vector_indexing_suite.hpp"

#include "BoostMetaGraph.h"

//using namespace boost::python;




/////////////////////EXPOSURE CLASSES//////////////////////
//
//struct TestClass
//{
//	boost::python::list GetTestVecD();
//	void SetTestVecD(boost::python::list& inList);
//
//private:
//	std::vector<float> myDVec;
//};
//
//struct PyRootAttributes
//{
//	PyRootAttributes(Roots::RootAttributes *srcAttributes);
//	PyRootAttributes();
//
//	//reload the public member from the source attributes
//	void reload();
//
//	//list of floats describing the root attributes
//	boost::python::list data;
//private:
//	//the source attribute object
//	Roots::RootAttributes *mSrcAttributes;
//};
//
//struct PyPoint3d
//{
//	PyPoint3d(Point3d p);
//	PyPoint3d();
//	float x, y, z;
//};
//
//
//
//struct PySkeleton
//{
//	PySkeleton();
//
//	PySkeleton(Roots::Skeleton *srcSkeleton);
//	
//	PySkeleton(std::string skelFile);
//
//	//reload the public members from the source skeleton
//	void reload();
//
//	void MoveCenterTo(float cx, float cy, float cz);
//	void ResetCenter();
//
//	//list of pypoint3d vertices
//	boost::python::list mVertices;
//	boost::python::list mVertArray;
//	//list of integer pairs eg [p11, p12, p21, p22, p31, p32, p41, p42]
//	//corresponding to vertex indices in above list
//	boost::python::list mEdges;
//	//list of edge attributes eg [a1, a2, a3] corresponding to each pair in the above list
//	boost::python::list mAttributes;
//	
//	float mSphereR;
//	PyPoint3d mSphereCenter;
//
//private:
//	//the source skeleton
//	Roots::Skeleton mSrcSkeleton;
//	std::vector<float> originalCenter;
//};
//
//struct PyMetaNode
//{
//	PyMetaNode(Roots::MetaNode srcNode);
//	PyMetaNode();
//
//	int mVertex;
//	int mIdx;
//
//	//list of neighbors by metanode id
//	boost::python::list mNeighbors;
//};
//
//struct PyMetaEdge
//{
//	PyMetaEdge(Roots::MetaEdge *srcEdge);
//	PyMetaEdge();
//
//	void reload();
//
//	//hierarchy order of this edge
//	int mOrder;
//	//index of this meta edge in the pymetagraph list
//	int mIndex;
//	//index of the parent edge in the pymetagraph list
//	int mParentIndex;
//	
//	//index of the parent and child nodes in the pymetagraph list
//	int mParentNodeIndex, mChildNodeIndex;
//
//	//indices of the child edges in the pymetagraph list
//	boost::python::list mChildIndices;
//	//list of member vertices from the pyskeleton in the metagraph by index
//	boost::python::list mVertIndices;
//	//list containing copies of edge [p1v1 p1v2 p2v1 p2v2 ....] redundant with the skeleton
//	boost::python::list mEdges;
//	//list containing copies of edge attributes [a1, a2, a3 ...] redundant with the skeleton
//	boost::python::list mAttributes;
//
//private:
//	//the source MetaEdge
//	Roots::MetaEdge *mSrcEdge;
//};
//
//struct PyMetaGraph
//{
//	PyMetaGraph(std::string srcFile);
//	PyMetaGraph(Roots::MetaGraph *srcGraph);
//
//	void reload();
//
//	//public member for the src skeleton
//	PySkeleton mPySkel;
//
//	//1-d list of PyMetaEdges in the graph
//	boost::python::list mPyEdges;
//
//	//1d list of PyMetaNodes in the graph
//	boost::python::list mNodes;
//
//private:
//	//the source metagraph
//	Roots::MetaGraph mSrcGraph;
//};

//
///////////////////////WRAPPER CODE//////////////////////////
//


//BOOST_PYTHON_MODULE(RootsTool)
//{
//
//	class_<TestClass>("TestClass")
//		.def("GetTestVecD", &TestClass::GetTestVecD)
//		.def("SetTestVecD", &TestClass::SetTestVecD)
//		;
//	
//	class_<Point3d>("Point3d")
//		.def_readwrite("x", &Point3d::x)
//		.def_readwrite("y", &Point3d::y)
//		.def_readwrite("z", &Point3d::z)
//		.def_readwrite("id", &Point3d::id)
//		.def("mag", &Point3d::mag)
//		;
//
//	class_<PyRootAttributes>("PyRootAttributes")
//		.def_readwrite("data", &PyRootAttributes::data)
//		.def("reload", &PyRootAttributes::reload)
//		;
//
	//class_<PySkeleton>("PySkeleton", init<std::string>())
//		.def_readwrite("vertices", &PySkeleton::mVertices)
//		.def_readwrite("edges", &PySkeleton::mEdges)
//		.def_readwrite("attributes", &PySkeleton::mAttributes)
//		.def_readwrite("sphereCenter", &PySkeleton::mSphereCenter)
//		.def_readwrite("sphereR", &PySkeleton::mSphereR)
//		.def("reload", &PySkeleton::reload)
//		.def("moveCenterTo", &PySkeleton::MoveCenterTo)
//		.def("resetCenter", &PySkeleton::ResetCenter)
//		;
//
//	class_<PyMetaNode>("PyMetaNode")
//		.def_readwrite("vertex", &PyMetaNode::mVertex)
//		.def_readwrite("idx", &PyMetaNode::mIdx)
//		.def_readwrite("neighbors", &PyMetaNode::mNeighbors)
//		;
//
//	class_<PyMetaEdge>("PyMetaEdge")
//		.def_readwrite("order", &PyMetaEdge::mOrder)
//		.def_readwrite("idx", &PyMetaEdge::mIndex)
//		.def_readwrite("pIdx", &PyMetaEdge::mParentIndex)
//		.def_readwrite("pNodeIndex", &PyMetaEdge::mParentNodeIndex)
//		.def_readwrite("cNodeIndex", &PyMetaEdge::mChildNodeIndex)
//		.def_readwrite("children", &PyMetaEdge::mChildIndices)
//		.def_readwrite("vertices", &PyMetaEdge::mVertIndices)
//		.def_readwrite("edges", &PyMetaEdge::mEdges)
//		.def_readwrite("attributes", &PyMetaEdge::mAttributes)
//		.def("reload", &PyMetaEdge::reload)
//		;
//
//	class_<PyMetaGraph>("PyMetaGraph", init<std::string>())
//		.def_readwrite("skeleton", &PyMetaGraph::mPySkel)
//		.def_readwrite("meta_edges", &PyMetaGraph::mPyEdges)
//		.def_readwrite("nodes", &PyMetaGraph::mNodes)
//		.def("reload", &PyMetaGraph::reload)
//		;
//};