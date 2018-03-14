#pragma once

#include "boost/python.hpp"
#include "boost\python\suite\indexing\map_indexing_suite.hpp"
#include "boost\python\suite\indexing\vector_indexing_suite.hpp"

#include "BoostMetaGraph.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(RootsTool)
{


	class_<Point3d>("Point3d")
		.def_readwrite("x", &Point3d::x)
		.def_readwrite("y", &Point3d::y)
		.def_readwrite("z", &Point3d::z)
		.def_readwrite("id", &Point3d::id)
		.def("mag", &Point3d::mag)
		;

	class_<Roots::RootAttributes>("RootAttributes")
		.def_readwrite("thickness", &Roots::RootAttributes::thickness)
		.def_readwrite("width", &Roots::RootAttributes::width)
		.def_readwrite("length", &Roots::RootAttributes::length)
		.def_readwrite("v0id", &Roots::RootAttributes::v0id)
		.def_readwrite("v1id", &Roots::RootAttributes::v1id)
		.def_readwrite("euclidLength", &Roots::RootAttributes::euclidLength)
		;

	class_<PySkeleton>("Skeleton")
		.def_readwrite("vertices", &PySkeleton::mVertexList)
		.def_readwrite("edges", &PySkeleton::mEdgeList)
		.def_readwrite("center", &PySkeleton::center)
		.def_readwrite("radius", &PySkeleton::radius)
		.def("findBoundingSphere", &PySkeleton::findBoundingSphere)
		.def("reload", &PySkeleton::reload)
		;

	class_<Roots::MetaNode3d>("MetaNode3d")
		.def_readwrite("x", &Roots::MetaNode3d::x)
		.def_readwrite("y", &Roots::MetaNode3d::y)
		.def_readwrite("z", &Roots::MetaNode3d::z)
		.def_readwrite("id", &Roots::MetaNode3d::id)
		.def_readwrite("order", &Roots::MetaNode3d::order)
		.def_readwrite("component", &Roots::MetaNode3d::connectedComponent)
		.def_readwrite("size", &Roots::MetaNode3d::nodeSize)
		;

	class_<Roots::MetaEdge3d>("MetaEdge3d")
		.def_readwrite("node0", &Roots::MetaEdge3d::node0)
		.def_readwrite("node1", &Roots::MetaEdge3d::node1)
		.def_readwrite("thickness", &Roots::MetaEdge3d::avgThickness)
		.def_readwrite("component", &Roots::MetaEdge3d::connectedComponent)
		.def_readwrite("order", &Roots::MetaEdge3d::order)
		.def_readwrite("edges", &Roots::MetaEdge3d::edges)
		;

	class_<PyMetaGraph>("MetaGraph", init<std::string>())
		.def_readwrite("skeleton", &PyMetaGraph::mSkeleton)
		.def_readwrite("nodeLocations", &PyMetaGraph::mMetaNodeLocations)
		.def_readwrite("edgeConnections", &PyMetaGraph::mMetaEdgeConnections)
		.def_readwrite("componentNodeMap", &PyMetaGraph::componentNodeMap)
		.def_readwrite("componentEdgeMap", &PyMetaGraph::componentEdgeMap)
		.def_readwrite("componentNodesOfInterestMap", &PyMetaGraph::componentNodesOfInterestMap)
		.def("initializeFromSkeleton", &PyMetaGraph::initializeFromSkeleton)
		.def("labelComponents", &PyMetaGraph::labelComponents)
		.def("reload", &PyMetaGraph::reload)
		.def("joinOperation", &PyMetaGraph::joinOperation)
		.def("breakOperation", &PyMetaGraph::breakOperation)
		;
};