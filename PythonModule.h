#pragma once

#include "boost/python.hpp"
#include "boost\python\suite\indexing\map_indexing_suite.hpp"
#include "boost\python\suite\indexing\vector_indexing_suite.hpp"

#include "BoostMetaGraph.h"
#include "IssuesGL.h"
#include "Sphere.h"

using namespace boost::python;
using namespace drawing;
using namespace Roots;
BOOST_PYTHON_MODULE(RootsTool)
{


	class_<Point3d>("Point3d")
		.def("x", &Point3d::x)
		.def("y", &Point3d::y)
		.def("z", &Point3d::z)
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

	//class_<PySkeleton>("Skeleton")
	//	.def_readwrite("vertices", &PySkeleton::mVertexList)
	//	.def_readwrite("edges", &PySkeleton::mEdgeList)
	//	.def_readwrite("center", &PySkeleton::center)
	//	.def_readwrite("radius", &PySkeleton::radius)
	//	.def_readwrite("minThickness", &PySkeleton::minThickness)
	//	.def_readwrite("maxThickness", &PySkeleton::maxThickness)
	//	.def_readwrite("minWidth", &PySkeleton::minWidth)
	//	.def_readwrite("maxWidth", &PySkeleton::maxWidth)
	//	.def_readwrite("minRatio", &PySkeleton::minRatio)
	//	.def_readwrite("maxRatio", &PySkeleton::maxRatio)
	//	.def_readwrite("thicknessPercentiles", &PySkeleton::thicknessPercentiles)
	//	.def_readwrite("widthPercentiles", &PySkeleton::widthPercentiles)
	//	.def_readwrite("ratioPercentiles", &PySkeleton::ratioPercentiles)
	//	.def("findBoundingSphere", &PySkeleton::findBoundingSphere)
	//	.def("reload", &PySkeleton::reload)
	//	;

	class_<Roots::BSkeleton>("Skeleton")
		.def("findBoundingSphere", &BSkeleton::findBoundingSphere)
		.def_readwrite("center", &BSkeleton::mCenter)
		.def_readwrite("radius", &BSkeleton::mRadius)
		;
	/*class_<Roots::MetaNode3d>("MetaNode3d")
		.def("x", &Roots::MetaNode3d::x)
		.def("y", &Roots::MetaNode3d::y)
		.def("z", &Roots::MetaNode3d::z)
		.def_readwrite("id", &Roots::MetaNode3d::id)
		.def_readwrite("order", &Roots::MetaNode3d::order)
		.def_readwrite("degree", &Roots::MetaNode3d::degree)
		.def_readwrite("component", &Roots::MetaNode3d::connectedComponent)
		.def_readwrite("thickness", &Roots::MetaNode3d::nodeThickness)
		.def_readwrite("width", &Roots::MetaNode3d::nodeWidth)
		;

	class_<Roots::MetaEdge3d>("MetaEdge3d")
		.def_readwrite("node0", &Roots::MetaEdge3d::node0)
		.def_readwrite("node1", &Roots::MetaEdge3d::node1)
		.def_readwrite("thickness", &Roots::MetaEdge3d::avgThickness)
		.def_readwrite("width", &Roots::MetaEdge3d::avgWidth)
		.def_readwrite("component", &Roots::MetaEdge3d::connectedComponent)
		.def_readwrite("order", &Roots::MetaEdge3d::order)
		.def_readwrite("edges", &Roots::MetaEdge3d::edges)
		.def_readwrite("isBridge", &Roots::MetaEdge3d::isBridge)
		;

	class_<PyMetaGraph>("MetaGraph", init<std::string>())
		.def_readwrite("skeleton", &PyMetaGraph::mSkeleton)
		.def_readwrite("srcGraph", &PyMetaGraph::mGraph)
		.def_readwrite("nodeLocations", &PyMetaGraph::mMetaNodeLocations)
		.def_readwrite("edgeConnections", &PyMetaGraph::mMetaEdgeConnections)
		.def_readwrite("componentNodeMap", &PyMetaGraph::componentNodeMap)
		.def_readwrite("componentEdgeMap", &PyMetaGraph::componentEdgeMap)
		.def_readwrite("componentNodesOfInterestMap", &PyMetaGraph::componentNodesOfInterestMap)
		.def_readwrite("numEdgesToBreak", &PyMetaGraph::numEdgesToBreak)
		.def("initializeFromSkeleton", &PyMetaGraph::initializeFromSkeleton)
		.def("labelComponents", &PyMetaGraph::labelComponents)
		.def("reload", &PyMetaGraph::reload)
		.def("joinOperation", &PyMetaGraph::joinOperation)
		.def("breakOperation", &PyMetaGraph::breakOperation)
		.def("splitOperation", &PyMetaGraph::splitOperation)
		;*/


	class_<IssuesGL>("IssuesGL")
		.def("issuegl", &IssuesGL::issueGL)
		;


	class_<VBOSphere>("VBOSphere")
		.def("draw", &VBOSphere::draw)
		.def("fancyDraw", &VBOSphere::fancyDraw)
		;

	class_<BMetaGraph>("mgraph")
		.def("loadFromFile", &BMetaGraph::loadFromFile)
		.def("showEdges", &BMetaGraph::showEdges)
		.def("joinOperation", &BMetaGraph::JoinOperation)
		.def("breakOperation", &BMetaGraph::BreakOperation)
		.def("splitOperation", &BMetaGraph::SplitOperation)
		.def("promoteOperation", &BMetaGraph::PromoteOperation)
		.def("assignEdgeHeatmap", &BMetaGraph::assignEdgeHeatMap)
		.def("showEdges", &BMetaGraph::showEdges)
		.def("setEdgeScale", &BMetaGraph::setEdgeScale)
		.def("magnifyNonBridges", &BMetaGraph::magnifyNonBridges)
		.def("showOnlyNonBridges", &BMetaGraph::showOnlyNonBridges)
		.def("colorizeEdgesByThickness", &BMetaGraph::colorizeEdgesByThickness)
		.def("colorizeEdgesByWidth", &BMetaGraph::colorizeEdgesByWidth)
		.def("colorizeEdgesByRatio", &BMetaGraph::colorizeEdgesByRatio)
		.def("colorizeEdgesByComponent", &BMetaGraph::colorizeEdgesByComponent)
		.def("setEdgeColorFloor", &BMetaGraph::setEdgeColorFloor)
		.def("setEdgeColorCeiling", &BMetaGraph::setEdgeColorCeiling)
		.def("setEdgeSelectionColor", &BMetaGraph::setEdgeSelectionColor)
		.def("assignNodeHeatmap", &BMetaGraph::assignNodeHeatMap)
		.def("showEndpoints", &BMetaGraph::showEndpoints)
		.def("showJunctions", &BMetaGraph::showJunctions)
		.def("setEndpointScale", &BMetaGraph::setEndpointScale)
		.def("setJunctionScale", &BMetaGraph::setJunctionScale)
		.def("setConstantNodeColor", &BMetaGraph::setConstantNodeColor)
		.def("colorizeNodesByThickness", &BMetaGraph::colorizeNodesByThickness)
		.def("colorizeNodesByWidth", &BMetaGraph::colorizeNodesByWidth)
		.def("colorizeNodesByDegree", &BMetaGraph::colorizeNodesByDegree)
		.def("colorizeNodesByComponent", &BMetaGraph::colorizeNodesByComponent)
		.def("colorizeNodesByConstantColor", &BMetaGraph::colorizeNodesByConstantColor)
		.def("setNodeColorFloor", &BMetaGraph::setNodeColorFloor)
		.def("setNodeColorCeiling", &BMetaGraph::setNodeColorCeiling)
		.def("setDisplayOnlySelectedComponents", &BMetaGraph::setDisplayOnlySelectedComponents)
		.def("setComponent1", &BMetaGraph::setComponent1)
		.def("setComponent2", &BMetaGraph::setComponent2)
		.def("setShowBoundingBoxes", &BMetaGraph::setShowBoundingBoxes)
		.def("drawEdges", &BMetaGraph::drawEdges)
		.def("drawNodes", &BMetaGraph::drawNodes)
		.def("drawBoxes", &BMetaGraph::drawBoxes)
		.def("draw", &BMetaGraph::draw)
		.def("startRotation", &BMetaGraph::startRotation)
		.def("mouseMoved", &BMetaGraph::mouseMoved)
		.def("setZoom", &BMetaGraph::setZoom)
		.def("selectConnectionNode", &BMetaGraph::selectConnectionNode)
		.def("selectBreakEdge", &BMetaGraph::selectBreakEdge)
		.def("selectSplitEdge", &BMetaGraph::selectSplitEdge)
		.def("setSelectionColor", &BMetaGraph::setSelectionColor)
		.def("unselectAll", &BMetaGraph::unselectAll)
		.def("getComponentSizes", &BMetaGraph::getComponentSizes)
		.def("getNumEdgesToBreak", &BMetaGraph::GetNumEdgesToBreak)
		.def("getNumComponents", &BMetaGraph::getNumComponents)
		.def_readwrite("skeleton", &BMetaGraph::mSkeleton)
		;

};