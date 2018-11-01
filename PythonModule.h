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
		.def_readwrite("v0id", &Roots::RootAttributes::v0id)
		.def_readwrite("v1id", &Roots::RootAttributes::v1id)
		.def_readwrite("euclidLength", &Roots::RootAttributes::euclidLength)
		;


	class_<Roots::BSkeleton>("Skeleton")
		.def("findBoundingSphere", &BSkeleton::findBoundingSphere)
		.def_readwrite("center", &BSkeleton::mCenter)
		.def_readwrite("radius", &BSkeleton::mRadius)
		;



	class_<IssuesGL>("IssuesGL")
		.def("issuegl", &IssuesGL::issueGL)
		;


	class_<VBOSphere>("VBOSphere")
		.def("draw", &VBOSphere::draw)
		.def("fancyDraw", &VBOSphere::fancyDraw)
		;

	class_<BMetaGraph>("mgraph")
		.def("loadFromFile", &BMetaGraph::loadFromFile)
		.def("loadMeshFromFile", &BMetaGraph::loadMeshFromFile)
		.def("saveToFile", &BMetaGraph::saveToFile)
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
		.def("showMesh", &BMetaGraph::showMesh)
		.def("setMeshAlpha", &BMetaGraph::setMeshAlpha)
		.def("setMeshColor", &BMetaGraph::setMeshColor)
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
		.def("shiftEye", &BMetaGraph::shiftEye)
		.def("pickNewViewCenter", &BMetaGraph::pickNewViewCenter)
		.def("changeRotationSpeed", &BMetaGraph::changeRotationSpeed)
		.def("getComponentSizes", &BMetaGraph::getComponentSizes)
		.def("getNumEdgesToBreak", &BMetaGraph::GetNumEdgesToBreak)
		.def("getNumComponents", &BMetaGraph::getNumComponents)
		.def_readwrite("skeleton", &BMetaGraph::mSkeleton)
		;

};