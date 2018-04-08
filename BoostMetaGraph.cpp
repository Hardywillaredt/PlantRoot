#include "BoostMetaGraph.h"
#include "boost/graph/connected_components.hpp"
#include "boost/graph/kruskal_min_spanning_tree.hpp"
#include <iostream>
#include <fstream>
namespace
{
	typedef boost::graph_traits<Roots::BMetaGraph>::vertex_iterator vertIter;
	typedef boost::graph_traits<Roots::BMetaGraph>::edge_iterator edgeIter;
	using namespace boost;
	struct metaVertIter : std::pair<vertIter, vertIter>
	{
		metaVertIter(const std::pair<vertIter, vertIter> &other)
			:std::pair<vertIter, vertIter>(other)
		{
		}
		metaVertIter operator++()
		{
			++this->first;
			return *this;
		}
		metaVertIter operator--()
		{
			--this->second;
			return *this;
		}
	};

	struct metaEdgeIter : std::pair<edgeIter, edgeIter>
	{
		metaEdgeIter(const std::pair<edgeIter, edgeIter> &other)
			:std::pair<edgeIter, edgeIter>(other)
		{
		}
		metaEdgeIter operator++()
		{
			++this->first;
			return *this;
		}

		metaEdgeIter operator--()
		{
			--this->second;
			return *this;
		}
	};

	template <class T>
	boost::python::list toPythonList(std::vector<T> vector) {
		typename std::vector<T>::iterator iter;
		boost::python::list list;
		for (iter = vector.begin(); iter != vector.end(); ++iter) {
			list.append(*iter);
		}
		return list;
	}

	template <class T>
	std::vector<T> toStdVector(boost::python::list list)
	{
		std::vector<T> result = {};
		for (int i = 0; i < boost::python::len(list); ++i)
		{
			result.push_back(boost::python::extract<T>(list[i]));
		}
		return result;
	}
}
namespace Roots
{
	/////////////////////////////////////////////// BMetaNode /////////////////////////////////////
	BMetaNode::BMetaNode()
	{
		mSrcVert = 0;
		mSrcSkeleton = nullptr;
		connectedComponent = -1;
		nodeThickness = 0.0;
		nodeWidth = 0.0;
		hasGeom = false;
		x = 0;
		y = 0;
		z = 0;
	}
	BMetaNode::BMetaNode(SkelVert srcId, BSkeleton *skel)
	{
		mSrcVert = srcId;
		mSrcSkeleton = skel;
		//std::cout << "number of vertices in pointed skeleton " << boost::num_vertices(*skel) << std::endl;
		connectedComponent = -1;
		nodeThickness = 0.0;
		nodeWidth = 0.0;
		hasGeom = true;
		x = skel->operator[](srcId).x;
		y = skel->operator[](srcId).y;
		z = skel->operator[](srcId).z;
	}

	////////////////////////////////////////////// BMetaEdge ///////////////////////////////////////

	BMetaEdge::BMetaEdge()
	{
		mVertices = {};
		mEdges = {};
		averageThickness = 0.0;
		averageWidth = 0.0;
		isBridge = false;
		//mSrcSkeleton = nullptr;

	}

	BMetaEdge::BMetaEdge(std::vector<SkelVert> vertices, BSkeleton* srcSkeleton, std::vector<GLfloat> &edgeVertices, std::vector<GLVertexIndex> &edgeIndices)
	{
		mVertices = vertices;
		mEdges = {};
		float weightedAvgThickness = 0.0;
		float weightedAvgWidth = 0.0;
		float skelEdgeLength = 0.0;
		averageThickness = 0.0;
		averageWidth = 0.0;
		isBridge = false;
		for (int i = 0; i < mVertices.size() - 1; ++i)
		{
			SkelVert v0 = mVertices[i];
			SkelVert v1 = mVertices[i + 1];
			SkelEdge e;
			bool exists;
			boost::tie(e, exists) = boost::edge(v0, v1, *srcSkeleton);
			if (exists)
			{
				RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
				mEdges.push_back(skelRootAttributes);
				float thickness = skelRootAttributes[Thickness];
				float length = skelRootAttributes.euclidLength;
				weightedAvgThickness += thickness * length;
				weightedAvgWidth += skelRootAttributes[Width] * length;
				skelEdgeLength += length;
			}
			else
			{
				//this is okay - this just means that this edge does not exist at the skeleton
				//level, eg it is an added edge at the metagraph level
				//do nothing here
				float cumThickness = 0.0;
				float cumWidth = 0.0;
				int numEdges = 0;

				BSkeleton::adjacency_iterator adjIt, adjEnd;
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *srcSkeleton);
				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelVert leadVert = *adjIt;
					boost::tie(e, exists) = boost::edge(v0, leadVert, *srcSkeleton);
					RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
					cumThickness += skelRootAttributes[Thickness];
					cumWidth += skelRootAttributes[Width];
					++numEdges;
				}
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *srcSkeleton);
				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelVert leadVert = *adjIt;
					boost::tie(e, exists) = boost::edge(v1, leadVert, *srcSkeleton);
					RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
					cumThickness += skelRootAttributes[Thickness];
					cumWidth += skelRootAttributes[Width];
					++numEdges;
				}
				if (numEdges == 0)
				{
					numEdges = 1;
				}
				weightedAvgThickness += cumThickness / numEdges;
				weightedAvgThickness += cumWidth / numEdges;
			}
		}
		if (skelEdgeLength == 0)
		{
			skelEdgeLength = 1;
		}
		averageThickness = weightedAvgThickness / skelEdgeLength;
		
		if (isnan(averageThickness))
		{
			averageThickness = 1.0;
			std::cout << "Note : improper MetaEdge thickness, skeleton data may be missing thickness data on one or more roots" << std::endl;
		}

		averageWidth = weightedAvgWidth / skelEdgeLength;

		if (isnan(averageWidth))
		{
			averageWidth = 1.0;
			std::cout << "Note : improper MetaEdge width, skeleton data may be missing width data on one or more roots" << std::endl;
		}

		glEdgeVerticesStart = edgeVertices.size() / 3;
		glEdgeLinesStart = edgeIndices.size() / 2;

		for (SkelVert skelvert : mVertices)
		{
			Point3d p = srcSkeleton->operator[](skelvert);
			edgeVertices.push_back(p.x);
			edgeVertices.push_back(p.y);
			edgeVertices.push_back(p.z);
		}

		//mSrcSkeleton = srcSkeleton;

	}

	BMetaEdge BMetaEdge::join(BMetaEdge &other, BSkeleton *srcSkeleton)
	{
		int myMaxId = mVertices.size() - 1;
		int otherMaxId = other.mVertices.size() - 1;
		std::vector<SkelVert> joinedVertices = {};
		if (mVertices[myMaxId] != other.mVertices[0] &&
			mVertices[myMaxId] != other.mVertices[otherMaxId] &&
			mVertices[0] != other.mVertices[0] &&
			mVertices[0] != other.mVertices[otherMaxId])
		{
			//std::cout << "Meta edges are not successfully joined.  No endpoints match between either edge." << std::endl;
			return BMetaEdge();
		}
		//in the case that the last vertex in this edge is the starting vertex in the 
		//other edge, simply append all but the first vertex of that edge to the vertex list
		//of this
		if (mVertices[myMaxId] == other.mVertices[0])
		{
			for (int i = 0; i < mVertices.size(); ++i)
			{
				joinedVertices.push_back(mVertices[i]);
			}
			for (int i = 1; i < other.mVertices.size(); ++i)
			{
				joinedVertices.push_back(other.mVertices[i]);
			}
		}
		else if (mVertices[myMaxId] == other.mVertices[otherMaxId])
		{
			for (int i = 0; i < mVertices.size(); ++i)
			{
				joinedVertices.push_back(mVertices[i]);
			}
			for (int i = otherMaxId - 1; i >= 0; --i)
			{
				joinedVertices.push_back(other.mVertices[i]);
			}
		}
		else
		{
			//to avoid repetitive code, simply flip the of the ordering for the same result
			//if matches are not made on this end.
			return other.join(*this, srcSkeleton);
		}

		return BMetaEdge(joinedVertices, srcSkeleton);
	}

	float BMetaEdge::getAvgThickness()
	{
		return averageThickness;
		//float weightedAvgThickness = 0.0;
		//float skelEdgeLength = 0.0;
		//for (int i = 0; i < mVertices.size() - 1; ++i)
		//{
		//	SkelVert v0 = mVertices[i];
		//	SkelVert v1 = mVertices[i + 1];
		//	SkelEdge e;
		//	bool exists;
		//	boost::tie(e, exists) = boost::edge(v0, v1, *mSrcSkeleton);
		//	if (exists)
		//	{
		//		RootAttributes skelRootAttributes = mSrcSkeleton->operator[](e);
		//		float thickness = skelRootAttributes[RootAttributes::Thickness];
		//		float length = skelRootAttributes.euclidLength;
		//		weightedAvgThickness += thickness * length;
		//		skelEdgeLength += length;
		//	}
		//	else
		//	{
		//		//this is okay - this just means that this edge does not exist at the skeleton
		//		//level, eg it is an added edge at the metagraph level
		//		//do nothing here
		//	}
		//}
		//return weightedAvgThickness / skelEdgeLength;
	}

	SkelVert BMetaEdge::start()
	{
		return mVertices[0];
	}
	SkelVert BMetaEdge::end()
	{
		return mVertices[mVertices.size() - 1];
	}

	///////////////////////////////////////// MetaNode3d //////////////////////////////////////////

	MetaNode3d::MetaNode3d()
	{
		order = -1;
		connectedComponent = -1;
		nodeThickness = -1;
		nodeWidth = -1;
		degree = 0;
	}

	MetaNode3d::MetaNode3d(BMetaNode src, Roots::BSkeleton *srcSkeleton, int aOrder, int aDegree)
		:Point3d(srcSkeleton->operator[](src.mSrcVert))
	{
		if (src.hasGeom)
		{
			x = src.x;
			y = src.y;
			z = src.z;
		}
		order = aOrder;
		connectedComponent = src.connectedComponent;
		nodeThickness = src.nodeThickness;
		nodeWidth = src.nodeWidth;
		degree = aDegree;
	}



	///////////////////////////////////////// MetaEdge3d //////////////////////////////////////////

	MetaEdge3d::MetaEdge3d()
	{
		node0 = 0;
		node1 = 0;
		avgThickness = 0;
		avgWidth = 0;
		connectedComponent = -1;
		isBridge = false;
		//edges = boost::python::list();

	}

	MetaEdge3d::MetaEdge3d(int n0, int n1, float thickness, float width, int component, int aOrder, std::vector<RootAttributes> aEdges, bool aIsBridge)
	{
		node0 = n0;
		node1 = n1;
		avgThickness = thickness;
		avgWidth = width;
		connectedComponent = component;
		order = aOrder;
		edges = toPythonList<RootAttributes>(aEdges);
		isBridge = aIsBridge;
	}


	///////////////////////////////////////// BMetaGraph //////////////////////////////////////////

	BMetaGraph::BMetaGraph()
	{
		mSkeleton = BSkeleton();
		vertNodeMap = std::map<SkelVert, MetaV>();
		mComponentMap = {};
		pythonUpToDate = false;
		mDuplicateEdgeMap = std::map<MetaE, MetaE>();
		mDuplicateNodeMap = std::map<MetaV, MetaV>();
	}

	BMetaGraph::BMetaGraph(std::string filename)
	{
		
		mSkeleton = BSkeleton();
		
		vertNodeMap = std::map<SkelVert, MetaV>();
		mComponentMap = {};
		mDuplicateEdgeMap = std::map<MetaE, MetaE>();
		mDuplicateNodeMap = std::map<MetaV, MetaV>();
		
		pythonUpToDate = false;
		
		loadFromFile(filename);
		//std::cout << "Finished loading file" << std::endl;
	}


	void BMetaGraph::loadFromFile(std::string filename)
	{
		std::ifstream filestream;
		filestream.open(filename);
		std::vector<std::string> lines = {};
		std::string line;
		while (!filestream.eof())
		{
			std::getline(filestream, line);
			lines.push_back(line);
		}
		filestream.close();
		//for each(std::string line in lines)
		//{
		//	std::cout << line << std::endl;
		//}
		loadFromLines(lines, 0);
	}

	//void BMetaGraph::loadFromFile(char *filename)
	//{
	//	std::string wrappedFile = std::string(filename);
	//	loadFromFileString(wrappedFile);
	//}

	int BMetaGraph::loadFromLines(std::vector<std::string> lines, int startingLine)
	{
		int result = loadSkeletonFromLines(lines, startingLine);
		//std::cout << "sucessfully loaded skeleton from file" << std::endl;
		pythonUpToDate = false;
		//std::cout << "successfully updated python " << std::endl;
		mSkeleton.findBoundingSphere();
		initializeFromSkeleton();
		return result;
	}

	void BMetaGraph::writeToStream(std::ostream & out)
	{
		mSkeleton.writeToStream(out);
	}

	void BMetaGraph::writeToFile(std::string filename)
	{
		std::ofstream filestream;
		filestream.open(filename);
		writeToStream(filestream);
	}

	int BMetaGraph::loadSkeletonFromLines(std::vector<std::string> lines, int & startingLine)
	{
		return mSkeleton.loadFromLines(lines, startingLine);
	}

	

	void BMetaGraph::initializeFromSkeleton()
	{
		std::vector<SkelVert> nodesOfInterest = mSkeleton.GetNotableVertices();
		std::vector<bool> visitedNodes = std::vector<bool>(mSkeleton.m_vertices.size(), false);

		for each(SkelVert node in nodesOfInterest)
		{
			MetaV startV = addNode(node, &mSkeleton);
			BSkeleton::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(node, mSkeleton);

			for (; adjIt != adjEnd; ++adjIt)
			{
				SkelVert leadVert = *adjIt;
				buildMetaEdge(startV, leadVert, visitedNodes);
			}
		}



		pythonUpToDate = false;
	}

	MetaE BMetaGraph::buildMetaEdge(MetaV startV, SkelVert &lead, std::vector<bool> &visitedSkelVerts)
	{
		//if the metaedge leading out of this node towards the lead has already been visited, then
		//we do not need to create a new edge
		if (visitedSkelVerts[lead])
		{
			return MetaE();
		}

		std::vector<SkelVert> edgeVerts = {};
		SkelVert startSkelVert = operator[](startV).mSrcVert;
		edgeVerts.push_back(startSkelVert);
		visitedSkelVerts[startSkelVert] = true;
		SkelVert last = startSkelVert;
		while (boost::degree(lead, mSkeleton) == 2 && lead != startSkelVert)
		{
			edgeVerts.push_back(lead);
			visitedSkelVerts[lead] = true;

			BSkeleton::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(lead, mSkeleton);
			if (*adjIt != last)
			{
				last = lead;
				lead = *adjIt;
				continue;
			}
			++adjIt;
			if (*adjIt != last)
			{
				last = lead;
				lead = *adjIt;
				continue;
			}

			break;
		}

		edgeVerts.push_back(lead);
		MetaV endV = addNode(lead, &mSkeleton);
		return addEdge(startV, endV, edgeVerts, &mSkeleton);
	}

	MetaV BMetaGraph::addNode(SkelVert srcVert, BSkeleton *srcSkel)
	{
		if (vertNodeMap.count(srcVert) == 0)
		{
			BMetaNode node = BMetaNode(srcVert, srcSkel);
			MetaV v = boost::add_vertex(*this);
			operator[](v) = node;
			vertNodeMap[srcVert] = v;
		}

		pythonUpToDate = false;
		return vertNodeMap[srcVert];
	}

	MetaE BMetaGraph::addEdge(MetaV v0, MetaV v1, std::vector<SkelVert> skelVerts, BSkeleton *srcSkel)
	{
		BMetaEdge edge = BMetaEdge(skelVerts, srcSkel);
		MetaE e;
		bool success;
		boost::tie(e, success) = boost::add_edge(v0, v1, *this);

		operator[](e) = edge;

		operator[](v0).nodeThickness = std::max(operator[](v0).nodeThickness, edge.averageThickness);
		operator[](v1).nodeThickness = std::max(operator[](v1).nodeThickness, edge.averageThickness);

		operator[](v0).nodeWidth = std::max(operator[](v0).nodeWidth, edge.averageWidth);
		operator[](v1).nodeWidth = std::max(operator[](v1).nodeWidth, edge.averageWidth);

		pythonUpToDate = false;
		return e;
	}

	void BMetaGraph::duplicateEdge(MetaE e0, BSkeleton *srcSkel, MetaV &dupV0, MetaV &dupV1, MetaE &dupE)
	{
		BMetaEdge edge = operator[](e0);

		std::vector<SkelVert> addedVerts = {};
		
		//iterate over source vertices, and add copies of their information to the skeleton
		if (edge.mVertices.size() > 0)
		{
			for each(SkelVert v in edge.mVertices)
			{
				addedVerts.push_back(srcSkel->addVertex(srcSkel->getVertData(v)));
			}

			bool exists = false;
			std::vector<SkelEdge> srcEdges = {};
			//find all the edge data associated with the original set of vertices
			for (int i = 0; i < edge.mVertices.size() - 1; ++i)
			{
				SkelEdge e = srcSkel->getEdge(edge.mVertices[i], edge.mVertices[i + 1], exists);
				if (exists)
				{
					srcEdges.push_back(e);
				}
			}

			//add edges between all of the new added vertices containing the info from the original edge data
			for (int i = 0; i < edge.mVertices.size() - 1; ++i)
			{
				srcSkel->addEdge(addedVerts[i], addedVerts[i + 1], srcSkel->getEdgeData(srcEdges[i]));
			}
		}


		//add metanodes and edges for the augmented elements in the skeleton
		dupV0 = addNode(addedVerts[0], srcSkel);
		dupV1 = addNode(addedVerts[addedVerts.size() - 1], srcSkel);
		MetaV srcNode0 = vertNodeMap[edge.mVertices[0]];
		MetaV srcNode1 = vertNodeMap[edge.mVertices[edge.mVertices.size() - 1]];
		dupE = addEdge(dupV0, dupV1, addedVerts, srcSkel);

		
		mDuplicateEdgeMap[dupE] = e0;
		mDuplicateNodeMap[dupV0] = srcNode0;
		mDuplicateNodeMap[dupV1] = srcNode1;
		return;
	}

	void BMetaGraph::removeNode(MetaV nodeToRemove)
	{
		BMetaNode toRemove = operator[](nodeToRemove);
		vertNodeMap.erase(toRemove.mSrcVert);
		boost::remove_vertex(nodeToRemove, *this);
		updateVertNodeMap();
		pythonUpToDate = false;
	}

	void BMetaGraph::removeEdge(MetaE edgeToRemove)
	{
		
		MetaV v0 = boost::source(edgeToRemove, *this);
		MetaV v1 = boost::target(edgeToRemove, *this);

		boost::remove_edge(edgeToRemove, *this);


		BMetaGraph::adjacency_iterator adjIt, adjEnd;

		BMetaGraph::out_edge_iterator edgeIt, edgeEnd;

		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *this);

		float maxThickness = -1;
		float maxWidth = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v0, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
			maxWidth = std::max(maxWidth, operator[](adjEdge).averageWidth);
			++adjIt;
		}
		if (maxThickness <= 0)
		{
			maxThickness = 1;
		}
		if (maxWidth <= 0)
		{
			maxWidth = 1;
		}
		operator[](v0).nodeThickness = maxThickness;
		operator[](v0).nodeWidth = maxWidth;


		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *this);

		maxThickness = -1;
		maxWidth = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v1, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
			maxWidth = std::max(maxWidth, operator[](adjEdge).averageWidth);
			++adjIt;
		}
		if (maxThickness <= 0)
		{
			maxThickness = 1;
		}
		if (maxWidth <= 0)
		{
			maxWidth = 1;
		}
		operator[](v1).nodeThickness = maxThickness / 2.0;

		

		pythonUpToDate = false;
	}

	void BMetaGraph::updateVertNodeMap()
	{
		vertNodeMap.clear();
		
		metaVertIter mvi = boost::vertices(*this);
		for (; mvi.first != mvi.second; ++mvi)
		{
			MetaV  node = *mvi.first;
			SkelVert vert = operator[](node).mSrcVert;
			vertNodeMap[vert] = node;
		}
	}


	void BMetaGraph::findAndLabelConnectedComponents()
	{
		std::cout << "Computing connected components" << std::endl;

		mComponentMap.resize(vertNodeMap.size(), 0);
		int numVertices = boost::num_vertices(*this);

		int numComponents = boost::connected_components(*this, &mComponentMap[0]);
		//std::cout << "Number of components " << numComponents << std::endl;
		std::map<int, int> componentSizeMap = std::map<int, int>();
		for (int i = 0; i < numComponents; ++i)
		{
			componentSizeMap[i] = 0;
		}
		int i = 0;
		for each(int vertComponent in mComponentMap)
		{
			componentSizeMap[vertComponent] += 1;
			++i;
		}
		std::vector<int> allSizes = {};
		for each(std::pair<int, int> componentSizePair in componentSizeMap)
		{
			allSizes.push_back(componentSizePair.second);
		}
		std::sort(allSizes.begin(), allSizes.end());

		std::map<int, int> componentPriorityMap = std::map<int, int>();
		int priority = 0;
		for (int i = allSizes.size() -1 ; i >= 0; --i)
		{
			for each(std::pair<int, int> componentSizePair in componentSizeMap)
			{
				if (componentSizePair.second == allSizes[i] && 
					componentPriorityMap.count(componentSizePair.first) == 0)
				{
					componentPriorityMap[componentSizePair.first] = priority;
					++priority;
				}
			}
		}
		for (int i = 0; i < m_vertices.size(); ++i)
		{
			mComponentMap[i] = componentPriorityMap[mComponentMap[i]];
		}
		for (MetaV node = 0; node < m_vertices.size(); ++node)
		{
			operator[](node).connectedComponent = mComponentMap[node];
		}
	}

	void BMetaGraph::ConnectComponents(int compOne, int compTwo)
	{
		int lowerComponent = std::min(compOne, compTwo);
		int higherComponent = std::max(compOne, compTwo);
		metaVertIter mvi = boost::vertices(*this);

		for (; mvi.first != mvi.second; ++mvi)
		{
			int nodeComponent = operator[](*mvi.first).connectedComponent;
			if (nodeComponent == compOne || nodeComponent == compTwo)
			{
				operator[](*mvi.first).connectedComponent = lowerComponent;
			}
			if (nodeComponent > higherComponent)
			{
				operator[](*mvi.first).connectedComponent--;
			}
		}
	}

	void BMetaGraph::FixUpComponents()
	{
		metaVertIter mvi = boost::vertices(*this);

		for (; mvi.first != mvi.second; ++mvi)
		{
			int nodeComponent = operator[](*mvi.first).connectedComponent;
			BMetaGraph::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(*mvi.first, *this);
			for (; adjIt != adjEnd; ++adjIt)
			{
				int neighborComponent = operator[](*adjIt).connectedComponent;
				if (neighborComponent != nodeComponent)
				{
					ConnectComponents(nodeComponent, neighborComponent);
				}
			}
		}
	}

	void BMetaGraph::FindMinimumSpanningTree()
	{
		mMinimumSpanningTree = {};

		kruskal_minimum_spanning_tree(*this, std::back_inserter(mMinimumSpanningTree),
			weight_map(get(&BMetaEdge::averageThickness, *this)));
	}

	int BMetaGraph::GetNumEdgesToBreak()
	{
		return m_edges.size() - mMinimumSpanningTree.size();
	}

	void BMetaGraph::JoinOperation(Point3d p1picked, Point3d p2picked)
	{
		
		//check if both picked points are in the node map
		if (vertNodeMap.count(p1picked.id) == 0 ||
			vertNodeMap.count(p2picked.id) == 0)
		{
			//std::cout << "Neither point that has been picked is a MetaNode.  Fear ye who enter here" << std::endl;
			return;
		}
		
		//create an edge between these nodes and add it to the metagraph
		std::vector<SkelVert> bridgeVerts = {};
		bridgeVerts.push_back(p1picked.id);
		bridgeVerts.push_back(p2picked.id);
		MetaV node1 = vertNodeMap[p1picked.id];
		MetaV node2 = vertNodeMap[p2picked.id];
		MetaE newMetaEdge = addEdge(node1, node2, bridgeVerts, &mSkeleton);

		BSkeleton::adjacency_iterator adjIt, adjEnd;
		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(p1picked.id, mSkeleton);
		float avgThick1=0.0, avgWidth1=0.0;
		int count1 = 0;
		for (; adjIt != adjEnd; ++adjIt)
		{
			SkelEdge e;
			bool exists;
			boost::tie(e, exists) = boost::edge(p1picked.id, *adjIt, mSkeleton);
			if (exists)
			{
				avgThick1 += mSkeleton[e].thickness;
				avgWidth1 += mSkeleton[e].width;
			}
			++count1;
		}
		avgThick1 /= count1;
		avgWidth1 /= count1;


		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(p2picked.id, mSkeleton);
		float avgThick2 = 0.0, avgWidth2 = 0.0;
		int count2 = 0;
		for (; adjIt != adjEnd; ++adjIt)
		{
			SkelEdge e;
			bool exists;
			boost::tie(e, exists) = boost::edge(p2picked.id, *adjIt, mSkeleton);
			if (exists)
			{
				avgThick2 += mSkeleton[e].thickness;
				avgWidth2 += mSkeleton[e].width;
			}
			++count2;
		}
		avgThick2 /= count2;
		avgWidth2 /= count2;

		RootAttributes newAttributes = RootAttributes((avgThick1 + avgThick2) / 2.0, (avgWidth1 + avgWidth2) / 2.0, p1picked, p2picked);



		SkelEdge newSkelEdge = mSkeleton.addEdge(p1picked.id, p1picked.id, newAttributes);

		operator[](newMetaEdge).mEdges.push_back(newAttributes);
		ConnectComponents(operator[](node1).connectedComponent, operator[](node2).connectedComponent);

		

		//join the components of the two nodes that have been picked
		//int componentOne = operator[](node1).connectedComponent;
		//int componentTwo = operator[](node2).connectedComponent;
		//if (componentOne != componentTwo)
		//{
		//	ConnectComponents(componentOne, componentTwo);
		//}


		//attempt to bridge both of the adjoining nodes (eg. if they are degree 2, then remove them
		//from the metagraph to be replaced by a single edge
		BridgeNode(node1);
		BridgeNode(node2);

		pythonUpToDate = false;
		return;
	}

	void BMetaGraph::BreakOperation(MetaEdge3d toBreak) 
	{
		MetaE e;
		bool exists;
		boost::tie(e, exists) = boost::edge(toBreak.node0, toBreak.node1, *this);
		
		if (exists)
		{
			std::cout << "Matching edge found to break" << std::endl;
			removeEdge(e);
		}
		else
		{
			std::cout << "No matching edge found to break" << std::endl;
		}
		BridgeNode(toBreak.node0);
		BridgeNode(toBreak.node1);

		pythonUpToDate = false;
		return;
	}

	void BMetaGraph::SplitOperation(MetaEdge3d toSplit, std::vector<MetaEdge3d> connectedEdges)
	{
		MetaE dupE;
		MetaV dupV0, dupV1;
		MetaE srcE;
		bool exists = false;
		boost::tie(srcE, exists) = boost::edge(toSplit.node0, toSplit.node1, *this);
		if (exists)
		{
			duplicateEdge(srcE, &mSkeleton, dupV0, dupV1, dupE);
			for (int i = 0; i < connectedEdges.size(); ++i)
			{
				MetaEdge3d srcEdge =connectedEdges[i];
				int srcNode, oldTarget, newTarget;
				//if the srcEdge node0 is equal to either of the nodes on the splitting edge, then 
				//we will break the connection between srcEdge node0 and node1, and create a new edge
				//between node0 and the new target node (determined by mapping of the duplicated nodes
				//to the old target node of the srcEdge
				if (srcEdge.node0 == toSplit.node0 || srcEdge.node0 == toSplit.node1)
				{
					srcNode = srcEdge.node1;
					oldTarget = srcEdge.node0;
					newTarget = mDuplicateNodeMap[dupV0] == oldTarget ? dupV0 : dupV1;
				}
				else if (srcEdge.node1 == toSplit.node0 || srcEdge.node1 == toSplit.node1)
				{
					srcNode = srcEdge.node0;
					oldTarget = srcEdge.node1;
					newTarget = mDuplicateNodeMap[dupV0] == oldTarget ? dupV0 : dupV1;
				}

				MetaE removeE;
				boost::tie(removeE, exists) = boost::edge(srcNode, oldTarget, *this);

				if (exists)
				{
					std::vector<SkelVert> verts = operator[](removeE).mVertices;
					removeEdge(removeE);
					addEdge(srcNode, newTarget, verts, &mSkeleton);
				}
				BMetaNode &newNode = operator[](newTarget);
				BMetaNode &goTowardsNode = operator[](srcNode);
				float deltaX = goTowardsNode.x - newNode.x;
				float deltaY = goTowardsNode.y - newNode.y;
				float deltaZ = goTowardsNode.z - newNode.z;
				float deltaNorm = sqrt(deltaX*deltaX + deltaY * deltaY + deltaZ * deltaZ);
				float xProp = deltaX / deltaNorm;
				float yProp = deltaY / deltaNorm;
				float zProp = deltaZ / deltaNorm;

				float xDist = std::min(abs(xProp *toSplit.avgThickness * 2.0), abs(deltaX / 3.0));
				float yDist = std::min(abs(yProp *toSplit.avgThickness * 2.0), abs(deltaY / 3.0));
				float zDist = std::min(abs(zProp *toSplit.avgThickness * 2.0), abs(deltaZ / 3.0));

				newNode.x += xProp < 0 ? -xDist : xDist;
				newNode.y += yProp < 0 ? -yDist : yDist;
				newNode.z += zProp < 0 ? -zDist : zDist;
				
			}

			//attempt to bridge the nodes at either end of both sides of the split
			//eg. if both the duplicated edge and the src edge will now have a single edge
			//connected at 1 end, then join those two edges together
			BridgeNode(toSplit.node0);
			BridgeNode(toSplit.node1);

			BridgeNode(dupV0);
			BridgeNode(dupV1);

			pythonUpToDate = false;
		}
		return;
	}

	void BMetaGraph::BridgeNode(MetaV nodeToBridge)
	{
		std::cout << "Bridging disabled" << std::endl;
		return;
		if (boost::degree(nodeToBridge, *this) == 2)
		{
			BMetaGraph::out_edge_iterator outIt, outEnd;
			boost::tie(outIt, outEnd) = boost::out_edges(nodeToBridge, *this);
			MetaE e1 = *outIt;
			MetaV v1 = outIt->m_target;
			++outIt;
			MetaE e2 = *outIt;
			MetaV v2 = outIt->m_target;
			BMetaEdge joinedEdge = operator[](e1).join(operator[](e2), &mSkeleton);
			addEdge(v1, v2, joinedEdge.mVertices, &mSkeleton);
			removeEdge(e1);
			removeEdge(e2);
			removeNode(nodeToBridge);

		}
	}

	int BMetaGraph::GetEdgeComponent(MetaE edge)
	{
		BMetaEdge src = operator[](edge);
		//ensure that both endpoint vertices are mapped as metanodes

		if (vertNodeMap.count(src.start()) == 0 ||
			vertNodeMap.count(src.end()) == 0)
		{
			//std::cout << "The endpoints of this metaedge are not mapped as metanodes" << std::endl;
			return -1;
		}

		MetaV node1 = vertNodeMap[src.start()];
		MetaV node2 = vertNodeMap[src.end()];

		int component1 = operator[](node1).connectedComponent;
		int component2 = operator[](node2).connectedComponent;

		if (component1 != component2)
		{
			ConnectComponents(component1, component2);
		}

		return operator[](node1).connectedComponent;
	}


	void BMetaGraph::bridgeUtil(int u, bool visited[], int disc[],
		int low[], int parent[], std::vector<std::pair<MetaV, MetaV>> &bridgeNodes)
	{
		// A static variable is used for simplicity, we can 
		// avoid use of static variable by passing a pointer.
		static int time = 0;

		// Mark the current node as visited
		visited[u] = true;

		// Initialize discovery time and low value
		disc[u] = low[u] = ++time;

		// Go through all vertices aadjacent to this
		BMetaGraph::adjacency_iterator adjIt, adjEnd;
		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(u, *this);
		
		for (; adjIt != adjEnd; ++adjIt)
		{
			int v = *adjIt;

			if (!visited[v])
			{
				parent[v] = u;
				bridgeUtil(v, visited, disc, low, parent, bridgeNodes);

				low[u] = std::min(low[u], low[v]);

				if (low[v] > disc[u])
				{
					bridgeNodes.push_back(std::make_pair(u, v));
					//std::cout << u << " " << v << std::endl;
				}
			}
			else if (v != parent[u])
			{
				low[u] = std::min(low[u], disc[v]);
			}
		}
	}

	// DFS based function to find all bridges. It uses recursive 
	// function bridgeUtil()
	void BMetaGraph::bridge()
	{
		// Mark all the vertices as not visited
		bool *visited = new bool[m_vertices.size()];
		int *disc = new int[m_vertices.size()];
		int *low = new int[m_vertices.size()];
		int *parent = new int[m_vertices.size()];

		boost::python::list bridgeEdges = boost::python::list();

		// Initialize parent and visited arrays
		for (int i = 0; i < m_vertices.size(); i++)
		{
			parent[i] = -1;
			visited[i] = false;
		}

		std::vector<std::pair<MetaV, MetaV>> bridgeNodes = {};

		std::cout << "==================Bridges exist between the following vertices==================" << std::endl;
		// Call the recursive helper function to find Bridges
		// in DFS tree rooted with vertex 'i'
		for (int i = 0; i < m_vertices.size(); i++)
			if (visited[i] == false)
				bridgeUtil(i, visited, disc, low, parent, bridgeNodes);

		
		for each(auto nodepair in bridgeNodes)
		{
			MetaE e;
			bool exists;
			boost::tie(e, exists) = boost::edge(nodepair.first, nodepair.second, *this);

			if (exists)
			{
				operator[](e).isBridge = true;
			}
			
		}

		delete[] visited;
		delete[] disc;
		delete[] low;
		delete[] parent;
		std::cout << "==================End Bridges==================" << std::endl;
	}

	
}


PyMetaGraph::PyMetaGraph(std::string filename)
{
	mGraph = Roots::BMetaGraph(filename);
	mSkeleton = PySkeleton(&mGraph.mSkeleton);
	mGraph.findAndLabelConnectedComponents();
	reload();
}

void PyMetaGraph::initializeFromSkeleton()
{
	mGraph.initializeFromSkeleton();
	mGraph.findAndLabelConnectedComponents();
	reload();
}

void PyMetaGraph::labelComponents()
{
	mGraph.findAndLabelConnectedComponents();
	reload();
}

void PyMetaGraph::joinOperation(int v0, int v1)
{
	//std::cout << "Attempting to join vertices " << v0 << " " << v1 << std::endl;
	Point3d p1 = mGraph.mSkeleton[v0];
	Point3d p2 = mGraph.mSkeleton[v1];
	mGraph.JoinOperation(p1, p2);
	//mSkeleton = PySkeleton(&mGraph.mSkeleton);
	reload();

}

void PyMetaGraph::breakOperation(Roots::MetaEdge3d edge)
{
	//std::cout << "Attempting to break edge between " << edge.node0 << " " << edge.node1 << std::endl;
	std::cout << "Breaking edge between " << edge.node0 << " and " << edge.node1 << ". This is the " << edge.order << "th edge." << std::endl;
	mGraph.BreakOperation(edge);
	mGraph.findAndLabelConnectedComponents();
	reload();
}

void PyMetaGraph::splitOperation(Roots::MetaEdge3d toSplit, boost::python::list connectedEdgeSet)
{
	std::vector<Roots::MetaEdge3d> edgeSet = toStdVector<Roots::MetaEdge3d>(connectedEdgeSet);
	mGraph.SplitOperation(toSplit, edgeSet);
	this->mSkeleton = PySkeleton(&mGraph.mSkeleton);
	mGraph.findAndLabelConnectedComponents();
	reload();
}

void PyMetaGraph::reload()
{
	if (!mGraph.pythonUpToDate)
	{
		//mGraph.findAndLabelConnectedComponents();
		//std::cout << "Entering boost python update loop " << std::endl;
		mGraph.FindMinimumSpanningTree();
		//std::cout << "Num edges to break " << mGraph.GetNumEdgesToBreak() << std::endl;
		numEdgesToBreak = mGraph.GetNumEdgesToBreak();
		mGraph.bridge();
		mMetaNodeLocations = boost::python::list();
		mMetaEdgeConnections = boost::python::list();
		std::vector<Roots::MetaNode3d> metaNodeLocs = {};
		std::vector<Roots::MetaEdge3d> metaEdgeCons = {};
		componentNodeMap = boost::python::dict();
		componentEdgeMap = boost::python::dict();
		componentNodesOfInterestMap = boost::python::dict();

		//stand in std format maps to edit in place before saving to python equivalents
		std::map<int, boost::python::list> tempComponentNodeMap;
		std::map<int, boost::python::list> tempComponentNodesOfInterestMap;
		std::map<int, boost::python::list> tempComponentEdgeMap;
		

		metaVertIter mvi = boost::vertices(mGraph);
		int vertI = 0;
		for (; mvi.first != mvi.second; ++mvi, ++vertI)
		{
			//std::cout << mGraph.m_vertices.size() << std::endl;
			//std::cout << "Trying to find vertex to add : ";
			//std::cout << *mvi.first << std::endl;

			//std::cout << "Getting meta node " << std::endl;
			//Roots::BMetaNode curNode = mGraph[*mvi.first];
			//std::cout << curNode.mSrcVert << std::endl;
			//std::cout << "Successfully accessed meta node " << std::endl;

			int nodeDegree = boost::degree(*mvi.first, mGraph);

			Roots::MetaNode3d toAdd = Roots::MetaNode3d(mGraph[*mvi.first], &mGraph.mSkeleton, vertI, nodeDegree);
			//std::cout << toAdd << std::endl;
			metaNodeLocs.push_back(toAdd);

			if (tempComponentNodeMap.count(toAdd.connectedComponent) == 0)
			{
				boost::python::list toStart = boost::python::list();
				tempComponentNodeMap[toAdd.connectedComponent] = toStart;
				tempComponentNodesOfInterestMap[toAdd.connectedComponent] = boost::python::list();
			}
			//std::cout << "Finding node degree for " << *mvi.first << "th node.  It is : ";
			
			//std::cout << nodeDegree << std::endl;
			if (nodeDegree < 2)
			{
				tempComponentNodesOfInterestMap[toAdd.connectedComponent].append<int>(vertI);
			}
			//std::cout << "Adding node to tempConnectedComponentMap " << std::endl;
			tempComponentNodeMap[toAdd.connectedComponent].append<int>(vertI);
		}
		//std::cout << "completed mapping nodes and components " << std::endl;

		//std::cout << "Beginning to map edges and components" << std::endl;
		metaEdgeIter mei = boost::edges(mGraph);
		int edgeI = 0;
		for (; mei.first != mei.second; ++mei, ++edgeI)
		{
			MetaE edge = *mei.first;
			Roots::BMetaEdge srcEdge = mGraph[edge];
			Roots::MetaEdge3d toAdd = Roots::MetaEdge3d(edge.m_source, edge.m_target,
				srcEdge.getAvgThickness(), srcEdge.averageWidth, mGraph.GetEdgeComponent(edge), edgeI, srcEdge.mEdges, srcEdge.isBridge);
			metaEdgeCons.push_back(toAdd);
			//mMetaEdgeConnections.append<MetaEdge3d>(toAdd);

			if (tempComponentEdgeMap.count(toAdd.connectedComponent) == 0)
			{
				boost::python::list toStart = boost::python::list();
				tempComponentEdgeMap[toAdd.connectedComponent] = toStart;
			}
			tempComponentEdgeMap[toAdd.connectedComponent].append<int>(edgeI);
		}

		for each(std::pair<int, boost::python::list> nodeMap in tempComponentNodeMap)
		{
			componentNodeMap[nodeMap.first] = nodeMap.second;
		}

		for each(std::pair<int, boost::python::list> nodeMap in tempComponentNodesOfInterestMap)
		{
			componentNodesOfInterestMap[nodeMap.first] = nodeMap.second;
		}

		for each(std::pair<int, boost::python::list> edgeMap in tempComponentEdgeMap)
		{
			componentEdgeMap[edgeMap.first] = edgeMap.second;
		}

		mMetaNodeLocations = toPythonList<Roots::MetaNode3d>(metaNodeLocs);

		mMetaEdgeConnections = toPythonList<Roots::MetaEdge3d>(metaEdgeCons);

		mGraph.pythonUpToDate = true;


	}
	//std::cout << "completed boost python update loop " << std::endl;
}