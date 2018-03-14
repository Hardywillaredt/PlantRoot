#include "BoostMetaGraph.h"
#include "boost/graph/connected_components.hpp"
#include <iostream>
#include <fstream>
namespace
{
	typedef boost::graph_traits<Roots::BMetaGraph>::vertex_iterator vertIter;
	typedef boost::graph_traits<Roots::BMetaGraph>::edge_iterator edgeIter;

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
}
namespace Roots
{
	/////////////////////////////////////////////// BMetaNode /////////////////////////////////////
	BMetaNode::BMetaNode()
	{
		mSrcVert = 0;
		mSrcSkeleton = nullptr;
		connectedComponent = -1;
		nodeSize = 0.0;
	}
	BMetaNode::BMetaNode(SkelVert srcId, BSkeleton *skel)
	{
		mSrcVert = srcId;
		mSrcSkeleton = skel;
		//std::cout << "number of vertices in pointed skeleton " << boost::num_vertices(*skel) << std::endl;
		connectedComponent = -1;
		nodeSize = 0.0;
	}

	////////////////////////////////////////////// BMetaEdge ///////////////////////////////////////

	BMetaEdge::BMetaEdge()
	{
		mVertices = {};
		mEdges = {};
		averageThickness = 0.0;
		//mSrcSkeleton = nullptr;

	}

	BMetaEdge::BMetaEdge(std::vector<SkelVert> vertices, BSkeleton* srcSkeleton)
	{
		mVertices = vertices;
		mEdges = {};
		float weightedAvgThickness = 0.0;
		float skelEdgeLength = 0.0;
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
				skelEdgeLength += length;
			}
			else
			{
				//this is okay - this just means that this edge does not exist at the skeleton
				//level, eg it is an added edge at the metagraph level
				//do nothing here
				float cumThickness = 0.0;
				int numEdges = 0;

				BSkeleton::adjacency_iterator adjIt, adjEnd;
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *srcSkeleton);
				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelVert leadVert = *adjIt;
					boost::tie(e, exists) = boost::edge(v0, leadVert, *srcSkeleton);
					RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
					cumThickness += skelRootAttributes[Thickness];
					++numEdges;
				}
				boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *srcSkeleton);
				for (; adjIt != adjEnd; ++adjIt)
				{
					SkelVert leadVert = *adjIt;
					boost::tie(e, exists) = boost::edge(v1, leadVert, *srcSkeleton);
					RootAttributes skelRootAttributes = srcSkeleton->operator[](e);
					cumThickness += skelRootAttributes[Thickness];
					++numEdges;
				}
				skelEdgeLength = 1.0;
				weightedAvgThickness = cumThickness / numEdges;
			}
		}
		averageThickness = weightedAvgThickness / skelEdgeLength;
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
		nodeSize = -1;
	}

	MetaNode3d::MetaNode3d(BMetaNode src, Roots::BSkeleton *srcSkeleton, int aOrder)
		:Point3d(srcSkeleton->operator[](src.mSrcVert))
	{
		order = aOrder;
		connectedComponent = src.connectedComponent;
		nodeSize = src.nodeSize;
	}



	///////////////////////////////////////// MetaEdge3d //////////////////////////////////////////

	MetaEdge3d::MetaEdge3d()
	{
		node0 = 0;
		node1 = 0;
		avgThickness = 0;
		connectedComponent = -1;
		edges = boost::python::list();

	}

	MetaEdge3d::MetaEdge3d(int n0, int n1, float thickness, int component, int aOrder, std::vector<RootAttributes> aEdges)
	{
		node0 = n0;
		node1 = n1;
		avgThickness = thickness;
		connectedComponent = component;
		order = aOrder;
		edges = toPythonList<RootAttributes>(aEdges);
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
		findAndLabelConnectedComponents();
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

		while (boost::degree(lead, mSkeleton) == 2)
		{
			edgeVerts.push_back(lead);
			visitedSkelVerts[lead] = true;

			BSkeleton::adjacency_iterator adjIt, adjEnd;
			boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(lead, mSkeleton);

			for (; adjIt != adjEnd; ++adjIt)
			{
				if (!visitedSkelVerts[*adjIt])
				{
					lead = *adjIt;
				}
			}
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

		operator[](v0).nodeSize = std::max(operator[](v0).nodeSize, edge.averageThickness);
		operator[](v1).nodeSize = std::max(operator[](v1).nodeSize, edge.averageThickness);

		pythonUpToDate = false;
		return e;
	}

	void BMetaGraph::duplicateEdge(MetaE e0, BSkeleton *srcSkel, MetaV &dupV0, MetaV &dupV1, MetaE &dupE)
	{
		BMetaEdge edge = operator[](e0);

		std::vector<SkelVert> addedVerts = {};
		
		//iterate over source vertices, and add copies of their information to the skeleton
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

		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v0, *this);

		float maxThickness = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v0, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
		}
		operator[](v0).nodeSize = maxThickness;


		boost::tie(adjIt, adjEnd) = boost::adjacent_vertices(v1, *this);

		maxThickness = -1;
		while (adjIt != adjEnd)
		{
			bool exists;
			MetaE adjEdge;
			boost::tie(adjEdge, exists) = boost::edge(v1, *adjIt, *this);

			maxThickness = std::max(maxThickness, operator[](adjEdge).averageThickness);
		}
		operator[](v1).nodeSize = maxThickness;
		
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
		//std::cout << "Computing connected components" << std::endl;

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
		metaVertIter mvi = boost::vertices(*this);

		for (; mvi.first != mvi.second; ++mvi)
		{
			int nodeComponent = operator[](*mvi.first).connectedComponent;
			if (nodeComponent == compOne || nodeComponent == compTwo)
			{
				operator[](*mvi.first).connectedComponent = lowerComponent;
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
		addEdge(node1, node2, bridgeVerts, &mSkeleton);

		//join the components of the two nodes that have been picked
		//int componentOne = operator[](node1).connectedComponent;
		//int componentTwo = operator[](node2).connectedComponent;
		//if (componentOne != componentTwo)
		//{
		//	ConnectComponents(componentOne, componentTwo);
		//}

		//findAndLabelConnectedComponents();

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
			removeEdge(e);
		}
		BridgeNode(toBreak.node0);
		BridgeNode(toBreak.node1);

		pythonUpToDate = false;
		return;
	}

	void BMetaGraph::SplitOperation(MetaEdge3d toSplit, boost::python::list connectedEdges)
	{
		MetaE dupE;
		MetaV dupV0, dupV1;
		MetaE srcE;
		bool exists = false;
		boost::tie(srcE, exists) = boost::edge(toSplit.node0, toSplit.node1, *this);
		if (exists)
		{
			duplicateEdge(srcE, &mSkeleton, dupV0, dupV1, dupE);
			for (int i = 0; i < boost::python::len(connectedEdges); ++i)
			{
				MetaEdge3d srcEdge = boost::python::extract<MetaEdge3d>(connectedEdges[i]);
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
				
			}

			//attempt to bridge the nodes at either end of both sides of the split
			//eg. if both the duplicated edge and the src edge will now have a single edge
			//connected at 1 end, then join those two edges together
			BridgeNode(toSplit.node0);
			BridgeNode(toSplit.node1);

			BridgeNode(dupV0);
			BridgeNode(dupV1);

			findAndLabelConnectedComponents();
			pythonUpToDate = false;
		}
		return;
	}

	void BMetaGraph::BridgeNode(MetaV nodeToBridge)
	{
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
	
}


PyMetaGraph::PyMetaGraph(std::string filename)
{
	mGraph = Roots::BMetaGraph(filename);
	mSkeleton = PySkeleton(&mGraph.mSkeleton);
}

void PyMetaGraph::initializeFromSkeleton()
{
	mGraph.initializeFromSkeleton();
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
	reload();

}

void PyMetaGraph::breakOperation(Roots::MetaEdge3d edge)
{
	//std::cout << "Attempting to break edge between " << edge.node0 << " " << edge.node1 << std::endl;
	mGraph.BreakOperation(edge);
	reload();
}

void PyMetaGraph::reload()
{
	if (!mGraph.pythonUpToDate)
	{
		mGraph.findAndLabelConnectedComponents();
		//std::cout << "Entering boost python update loop " << std::endl;
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

			Roots::MetaNode3d toAdd = Roots::MetaNode3d(mGraph[*mvi.first], &mGraph.mSkeleton, vertI);
			//std::cout << toAdd << std::endl;
			metaNodeLocs.push_back(toAdd);

			if (tempComponentNodeMap.count(toAdd.connectedComponent) == 0)
			{
				boost::python::list toStart = boost::python::list();
				tempComponentNodeMap[toAdd.connectedComponent] = toStart;
				tempComponentNodesOfInterestMap[toAdd.connectedComponent] = boost::python::list();
			}
			//std::cout << "Finding node degree for " << *mvi.first << "th node.  It is : ";
			int nodeDegree = boost::degree(*mvi.first, mGraph);
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
				srcEdge.getAvgThickness(), mGraph.GetEdgeComponent(edge), edgeI, srcEdge.mEdges);
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