#pragma once

#include "BoostSkeleton.h"
#include "boost/python.hpp"

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>

#include "Sphere.h"

typedef GLuint GLLineIndex;
typedef GLuint GLVertexIndex;


enum ColorizationOptions
{
	ByThickness = 0,
	ByWidth = 1,
	ByRatio = 2,
	ByDegree = 3,
	ByComponent = 4
};

enum OperationModes
{
	None = 0,
	Connection = 1,
	AddNodes = 2,
	Break = 3,
	Split = 4
};

struct NodeVisualizationOptions
{
	bool show;
	bool highlightEndpoints;
	bool highlightJunctions;
	ColorizationOptions colorization;
	boost::python::list heatmap;
	float baseScale;
	float highlightScale;
};

struct EdgeVisualizationOptions
{
	bool show;
	ColorizationOptions colorization;
	boost::python::list heatmap;
	float baseScale;
	float highlightScale;
};

namespace Roots
{
	struct BMetaNode
	{
		SkelVert mSrcVert;
		BSkeleton* mSrcSkeleton;
		bool hasGeom;
		float x, y, z;
		int connectedComponent;
		float nodeThickness;
		float nodeWidth;

		BMetaNode();
		BMetaNode(SkelVert srcId, BSkeleton *skel);

	};

	struct BMetaEdge
	{
		std::vector<SkelVert> mVertices;
		std::vector<RootAttributes> mEdges;
		//BSkeleton *mSrcSkeleton;
		float averageThickness, averageWidth;
		bool isBridge;

		int node0, node1;
		int id;
		int component;

		GLVertexIndex glEdgeVerticesStart, glEdgeVerticesEnd;
		GLLineIndex glEdgeLinesStart, glEdgeLinesEnd;
		bool isPartOfLinesList;



		BMetaEdge();
		BMetaEdge(std::vector<SkelVert> vertices, BSkeleton* srcSkeleton, std::vector<GLfloat> &edgeVertices, std::vector<GLVertexIndex> &edgeIndices);
		BMetaEdge join(BMetaEdge &other, BSkeleton* srcSkeleton);
		float getAvgThickness();
		SkelVert start();
		SkelVert end();
	};

	struct MetaNode3d : public Point3d
	{
		int connectedComponent;
		int order;
		float nodeThickness;
		float nodeWidth;
		int degree;
		MetaNode3d();
		MetaNode3d(BMetaNode srcNode, BSkeleton *srcSkel, int aOrder, int aDegree);
	};

	struct MetaEdge3d
	{
		int node0, node1;
		float avgThickness, avgWidth;
		int connectedComponent;
		int order;
		boost::python::list edges;
		bool isBridge;

		MetaEdge3d();
		MetaEdge3d(int n0, int n1, float thickness, float width, int connectedComponent, int aOrder, std::vector<RootAttributes> aEdges, bool aIsBridge);
	};
}

namespace
{
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Roots::BMetaNode, Roots::BMetaEdge, boost::no_property, boost::listS> BoostMetaGraph;
	typedef boost::graph_traits<BoostMetaGraph>::vertex_descriptor MetaV;
	typedef boost::graph_traits<BoostMetaGraph>::edge_descriptor MetaE;
}

namespace Roots{
	struct BMetaGraph : public BoostMetaGraph
	{
		/*vvvvvvvvvvvvvvvvvvvvv MEMBER VARIABLES vvvvvvvvvvvvvvvvvvvvv*/
		//the skeleton from which this metagraph is constructed
		BSkeleton mSkeleton;

		//map between skeleton vertices and their corresponding metanodes
		std::map<SkelVert, MetaV> vertNodeMap;

		//vector correlating each edge to a connected component
		std::vector<int> mComponentMap;
		std::vector<MetaE> mMinimumSpanningTree;

		//map backwards from the duplicated element to the src element
		//these maps are only to be used for writing out the files
		std::map<MetaE, MetaE> mDuplicateEdgeMap;
		std::map<MetaV, MetaV> mDuplicateNodeMap;


		//visualization and drawing options
		EdgeVisualizationOptions edgeOptions;
		NodeVisualizationOptions nodeOptions;
		drawing::VBOSphere drawSphere;

		//vbo objects
		std::vector<GLfloat> edgeVertices;
		int numEdges;
		std::vector<GLfloat> edgeColors;
		std::vector<GLuint> highlightEdges;
		std::vector<GLuint> baseEdges;
		

		void drawEdges();
		void drawNodes();

		bool pythonUpToDate;

		/*vvvvvvvvvvvvvvvvvvvvv CONSTRUCTORS vvvvvvvvvvvvvvvvvvvvv*/

		/*default constructor.  Only constructor defined.  All operations such as loading or
		adding edges should be done with appropriate member functions*/
		BMetaGraph();
		BMetaGraph(std::string filename);

		/*vvvvvvvvvvvvvvvvvvvvv LOADERS AND SAVERS vvvvvvvvvvvvvvvvvvvvv*/

		void loadFromFile(std::string filename);
		//void loadFromFile(boost::python::str filename);
		/*
		External load operation, this is the only one that should be used outside this class.
		Aliases loading operations for different file structures internally.
		
		Parameters: 
		lines - Set of line strings defining the skeleton
		startingLine - Index of first line from provided set to load from

		returns - line number after the last line consumed for the loading process
		*/
		int loadFromLines(std::vector<std::string> lines, int startingLine);

		void writeToFile(std::string filename);
		/*
		External save operation, only save functionality for this class.
		Inserts this metagraph as a string with a ply-style header into the provided ostream.

		Parameters:
		out - ostream to push this skelton into
		*/
		void writeToStream(std::ostream &out);

	private:
		/*
		Load just a skeleton, no metagraph information
		*/
		int loadSkeletonFromLines(std::vector<std::string> lines, int &startingLine);



	public:
		/*
		Initializes the metagraph from the structure of the underlying skeleton.  Creates metanodes
		for all vertices of the skeleton which are not of degree two, and creates connections
		between those nodes based on the connections between these 'interesting' vertices in the skeleton.

		Builds the underyling lists of member vertices for each of the meta edges
		*/
		void initializeFromSkeleton();

	private:
		/*
		Build a metaedge by iterating from the starting metavert in the direction of the lead skelvert.

		Stop upon reaching a skelvert of degree != 2, add that vert as a metanode, and create a
		MetaEdge between the start vert and it.  Both the start vert and the end vert need to be in
		the metaedge
		*/
		MetaE buildMetaEdge(MetaV start, SkelVert &lead, std::vector<bool> &skelVertsVisited);

	public:

		/*
		Only call that should be made to add vertices
		*/
		MetaV addNode(SkelVert srcVert, BSkeleton *srcSkel);
		/*
		Only call that should be made to add edges
		*/
		MetaE addEdge(MetaV v0, MetaV v1, std::vector<SkelVert> skelVerts, BSkeleton *srcSkel);

		/*
		Only call that should be made to duplicate edges
		*/
		void duplicateEdge(MetaE e0, BSkeleton *srcSkel, MetaV &dupV0, MetaV &dupV1, MetaE &dupE);

		void removeNode(MetaV nodeToRemove);

		void removeEdge(MetaE edgeToRemove);

	private:
		void updateVertNodeMap();

	public:
		/*
		Finds all of the connected components of this graph, and stores the components in the
		component map sorted by their size (eg. component 0 is the largest component)
		*/
		void findAndLabelConnectedComponents();

		/*
		Iterates over all vertices, and on each vertex belonging to one of the provided components,
		stores their mapped component id as the lower id between the two components 
		(thus the largest components will tend to have the smallest id's)
		*/
		void ConnectComponents(int compOne, int compOther);

		/*
		Iterates over all vertices, and at each vertex iterates over its neighbors.  If one of its
		neighbors belongs to a different component, run the connect components routine for their
		two componentIDs
		*/
		void FixUpComponents();

		

		void FindMinimumSpanningTree();

		int GetNumEdgesToBreak();

		void bridgeUtil(int u, bool visited[], int disc[],
			int low[], int parent[], std::vector<std::pair<MetaV, MetaV>> &bridgeNodes);

		void bridge();
		/*
		This method should support the first operation for this tool - making new connections 
		between isolated components.  Provided two meta nodes, create a new connection between them.
		Most likely making this new connection will result in both nodes being degree two, so the
		metaedge between them, and the meta edges leading to the rest of their connected components
		should be merged, and then these nodes should be deleted, and a new edge created between
		the remaining nodes
		*/
		void JoinOperation(Point3d p1picked, Point3d p2picked);

		/*
		This method will also support the first operation.  It will help to pick out the nodes and 
		edges between them belonging to a single connected component, returned as an adjacency list
		between points of the same id.
		*/

		void BreakOperation(MetaEdge3d toBreak);

		/*
		This method will support the splitting operation.  Specifically, one central edge will be selected, toSplit.
		Then, a group of edges connected at the source node and target nodes of that edge of that edge will be selected.
		This set of connected edges edges at the source and target nodes will then be disconnected from the
		source edge, and reconnected to a duplicate of the central edge.
		*/
		void SplitOperation(MetaEdge3d toSplit, std::vector<MetaEdge3d> connectedEdgeSet);
	
	private:
		/*
		Determines if the provided metanode is of degree 2, and if so, joines the edges leading
		into this node, creates a new edge in the graph bridging its two adjacent nodes, and 
		removing this node and its incident edges from the graph.
		*/
		void BridgeNode(MetaV nodeToBridge);
	public:
		/*
		Use the vert node map to determine the connected component to which an edge belongs
		*/
		int GetEdgeComponent(MetaE edge);
	};



}

struct PyMetaGraph
{
	PySkeleton mSkeleton;
	//list of MetaNode3d values in the same order as the vertices of this metagraph
	boost::python::list mMetaNodeLocations;
	//list of MetaEdge3d values providing the connections between the metanodes in the previous list 
	boost::python::list mMetaEdgeConnections;
	//connected component maps for vertices and edges, each maps between the connected components
	//and a list of the edges/vertices which map to it
	boost::python::dict componentNodeMap;
	boost::python::dict componentEdgeMap;
	boost::python::dict componentNodesOfInterestMap;
	int numEdgesToBreak;
	
	PyMetaGraph(std::string filename);

	boost::python::list getComponentNodes(int component);
	boost::python::list getComponentEdges(int component);
	

	void initializeFromSkeleton();
	void reload();
	void labelComponents();
	void joinOperation(int v0, int v1);
	void breakOperation(Roots::MetaEdge3d edge);
	void splitOperation(Roots::MetaEdge3d toSplit, boost::python::list connectedEdgeSet);


private:
	Roots::BMetaGraph mGraph;
};