#pragma once

#include "BoostSkeleton.h"
#include "boost/python.hpp"
#include "arcball.h"
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
	ByComponent = 4,
	Other = 5
};

enum OperationMode
{
	None = 0,
	Connection = 1,
	AddNodes = 2,
	Break = 3,
	Split = 4
};


#define nodeSelectionScaling 2.5
#define nonBridgeScaling 3.0
#define edgeSelectionScaling 5.0
#define useArcball false


namespace Roots
{
	struct ColorTable
	{
		static std::vector<std::vector<float>> colors;
		ColorTable()
		{
			return;
		}
		static std::vector<float> getComponentColor(int component)
		{
			if (component > ColorTable::colors.size())
			{
				colors.push_back({});
				for (int c = 0; c < 3; ++c)
				{
					float randomC = ((float)rand() / RAND_MAX);
					colors.back().push_back(randomC);
				}
				colors.back().push_back(1.0);
			}
			return colors[component];
		}
	};

	struct EdgeVisualizationOptions
	{
		bool show;
		float scale;
		bool magnifyNonBridges;
		bool showOnlyNonBridges;

		//scale of 0-1 between min and max for whatever attribute is being colorized
		//lower/higher values are flooded
		float minColorCutoff, maxColorCutoff;
		ColorizationOptions colorization;
		std::vector<std::vector<GLfloat>> heatmap;
		float flatSelectionColor[4];

		EdgeVisualizationOptions();
	};

	struct NodeVisualizationOptions
	{
		bool showEndpoints;
		float endpointScale;
		bool showJunctions;
		float junctionScale;

		bool useConstantNodecolor;
		float constantNodeColor[4];

		float minColorCutoff, maxColorCutoff;
		ColorizationOptions colorization;
		std::vector<std::vector<GLfloat>> heatmap;



		NodeVisualizationOptions();
	};

	struct BMetaNode
	{
		SkelVert mSrcVert;
		BSkeleton* mSrcSkeleton;
		bool hasGeom;
		float p[3];
		int degree;
		int connectedComponent;
		float nodeThickness;
		float nodeWidth;
		float glThicknessColor[4];
		float glWidthColor[4];
		float glDegreeColor[4];
		float glComponentColor[4];
		float *currentColor;
		

		BMetaNode();
		BMetaNode(SkelVert srcId, BSkeleton *skel);
		void updateColors(NodeVisualizationOptions options);
		void updateComponentColor();
		float &operator[](size_t i)
		{
			return p[i];
		}
		float x()
		{
			return p[0];
		}
		float y()
		{
			return p[1];
		}
		float z()
		{
			return p[2];
		}
	};

	struct BMetaEdge
	{
		std::vector<SkelVert> mVertices;
		std::vector<RootAttributes> mEdges;
		//BSkeleton *mSrcSkeleton;
		float averageThickness, averageWidth, mLength;
		bool isBridge;
		bool isSelected;

		int connectedComponent;
		

		//whenever we remove a meta edge, we will just completely rebuild the indices lists
		//this is the list of this edges indices for each line referencing the vertex indices in 
		//a glvertices list
		std::vector<GLuint> indicesList;
		//GLfloat *glVertices;

		std::vector<GLfloat> localThicknessColors, localWidthColors, localRatioColors, localComponentColors;


		BMetaEdge();
		BMetaEdge(std::vector<SkelVert> vertices, BSkeleton* srcSkeleton, std::vector<GLfloat> &glVertices);
		BMetaEdge join(BMetaEdge &other, BSkeleton* srcSkeleton);
		std::pair<BMetaEdge, BMetaEdge> split(SkelVert toSplitOn, BSkeleton *srcSkeleton);
		float getAvgThickness();
		SkelVert start();
		SkelVert end();
		void addToIndicesList(std::vector<GLuint> &edgeIndices);

		void updateColors(EdgeVisualizationOptions options, std::vector<std::vector<GLfloat>> &edgeColors);
		void updateComponentColor(std::vector<std::vector<GLfloat>> &edgeColors);
		void select(GLfloat *selectionColor, std::vector<std::vector<GLfloat>> &edgeColors);
		void unselect(std::vector<std::vector<GLfloat>> &edgeColors);
		void updateGraphColors(std::vector<std::vector<GLfloat>> &edgeColors);
	};

	struct BoundingBox
	{
		std::vector<std::vector<float>> corners;
		bool hasPoints;
		float minx, maxx, miny, maxy, minz, maxz;

		BoundingBox();

		void addPoint(float *p);
		void draw(std::vector<GLfloat> componentColor, float lineWidth);
	};

	//struct MetaNode3d : public Point3d
	//{
	//	int connectedComponent;
	//	int order;
	//	float nodeThickness;
	//	float nodeWidth;
	//	int degree;
	//	MetaNode3d();
	//	MetaNode3d(BMetaNode srcNode, BSkeleton *srcSkel, int aOrder, int aDegree);
	//};

	//struct MetaEdge3d
	//{
	//	int node0, node1;
	//	float avgThickness, avgWidth;
	//	int connectedComponent;
	//	int order;
	//	boost::python::list edges;
	//	bool isBridge;

	//	MetaEdge3d();
	//	MetaEdge3d(int n0, int n1, float thickness, float width, int connectedComponent, int aOrder, std::vector<RootAttributes> aEdges, bool aIsBridge);
	//};
}

namespace
{
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Roots::BMetaNode, Roots::BMetaEdge, boost::no_property, boost::listS> BoostMetaGraph;
	typedef boost::graph_traits<BoostMetaGraph>::vertex_descriptor MetaV;
	typedef boost::graph_traits<BoostMetaGraph>::edge_descriptor MetaE;
}

namespace Roots{
	struct MinMaxStruct
	{
		static float minThickness, maxThickness, minWidth, maxWidth, minRatio, maxRatio, minComponent, maxComponent, minDegree, maxDegree;
	};
	



	struct BMetaGraph : public BoostMetaGraph
	{
		/*vvvvvvvvvvvvvvvvvvvvv MEMBER VARIABLES vvvvvvvvvvvvvvvvvvvvv*/
		//the skeleton from which this metagraph is constructed
		BSkeleton mSkeleton;

		//map between skeleton vertices and their corresponding metanodes
		std::map<SkelVert, MetaV> vertNodeMap;

		//vector correlating each edge to a connected component
		std::vector<int> mComponentMap;
		std::vector<float> mComponentSizeMap;
		std::vector<MetaE> mMinimumSpanningTree;
		std::vector<BoundingBox> componentBounds;
		int numComponents;

		//map backwards from the duplicated element to the src element
		//these maps are only to be used for writing out the files
		std::map<MetaE, MetaE> mDuplicateEdgeMap;
		std::map<MetaV, MetaV> mDuplicateNodeMap;
		bool isLoaded;


		//visualization and drawing options
		EdgeVisualizationOptions edgeOptions;
		NodeVisualizationOptions nodeOptions;
		OperationMode mode;
		static drawing::VBOSphere drawSphere;
		static drawing::VBOCube drawCube;

		//vbo objects
		std::vector<GLfloat> edgeVertices;
		std::vector<std::vector<GLfloat>> edgeColors;
		//std::vector<GLuint> nonBridgeEdgeIndices;
		//std::vector<GLuint> bridgeEdgeIndices;
		//std::vector<GLuint> baseEdgeIndices;
		//std::vector<std::vector<GLuint>> componentIndices;
		GLfloat selectionColor[4];
		
		bool onlyDisplaySelectedComponents;
		int selectedComponent1, selectedComponent2;
		bool showBoundingBoxes;

		MetaV selectNode1, selectNode2;
		bool selectNode1Valid, selectNode2Valid;

		std::pair<MetaV, MetaV> breakEdge;
		bool breakEdgeValid;

		std::pair<MetaV, MetaV> splitEdge;
		bool splitEdgeValid;
		std::vector<std::pair<MetaV, MetaV>> splitNeighbors;

		

		//this should not be used anymore
		bool pythonUpToDate;



		int getNumGLEdges();
		int getNumGLVertices();

		
		//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvSettingsvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

		
		void assignEdgeHeatMap(boost::python::list heatmap);
		void showEdges(bool doShow);
		void setEdgeScale(float scale);
		void magnifyNonBridges(bool doMagnify);
		void showOnlyNonBridges(bool showOnly);
		void colorizeEdgesByThickness();
		void colorizeEdgesByWidth();
		void colorizeEdgesByRatio();
		void colorizeEdgesByComponent();
		void setEdgeColorFloor(float val);
		void setEdgeColorCeiling(float val);
		void setEdgeSelectionColor(float r, float g, float b);

		void assignNodeHeatMap(boost::python::list heatmap);
		void showEndpoints(bool doShow);
		void showJunctions(bool doShow);
		void setEndpointScale(float scale);
		void setJunctionScale(float scale);
		void setConstantNodeColor(float r, float g, float b);
		void colorizeNodesByThickness();
		void colorizeNodesByWidth();
		void colorizeNodesByDegree();
		void colorizeNodesByComponent();
		void colorizeNodesByConstantColor();
		void setNodeColorFloor(float val);
		void setNodeColorCeiling(float val);

		void setDisplayOnlySelectedComponents(bool showOnlySelected);
		void setComponent1(int component);
		void setComponent2(int component);
		void setShowBoundingBoxes(bool doShow);

		void drawEdges();
		void drawNodes();
		void drawBoxes();
		void draw();

		void startRotation(int mx, int my);
		void mouseMoved(int mx, int my);
		void setZoom(float rad, float eyex, float eyey, float eyez, float upx, float upy, float upz);

		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Settings^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

		//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvInteractionsvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

		void selectConnectionNode(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ);

		//void selectMetaNode(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ);
		//void selectPromotionNode(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ);
		void selectBreakEdge(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ);
		void selectSplitEdge(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ);
		void setSelectionColor(float r, float g, float b);
		void unselectAll();

	private:
		MetaV getFirstMetaNodeHit(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ, bool &isValid);
		SkelVert getFirstSkelNodeHit(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ, bool &isValid);
		std::pair<MetaV, MetaV> getFirstMetaEdgeHit(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ, bool &isValid);
		void privateSelectConnectionNode(SkelVert nodeVert, int selectComponent);
		void unselectAllEdges();

		void unselectEdge(std::pair<MetaV, MetaV> toUnselect);
		void selectEdge(std::pair<MetaV, MetaV> toSelect);
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Interactions^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	public:





		boost::python::list getComponentSizes();

		int getNumComponents()
		{
			return numComponents;
		}

		/*vvvvvvvvvvvvvvvvvvvvv CONSTRUCTORS vvvvvvvvvvvvvvvvvvvvv*/

		/*default constructor.  Only constructor defined.  All operations such as loading or
		adding edges should be done with appropriate member functions*/
		BMetaGraph();
		void Initialize();
		//BMetaGraph(std::string filename);

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
		MetaE buildMetaEdge(MetaV start, SkelVert &lead, std::vector<bool> &skelVertsVisited, bool isLoading);

	public:

		/*
		Only call that should be made to add vertices
		*/
		MetaV addNode(SkelVert srcVert, BSkeleton *srcSkel, bool inheretEdgeParams=false);
		/*
		Only call that should be made to add edges
		*/
		MetaE addEdge(MetaV v0, MetaV v1, std::vector<SkelVert> skelVerts, BSkeleton *srcSkel, bool isLoading, bool isJoining=false);

		/*
		Only call that should be made to duplicate edges
		*/
		void duplicateEdge(MetaE e0, BSkeleton *srcSkel, MetaV &dupV0, MetaV &dupV1, MetaE &dupE);

		void removeNode(MetaV nodeToRemove);

		void removeEdge(MetaE edgeToRemove);

		void removeEdgeNoBridging(MetaE edgeToRemove);

	private:
		void updateVertNodeMap();

		void updateNodeDegrees();

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

		void findBridges();
		/*
		This method should support the first operation for this tool - making new connections 
		between isolated components.  Provided two meta nodes, create a new connection between them.
		Most likely making this new connection will result in both nodes being degree two, so the
		metaedge between them, and the meta edges leading to the rest of their connected components
		should be merged, and then these nodes should be deleted, and a new edge created between
		the remaining nodes
		*/
		void JoinOperation();
		//void JoinOperation(Point3d p1picked, Point3d p2picked);

		/*
		This method will also support the first operation.  It will help to pick out the nodes and 
		edges between them belonging to a single connected component, returned as an adjacency list
		between points of the same id.
		*/
		void BreakOperation();
		//void BreakOperation(MetaEdge3d toBreak);

		/*
		This method will support the splitting operation.  Specifically, one central edge will be selected, toSplit.
		Then, a group of edges connected at the source node and target nodes of that edge of that edge will be selected.
		This set of connected edges edges at the source and target nodes will then be disconnected from the
		source edge, and reconnected to a duplicate of the central edge.
		*/
		void SplitOperation();
		//void SplitOperation(MetaEdge3d toSplit, std::vector<MetaEdge3d> connectedEdgeSet);


		void PromoteOperation(SkelVert toPromote);
	
	private:
		/*
		Determines if the provided metanode is of degree 2, and if so, joines the edges leading
		into this node, creates a new edge in the graph bridging its two adjacent nodes, and 
		removing this node and its incident edges from the graph.
		*/
		bool BridgeNode(SkelVert vertToBridge);

		//void BridgeNodeSweeper();
	public:
		/*
		Use the vert node map to determine the connected component to which an edge belongs
		*/
		int GetEdgeComponent(MetaE edge);
	};

}

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


//struct PyMetaGraph
//{
//	PySkeleton mSkeleton;
//	//list of MetaNode3d values in the same order as the vertices of this metagraph
//	boost::python::list mMetaNodeLocations;
//	//list of MetaEdge3d values providing the connections between the metanodes in the previous list 
//	boost::python::list mMetaEdgeConnections;
//	//connected component maps for vertices and edges, each maps between the connected components
//	//and a list of the edges/vertices which map to it
//	boost::python::dict componentNodeMap;
//	boost::python::dict componentEdgeMap;
//	boost::python::dict componentNodesOfInterestMap;
//	int numEdgesToBreak;
//	
//	PyMetaGraph(std::string filename);
//
//	boost::python::list getComponentNodes(int component);
//	boost::python::list getComponentEdges(int component);
//	
//
//	void initializeFromSkeleton();
//	void reload();
//	void labelComponents();
//	void joinOperation(int v0, int v1);
//	void breakOperation(Roots::MetaEdge3d edge);
//	void splitOperation(Roots::MetaEdge3d toSplit, boost::python::list connectedEdgeSet);
//	Roots::BMetaGraph mGraph;
//
//	
//};