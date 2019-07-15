#pragma once
#include "BoostSkeleton.h"
#include "boost/python.hpp"
#include "arcball.h"
#include <GLFW/glfw3.h>
#include "boost/graph/dijkstra_shortest_paths.hpp"
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "Sphere.h"
#include "Mesh.h"
#include "Sorghum.h"

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
#define useArcball true


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
			
			while (ColorTable::colors.size() <= component)
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
		float stemScale;
		float branchScale;
		//scale of 0-1 between min and max for whatever attribute is being colorized
		//lower/higher values are flooded
		float minColorCutoff, maxColorCutoff;
		ColorizationOptions colorization;
		std::vector<std::vector<GLfloat>> heatmap;
		float flatSelectionColor[4];
		float traitSelectionColor[4];
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

		int connectedPrimaryNode;
		

		BMetaNode();
		BMetaNode(SkelVert srcId, BSkeleton *skel);
		void updateColors(NodeVisualizationOptions options, bool isUpdateComponentColor = 1);
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
		float averageThickness, averageWidth, mLength;
		bool isBridge;
		bool isSelected;
		int instanceId;
		static int instanceCounter;
		float curvature;
		int connectedComponent;

		int connectedPrimaryNode;
		

		//whenever we remove a meta edge, we will just completely rebuild the indices lists
		//this is the list of this edges indices for each line referencing the vertex indices in 
		//a glvertices list
		std::vector<GLuint> indicesList;
		//GLfloat *glVertices;

		std::vector<GLfloat> localThicknessColors, localWidthColors, localRatioColors, localComponentColors;


		BMetaEdge();
		BMetaEdge(std::vector<SkelVert> vertices, BSkeleton* srcSkeleton);
		BMetaEdge join(BMetaEdge &other, BSkeleton* srcSkeleton);
		std::pair<BMetaEdge, BMetaEdge> split(SkelVert toSplitOn, BSkeleton *srcSkeleton);
		float getAvgThickness();
		SkelVert start();
		SkelVert end();
		void addToIndicesList(std::vector<GLuint> &edgeIndices);

		void updateColors(EdgeVisualizationOptions options, std::vector<std::vector<GLfloat>> &vertexColors, BSkeleton* srcSkeleton);
		void updateComponentColor(std::vector<std::vector<GLfloat>> &vertexColors);
		void select(GLfloat *selectionColor, std::vector<std::vector<GLfloat>> &vertexColors);
		void unselect(std::vector<std::vector<GLfloat>> &vertexColors);
		void highLight(GLfloat *selectionColor, std::vector<std::vector<GLfloat>> &vertexColors);
		void unhighLigh(std::vector<std::vector<GLfloat>> &vertexColors);
		void updateGraphColors(std::vector<std::vector<GLfloat>> &vertexColors);
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

	struct MetaFace
	{
		std::set<int> faceIndices;
		std::vector<GLuint> vertices;
		Point3d center;

		MetaFace(std::set<int> memberFaces, std::vector<Face>& skelFaces);

		static std::vector<MetaFace> findMetaFaces(std::vector<Face> &allFaces);
	};
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
		Sorghum sorghum;
		BSkeleton mSkeleton;
		
		//map between skeleton vertices and their corresponding metanodes
		std::map<SkelVert, MetaV> vertNodeMap;

		//vector correlating each edge to a connected component
		std::vector<int> mComponentMap;
		std::vector<float> mComponentSizeMap;
		std::vector<MetaE> mMinimumSpanningTree;
		std::vector<BoundingBox> componentBounds;
		std::vector<MetaFace> faces;
		int numComponents;
		

		//map backwards from the duplicated element to the src element
		//these maps are only to be used for writing out the files
		std::map<MetaE, MetaE> mDuplicateEdgeMap;
		std::map<MetaV, MetaV> mDuplicateNodeMap;
		bool isLoaded;


		//visualization and drawing options
		bool displayFaces;
		EdgeVisualizationOptions edgeOptions;
		NodeVisualizationOptions nodeOptions;
		OperationMode mode;
		static drawing::VBOSphere drawSphere;
		static drawing::VBOCube drawCube;

		//graphics objects and values
		std::vector<std::vector<GLfloat>> vertexColors;
		std::vector<GLfloat> stemColors;
		std::vector<GLuint> bridgeVBO;
		std::vector<GLuint> nonBridgeVBO;
		std::vector<GLuint> selectionVBO;

		std::vector<GLuint> selectedSegmentVBO;
		std::vector<GLuint> primaryBranchesVBO;
		std::vector<GLuint> stemVBO;

		std::vector<GLuint> autoStemVBO;
		std::vector<GLuint> sorghumBranchVBO;
		std::vector<GLuint> testVBO; // use for debug
		GLfloat selectionColor[4];
		float eyeShiftX, eyeShiftY;
		Mesh alphaMesh;
		Point3d viewCenter;

		bool displayMesh;

		float projectionTransform[16];
		float modelViewTransform[16];

		//editing options
		bool onlyDisplaySelectedComponents;
		int selectedComponent1, selectedComponent2;
		bool showBoundingBoxes;

		MetaV selectNode1, selectNode2;
		bool selectNode1Valid, selectNode2Valid;

		MetaE breakEdge;
		bool breakEdgeValid;

		MetaE splitEdge;
		bool splitEdgeValid;
		std::vector<MetaE> splitNeighbors;

		MetaE removeComponentEdge;
		bool removeComponentEdgeValid;
		
		int getNumGLEdges();
		int getNumGLVertices();

		// find traits
		MetaV nodeToView;
		bool viewNodeInfoValid;

		bool showSuggestedStem;
		bool showStem;
		bool showBranch;
		bool showTest = false;
		std::vector<MetaE> StemPath;
		std::vector<MetaV> StemPath_node;
		std::vector<SkelEdge> auto_stem;

		MetaV selectStemStart, selectStemEnd;
		bool selectStemStartValid, selectStemEndValid;
		bool stemSelected;
		float rootClipPlaneNormal[2][4] = {0};
		
		bool showSuggestedNode;
		std::vector<BoundingBox> auto_node;
		bool showPrimaryNodes;
		std::vector< std::pair<int, MetaV> > PrimaryNodes;
		bool selectPrimaryNodesValid;
		struct PrimaryNodesPredecessorsStruct
		{
			int nodeID;
			std::vector<MetaV> predecessors;
		};
		std::vector<PrimaryNodesPredecessorsStruct> NodesPredecessors;

		int CurrentPrimaryNode;
		std::vector<MetaV> CurrentPredecessors;
		GLfloat selectedPrimaryNodeColor[4] = { 0.0, 1.0, 0.0, 1.0 };
		GLfloat YELLOW[4] = { 1.0, 1.0, 0.0, 1.0 };
		GLfloat RED[4] = { 1.0, 0.0, 0.0, 1.0 };
		GLfloat GREEN[4] = { 0.0, 1.0, 0.0, 1.0 };
		GLfloat BLUE[4] = { 0.0, 0.0, 1.0, 1.0 };
		std::vector<std::vector<GLfloat>> randomColorLoopUpTable{
			{0.824, 0.961, 0.235, 1.0},
			{0.569, 0.118, 0.706, 1.0},
			{0.275, 0.941, 0.941, 1.0},
			{0.941, 0.196, 0.902, 1.0},
			{1.0, 0.882, 0.098, 1.0},
			{0.98, 0.745, 0.745, 1.0},
			{0.0, 0.51, 0.784, 1.0},
			{0.51, 0.961, 0.188, 1.0},
			{0.0, 0.502, 0.502, 1.0},
			{0.902, 0.745, 1.0, 1.0},
			{0.667f, 0.431f, 0.157f, 1.0f},
			{1.0f, 0.98f, 0.784f, 1.0f},
			{0.502f, 0.0f, 0.0f, 1.0f},
			{0.667f, 1.0f, 0.765f, 1.0f},
			{0.502f, 0.502f, 0.0f, 1.0f},
			{1.0f, 0.843f, 0.706f, 1.0f},
			{0.0f, 0.0f, 0.502f, 1.0f},
			{0.502f, 0.502f, 0.502f, 1.0f},
			{0.902f, 0.098f, 0.294f, 1.0f},
			{0.235f, 0.706f, 0.294f, 1.0f},
			{1.0f, 0.882f, 0.098f, 1.0f},
			{0.0f, 0.51f, 0.784f, 1.0f},
			{0.961f, 0.51f, 0.188f, 1.0f},
			{0.569f, 0.118f, 0.706f, 1.0f},
			{0.275f, 0.941f, 0.941f, 1.0f},
			{0.941f, 0.196f, 0.902f, 1.0f},
			{0.824f, 0.961f, 0.235f, 1.0f},
			{0.98f, 0.745f, 0.745f, 1.0f},
			{0.0f, 0.502f, 0.502f, 1.0f},
			{0.902f, 0.745f, 1.0f, 1.0f},
			{0.667f, 0.431f, 0.157f, 1.0f},
			{1.0f, 0.98f, 0.784f, 1.0f},
			{0.502f, 0.0f, 0.0f, 1.0f},
			{0.667f, 1.0f, 0.765f, 1.0f},
			{0.502f, 0.502f, 0.0f, 1.0f},
			{1.0f, 0.843f, 0.706f, 1.0f},
			{0.0f, 0.0f, 0.502f, 1.0f},
			{0.502f, 0.502f, 0.502f, 1.0f} };
		
		
		bool PrimaryBranchSelectionValid;
		MetaV PrimaryBranchSelection;
		struct PrimaryBranchStruct
		{
			int primaryNodeIndex;
			MetaV branchEndNode;
			std::vector<MetaE> metaEdges;
		};std::vector<PrimaryBranchStruct> PrimaryBranchesObj;
		bool showConfirmedPrimaryBranches;
		bool showTraitsOnly;

		bool isTraitsLoad;
		bool isRandomColorizePrimaryNodes;
		bool showOnlyBranchesOfCurrentPrimaryNode;
		bool showSelectedSegment;
		std::vector<MetaV> TraitsNodes;

		MetaV selectSegmentPoint1, selectSegmentPoint2;
		bool selectSegmentPoint1Valid, selectSegmentPoint2Valid;
		std::vector<MetaE> SegmentMetaEdges;
		std::vector<MetaV> SegmentMetaNodes;
		float SegmentHorizontalRadius;
		std::vector<MetaE> SegmentPath;
		std::vector<MetaV> SegmentPath_node;
		struct SegmentNodesDistanceStruct
		{
			MetaV node;
			std::vector<float> distance;
		}; // distance map of each mode on the path of selection
		std::vector<SegmentNodesDistanceStruct> SegmentNodesDistances;

		//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvSettingsvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//

		
		void assignEdgeHeatMap(boost::python::list heatmap);
		void showEdges(bool doShow);
		void setEdgeScale(float scale);
		void setStemScale(float scale);
		void setBranchScale(float scale);
		void magnifyNonBridges(bool doMagnify);
		void showOnlyNonBridges(bool showOnly);
		void colorizeEdgesByThickness();
		void colorizeEdgesByWidth();
		void colorizeEdgesByRatio();
		void colorizeEdgesByComponent();
		void setEdgeColorFloor(float val);
		void setEdgeColorCeiling(float val);
		void setEdgeSelectionColor(float r, float g, float b);
		void setCurrentPrimaryNodeSelectionColor(float r, float g, float b);
		void showMesh(bool doShow);
		void setMeshAlpha(float alpha);
		void setMeshColor(float red, float green, float blue);

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

		void setDisplayStem(bool doShowStem);
		void setDisplayBranch(bool doShowBranch);
		void setDisplaySuggestedStem(bool doShow);
		void setDisplayPrimaryNodes(bool doShowPriamryStem);
		void setDisplaySuggestedNode(bool doShow);
		void setCurrentPrimaryNode(int layer);
		void setDisplayConfirmedPrimaryBranches(bool doShowConfirmedBranches);
		void setRandomColorizePrimaryNodes(bool doShow);
		void setDisplayOnlyBranchesOfCurrentPrimaryNode(bool doShow);
		void setDisplayTraitsOnly(bool doShow);
		void setDisplaySelectedSegment(bool doShow);
		void setSegmentHorizontalSliderRadius(float val);
		
		void drawEdges();
		void edgePickRender();
		void drawNodes();
		void nodePickRender();
		void vertPickRender();
		void drawBoxes();
		void draw();

		void startRotation(int mx, int my);
		void mouseMoved(int mx, int my, float zoom);
		void setZoom(float rad, float eyex, float eyey, float eyez, float upx, float upy, float upz);
		void writeToBinary(std::string filename);

	private:
		float DistanceFromPointToLine(MetaV point, MetaV line1, MetaV line2);
		float MinDistanceToSelectedSegment(MetaV point);
		float* Substract3DPoint(float *point0, float *point1);

		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Settings^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

		//vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvInteractionsvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv//
	public:
		void selectConnectionNode(int mouseX, int mouseY);
		void selectBreakEdge(int mouseX, int mouseY);
		void selectSplitEdge(int mouseX, int mouseY);
		void selectRemoveComponentEdge(int mouseX, int mouseY);
		void selectStemStartEnd(int mouseX, int mouseY);
		void selectStemPrimaryNode(int mouseX, int mouseY);
		void selectPrimaryBranches(int mouseX, int mouseY);
		void selectSegmentPointAction(int mouseX, int mouseY);
		void viewNodeInfoAction(int mouseX, int mouseY);

		void setSelectionColor(float r, float g, float b);
		void unselectAll();
		void shiftEye(float xShift, float yShift);

		bool pickNewViewCenter(int mouseX, int mouseY);
		void changeRotationSpeed(bool increase);
		void runSorghumAlgorithm();


	private:
		MetaV getFirstMetaNodeHit(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ, bool &isValid);
		SkelVert getFirstSkelNodeHit(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ, bool &isValid);
		std::pair<MetaV, MetaV> getFirstMetaEdgeHit(float eyeX, float eyeY, float eyeZ, float lookX, float lookY, float lookZ, bool &isValid);

		MetaE selectEdgeByRender(int mouseX, int mouseY, bool &isValid);
		MetaV selectNodeByRender(int mouseX, int mouseY, bool &isValid);
		SkelVert selectVertByRender(int mouseX, int mouseY, bool &isValid);

		void privateSelectConnectionNode(SkelVert nodeVert, int selectComponent);
		void privateSelectStemStartEnd(SkelVert nodeVert, int selectComponent);
		void privateSelectSegmentPointAction(SkelVert nodeVert, int selectComponent);
		void privateSelectPrimaryBranches(MetaV target, int selectComponent);
		std::vector<MetaE> privateShortestPath(MetaV source, MetaV target);
		void privateSelectNodeToView(SkelVert nodeVert, int selectComponent);
		void unselectAllTraits();
		void unselectAllEdges();

		void unselectEdge(std::pair<MetaV, MetaV> toUnselect);
		void unselectEdge(MetaE toUnselect);
		void selectEdge(std::pair<MetaV, MetaV> toSelect);
		void selectEdge(MetaE toSelect);
		void highLightEdge(MetaE edge);
		void unhighLightEdge(MetaE edge);
		void colorEdge(MetaE edge, GLfloat *color);
		void uncolorEdge(MetaE edge);
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^Interactions^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^//

	public:

		boost::python::list getComponentSizes();

		int getNumComponents()
		{
			return numComponents;
		}

		int getNumPrimaryNodes()
		{
			return PrimaryNodes.size();
		}

		bool isDisplaySelectedSegment();

		/*vvvvvvvvvvvvvvvvvvvvv CONSTRUCTORS vvvvvvvvvvvvvvvvvvvvv*/

		/*default constructor.  Only constructor defined.  All operations such as loading or
		adding edges should be done with appropriate member functions*/
		BMetaGraph();
		void Initialize();
		void InitializeTraitParameters();

		/*vvvvvvvvvvvvvvvvvvvvv LOADERS AND SAVERS vvvvvvvvvvvvvvvvvvvvv*/

		void loadFromFile(std::string filename);
		void loadMeshFromFile(std::string filename);
		//void loadFromFile(boost::python::str filename);
		/*
		External load operation, this is the only one that should be used outside this class.
		Aliases loading operations for different file structures internally.
		
		Parameters: 
		lines - Set of line strings defining the skeleton
		startingLine - Index of first line from provided set to load from

		returns - line number after the last line consumed for the loading process
		*/
		int loadFromLines(std::vector<std::string> &lines, int startingLine);

		void saveToFile(std::string filename);
		/*
		External save operation, only save functionality for this class.
		Inserts this metagraph as a string with a ply-style header into the provided ostream.

		Parameters:
		out - ostream to push this skelton into
		*/
		void writeToStream(std::ostream &out);

		/*
		Traits 
		save operation
		*/
		void loadTraitsFromFile(std::string filename);
		int loadTraitsFromLines(std::vector<std::string> &lines, int startingLine);

		void saveTraitsToFile(std::string filename);


		int loadFromLines_Binary(FILE *fp);
		//void loadFromFile_Binary(std::string filename);

		int loadSkeletonFromLines_Binary(FILE *fp);

	private:
		/*
		Load just a skeleton, no metagraph information
		*/
		int loadSkeletonFromLines(std::vector<std::string> &lines, int &startingLine);

		int loadGraphFromLines(std::vector<std::string> &lines, int &metaLinesStart);

		int loadTraitsPly(std::vector<std::string> &lines, int startingLine);
		void loadStem(std::vector<std::string> &lines, int &lineOn);
		void loadPrimaryNodes(std::vector<std::string> &lines, int &lineOn, int numPrimaryNodes);
		void loadPrimaryBranches(std::vector<std::string> &lines, int &lineOn, int numPrimaryBranches);

		bool checkHeaderForMetaInfo(std::vector<std::string> &lines);

		void buildEdgeVBOs();

		


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
		MetaV addNode(SkelVert srcVert, BSkeleton *srcSkel);
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
		void setSorghumBranchParameters(int minBranchSize, int maxBranchSize, float radiusTolerance, float tipAngleThresh, float tortuosityThresh);
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
		void PromoteOperation(SkelVert toPromote);
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

		/*
		RemoveComponentOperation() will support the component removal operation. User can select a edge, then remove entire component
		*/
		void RemoveComponentOperation();

		/*
		This method is the first step to build root hierarchy. User click on two end point of the stem.
		*/
		void SelectStemOperation();
		void SelectStemPrimaryNodeOperation();
		void SelectPrimaryBranchesOperation();
		void RemovePrimaryBranchesOperation();
		void SelectSegmentPointOperation();

		/*
		Automatic find stem. 
		Require input: lower threshold
		The algorithm first find a set of edges which has radius larger than the input threshold.
		Then find the minima spanning tree out of the selected edges. Finally, find the longest path
		within the MST
		*/
		void FindStemOperation(float lowThreshold);
		void sorghumBranchOperation();
		float getEdgeEuclidLength(SkelEdge srcId, BSkeleton *skel);
		float getVertThickness(SkelVert srcId, BSkeleton *skel);
		float getVertWidth(SkelVert srcId, BSkeleton *skel);

		/*
		Automatically find priamry nodes
		*/
		void FindPrimaryNodeOperation(float look_distance, float kernel_bandwidth);
		std::vector<float> neighbourhoodPoints(std::vector<float> positions, float x_centroid, float distance = 5);
		float gaussian_kernal(float distance, int bandwidth);

	private:

		/*
		Determines if the provided metanode is of degree 2, and if so, joines the edges leading
		into this node, creates a new edge in the graph bridging its two adjacent nodes, and 
		removing this node and its incident edges from the graph.
		*/
		bool BridgeNode(SkelVert vertToBridge);
		
		void SetSelectSegmentPointOperation();
		//void BridgeNodeSweeper();
	public:
		/*
		Use the vert node map to determine the connected component to which an edge belongs
		*/
		int GetEdgeComponent(MetaE edge);
	};

	struct DisjointSets
	{
		int *parent, *rnk;
		int n;

		// Constructor. 
		DisjointSets(int n);
		// Find the parent of a node 'u' 
		// Path Compression 
		int find(int u);
		// Union by rank 
		void merge(int x, int y);
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