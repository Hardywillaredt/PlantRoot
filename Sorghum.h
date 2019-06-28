#pragma once
#define USE_MATH_DEFINES
#include <GLFW/glfw3.h>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <math.h>
#include <algorithm>
#include <fstream>
#include <ostream>
#include <algorithm>
#include <queue>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <cmath>

using namespace std;

struct Internal {
	struct Vertex { float x, y, z; /*unsigned int id; */};
	struct Edge { int v1, v2; /*unsigned int id;*/ };
	struct Face {
		unsigned char nvts;
		int verts[3];
		float s;
	};
};
struct VertexWithMsure
{
	float width, radius;
	float x, y, z;
	float compType;
	int inDiameter;
	//unsigned int id;

};
typedef std::pair<int, int> E;

class Sorghum {
public:
	std::vector<VertexWithMsure> vts;
	std::vector<Internal::Edge> edges;
	std::vector<Internal::Face> faces;
	float compThreshSize;
	int branchSizeMin;
	int curvatureWindow;
	int maxBranchLength;
	float radiusTol;
	float emergenceAngleThresh;
	float tipAngleThresh;
	float tortuosityThresh;
	float emergenceLowerRadius;
	float lowerStemThresh;
	int emergeWindow;
	string inFile;
	float upperRadiiThresh;
	float lowerRadiiThresh ;
	Sorghum() {
		compThreshSize = 200.0;
		branchSizeMin = 200;
		curvatureWindow = 30;
		maxBranchLength = 900;
		radiusTol = 1.2;
		emergenceAngleThresh = 360.0;
		tipAngleThresh = 120.0;
		tortuosityThresh = 2.5;
		emergenceLowerRadius = 2.0;
		lowerStemThresh = 1.0;
		emergeWindow = 30;
		upperRadiiThresh = 1.0;
		lowerRadiiThresh = 0.1;
	}
	void fillComponentIterative(std::map<int, int> &visitedMap, int i, std::map<int, std::vector<int> > vertEdgeMapping, std::vector<Internal::Edge> edges, std::vector<int> vertComps, std::vector<int> &newEdgeComp, std::vector<int> &newVertComp, float lowerRadiiThresh, vector<VertexWithMsure> vts);
	float euclideanDistance(VertexWithMsure v1, VertexWithMsure v2);
	float getLengthOfComponent(std::vector<int> compEdges, std::vector<Internal::Edge> edges, std::vector<VertexWithMsure> vertices);
	float gaussianFactor(float val, float factor);
	float getRangeScore(vector< vector<E> >& mstVertexEdges, int seedVertex, vector<VertexWithMsure>& allVertices, map<int, bool>& vertexInDiameter, int range);
	double dotProd(vector<float> v1, vector<float> v2);
	vector<float> vecDiff(VertexWithMsure v1, VertexWithMsure v2);
	float getMagnitude(vector<float> vec);
	float angleBetweenVector(vector<float> v1, vector<float> v2);
	vector<float> crossProduct(vector<float> v1, vector<float> v2);
	void normalize(vector<float>& vec);
	vector<vector<int> > buildBranchesFromVertexForwardBFS(int& vertexSeed, vector<bool>& inDiameter, std::vector<Internal::Edge>& edges, vector< std::vector<int> >& vertIndexEdgeIndices, vector<VertexWithMsure>& stemVertices, vector<int>& visited, int& branchSizeMin, int& curvatureWindow, int& maxBranchLen,
		vector<float>& stemCentroid, vector<int>& loop, float& radiusTolerance, float& emergenceAngleThresh, float& tipAngleThresh, float& tortuosityThresh, vector<float>& stemVec, vector<float>& emergenceAngles, vector<float>& tipAngles, vector<float>& branchLengths, vector<float>& tipLengths, float& emergenceLowerRadius, int& emergeWindow);
	int sorghumAlgorithm(std::vector<GLuint> &sorghumBranchVBO);
	void getFilename(string filename);
};
