#pragma once

#include "Geometry.h"
#include "RootAttributes.h"


#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/config.hpp"
#include "boost/graph/undirected_graph.hpp"
#include "boost/graph/lookup_edge.hpp"
#include <stdlib.h>
#include <GLFW/glfw3.h>
#include <stdio.h>

namespace
{
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Point3d, Roots::RootAttributes, boost::no_property, boost::listS> BoostSkeleton;
	typedef boost::graph_traits<BoostSkeleton>::vertex_descriptor SkelVert;
	typedef boost::graph_traits<BoostSkeleton>::edge_descriptor SkelEdge;
}

typedef boost::graph_traits<BoostSkeleton>::vertex_iterator skelVIter;
typedef boost::graph_traits<BoostSkeleton>::edge_iterator skelEIter;

struct skelVertIter : std::pair<skelVIter, skelVIter>
{
	skelVertIter(const std::pair<skelVIter, skelVIter> &other)
		:std::pair<skelVIter, skelVIter>(other)
	{
	}
	skelVertIter operator++()
	{
		++this->first;
		return *this;
	}
	skelVertIter operator--()
	{
		--this->second;
		return *this;
	}
};

struct skelEdgeIter : std::pair<skelEIter, skelEIter>
{
	skelEdgeIter(const std::pair<skelEIter, skelEIter> &other)
		:std::pair<skelEIter, skelEIter>(other)
	{
	}
	skelEdgeIter operator++()
	{
		++this->first;
		return *this;
	}

	skelEdgeIter operator--()
	{
		--this->second;
		return *this;
	}
};



namespace Roots
{

	struct BSkeleton : public BoostSkeleton
	{
		/*vvvvvvvvvvvvvvvvvvvvv MEMBER VARIABLES vvvvvvvvvvvvvvvvvvvvv*/
		
		/*3-axis bounds of the graph*/
		float mMinX, mMaxX, mMinY, mMaxY, mMinZ, mMaxZ;
		/*indicator whether 3-axis bounds have been determined*/
		bool mBoundsFound;
		/*Geometric center of the graph defined as the 3-axis midpoint of the bounds*/
		Point3d mCenter;
		/*Geometric radius of the bounding sphere defined by the */
		float mRadius;

		

		static std::string beginSkeletonString;
		static std::string endSkeletonString;
		static std::string vertexString;
		static std::string edgeString;
		static std::string endHeaderString;
		static std::string beginPlyString;

		/*vvvvvvvvvvvvvvvvvvvvv CONSTRUCTORS vvvvvvvvvvvvvvvvvvvvv*/
		
		/*default constructor.  Only constructor defined.  All operations such as loading or 
		adding edges should be done with appropriate member functions*/
		BSkeleton();

		void Initialize();

		/*vvvvvvvvvvvvvvvvvvvvv LOADERS AND SAVERS vvvvvvvvvvvvvvvvvvvvv*/

		/*
		External load operation, this is the only one that should be used outside this class.
		Aliases loading operations for different file structures internally.
		
		Parameters: 
		lines - Set of line strings defining the skeleton
		startingLine - Index of first line from provided set to load from

		returns - line number after the last line consumed for the loading process
		*/
		int loadFromLines(std::vector<std::string> lines, int startingLine, std::vector<GLfloat> &glVertices);

		/*
		External save operation, only save functionality for this class.
		Inserts this skeleton as a string with a ply-style header into the provided ostream.

		Parameters:
		out - ostream to push this skelton into
		*/
		void writeToStream(std::ostream &out);

	private:
		/*
		Internal load operation for Wenzhen style skeletons.

		Parameters and return same as LoadFromLines
		*/
		int loadWenzhenLines(std::vector<std::string> lines, std::vector<GLfloat> &glVertices, int startingLine = 0);

		int loadDanPly(std::vector<std::string> lines, std::vector<GLfloat> &glVertices, int startingLine = 0);

		/*
		Internal load operation for Ply-style skeletons.

		Parameters and returns same as LoadFromLines
		*/
		int loadPlyStyleLines(std::vector<std::string> lines, std::vector<GLfloat> &glVertices, int startingLing = 0);

		/*Simple callout function to load all vertices described by the provided set of lines.*/
		void loadVertices(std::vector<std::string> lines, int &lineOn, int numVertices);
		
		/*Simple callout function to load all edges described by the provided set of lines.*/
		void loadEdges(std::vector<std::string> lines, int &lineOn, int numEdges, std::vector<GLfloat> &glVertices);


	public:
		/*
		alias that both adds the connection between v0 and v1 in the boost graph,
		and the decoration edge_descriptor for the root attributes
		*/
		SkelEdge addEdge(int v0, int v1, RootAttributes attributes, std::vector<GLfloat> &glVertices);

		/*
		alias that both adds a new vertex to the boost graph, and adds the decoration
		vertex_descriptor for the vertex 3d location
		*/
		SkelVert addVertex(Point3d pointLocation);


		/*vvvvvvvvvvvvvvvvvvvvv GEOMETRIC MANIPULATIONS vvvvvvvvvvvvvvvvvvvvv*/
	public:
		
		/*
		Finds the euclidean bounds of the skeleton and stores them in the class variables
		*/
		void findBounds();

		/*
		Finds the bounding sphere (if it hasn't been determined) and stores the results in the
		class parameters
		*/
		void findBoundingSphere();

		/*
		Returns a copy of this skeleton graph, but with the center of its bounding sphere at the 
		provided point.

		Parameters:
		newCenter - the desired centerpoint of the result graph
		*/
		BSkeleton recenterSkeleton(Point3d newCenter, std::vector<float> &vertices);

	


		/*vvvvvvvvvvvvvvvvvvvvv GENERIC OPERATIONS vvvvvvvvvvvvvvvvvvvvv*/

		/*
		Finds the edge_iterator for the provided vertex pair (aliases boost::edge(v0, v1, skeleton))
		*/
		SkelEdge getEdge(int v0, int v1, bool &exists);

		RootAttributes& getEdgeData(int v0, int v1, bool &exists);

		RootAttributes& getEdgeData(SkelEdge e);

		Point3d& getVertData(SkelVert v);

		/*
		Finds all vertices of degree != 2
		*/
		std::vector<SkelVert> GetNotableVertices();

		void normalizeEdgeAttributes();
	};
}

struct PySkeleton
{
	float radius;
	Point3d center;
	float minThickness, maxThickness, minWidth, maxWidth, minRatio, maxRatio;

	boost::python::list thicknessPercentiles, widthPercentiles, ratioPercentiles;

	boost::python::list mVertexList;
	boost::python::list mEdgeList;


	PySkeleton();
	PySkeleton(Roots::BSkeleton *srcSkel);
	void findBoundingSphere();
	void reload();

private:
	Roots::BSkeleton *mSkeleton;
};