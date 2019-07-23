#pragma once

#include "Geometry.h"
#include "RootAttributes.h"
#include <map>
#include "Face.h"
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
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Point3d, Roots::RootAttributes, boost::no_property, boost::vecS> BoostSkeleton;
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

enum ParsingOrder
{
	X = 0,
	Y = 1,
	Z = 2,
	Thickness = 3,
	Width = 4,
	ParsingCount = 5
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
		Point3d originalCenter;
		/*Geometric radius of the bounding sphere defined by the */
		float mRadius;
		std::map<ParsingOrder, int> mVertexParseOrder;
		std::map<int, ParsingOrder> mVertexWriteOrder;

		std::vector<GLfloat> glVertices;
		std::vector<Face> faces;

		

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
		int loadFromLines(std::vector<std::string> &lines, int startingLine);

		/*
		External save operation, only save functionality for this class.
		Inserts this skeleton as a string with a ply-style header into the provided ostream.

		Parameters:
		out - ostream to push this skelton into
		*/
		void writeToStream(std::ostream &out);

		void writeDanPly(std::ostream &out);

		void writeToBinary(std::string filename);

		void loadVertices_binary(FILE *fp);

		void loadEdges_binary(FILE *fp);

		int loadFromLines_Binary(FILE *fp);



	private:
		/*
		Internal load operation for Wenzhen style skeletons.

		Parameters and return same as LoadFromLines
		*/
		int loadWenzhenLines(std::vector<std::string> &lines, int startingLine = 0);

		int loadDanPly(std::vector<std::string> &lines, int startingLine = 0);

		/*
		Internal load operation for Ply-style skeletons.

		Parameters and returns same as LoadFromLines
		*/
		int loadPlyStyleLines(std::vector<std::string> &lines, int startingLing = 0);

		/*Simple callout function to load all vertices described by the provided set of lines.*/
		void loadVertices(std::vector<std::string> &lines, int &lineOn, int numVertices);
		
		/*Simple callout function to load all edges described by the provided set of lines.*/
		void loadEdges(std::vector<std::string> &lines, int &lineOn, int numEdges);

		void loadFaces(std::vector<std::string> &lines, int &lineOn, int numFaces);


	public:
		/*
		alias that both adds the connection between v0 and v1 in the boost graph,
		and the decoration edge_descriptor for the root attributes
		*/
		inline SkelEdge addEdge(int v0, int v1) {

			SkelEdge e;
			bool edgeAdded;
			boost::tie(e, edgeAdded) = boost::add_edge(v0, v1, *this);
			if (edgeAdded)
			{
				RootAttributes ra = RootAttributes();
				ra.euclidLength = (operator[](v0) - operator[](v1)).mag();
				ra.v0id = v0;
				ra.v1id = v1;
				operator[](e) = ra;

			}
			return e;
		
		}

		/*
		alias that both adds a new vertex to the boost graph, and adds the decoration
		vertex_descriptor for the vertex 3d location
		*/
		inline SkelVert addVertex(Point3d pointLocation) {

			SkelVert v;
			v = boost::add_vertex(*this);
			operator[](v) = pointLocation;
			operator[](v).id = v;
			mBoundsFound = false;
			return v;
		}
		
		void updateGLVertices();

		/*vvvvvvvvvvvvvvvvvvvvv GEOMETRIC MANIPULATIONS vvvvvvvvvvvvvvvvvvvvv*/
	public:
		
		/*
		Finds the 
		bounds of the skeleton and stores them in the class variables
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
		void recenterSkeleton(Point3d newCenter);

	


		/*vvvvvvvvvvvvvvvvvvvvv GENERIC OPERATIONS vvvvvvvvvvvvvvvvvvvvv*/

		/*
		Finds the edge_iterator for the provided vertex pair (aliases boost::edge(v0, v1, skeleton))
		*/
		// inline for these four functions
		inline SkelEdge getEdge(int v0, int v1, bool &exists) {
			SkelEdge result;
			boost::tie(result, exists) = boost::edge(v0, v1, *this);
			return result;
		}

		 inline RootAttributes& getEdgeData(int v0, int v1, bool &exists)
		 {
			 SkelEdge temp = getEdge(v0, v1, exists);
			 if (exists)
			 {
				 return getEdgeData(temp);
			 }
			 else
			 {
				 return RootAttributes();
			 }
		 }

		 inline RootAttributes& getEdgeData(SkelEdge e)
		 {
			 return operator[](e);
		 }

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