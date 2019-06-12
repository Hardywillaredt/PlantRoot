#pragma once

#include "Geometry.h"
#include "RootAttributes.h"
#include "SkeletonEdge.h"


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <sstream>
#include <map>
#include <vector>

#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/config.hpp"
#include "boost/graph/undirected_graph.hpp"


namespace Roots
{

	typedef boost::adjacency_list<boost::listS, boost::vecS, boost::undirectedS, Point3d, RootAttributes> BoostSkeleton;
	typedef boost::graph_traits<BoostSkeleton>::vertex_descriptor SkelVert;
	typedef boost::graph_traits<BoostSkeleton>::edge_descriptor SkelEdge;

	typedef std::vector<SkeletonEdge> edgeList;
	typedef edgeList::iterator edgeIter;

	typedef std::vector<SkeletonEdge*> edgePtrList;
	typedef edgePtrList::iterator edgePtrIter;

	typedef std::vector<Point3d> vertList;
	typedef vertList::iterator vertIter;




	class Skeleton
	{
	public:
		//Serializable Implementation >> //
		//Skeleton(Json::Value json);
		//Json::Value ToJson();

		static std::string beginSkeletonString;
		static std::string endSkeletonString;
		static std::string vertexString;
		static std::string edgeString;
		static std::string endHeaderString ;

		friend std::ostream& operator<<(std::ostream& out, const Skeleton &skel);
		friend std::istream& operator>>(std::istream& in, Skeleton &skel);

		// << Serializable Implementation//

		Skeleton();
		Skeleton(std::string filename);
		Skeleton(vertList aVertices, std::vector<edgeList> aEdges);

		static Skeleton LoadFromTextFile(std::string filename);
		void LoadWenzhenStyleFile(std::string filename);
		void LoadWenzhenStyleFile(std::istream &in, int numVerts, int numEdges, int numFaces);
		//void LoadFromJson(Json::Value skelJson);
		//void LoadFromJsonFile(char *filename);
		

		vertList getVertices();
		std::vector<edgeList> getEdges();


		void AddEdge(SkeletonEdge toAdd);

		void RemoveEdges(edgeList toRemove, bool careAboutAttributes = false);
		void RemoveEdge(SkeletonEdge toRemove, bool careAboutAttributes = false);
		void RemoveEdge(int v0, int v1);

		SkeletonEdge* GetEdge(int v0, int v1);

		void GetBounds(float &leftX, float &rightX, float &botY, float &topY, float &backZ, float &frontZ);
		
		void GetBoundingSphere(float &cx, float &cy, float &cz, float &r);

		void ResetCenter();
		void MoveCenterTo(Point3d targetCenter);
		//static Skeleton FromJson();


	private:

		//vertex list, potentially sorted in some way, but for now, nada
		vertList mVerts;
		//list of edges sorted by the id of the first vertex, then the second
		std::vector<edgeList> mEdges;
		float mLeftX, mRightX, mBotY, mTopY, mBackZ, mFrontZ;
		bool boundsFound;
		
		void FindBounds();

		//this is just convenience information, no need to store or load it
		Point3d mOriginalCenter;
		Point3d mCurrentCenter;
		bool originalCenterSet;
	public:
		std::vector<std::vector<int>> mNeighbors;
		int mNumVertices;
		int mNumEdges;
	};

	/*
	This function will take a single root (edgePtrList) and order that list of edgePtrs from beginning to end.  If no edge with a degree 1 vertex is provided as
	the forward end, then the returned edge-list will be ordered by the lower degree 1 vertex id.  Will return the 
	*/
	edgePtrList tidyUpRoot(edgePtrList toTidy, bool &valid, int &end1, int &end2, SkeletonEdge *beginning = nullptr);

}
