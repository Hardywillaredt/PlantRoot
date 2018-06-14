#include "Skeleton.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>

namespace
{
	std::string edgeString = "Edges";
	std::string vertexString = "Vertices";
}

std::string Roots::Skeleton::beginSkeletonString = "begin_skeleton";
std::string Roots::Skeleton::endSkeletonString = "end_skeleton";
std::string Roots::Skeleton::vertexString = "vertex";
std::string Roots::Skeleton::edgeString = "edge";
std::string Roots::Skeleton::endHeaderString = "end_header";

namespace Roots
{

	//Skeleton::Skeleton(Json::Value json)
	//{
	//	int numVertices = json[vertexString].size();
	//	int numEdges = json[edgeString].size();

	//	mNumEdges = numEdges;
	//	mNumVertices = numVertices;

	//	mVerts = vertList(numVertices, Point3d());
	//	Json::Value vertJson = json[vertexString];
	//	for (int i = 0; i < numVertices; ++i)
	//	{

	//		mVerts[i] = Point3d(vertJson[0].asfloat(), vertJson[1].asfloat(), vertJson[2].asfloat());
	//	}

	//	Json::Value edgeJson = json[edgeString];
	//	mEdges = std::vector<edgePtrList>(numVertices, edgePtrList());
	//	for (int i = 0; i < numEdges; ++i)
	//	{
	//		SkeletonEdge edge = SkeletonEdge(edgeJson[i]);
	//		AddEdge(edge);
	//		//mEdges[edge.v0].push_back(new SkeletonEdge(edge));
	//	}

	//	for (int i = 0; i < mEdges.size(); ++i)
	//	{
	//		std::sort(mEdges[i].begin(), mEdges[i].end());
	//	}
	//}

	Skeleton::Skeleton()
		: mVerts(), mEdges(), mNumEdges(), mNumVertices(), mNeighbors() 
	{
		boundsFound = false;
		originalCenterSet = false;
	}


	Skeleton::Skeleton(vertList aVertices, std::vector<edgeList> aEdges)
		: mVerts(aVertices), mEdges(aEdges)
	{
		boundsFound = false;
		originalCenterSet = false;
		mNumEdges = 0;
		mNumVertices = mVerts.size();
		mNeighbors = std::vector<std::vector<int>>(mVerts.size());
		for (int i = 0; i < mEdges.size(); ++i)
		{
			mNumEdges += mEdges[i].size();
			for (int j = 0; j < mEdges[i].size(); ++j)
			{
				AddEdge(mEdges[i][j]);
			}
		}
		float cx, cy, cz, r;
		GetBoundingSphere(cx, cy, cz, r);
	}

	Skeleton::Skeleton(std::string filename)
	{
		boundsFound = false;
		originalCenterSet = false;
		std::ifstream in;
		in.open(filename);

		bool endHeaderReached = false;
		std::string line;
		int numEdges = 0, numVerts = 0;
		Log::out << "Istream operation for skeleton" << std::endl;
		std::getline(in, line);
		if (!boost::iequals(line, Skeleton::beginSkeletonString))
		{
			Log::out << "The skeleton file lacks the " << Skeleton::beginSkeletonString << "header. Should attempt to load wenzhen style file " << std::endl;
			std::vector<std::string> words = {};
			boost::split(words, line, boost::is_any_of(" "));
			int numVerts, numEdges, numFaces;
			if (words.size() < 3)
			{
				if (words.size() == 2)
				{
					numVerts = boost::lexical_cast<int>(words[0]);
					numEdges = boost::lexical_cast<int>(words[1]);
				}
				else
				{
					Log::out << "The skeleton header lacks sufficient info to parse.  Ending read attempt" << std::endl;
					return;
				}
			}
			else
			{
				numVerts = boost::lexical_cast<int>(words[0]);
				numEdges = boost::lexical_cast<int>(words[1]);
				numFaces = boost::lexical_cast<int>(words[2]);
			}


			LoadWenzhenStyleFile(in, numVerts, numEdges, numFaces);

			return;
		}
		Log::out << "This is not being recognized as a wenzhen style file" << std::endl;
		Log::out << "The first input line is " << line << std::endl;
		std::vector<std::vector<std::string>> headerInfo;
		while (!endHeaderReached && !in.eof())
		{
			std::getline(in, line);
			std::vector<std::string> words = {};
			boost::split(words, line, boost::is_any_of(" "));
			if (boost::iequals(words[0], Skeleton::endHeaderString))
			{
				endHeaderReached = true;
				break;
			}
			headerInfo.push_back(words);
		}
		if (in.eof())
		{
			Log::out << "End of file reached before end of header info, improper file" << std::endl;
			return;
		}
		for (int i = 0; i < headerInfo.size(); ++i)
		{
			std::vector<std::string> headerLine = headerInfo[i];
			if (boost::iequals(headerLine[0], vertexString))
			{
				//handle vertices
				numVerts = boost::lexical_cast<int>(headerLine[1]);
				mVerts = vertList(numVerts);
				for (int i = 0; i < numVerts; ++i)
				{
					in >> mVerts[i];
				}
			}
			else if (boost::iequals(headerLine[0], edgeString))
			{
				//handle edges
				numEdges = boost::lexical_cast<int>(headerLine[1]);
				edgeList edges = edgeList(numEdges);
				for (int i = 0; i < numEdges; ++i)
				{
					in >> edges[i];
				}
				for each(SkeletonEdge edge in edges)
				{
					AddEdge(edge);
				}
			}
			else
			{
				//idk mang
				//do nothing
			}
		}
		float cx, cy, cz, r;
		GetBoundingSphere(cx, cy, cz, r);
	
	}

	Skeleton Skeleton::LoadFromTextFile(std::string filename)
	{
		

		std::ifstream in;
		in.open(filename);

		Skeleton result = Skeleton();

		Log::WriteLine("beginning istream for file");

		in >> result;

		Log::WriteLine("End of istream for skeleton file");

		float cx, cy, cz, r;
		result.GetBoundingSphere(cx, cy, cz, r);
		return result;
	}

	void Skeleton::LoadWenzhenStyleFile(std::string filename)
	{
		std::ifstream in;
		in.open(filename);
		int numVertices, numEdges, numFaces;
		in >> numVertices >> numEdges >> numFaces;

		LoadWenzhenStyleFile(in, numVertices, numEdges, numFaces);

	}

	void Skeleton::LoadWenzhenStyleFile(std::istream &in, int numVerts, int numEdges, int numFaces)
	{

		mNumEdges = numEdges;
		mNumVertices = numVerts;
		Log::out << "Wenzhen file - vertices: " << mNumVertices << " edges: " << numEdges << std::endl;
		mVerts = vertList(numVerts, Point3d());
		
		float x, y, z;
		for (int i = 0; i < numVerts; ++i)
		{
			in >> x >> y >> z;
			mVerts[i] = Point3d(x, y, z);
		}
		Log::out << "Successfully loaded all " << numVerts << " vertices" << std::endl;

		mEdges = std::vector<edgeList>(numVerts, edgeList());
		mNeighbors = std::vector<std::vector<int>>(mNumVertices, std::vector<int>());
		int v0, v1;
		std::vector<float> attributes = std::vector<float>(NumAttributes);
		for (int i = 0; i < numEdges; ++i)
		{
			//std::cout << "Adding edge" << std::endl;
			in >> v0 >> v1;
			for (int att = 0; att < NumAttributes; ++att)
			{
				in >> attributes[att];
			}
			SkeletonEdge toAdd = SkeletonEdge(v0, v1, attributes);
			AddEdge(toAdd);
			//std::cout << "Edge added " << std::endl;
		}

		Log::out << "Successfully loaded all " << numEdges << " edges" << std::endl;

		for (int i = 0; i < mEdges.size(); ++i)
		{
			std::sort(mEdges[i].begin(), mEdges[i].end());
		}

		if (!in.eof())
		{
			Log::out << "Not all data read from file, that is slightly strange" << std::endl;
		}
		float cx, cy, cz, r;
		GetBoundingSphere(cx, cy, cz, r);
		/////////////////////////////////// TO DO ->  DEAL WITH FACES   ///////////////////////////
	}


	vertList Skeleton::getVertices()
	{
		return mVerts;
	}

	std::vector<edgeList> Skeleton::getEdges()
	{
		return mEdges;
	}


	void Skeleton::AddEdge(SkeletonEdge toAdd)
	{
		
		mNeighbors[toAdd.v0].push_back(toAdd.v1);
		mNeighbors[toAdd.v1].push_back(toAdd.v0);
		if (toAdd.v0 > mEdges.size())
		{
			return;
		}
		if (mEdges[toAdd.v0].size() == 0)
		{
			SkeletonEdge addedEdge = SkeletonEdge(toAdd);
			mEdges[toAdd.v0].push_back(addedEdge);
		}
		else
		{
			for (int i = 0; i < mEdges[toAdd.v0].size(); ++i)
			{
				if ((toAdd.v1 < mEdges[toAdd.v0][i].v1))
				{
					mEdges[toAdd.v0].insert(mEdges[toAdd.v0].begin() + i, SkeletonEdge(toAdd));
					break;
				}
			}
		}
		
	}


	void Skeleton::RemoveEdges(edgeList toRemove, bool careAboutAttributes)
	{
		for each(SkeletonEdge edge in toRemove)
		{
			RemoveEdge(edge, careAboutAttributes);
		}
	}

	void Skeleton::RemoveEdge(int v0, int v1)
	{
		RemoveEdge(SkeletonEdge(v0, v1), false);
	}

	void Skeleton::RemoveEdge(SkeletonEdge toRemove, bool careAboutAttributes)
	{
		int v0 = toRemove.v0;
		int v1 = toRemove.v1;
		for (std::vector<int>::iterator iter = mNeighbors[v0].begin(); iter < mNeighbors[v0].end(); ++iter)
		{
			if (iter[0] == v1)
			{
				mNeighbors[v0].erase(iter, iter + 1);
			}
		}
		for (std::vector<int>::iterator iter = mNeighbors[v1].begin(); iter < mNeighbors[v1].end(); ++iter)
		{
			if (iter[0] == v0)
			{
				mNeighbors[v1].erase(iter, iter + 1);
			}
			
		}

		bool edgeFound = false;
		if (careAboutAttributes)
		{
			for (int i = 0; i < mEdges[v0].size(); ++i)
			{
				if (mEdges[v0][i]==(toRemove))
				{
					mEdges[v0].erase(mEdges[v0].begin() + i, mEdges[v0].begin() + i + 1);
					edgeFound = true;
				}
			}
		}
		else
		{
			for (int i = 0; i < mEdges[v0].size(); ++i)
			{
				if (mEdges[v0][i].v1 == v1)
				{
					mEdges[v0].erase(mEdges[v0].begin() + i, mEdges[v0].begin() + i + 1);
					edgeFound = true;
				}
			}
		}
		if (edgeFound)
		{
			mNumEdges -= 1;
		}
	}


	SkeletonEdge* Skeleton::GetEdge(int vert0, int vert1)
	{
		int v0 = std::min(vert0, vert1);
		int v1 = std::max(vert0, vert1);

		for each(SkeletonEdge edge in mEdges[v0])
		{
			if(edge.v1 == v0)
			{
				return &edge;
			}
		}
		return nullptr;
	}


	void Skeleton::GetBounds(float &leftX, float &rightX, float &botY, float &topY, float &backZ, float &frontZ)
	{
		FindBounds();
		leftX = mLeftX;
		rightX = mRightX;
		botY = mBotY;
		topY = mTopY;
		backZ = mBackZ;
		frontZ = mFrontZ;
		return;

	}

	void Skeleton::FindBounds()
	{
		if (boundsFound)
		{
			return;
		}
		else
		{
			float big = 99999999999;
			mLeftX = big;
			mRightX = -big;
			mBotY = big;
			mTopY = -big;
			mBackZ = big;
			mFrontZ = -big;

			float x, y, z;
			for each(Point3d p in mVerts)
			{
				x = p.x();
				y = p.y();
				z = p.z();

				mLeftX = std::min(x, mLeftX);
				mRightX = std::max(x, mRightX);

				mBotY = std::min(y, mBotY);
				mTopY = std::max(y, mTopY);

				mBackZ = std::min(z, mBackZ);
				mFrontZ = std::max(z, mFrontZ);

			}

			boundsFound = true;
		}
	}

	void Skeleton::GetBoundingSphere(float &cx, float &cy, float &cz, float &r)
	{
		FindBounds();
		std::vector<Point3d> corners = {};
		corners.push_back(Point3d(mLeftX, mBotY, mBackZ));
		corners.push_back(Point3d(mLeftX, mBotY, mFrontZ));
		corners.push_back(Point3d(mLeftX, mTopY, mBackZ));
		corners.push_back(Point3d(mLeftX, mTopY, mFrontZ));
		corners.push_back(Point3d(mRightX, mBotY, mBackZ));
		corners.push_back(Point3d(mRightX, mBotY, mFrontZ));
		corners.push_back(Point3d(mRightX, mTopY, mBackZ));
		corners.push_back(Point3d(mRightX, mTopY, mFrontZ));

		float maxDist = 0.0;
		float maxCX=0, maxCY=0, maxCZ=0, maxR=0;
		for (int i = 0; i < corners.size(); ++i)
		{
			for (int k = i + 1; k < corners.size(); ++k)
			{
				Point3d dif = corners[k] - corners[i];
				float dist = dif.mag();
				if (dist > maxDist)
				{
					maxR = dist / 2;
					Point3d sum = corners[k] + corners[i];
					Point3d centerpoint = sum / 2;
					//std::cout << "Point1 " << corners[k] << "Point2 " << corners[i] << "Center " << centerpoint;
					maxCX = centerpoint.x();
					maxCY = centerpoint.y();
					maxCZ = centerpoint.z();
				}
			}
		}
		cx = maxCX;
		cy = maxCY;
		cz = maxCZ;
		//pad the sphere a little bit for giggles
		r = maxR * 1.1;
		//std::cout << "Skeleton bounding sphere " << std::endl;
		//std::cout << "cx: " << cx << " cy: " << cy << " cz: " << cz << " r: " << r << std::endl;
		if (!originalCenterSet)
		{
			mOriginalCenter = Point3d(cx, cy, cz);
			mCurrentCenter = Point3d(cx, cy, cz);
			originalCenterSet = true;
		}
	}

	void Skeleton::MoveCenterTo(Point3d targetCenter)
	{
		Point3d offset = targetCenter - mCurrentCenter;
		for (int i = 0; i < mVerts.size(); ++i)
		{
			mVerts[i] = mVerts[i] + offset;
		}
		mCurrentCenter = targetCenter;
		boundsFound = false;
		FindBounds();
	}

	void Skeleton::ResetCenter()
	{
		MoveCenterTo(mOriginalCenter);
	}

	//Json::Value Skeleton::ToJson()
	//{
	//	Json::Value skelJson;
	//	skelJson["Vertices"] = Json::Value(Json::ValueType::arrayValue);
	//	skelJson["Vertices"].resize(mVerts.size());

	//	for (int i = 0; i < mVerts.size(); ++i)
	//	{
	//		skelJson["Vertices"][i] = mVerts[i].ToJson();
	//	}

	//	skelJson["Edges"] = Json::Value(Json::ValueType::arrayValue);

	//	for (int i = 0; i < mEdges.size(); ++i)
	//	{
	//		for (int j = 0; j < mEdges[i].size(); ++i)
	//		{
	//			SkeletonEdge *skelEdge = mEdges[i][j];
	//			skelJson["Edges"].append(skelEdge->ToJson());
	//		}
	//	}

	//	return skelJson;
	//}





	edgePtrList tidyUpRoot(edgePtrList toTidy, bool &valid, int &end1, int &end2, SkeletonEdge *beginning)
	{
		std::sort(toTidy.begin(), toTidy.end(), edgePtrLessThan);

		std::map<int, std::vector<SkeletonEdge*>> vertEdges;
		for (int i = 0; i < toTidy.size(); ++i)
		{
			SkeletonEdge *edge = toTidy[0];
			if (!edge)
			{
				valid = false;
				return edgePtrList();
			}

			if (vertEdges.count(edge->v0) == 0)
			{
				vertEdges[edge->v0] = {edge};
			}
			else
			{
				vertEdges[edge->v0].push_back(edge);
			}
			if (vertEdges.count(edge->v1) == 0)
			{
				vertEdges[edge->v1] = { edge };
			}
			else
			{
				vertEdges[edge->v1].push_back(edge);
			}
		}

		std::vector<int> endpoints;
		for each(auto keyValPair in vertEdges)
		{
			if (keyValPair.second.size() == 1)
			{
				endpoints.push_back(keyValPair.first);
			}
			if (keyValPair.second.size() > 2)
			{
				valid = false;
				return edgePtrList();
			}
		}

		//if the number of endpoints (vertices degree 1 or less) is not identically 2, the root cannot be valid
		if (endpoints.size() != 2)
		{
			valid = false;
			return edgePtrList();
		}

		
		int nextVert = 0;

		if (beginning != nullptr)
		{
			if (endpoints[0] == beginning->v0 || endpoints[0] == beginning->v1)
			{
				nextVert = endpoints[0];
			}
			else if (endpoints[1] == beginning->v1 || endpoints[1] == beginning->v1)
			{
				nextVert = endpoints[1];
			}
			else
			{
				nextVert = std::min(endpoints[0], endpoints[1]);
			}
		}
		else
		{
			nextVert = std::min(endpoints[0], endpoints[1]);
		}

		end1 = nextVert;

		SkeletonEdge *priorEdge = nullptr;
		bool hasNext = true;
		
		edgePtrList result = edgePtrList();
		while (hasNext)
		{
			if (vertEdges[nextVert].size() == 1)
			{
				if (vertEdges[nextVert][0] != priorEdge)
				{
					priorEdge = vertEdges[nextVert][0];
					result.push_back(priorEdge);
				}
				else
				{
					hasNext = false;
				}
			}
			else
			{
				if (vertEdges[nextVert][0] != priorEdge)
				{
					priorEdge = vertEdges[nextVert][0];
					result.push_back(priorEdge);
				}
				else
				{
					priorEdge = vertEdges[nextVert][1];
					result.push_back(priorEdge);
				}
			}
			
			nextVert = priorEdge->v0 != nextVert ? priorEdge->v0 : priorEdge->v1;
		}

		end2 = nextVert;
		valid = true;
		return result;
	}

	std::ostream& operator<<(std::ostream& out, const Skeleton& skel)
	{
		out << Skeleton::beginSkeletonString << std::endl;
		out << Skeleton::vertexString << " " << skel.mVerts.size() << std::endl;
		int numEdges = 0;
		for (int i = 0; i < skel.mEdges.size(); ++i)
		{
			for (int j = 0; j < skel.mEdges[i].size(); ++j)
			{
				++numEdges;
			}
		}

		out << Skeleton::edgeString << " " << numEdges << std::endl;
		out << Skeleton::endHeaderString << std::endl;

		for (int i = 0; i < skel.mVerts.size(); ++i)
		{
			out << skel.mVerts[i];
		}
		for (int i = 0; i < skel.mEdges.size(); ++i)
		{
			for (int j = 0; j < skel.mEdges[i].size(); ++j)
			{
				out << (skel.mEdges[i][j]);
			}
		}
		out << Skeleton::endSkeletonString << std::endl;
		return out;
	}

	std::istream& operator>>(std::istream& in, Skeleton &skel)
	{
		skel.boundsFound = false;
		bool endHeaderReached = false;
		std::string line;
		int numEdges = 0, numVerts = 0;
		Log::out << "Istream operation for skeleton" << std::endl;
		std::getline(in, line);
		if (!boost::iequals(line, Skeleton::beginSkeletonString))
		{
			Log::out << "The skeleton file lacks the " << Skeleton::beginSkeletonString << "header. Should attempt to load wenzhen style file " << std::endl;
			std::vector<std::string> words = {};
			boost::split(words, line, boost::is_any_of(" "));
			int numVerts, numEdges, numFaces;
			if (words.size() < 3)
			{
				if (words.size() == 2)
				{
					numVerts = boost::lexical_cast<int>(words[0]);
					numEdges = boost::lexical_cast<int>(words[1]);
				}
				else
				{
					Log::out << "The skeleton header lacks sufficient info to parse.  Ending read attempt" << std::endl;
					return in;
				}
			}
			else
			{
				numVerts = boost::lexical_cast<int>(words[0]);
				numEdges = boost::lexical_cast<int>(words[1]);
				numFaces = boost::lexical_cast<int>(words[2]);
			}
			

			skel.LoadWenzhenStyleFile(in, numVerts, numEdges, numFaces);

			return in;
		}
		Log::out << "This is not being recognized as a wenzhen style file" << std::endl;
		Log::out << "The first input line is " << line << std::endl;
		std::vector<std::vector<std::string>> headerInfo;
		while (!endHeaderReached && !in.eof())
		{
			std::getline(in, line);
			std::vector<std::string> words = {};
			boost::split(words, line, boost::is_any_of(" "));
			if (boost::iequals(words[0], Skeleton::endHeaderString))
			{
				endHeaderReached = true;
				break;
			}
			headerInfo.push_back(words);
		}
		if (in.eof())
		{
			Log::out << "End of file reached before end of header info, improper file" << std::endl;
			return in;
		}
		for (int i = 0; i < headerInfo.size(); ++i)
		{
			std::vector<std::string> headerLine = headerInfo[i];
			if (boost::iequals(headerLine[0], vertexString))
			{
				//handle vertices
				numVerts = boost::lexical_cast<int>(headerLine[1]);
				skel.mVerts = vertList(numVerts);
				for (int i = 0; i < numVerts; ++i)
				{
					in >> skel.mVerts[i];
				}
			}
			else if (boost::iequals(headerLine[0], edgeString))
			{
				//handle edges
				numEdges = boost::lexical_cast<int>(headerLine[1]);
				edgeList edges = edgeList(numEdges);
				for (int i = 0; i < numEdges; ++i)
				{
					in >> edges[i];
				}
				for each(SkeletonEdge edge in edges)
				{
					skel.AddEdge(edge);
				}
			}
			else
			{
				//idk mang
				//do nothing
			}
		}
		float cx, cy, cz, r;
		skel.GetBoundingSphere(cx, cy, cz, r);
	}
	
}