#include "BoostSkeleton.h"

namespace Roots
{
	std::string BSkeleton::beginSkeletonString = "begin_skeleton";
	std::string BSkeleton::endSkeletonString = "end_skeleton";
	std::string BSkeleton::vertexString = "vertex";
	std::string BSkeleton::edgeString = "edge";
	std::string BSkeleton::endHeaderString = "end_header";
	std::string BSkeleton::beginPlyString = "ply";


	BSkeleton::BSkeleton()
	{
		Initialize();
	}

	void BSkeleton::Initialize()
	{
		mMinX = 0, mMaxX = 0, mMinY = 0, mMaxY = 0, mMinZ = 0, mMaxZ = 0;
		mBoundsFound = false;
		mCenter = Point3d();
		originalCenter = Point3d();
		mRadius = 0;
		glVertices = {};
		faces = {};
		clear();
	}

	int BSkeleton::loadFromLines(std::vector<std::string> &lines, int startingLine)
	{
		Initialize();
		mBoundsFound = false;
		bool endHeaderReached = false;
		std::string line;
		int numEdges = 0, numVerts = 0;
		line = lines[startingLine];
		boost::trim(line);

		int result = startingLine;
		std::cout << "Skeleton file type: " << line << std::endl;
		if (boost::iequals(line, beginSkeletonString))
		{
			std::cout << "The Skeleton file lacks the " << beginSkeletonString << "header.  Should attempt to load wenzhen style file" << std::endl;
			result = loadWenzhenLines(lines, startingLine);
		}
		else if (boost::iequals(line, beginPlyString))
		{
			result = loadDanPly(lines, startingLine);
		}
		else
		{
			std::cout << "This is recognized as a PLY style file.  Attempting to read..." << std::endl;
			result = loadPlyStyleLines(lines, startingLine);
		}
		
		findBounds();
		findBoundingSphere();
		return result;
	}

	void BSkeleton::writeToStream(std::ostream & out)
	{

		out << beginSkeletonString << std::endl;
		out << vertexString << " " << m_vertices.size() << std::endl;
		out << edgeString << " " << m_edges.size() << std::endl;
		out << endHeaderString << std::endl;
		skelVertIter vi = boost::vertices(*this);
		while (vi.first != vi.second)
		{
			out << getVertData(*vi.first) << std::endl;
			++vi;
		}

		skelEdgeIter ei = boost::edges(*this);
		while (ei.first != ei.second)
		{
			out << getEdgeData(*ei.first) << std::endl;
			++ei;
		}
	}

	int BSkeleton::loadWenzhenLines(std::vector<std::string> &lines, int startingLine)
	{
		std::cout << "Load Wenzhen lines" << std::endl;
		int lineOn = startingLine;
		
		int numEdges=0, numVertices=0, numFaces=0;
		std::vector<std::string> words = {};

		boost::split(words, lines[lineOn], boost::is_any_of(" "));
		++lineOn;
		if (words.size() < 2)
		{

		}
		if (words.size() >= 2)
		{
			numVertices = boost::lexical_cast<int>(words[0]);
			numEdges = boost::lexical_cast<int>(words[1]);
		}
		if (words.size() == 3)
		{
			numFaces = boost::lexical_cast<int>(words[2]);
		}
		
		loadVertices(lines, lineOn, numVertices);
		loadEdges(lines, lineOn, numEdges);
		if (lineOn >= lines.size())
		{
			std::cout << "Line on greater than lines size" << std::endl;
		}
		return lineOn;
	}

	int BSkeleton::loadDanPly(std::vector<std::string> &lines, int startingLine)
	{
		int lineOn = startingLine;

		int numEdges = 0, numVertices = 0, numFaces = 0;
		std::vector<std::string> words = {};

		bool headerEnded = false;

		while (!headerEnded)
		{
			boost::split(words, lines[lineOn], boost::is_any_of(" "));
			//vertices
			if (boost::iequals(words[0], "element") && boost::iequals(words[1], "vertex"))
			{
				numVertices = boost::lexical_cast<int>(words[2]);
				mVertexParseOrder = std::map<ParsingOrder, int>();
				mVertexWriteOrder = std::map<int, ParsingOrder>();
				++lineOn;
				boost::split(words, lines[lineOn], boost::is_any_of(" "));
				int propertyPos = 0;
				while (boost::iequals(words[0], "property"))
				{
					boost::split(words, lines[lineOn], boost::is_any_of(" "));
					
					if (!boost::iequals(words[0], "property"))
					{
						break;
					}
					if (boost::iequals(words[2], "bt2"))
					{
						mVertexParseOrder[ParsingOrder::Thickness] = propertyPos;
						mVertexWriteOrder[propertyPos] = ParsingOrder::Thickness;
					}
					else if (boost::iequals(words[2], "radius"))
					{
						mVertexParseOrder[ParsingOrder::Width] = propertyPos;
						mVertexWriteOrder[propertyPos] = ParsingOrder::Width;
					}
					else if (boost::iequals(words[2], "x"))
					{
						mVertexParseOrder[ParsingOrder::X] = propertyPos;
						mVertexWriteOrder[propertyPos] = ParsingOrder::X;
					}
					else if (boost::iequals(words[2], "y"))
					{
						mVertexParseOrder[ParsingOrder::Y] = propertyPos;
						mVertexWriteOrder[propertyPos] = ParsingOrder::Y;
					}
					else if (boost::iequals(words[2], "z"))
					{
						mVertexParseOrder[ParsingOrder::Z] = propertyPos;
						mVertexWriteOrder[propertyPos] = ParsingOrder::Z;
					}

					++lineOn;
					++propertyPos;
				}
			}
			if (boost::iequals(words[0], "element") && boost::iequals(words[1], "edge"))
			{
				numEdges = boost::lexical_cast<int>(words[2]);
			}
			if (boost::iequals(words[0], "element") && boost::iequals(words[1], "face"))
			{
				numFaces = boost::lexical_cast<int>(words[2]);
			}
			if (boost::iequals(words[0], "end_header"))
			{
				headerEnded = true;
			}
			++lineOn;
		}
		std::cout << "loading vertices " << std::endl;
		loadVertices(lines, lineOn, numVertices);
		std::cout << "loading edges " << std::endl;
		loadEdges(lines, lineOn, numEdges);
		std::cout << "loading faces " << std::endl;
		loadFaces(lines, lineOn, numFaces);
		std::cout << "Finished loading skeleton " << std::endl;
		

		return lineOn;
	}

	void BSkeleton::writeDanPly(std::ostream &os)
	{

		skelVertIter svi = boost::vertices(*this);

		for (; svi.first != svi.second; ++svi)
		{
			Point3d *v = &operator[](*svi.first);
			for (int propertyPos = 0; propertyPos < ParsingOrder::ParsingCount; ++propertyPos)
			{
				os << v->operator[](mVertexWriteOrder[propertyPos]) << " ";
			}
			os << std::endl;
		}

		skelEdgeIter sei = boost::edges(*this);

		for (; sei.first != sei.second; ++sei)
		{
			os << sei.first->m_source << " " << sei.first->m_target << std::endl;
		}

		for (int i = 0; i < faces.size(); ++i)
		{
			os << faces[i] << std::endl;
		}
		return;
	}

	int BSkeleton::loadPlyStyleLines(std::vector<std::string> &lines, int startingLine)
	{
		std::cout << "Load ply style lines" << std::endl;
		std::vector<std::vector<std::string>> headerInfo;
		bool endHeaderReached = false;
		int lineOn = startingLine;
		while (!endHeaderReached && lineOn < lines.size())
		{
			std::vector<std::string> words = {};
			boost::split(words, lines[lineOn], boost::is_any_of(" "));
			++lineOn;
			if (boost::iequals(words[0], endHeaderString))
			{
				endHeaderReached = true;
			}
			headerInfo.push_back(words);

		}
		if (lineOn >= lines.size())
		{
			return lineOn;
		}

		for (int headerLine = 0; headerLine < headerInfo.size(); ++headerLine)
		{
			std::vector<std::string> headerWords = headerInfo[headerLine];

			if (boost::iequals(headerWords[0], vertexString))
			{
				//handle vertices
				int numVerts = boost::lexical_cast<int>(headerWords[1]);
				loadVertices(lines, lineOn, numVerts);
				
			}
			else if (boost::iequals(headerWords[0], edgeString))
			{
				//handle edges
				int numEdges = boost::lexical_cast<int>(headerWords[1]);
				loadEdges(lines, lineOn, numEdges);
			}
			else
			{
			}
		}
		return lineOn;
	}

	void BSkeleton::loadVertices(std::vector<std::string> &lines, int &lineOn, int numVertices)
	{
		std::vector<std::string> words;
		float vertData[5];
		for (int i = 0; i < numVertices; ++i, ++lineOn)
		{
			words = {};
			boost::split(words, lines[lineOn], boost::is_any_of(" "));
			for (int parsingOrder = 0; parsingOrder < ParsingOrder::ParsingCount; ++parsingOrder)
			{
				vertData[parsingOrder] = boost::lexical_cast<float>(words[mVertexParseOrder[(ParsingOrder)parsingOrder]]);
			}
			addVertex(Point3d(vertData[0], vertData[1], vertData[2], vertData[3], vertData[4]));
		}
		updateGLVertices();
	}

	void BSkeleton::loadEdges(std::vector<std::string> &lines, int &lineOn, int numEdges)
	{
		std::cout << "Beginning to load edges " << std::endl;
		std::cout << "first line is : " << lines[lineOn] << std::endl;
		int v0, v1;
		std::vector<std::string> words;
		for (int i = 0; i < numEdges && lineOn < lines.size(); ++i, ++lineOn)
		{
			words = {};
			boost::split(words, lines[lineOn], boost::is_any_of(" "));
			v0 = boost::lexical_cast<int>(words[0]);
			v1 = boost::lexical_cast<int>(words[1]);

			if (m_vertices.size() > v1 && m_vertices.size() > v0)
			{
				addEdge(v0, v1);
			}
		}
	}

	void BSkeleton::loadFaces(std::vector<std::string> &lines, int &lineOn, int numFaces)
	{
		faces.resize(numFaces);

		std::vector<GLuint> faceVerts = {0, 0, 0};
		std::vector<std::string> words;
		int i = 0;
		for (;i < numFaces && lineOn < lines.size(); ++i, ++lineOn)
		{
			words = {};
			boost::split(words, lines[lineOn], boost::is_any_of(" "));
			for (int v = 0; v < 3; ++v)
			{
				faceVerts[v] = (GLuint)boost::lexical_cast<int>(words[v + 1]);
			}
			faces.push_back(Face(faceVerts, &glVertices[0]));
		}
		if (i < numFaces)
		{
			std::cout << "Not all faces loaded?  " << std::endl << "First face " << faces[0] << std::endl;
			std::cout << "Last face " << faces[i] << std::endl;
		}
	}

	SkelEdge BSkeleton::addEdge(int v0, int v1)
	{
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
	SkelVert BSkeleton::addVertex(Point3d pointLocation)
	{
		SkelVert v;
		v = boost::add_vertex(*this);
		operator[](v) = pointLocation;
		operator[](v).id = v;
		mBoundsFound = false;
		return v;
	}

	void BSkeleton::updateGLVertices()
	{
		glVertices.resize(m_vertices.size() * 3);

		for (SkelVert i = 0; i < m_vertices.size(); ++i)
		{
			for (int dim = 0; dim < 3; ++dim)
			{
				glVertices[i * 3 + dim] = operator[](i)[dim];
			}
		}
	}

	void BSkeleton::findBounds()
	{
		if (mBoundsFound)
		{
			return;
		}
		else
		{
			float big = 99999999999;
			mMinX = big;
			mMaxX = -big;
			mMinY = big;
			mMaxY = -big;
			mMinZ = big;
			mMaxZ = -big;

			Point3d p;
			skelVertIter it = boost::vertices(*this);
			while (it.first != it.second)
			{
				p = getVertData(*(it.first));
				mMinX = std::min(p[0], mMinX);
				mMinY = std::min(p[1], mMinY);
				mMinZ = std::min(p[2], mMinZ);
				mMaxX = std::max(p[0], mMaxX);
				mMaxY = std::max(p[1], mMaxY);
				mMaxZ = std::max(p[2], mMaxZ);
				++it;
			}
			mBoundsFound = true;
		}

	}
	void BSkeleton::findBoundingSphere()
	{
		findBounds();
		Point3d minPoint = Point3d(mMinX, mMinY, mMinZ, 0, 0);
		Point3d maxPoint = Point3d(mMaxX, mMaxY, mMaxZ, 0, 0);

		Point3d sum = Point3d();
		skelVertIter vIter = boost::vertices(*this);
		for(; vIter.first != vIter.second; ++vIter)
		{
			sum = sum + operator[](*vIter.first);
		}

		sum = sum / m_vertices.size();
		mCenter = sum;
		originalCenter = mCenter;
		mRadius = (maxPoint - mCenter).mag();
	}

	void BSkeleton::recenterSkeleton(Point3d newCenter)
	{
		findBoundingSphere();
		Point3d offset = newCenter - mCenter;
		mCenter = newCenter;

		skelVertIter svi = boost::vertices(*this);
		while (svi.first != svi.second)
		{
			operator[](*svi.first) += offset;
			++svi;
		}
		updateGLVertices();
		return;
	}

	SkelEdge BSkeleton::getEdge(int v0, int v1, bool &exists)
	{
		SkelEdge result;
		boost::tie(result, exists) = boost::edge(v0, v1, *this);
		return result;
	}

	RootAttributes& BSkeleton::getEdgeData(int v0, int v1, bool &exists)
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

	RootAttributes& BSkeleton::getEdgeData(SkelEdge e)
	{
		return operator[](e);
	}
	Point3d& BSkeleton::getVertData(SkelVert v)
	{
		return operator[](v);
	}

	std::vector<SkelVert> BSkeleton::GetNotableVertices()
	{
		std::vector<SkelVert> result = {};
		skelVertIter vi = boost::vertices(*this);
		while (vi.first != vi.second)
		{
			int deg = boost::degree(*vi.first, *this);
			if (deg != 2 && deg != 0)
			{
				result.push_back(*vi.first);
			}
			++vi;
		}
		return result;
	}
}

PySkeleton::PySkeleton()
{
	mSkeleton = nullptr;
	minThickness = 0;
	maxThickness = 0;
	minWidth = 0;
	maxWidth = 0;
	minRatio = 0;
	maxRatio = 0;
	thicknessPercentiles = boost::python::list();
	widthPercentiles = boost::python::list();
	ratioPercentiles = boost::python::list();

}

PySkeleton::PySkeleton(Roots::BSkeleton *srcSkel)
{
	mSkeleton = srcSkel;
	reload();
}

void PySkeleton::findBoundingSphere()
{
	mSkeleton->findBoundingSphere();
}

void PySkeleton::reload()
{
	radius = 0;
	center = Point3d();
	mVertexList = boost::python::list();
	mEdgeList = boost::python::list();

	thicknessPercentiles = boost::python::list();
	widthPercentiles = boost::python::list();
	ratioPercentiles = boost::python::list();

	minThickness = 100000;
	maxThickness = 0;
	minWidth = 100000;
	maxWidth = 0;
	minRatio = 100000;
	maxRatio = 0;

	std::vector<float> thicknessList = {};
	std::vector<float> widthList = {};
	std::vector<float> ratioList = {};

	if (mSkeleton == nullptr)
	{
		return;
	}
	
	for (skelVertIter svi = boost::vertices(*mSkeleton); svi.first != svi.second; ++svi)
	{
		Point3d vertexAdded = mSkeleton->operator[](*svi.first);
		mVertexList.append<Point3d>(vertexAdded);
	}
	for (skelEdgeIter sei = boost::edges(*mSkeleton); sei.first != sei.second; ++sei)
	{
	}

	std::sort(thicknessList.begin(), thicknessList.end());
	std::sort(widthList.begin(), widthList.end());
	std::sort(ratioList.begin(), ratioList.end());

	for (int i = 0; i < 11; ++i)
	{
		float fi = i;
		int idx = (fi / 10.0) * thicknessList.size();
		idx -= 1;
		idx = std::max(idx, 0);
		thicknessPercentiles.append<float>(thicknessList[idx]);
		widthPercentiles.append<float>(widthList[idx]);
		ratioPercentiles.append<float>(ratioList[idx]);
	}


	std::cout << "min thickness " << minThickness << " max thickness " << maxThickness << " min width " << minWidth << " max width " << maxWidth << " min ratio " << minRatio << " max ratio " << maxRatio << std::endl;
	radius = mSkeleton->mRadius;
	center = mSkeleton->mCenter;
}