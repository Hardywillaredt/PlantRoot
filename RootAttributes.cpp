#include "RootAttributes.h"
#include <algorithm>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

namespace Roots
{
	RootAttributes::RootAttributes()
	{
		euclidLength = 0.0;
		v0id = 0;
		v1id = 0;
	}

	RootAttributes::RootAttributes(float * aData, Point3d v0, Point3d v1)
	{
		Point3d dif = v0 - v1;
		euclidLength = dif.mag();
		v0id = v0.id;
		v1id = v1.id;
	}

	RootAttributes::RootAttributes(std::vector<float> aData, Point3d v0, Point3d v1)
	{
		Point3d dif = v0 - v1;
		euclidLength = dif.mag();
		v0id = v0.id;
		v1id = v1.id;

	}

	RootAttributes::RootAttributes(float aThickness, float aWidth, Point3d v0, Point3d v1)
	{

		Point3d dif = v0 - v1;
		euclidLength = dif.mag();
		v0id = v0.id;
		v1id = v1.id;
	}

	RootAttributes::RootAttributes(int id0, int id1, Point3d &v0, Point3d &v1)
	{
		v0id = id0;
		v1id = id1;
		euclidLength = (v0 - v1).mag();
	}


	bool RootAttributes::operator==(RootAttributes& second)
	{
		bool same = true;
		if (v0id != second.v0id || v1id != second.v1id)
		{
			same = false;
		}
		return same;
	}



	std::ostream& operator<<(std::ostream &out, const RootAttributes &attribs)
	{
		out << attribs.v0id << " " << attribs.v1id << " ";
		
		return out;
	}

	std::istream& operator>>(std::istream &in, RootAttributes &attribs)
	{
		std::string line;
		std::getline(in, line);
		std::vector<std::string> words;
		boost::split(words, line, boost::is_any_of(" "));
		return in;
	}
}

	


