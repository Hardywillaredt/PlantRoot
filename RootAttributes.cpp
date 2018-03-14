#include "RootAttributes.h"
#include <algorithm>


#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

namespace Roots
{
	RootAttributes::RootAttributes()
	{
		//data = boost::python::list();
		thickness = 0;
		width = 0;
		length = 0;
		/*for (int i = 0; i < NumAttributes; ++i)
		{
			data.append<float>(0);
		}*/
		euclidLength = 0.0;
		v0id = 0;
		v1id = 0;
	}

	RootAttributes::RootAttributes(float * aData, Point3d v0, Point3d v1)
	{
		//data = boost::python::list();
		thickness = aData[Thickness];
		width = aData[Width];
		length = aData[Length];
		/*for (int i = 0; i < NumAttributes; ++i)
		{
			data.append<float>(aData[i]);
		}*/
		Point3d dif = v0 - v1;
		euclidLength = dif.mag();
		v0id = v0.id;
		v1id = v1.id;
	}

	RootAttributes::RootAttributes(std::vector<float> aData, Point3d v0, Point3d v1)
	{
		//data = boost::python::list();
		/*for (int i = 0; i < NumAttributes; ++i)
		{
			data.append<float>(aData[i]);
		}*/
		thickness = aData[Thickness];
		width = aData[Width];
		length = aData[Length];
		Point3d dif = v0 - v1;
		euclidLength = dif.mag();
		v0id = v0.id;
		v1id = v1.id;
	}

	RootAttributes::RootAttributes(float aThickness, float aWidth, float aLength, Point3d v0, Point3d v1)
	{
		/*data = boost::python::list();
		for (int i = 0; i < NumAttributes; ++i)
		{
			data.append<float>(0.0);
		}*/
		thickness = aThickness;
		width = aWidth;
		length = aLength;
		/*data[Thickness] = aThickness;
		data[Width] = aWidth;
		data[Length] = aLength;*/
		Point3d dif = v0 - v1;
		euclidLength = dif.mag();
		v0id = v0.id;
		v1id = v1.id;
	}

	//RootAttributes::RootAttributes(Json::Value rootJson)
	//{
	//	data = new float[NumAttributes];
	//	for (int i = 0; i < std::min((int)NumAttributes, (int)rootJson.size()); ++i)
	//	{
	//		data[i] = rootJson[i].asfloat();
	//	}
	//}

	//Json::Value RootAttributes::ToJson()
	//{
	//	Json::Value result = Json::Value(Json::ValueType::arrayValue);
	//	result.resize(NumAttributes);
	//	for (int i = 0; i < NumAttributes; ++i)
	//	{
	//		result[i] = data[i];
	//	}
	//	return result;
	//}

	bool RootAttributes::operator==(RootAttributes& second)
	{
		bool same = false;
		if (thickness == second.thickness && width == second.width && length == second.length)
		{
			same = true;
		}
		/*for (int i = 0; i < NumAttributes; ++i)
		{
			same = same && data[i] == second[i];
		}*/
		return same;
	}

	float RootAttributes::operator[](const int index)
	{
		switch (index)
		{
		case Thickness:
			return thickness;
		case Length:
			return length;
		case Width:
			return width;
		default:
			return -1;
			break;
		}
		
	}


	std::ostream& operator<<(std::ostream &out, const RootAttributes &attribs)
	{
		out << attribs.v0id << " " << attribs.v1id << " ";
		out << attribs.thickness << " " << attribs.width << " " << attribs.length;
		
		return out;
	}

	std::istream& operator>>(std::istream &in, RootAttributes &attribs)
	{
		std::string line;
		std::getline(in, line);
		std::vector<std::string> words;
		boost::split(words, line, boost::is_any_of(" "));
		attribs.thickness = boost::lexical_cast<float>(words[Thickness]);
		attribs.width = boost::lexical_cast<float>(words[Width]);
		attribs.length = boost::lexical_cast<float>(words[Length]);

		return in;
	}
}

	


