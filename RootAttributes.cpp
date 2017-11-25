#include "RootAttributes.h"
#include <algorithm>


#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"

namespace Roots
{
	RootAttributes::RootAttributes()
	{
		data = new double[NumAttributes];
		for (int i = 0; i < NumAttributes; ++i)
		{
			data[i] = 0;
		}
	}

	RootAttributes::RootAttributes(double * aData)
	{
		data = new double[NumAttributes];
		for (int i = 0; i < NumAttributes; ++i)
		{
			data[i] = aData[i];
		}
	}

	RootAttributes::RootAttributes(double aThickness, double aWidth, double aLength)
	{
		data = new double[NumAttributes];
		data[Thickness] = aThickness;
		data[Width] = aWidth;
		data[Length] = aLength;
	}

	//RootAttributes::RootAttributes(Json::Value rootJson)
	//{
	//	data = new double[NumAttributes];
	//	for (int i = 0; i < std::min((int)NumAttributes, (int)rootJson.size()); ++i)
	//	{
	//		data[i] = rootJson[i].asDouble();
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

	


	RootAttributes::~RootAttributes()
	{
		delete[] data;
	}

	double* RootAttributes::getData()
	{
		return data;
	}

	bool RootAttributes::operator==(RootAttributes& second)
	{
		bool same = true;
		for (int i = 0; i < NumAttributes; ++i)
		{
			same = same && data[i] == second[i];
		}
		return same;
	}

	double RootAttributes::operator[](const int index)
	{
		return data[index];
	}


	std::ostream& operator<<(std::ostream &out, const RootAttributes &attribs)
	{
		for (int i = 0; i < RootAttributes::NumAttributes; ++i)
		{
			out << attribs.data[i] << " ";
		}
		out << std::endl;
		return out;
	}

	std::istream& operator>>(std::istream &in, RootAttributes &attribs)
	{
		std::string line;
		std::getline(in, line);
		std::vector<std::string> words;
		boost::split(words, line, boost::is_any_of(" "));
		for (int i = 0; i < RootAttributes::NumAttributes; ++i)
		{
			attribs.data[i] = boost::lexical_cast<double>(words[i]);
		}

		return in;
	}
}

	


