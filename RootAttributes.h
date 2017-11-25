#pragma once

#include <iostream>
#include <vector>
#include <sstream>

namespace Roots
{

	struct RootAttributes
	{

	public:
		enum Attributes
		{
			Thickness = 0,
			Width,
			Length,


			NumAttributes
		};		


		friend std::ostream& operator<<(std::ostream &out, const RootAttributes &attribs);
		friend std::istream& operator>>(std::istream &in, RootAttributes &attribs);

		RootAttributes();
		RootAttributes(double *attributeData);
		RootAttributes(double aThickness, double aWidth, double aLength);
		

		~RootAttributes();
		double *getData();

		
		

		bool operator==(RootAttributes& second);

		double operator[](const int index);

	private:
		double *data;
	
	};
}
