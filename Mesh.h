#pragma once

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include "arcball.h"
#include <vector>
#include "boost/algorithm/string.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/config.hpp"

class Mesh
{
public:
	std::vector<GLfloat> vertices;
	std::vector<GLuint> faces;
	GLfloat mColor[4];
	GLfloat f_factor;
	Mesh();
	void render();
	void loadFromOff(std::string filename, float *shift);
	void setAlpha(float alpha);
	void setColor(float red, float green, float blue);
	void recenter(float *shift);

};