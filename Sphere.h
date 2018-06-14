#pragma once

#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>

#define xfloat .525731112119133606 
#define zfloat .850650808352039932
namespace drawing
{
	static GLfloat vdata[12][3] = {
		{ -xfloat, 0.0, zfloat },{ xfloat, 0.0, zfloat },{ -xfloat, 0.0, -zfloat },{ xfloat, 0.0, -zfloat },
		{ 0.0, zfloat, xfloat },{ 0.0, zfloat, -xfloat },{ 0.0, -zfloat, xfloat },{ 0.0, -zfloat, -xfloat },
		{ zfloat, xfloat, 0.0 },{ -zfloat, xfloat, 0.0 },{ zfloat, -xfloat, 0.0 },{ -zfloat, -xfloat, 0.0 }
	};
	static GLuint tindices[20][3] = {
		{ 0,4,1 },{ 0,9,4 },{ 9,5,4 },{ 4,5,8 },{ 4,8,1 },
		{ 8,10,1 },{ 8,3,10 },{ 5,3,8 },{ 5,2,3 },{ 2,7,3 },
		{ 7,10,3 },{ 7,6,10 },{ 7,11,6 },{ 11,0,6 },{ 0,1,6 },
		{ 6,1,10 },{ 9,0,11 },{ 9,11,2 },{ 9,2,5 },{ 7,2,11 } };

	void normalize(GLfloat *a);


	struct VBOSphere
	{
		double rad;
		int numVertices;
		int subdivisions;
		int numIndices;
		std::vector<GLfloat> vertices;
		std::vector<GLfloat> normals;
		std::vector<GLuint> indices;


		VBOSphere(double radius, int numSubdivisions);
		VBOSphere();

		void draw();
		void init(double radius, int numSubdivisions);
		void fancyDraw(float red, float green, float blue, float x, float y, float z, float scale);
		void fancierDraw(GLfloat *color, float x, float y, float z, float scale);
	};

	struct VBOCube
	{
		double scale;
		std::vector<GLfloat> vertices;
		std::vector<GLuint> indices;

		VBOCube(double radius);
		VBOCube();

		void init(float radius);
		void fancierDraw(GLfloat *color, float x, float y, float z, float scale);

	};


	void CreateSubdividedIcoFace(int side, int subdivisions, float radius, GLfloat *vertexStart, GLfloat *normalStart, GLuint *indexStart, int indexOffset, std::vector<GLuint> &indices);

	int GetNumVerticesForFace(int subdivisions);
	int GetNumIndicesForFace(int subdivisions);

	int GetVertexBotRowCount(int subdivisions);
	
	int GetVertexRowCount(int row, int subdivisions);

	int GetVertexRowStartI(int row, int subdivisions);

}

