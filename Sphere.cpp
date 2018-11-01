#include "Sphere.h"
#include <iostream>
namespace drawing
{


	void normalize(GLfloat *a) 
	{
		GLfloat d = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
		a[0] /= d; a[1] /= d; a[2] /= d;
	}

	
	int GetVertexBotRowCount(int subdivisions)
	{
		return pow(2, subdivisions)+1;
	}

	int GetVertexRowCount(int row, int subdivisions)
	{
		int botRowCount = GetVertexBotRowCount(subdivisions);
		return botRowCount - row;
	}

	int GetVertexRowStartI(int row, int subdivisions)
	{
		int curI = 0;
		for (int i = 0; i < row; ++i)
		{
			curI += GetVertexRowCount(i, subdivisions);
		}
		return curI;
	}
}

void drawing::VBOSphere::init(double radius, int numSubdivisions)
{
	rad = radius;
	subdivisions = numSubdivisions;
	
	int numVerticesPerFace = GetNumVerticesForFace(subdivisions);
	int numIndicesPerFace = GetNumIndicesForFace(subdivisions);
	numVertices = 20 * numVerticesPerFace;
	numIndices = 20 * numIndicesPerFace;


	vertices.resize(numVertices * 3);
	normals.resize(numVertices * 3);
	indices.resize(numIndices);

	for (int side = 0; side < 20; ++side)
	{
		int vertexStart = numVerticesPerFace * side * 3;
		int indexStart = numIndicesPerFace * side;
		CreateSubdividedIcoFace(side, subdivisions, radius, &vertices[vertexStart], &normals[vertexStart], &indices[0], indexStart, indices);
	}
	
}

drawing::VBOSphere::VBOSphere(double radius, int numSubdivisions)
{
	init(radius, numSubdivisions);
}

drawing::VBOSphere::VBOSphere()
{
	init(1.0, 1);
}


void drawing::VBOSphere::draw()
{
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);
	glNormalPointer(GL_FLOAT, 0, &normals[0]);

	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void drawing::VBOSphere::fancyDraw(float red, float green, float blue, float x, float y, float z, float scale)
{
	GLfloat color[4] = { red, green, blue, 1.0 };
	fancierDraw(color, x, y, z, scale);
}

void drawing::VBOSphere::fancierDraw(GLfloat* color, float x, float y, float z, float scale)
{
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
	glTranslatef(x, y, z);
	glScalef(scale, scale, scale);
	draw();
	glScalef(1.0 / scale, 1.0 / scale, 1.0 / scale);
	glTranslatef(-x, -y, -z);
}

void drawing::VBOSphere::pickDraw(GLubyte *color, float x, float y, float z, float scale)
{
	glColor3ub(color[0], color[1], color[2]);
	glTranslatef(x, y, z);
	glScalef(scale, scale, scale);
	
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);

	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, &indices[0]);

	glDisableClientState(GL_VERTEX_ARRAY);

	glScalef(1.0 / scale, 1.0 / scale, 1.0 / scale);
	glTranslatef(-x, -y, -z);
	
}

drawing::VBOCube::VBOCube(double radius)
{
	init(radius);
}

drawing::VBOCube::VBOCube()
{
	init(1.0);
}

void drawing::VBOCube::init(float radius)
{
	vertices = 
	{
		-radius, -radius, -radius,
		-radius, -radius, radius,
		 -radius, radius, radius ,
		-radius, radius, -radius,
		 radius, -radius, -radius ,
		 radius, -radius, radius ,
		 radius, radius, radius ,
		 radius, radius, -radius 
	};

	indices =
	{
		0, 1, 2, 3,
		4, 5, 6, 7,
		0, 1, 5, 4,
		1, 2, 6, 5,
		2, 3, 7, 6,
		3, 0, 4, 7
	};
}

void drawing::VBOCube::fancierDraw(GLfloat* color, float x, float y, float z, float scale)
{
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
	glTranslatef(x, y, z);
	glScalef(scale, scale, scale);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);
	glNormalPointer(GL_FLOAT, 0, &vertices[0]);

	glDrawElements(GL_QUADS, indices.size(), GL_UNSIGNED_INT, &indices[0]);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glScalef(1.0 / scale, 1.0 / scale, 1.0 / scale);
	glTranslatef(-x, -y, -z);
}

void drawing::VBOCube::pickDraw(GLubyte* color, float x, float y, float z, float scale)
{
	glColor3ub(color[0], color[1], color[2]);
	glTranslatef(x, y, z);
	glScalef(scale, scale, scale);

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);

	glDrawElements(GL_QUADS, indices.size(), GL_UNSIGNED_INT, &indices[0]);

	glDisableClientState(GL_VERTEX_ARRAY);

	glScalef(1.0 / scale, 1.0 / scale, 1.0 / scale);
	glTranslatef(-x, -y, -z);
}

int drawing::GetNumVerticesForFace(int subdivisions)
{
	return pow(2, 2 * subdivisions - 1) + pow(2, subdivisions) + pow(2, subdivisions - 1) + 1;
}

int drawing::GetNumIndicesForFace(int subdivisions)
{
	return pow(4, subdivisions) * 3;
}

void drawing::CreateSubdividedIcoFace(int side, int subdivisions, float radius, GLfloat *vertexStart, GLfloat *normalStart, GLuint *indexStart, int indexOffset, std::vector<GLuint> &indices)
{
	GLfloat *a, *b, *c;
	a = vdata[tindices[side][0]];
	b = vdata[tindices[side][1]];
	c = vdata[tindices[side][2]];

	GLfloat acTop[3];
	GLfloat bcTop[3];
	GLfloat acBot[3];
	GLfloat bcBot[3];

	GLuint triangleRowCount = (GLuint)pow(2, subdivisions);
	float fTriangleRowCount = triangleRowCount;
	GLuint localIndex = 0;
	GLuint globalVertexIndex = GetNumVerticesForFace(subdivisions) * side;
	indexStart = indexStart + indexOffset;
	for (int triangleRow = 0; triangleRow < triangleRowCount; ++triangleRow)
	{
		float fTriangleRow = triangleRow;
		
		float abFracBot = (fTriangleRowCount - fTriangleRow) / fTriangleRowCount;
		float cFracBot = 1.0 - abFracBot;

		float abFracTop = (fTriangleRowCount - (fTriangleRow + 1)) / fTriangleRowCount;
		float cFracTop = 1.0 - abFracTop;

		for (int i = 0; i < 3; ++i)
		{
			acBot[i] = a[i] * abFracBot + c[i] * cFracBot;
			bcBot[i] = b[i] * abFracBot + c[i] * cFracBot;

			acTop[i] = a[i] * abFracTop + c[i] * cFracTop;
			bcTop[i] = b[i] * abFracTop + c[i] * cFracTop;
		}

		//this triangle rows bot vertices
		GLuint botRowVertexCount = GetVertexRowCount(triangleRow, subdivisions);
		GLuint botRowVertexStartI = GetVertexRowStartI(triangleRow, subdivisions);
		GLuint botRowVertexEndI = botRowVertexStartI + botRowVertexCount;

		//this triangle rows top vertices
		GLuint topRowVertexCount = GetVertexRowCount(triangleRow + 1, subdivisions);
		GLuint topRowVertexStartI = GetVertexRowStartI(triangleRow + 1, subdivisions);
		GLuint topRowVertexEndI = topRowVertexStartI + topRowVertexCount;
		if (topRowVertexEndI == 45)
		{
			int p = 3;
		}

		GLfloat *toNormalize;

		for (int botRowVertexI = botRowVertexStartI; botRowVertexI < botRowVertexEndI; ++botRowVertexI)
		{
			float aFrac = ((float)(botRowVertexEndI - botRowVertexI - 1)) / (botRowVertexCount - 1);
			float bFrac = 1.0 - aFrac;
			for (int dim = 0; dim < 3; ++dim)
			{
				vertexStart[botRowVertexI * 3 + dim] = acBot[dim] * aFrac + bcBot[dim] * bFrac;
			}
			toNormalize = vertexStart + botRowVertexI * 3;
			normalize(toNormalize);
			for (int dim = 0; dim < 3; ++dim)
			{
				normalStart[botRowVertexI * 3 + dim] = vertexStart[botRowVertexI * 3 + dim];
				vertexStart[botRowVertexI * 3 + dim] *= radius;
			}
		}
		for (int topRowVertexI = topRowVertexStartI; topRowVertexI < topRowVertexEndI; ++topRowVertexI)
		{
			float aFrac, bFrac;
			if (topRowVertexCount == 1)
			{
				aFrac = 1.0;
				bFrac = 0.0;
			}
			else
			{
				aFrac = ((float)(topRowVertexEndI - topRowVertexI - 1)) / (topRowVertexCount - 1);
				bFrac = 1.0 - aFrac;
			}
			for (int dim = 0; dim < 3; ++dim)
			{
				vertexStart[topRowVertexI * 3 + dim] = acTop[dim] * aFrac + bcTop[dim] * bFrac;
			}
			toNormalize = vertexStart + topRowVertexI * 3;
			normalize(toNormalize);
			for (int dim = 0; dim < 3; ++dim)
			{
				normalStart[topRowVertexI * 3 + dim] = vertexStart[topRowVertexI * 3 + dim];
				vertexStart[topRowVertexI* 3 + dim] *= radius;
			}
		}
		
		GLuint botRowLeft = botRowVertexStartI;
		GLuint topRowLeft = topRowVertexStartI;

		bool isBotTriangle = true;
		
		while (botRowLeft < botRowVertexEndI - 1)
		{
			if (isBotTriangle)
			{
				indexStart[localIndex] = botRowLeft + globalVertexIndex;
				++localIndex;
				++botRowLeft;
				indexStart[localIndex] = botRowLeft + globalVertexIndex;
				++localIndex;
				indexStart[localIndex] = topRowLeft + globalVertexIndex;
				++localIndex;
			}
			else
			{
				indexStart[localIndex] = topRowLeft + globalVertexIndex;
				++localIndex;
				++topRowLeft;
				indexStart[localIndex] = botRowLeft + globalVertexIndex;
				++localIndex;
				indexStart[localIndex] = topRowLeft + globalVertexIndex;
				++localIndex;
			}
			isBotTriangle = !isBotTriangle;
		}

		
	}
}