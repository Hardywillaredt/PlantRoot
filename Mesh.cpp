#include "Mesh.h"

#include <fstream>
#include <stdio.h>
#include <iostream>

Mesh::Mesh()
{
	vertices = {};
	faces = {};
	f_factor = 0.75;
	mColor[0] = 0.0;
	mColor[1] = 0.0;
	mColor[2] = 1.0;
	mColor[3] = 0.3;
}

void Mesh::render()
{
	if (vertices.size() <= 0)
	{
		return;
	}
	float alpha = mColor[3];
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, &vertices[0]);

	// render alpha = 0, to prime the depth buffer
	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LESS);
	glColor4f(mColor[0], mColor[1], mColor[2], 0);
	glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_INT, &faces[0]);

	// render with alpha = f*alpha
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glDepthFunc(GL_ALWAYS);
	glColor4f(mColor[0], mColor[1], mColor[2], f_factor*alpha);
	glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_INT, &faces[0]);


	// render with alpha = (alpha-f*alpha)/(1.0-f*alpha)
	glEnable(GL_CULL_FACE);
	glCullFace(GL_FRONT);
	glDepthFunc(GL_LEQUAL);
	glColor4f(mColor[0], mColor[1], mColor[2], (alpha - f_factor*alpha) / (1.0 - f_factor*alpha));
	glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_INT, &faces[0]);

	// render with alpha = f*alpha
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glDepthFunc(GL_ALWAYS);
	glColor4f(mColor[0], mColor[1], mColor[2], alpha*f_factor);
	glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_INT, &faces[0]);

	// render with alpha = (alpha-f*alpha)/(1.0-f*alpha)
	glDisable(GL_CULL_FACE);
	glDepthFunc(GL_LEQUAL);
	glColor4f(mColor[0], mColor[1], mColor[2], (alpha - f_factor*alpha) / (1.0 - f_factor*alpha));
	glDrawElements(GL_TRIANGLES, faces.size(), GL_UNSIGNED_INT, &faces[0]);


	glEnable(GL_CULL_FACE);
	glDisableClientState(GL_VERTEX_ARRAY);
	

	


}

void Mesh::loadFromOff(std::string filename, float *shift)
{
	if (filename.length() < 4)
	{
		std::cout << "Invalid mesh file " << filename << std::endl;
		return;
	}
	std::cout << "Shift x : " << shift[0] << " y : " << shift[1] << " z : " << shift[2] << std::endl;
	std::cout << "Loading mesh from file " << filename << " ..." << std::endl;
	std::ifstream filestream;
	filestream.open(filename);
	std::vector<std::string> lines = {};
	std::string line;
	std::vector<std::string> words = {};
	int numVertices = 0;
	int numFaces = 0;
	int loadedVerts = 0;
	int loadedFaces = 0;
	bool loadingGemoetry = false;
	while (!filestream.eof())
	{
		std::getline(filestream, line);
		
		boost::split(words, line, boost::is_any_of(" "));
		if (words.size() < 3)
		{
			continue;
		}
		if (!loadingGemoetry)
		{
			numVertices = boost::lexical_cast<int>(words[0]);
			numFaces = boost::lexical_cast<int>(words[1]);
			std::cout << "Num vertices " << numVertices << " num faces " << numFaces << std::endl;
			loadingGemoetry = true;
			vertices.resize(numVertices * 3);
			faces.resize(numFaces * 3);
		}
		else 
		{
			if (loadedVerts < numVertices)
			{
				for (int i = 0; i < 3; ++i)
				{
					vertices[loadedVerts * 3 + i] = boost::lexical_cast<float>(words[i]) + shift[i];
				}
				++loadedVerts;
			}
			else
			{
				for (int i = 0; i < 3; ++i)
				{
					//because the format of OFF faces is "NUM_VERTS IN FACE vert1 vert2 vert3 ..." we have to offset by 1
					faces[loadedFaces * 3 + i] = boost::lexical_cast<int>(words[i + 1]);
				}
				++loadedFaces;
			}
		}
	}
	filestream.close();
	std::cout << "...Finished loading mesh" << std::endl;
	return;

}

void Mesh::setAlpha(float f)
{
	mColor[3] = f;
}

void Mesh::setColor(float red, float green, float blue)
{
	mColor[0] = red;
	mColor[1] = green;
	mColor[2] = blue;
}

void Mesh::recenter(float *shift)
{
	for (int v = 0; v < vertices.size() / 3; ++v)
	{
		for (int i = 0; i < 3; ++i)
		{
			vertices[v * 3 + i] += shift[i];
		}
	}
}