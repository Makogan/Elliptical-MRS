//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
/*
*	Author:	Camilo Talero
*
*
*	Version: 0.0.1
*
*	Header to define variables, structure definitoins, include libraries... 
*   Shared among all rendering functions.
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Includes and macros
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#define GLEW_DYNAMIC
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <time.h>
#include <queue>
#include <algorithm>    // std::sort

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtx/transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <ft2build.h>
#include FT_FREETYPE_H

#include "Camera.hpp"

using namespace std;

vector<dvec3> subdivision(vector<dvec3> verts, dvec3(*interp)(dvec3,dvec3,double));

dvec3 mlerp(dvec3 p1, dvec3 p2, double t);
dvec3 slerp(dvec3 p1, dvec3 p2, double t);
dvec3 splerp(dvec3 p1, dvec3 p2, double t);
dvec3 plerp(dvec3 p1, dvec3 p2, double u);
dvec3 ilerp(dvec3 p1, dvec3 p2, double t);
dvec3 glerp(dvec3 p1, dvec3 p2, double t);

void elliptical_D_Decomposition(vector<dvec3> fine, vector<double> w, 
	vector<dvec3> *coarse, vector<dvec4> *details, 
	dvec3(*interp)(dvec3,dvec3,double));
void elliptical_D_reconstruction(vector<dvec3> *fine, vector<double> w, 
	vector<dvec3> coarse, vector<dvec4> &details, 
    dvec3(*interp)(dvec3,dvec3,double));
void elliptical_P_Decomposition(vector<dvec3> fine, vector<double> w, 
    vector<dvec3> *coarse, vector<dvec4> *details, 
    dvec3(*interp)(dvec3,dvec3,double));
void elliptical_P_reconstruction(vector<dvec3> *fine, vector<double> w, 
    vector<dvec3> coarse, vector<dvec4> &details, 
    dvec3(*interp)(dvec3,dvec3,double));
void dtof(vector<dvec3> ds, vector<vec3> &fs);
void ftod(vector<vec3> fs, vector<dvec3> &ds);

dvec3 ellipse_project(dvec3 p, float a, float b, float c);
dvec3 sphere_project(dvec3 p, float r);

extern vector<double> weights;
extern vector<dvec3> holder;

#define a_axis 0.8
#define b_axis 1.d
#define c_axis 2.d

extern dvec3 (*INTERP)(dvec3, dvec3, double);

using namespace std;
using namespace glm;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Structure definitions
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
struct Graph
{
    vector<vec3> nodes;
    vector<vector<uint>> graph;
    vector<double> lengths;

    Graph(vector<vec3> *vertices, vector<uint> *indices);
    void djikstra(uint start);
    double node_length(uint);

    void toString();

    Graph(){}
};

//TODO: verify the following data structures, maybe program shoud be a list to allow
//  for multiple programs
struct Shader
{
	string fileName;    //name of file
	GLuint shaderID;    //generated OpenGL shader ID
	GLuint type;        //shader type
	GLuint program;     //Associated OpenGL shading program
};

struct Geometry
{
	GLuint vertexArray;     //Vertex array associated with this information

    GLuint vertexBuffer;    //Vertex data ID of the VAO
    GLuint normalsBuffer;   //Normal data (vec3)
    GLuint uvBuffer;        
    GLuint elmentBuffer;    //Element data ID of the VAO this specifies the 
                            //sequence in wich the vertices and normals will be read 

    vector<vec3> vertices;  //Vertex data
    vector<vec3> normals;   //Normal data
	vector<uint> indices;   //Element data (sequence in which data will be read)
    vector<vec2> uvs;       //Texture data for this geometry 
                            //(the associated coordinates on the mesh)
};

struct Texture
{
	GLuint textureID;   //OpenGL generated ID for the texture
    GLuint target;      //OpenGL target (Usuallly 2D texture or rectangle) check OpenGL doc

    //Dimensions of the texture
	int width;         
	int height;
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Global Values
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//TODO: delete or modify the geometry list as other representations are probably better

extern Camera cam; //Main camera for perspective projection

extern vector<GLuint> programs; //Global list of shading programs
extern vector<Shader> shaders;  //Global List of shaders
extern vector<Geometry> shapes; ///GLobal list of geometry shapes temporary
extern vector<Texture> textures; //Temporary
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//========================================================================================
/*
*	List of function headers:
*/
//========================================================================================

void loadGeometryArrays(GLuint program, Geometry &g);
void render(GLuint program, Geometry &g, GLenum drawType);

void compileShader(GLuint &shader, string &filepath, GLenum shaderType);
void createShader(Shader &s, string file, GLenum type);
void deleteShader(Shader &s);

void createGeometry(Geometry &g);
void createGeometry(Geometry &g, vector<vec3> vertices, vector<uint> indices);
void createGeometry(Geometry &g, vector<vec3> vertices,  vector<vec3> normals, 
	vector<vec2> uvs, vector<uint> indices);
void deleteGeometry(Geometry &g);

bool createTexture(Texture &texture, const char* filename, GLuint target = GL_TEXTURE_2D);
void DestroyTexture(Texture &texture);

string loadSourceFile(string &filepath);

GLuint createShadingProgram(GLuint vertexShader, GLuint fragmentShader);
double calculateAverageLength(vector<vec3> vertices, Graph g);
//########################################################################################