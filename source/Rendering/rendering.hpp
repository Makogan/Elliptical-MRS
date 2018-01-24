//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
/*
*	Author:	Camilo Talero
*
*
*	Version: 0.0.1
*
*	Header to expose rendering functions to other compilation units.
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Includes and macros
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#pragma once

#include "wavefront-loader.hpp"
#include "rendering-commons.hpp"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

struct Graph
{
    vector<vec3> nodes;
    vector<vector<uint>> graph;
    vector<double> lengths;

    Graph(vector<vec3> *vertices, vector<uint> *indices);
    void djikstra(uint start);

    void toString();
};

//========================================================================================
/*
*	List of function headers:
*/
//========================================================================================

//Load shading information 
int loadViewProjMatrix(Camera &c, GLuint &program);
void loadColor(vec4 color, GLuint program);
int loadCamera(vec3 cameraPos, GLuint program);

//Major wrap/control functions
void render_loop(GLFWwindow*);
void end_rendering(GLFWwindow* window);

//########################################################################################