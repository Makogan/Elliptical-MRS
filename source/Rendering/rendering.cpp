//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
/*
*	Author:	Camilo Talero
*
*
*	Version: 0.0.1
*
*	File to define methods that abstract and simplify the use of opengl for rendering.
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Includes and macros
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "rendering.hpp"

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>
#include <glm/gtx/transform.hpp>

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Global Values
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

vector<GLuint> programs;//Global shading programs list
vector<Shader> shaders;//Global Shaders list 
//TODO: Maybe delete Geometry list and change it for a different data structure or 
//		Avoid it alltogether
vector<Geometry> shapes(10);//Global Shapes list Temporary!
vector<Texture> textures(2); //Temporary

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//========================================================================================
/*
*	Rendering Functions:
*/
//========================================================================================

//TODO: look into making other functions that also load geometry info into the shaders
//TODO: comment undocumented functions
void loadTexture(GLuint program, Texture &t)
{
	glUseProgram(program);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, t.textureID);
	GLint loc = glGetUniformLocation(program, "text");
	if(loc == GL_INVALID_VALUE || loc==GL_INVALID_OPERATION)
	{
		cerr << "Error returned when trying to find texture uniform."
			<< "\nuniform: text"
			<< "Error num: " << loc
			<< endl;
		return;
	}
	
	glUniform1i(loc,0);
}

/*
*	Function to load geometry information into an OpenGL shader program
*
*	Params:
*		program: the ID of the shader program into which we will looad the Geometry info
*		g:	the Geometry structure containing the information to load into teh program
*/
void loadGeometryArrays(GLuint program, Geometry &g)
{
	glUseProgram(program);//Make program the current shading program

	glBindVertexArray(g.vertexArray);//bind g's vertex array to be OpenGL's vertex arr

	glBindBuffer(GL_ARRAY_BUFFER, g.vertexBuffer);//set OpenGL's vertex buffer (vertex data)
	//Set the buffer data with, type, size of buffer in bytes, array pointer, drawing mode
	glBufferData(GL_ARRAY_BUFFER, g.vertices.size()*sizeof(vec3), 
		g.vertices.data(), GL_DYNAMIC_DRAW);

	//Wich layout in the shader to pass current information (layout 0 is vertex for us)
	glEnableVertexAttribArray(0);
	//Specify how to read the data: layout index, components per vertex attribute, 
	//Type of the data, normalization, byte offset between consecutive elements, 
	//offset of first component
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(vec3), (void*)0);

	//Check if there are normals to process for the geometry object
	if(g.normals.size()>0)
	{
		//Set attribute layout to 1 (layout 1 is normals for us)
		glEnableVertexAttribArray(1);
		//Set the buffer array
		glBindBuffer(GL_ARRAY_BUFFER, g.normalsBuffer);
		//As above specify how to read the information
		glVertexAttribPointer(1, 3, GL_FLOAT, GL_TRUE, sizeof(vec3), (void*)0);
		//Fill the buffer array with the data 
		glBufferData(GL_ARRAY_BUFFER, g.normals.size()*sizeof(vec3),
			g.normals.data(), GL_DYNAMIC_DRAW);
	}

	//Check if there are texture coordinates to process for the geometry object
	if(g.uvs.size()>0)
	{
		//Set attribute layout to 2 (layout 2 is uvs for us)
		glEnableVertexAttribArray(2);
		//Set the buffer array
		glBindBuffer(GL_ARRAY_BUFFER, g.uvBuffer);
		//As above specify how to read the information
		glVertexAttribPointer(2, 2, GL_FLOAT, GL_TRUE, sizeof(vec2), (void*)0);
		//Fill the buffer array with the data 
		glBufferData(GL_ARRAY_BUFFER, g.uvs.size()*sizeof(vec2),
			g.uvs.data(), GL_DYNAMIC_DRAW);
	}

	//Check if there are any indices specified in g
	if(g.indices.size()>0)
	{
		//Set the element buffer (non linear reading/loading of data information)
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g.elmentBuffer);
		//Indicate how to read the element data
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, g.indices.size()*sizeof(uint),
			g.indices.data(), GL_DYNAMIC_DRAW);
	}
}

/*
*	Fundamental rendering function, it assumes initialization of all global values 
*	used to organize the process of sending information to the GPU for rendering
*
*	Params:
*		program: OpenGL shading program used to render 
*		g: Geometry object containing mesh information
*		drawType: the drawing mode (Lines, triangles, points, strips...)
*/
void render(GLuint program, Geometry &g, GLenum drawType)
{
	//Set the OpenGL shading program
	glUseProgram(program);

	//Bind the VAO of the geometry
	glBindVertexArray(g.vertexArray);

	//Chose appropriate drawing method based on whether there are indices in the
	//geometry strucutre has explicit indices
	if(g.indices.size()>0)
		//non sequencial vertex drawing
		glDrawElements(drawType, g.indices.size(), GL_UNSIGNED_INT, (void*)0);

	else
		//sequential vertex drawing
		glDrawArrays(drawType, 0, g.vertices.size());
}

/*
*	Load the color uniform for the fragment shader (Assumes existence of a color uniform)
*
*	Params:
*		color: a representation of the RGBA values of a color
*		program: the OpenGL shading program to use
*/
void loadColor(vec4 color, GLuint program)
{
	glUseProgram(program); 	//Active program
	GLint loc = glGetUniformLocation(program, "color");
	glUniform4f(loc, color[0], color[1], color[2], color[3]);
}

/*
*	Load camera parameters into a shading program for perspective projection
*	
*	Params:
*		c: the camera struct used to create the projection matrix
*		program: program in wich to load the projection matrix
*	
*	return: 0 on error or 1 on success
*/
int loadViewProjMatrix(Camera &c, GLuint &program)
{
	glUseProgram(program);
	GLint loc = glGetUniformLocation(program, "view");
	if(loc == GL_INVALID_VALUE || loc==GL_INVALID_OPERATION)
	{
		cerr << "Error returned when trying to find view matrix."
			<< "\nuniform: view"
			<< "Error num: " << loc
			<< endl;
		return 0;
	}
	//Pass the calculated view matrix onto the shader
	glUniformMatrix4fv(loc, 1, GL_FALSE, value_ptr(c.getViewMatrix()));

	loc = glGetUniformLocation(program, "proj");
	if(loc == GL_INVALID_VALUE || loc==GL_INVALID_OPERATION)
	{

		cerr << "Error returned when trying to find projection matrix."
			<< "\nuniform: proj"
			<< "Error num: " << loc
			<< endl;
		return 0;
	}
	//Pass the calculated projection/perspective matrix onto the shader
	glUniformMatrix4fv(loc, 1, GL_FALSE, value_ptr(c.getPerspectiveMatrix()));

	return 1;
}


/*
* Load camera position into the current shader program
*
*	Params:
*		c: the camera struct used to create teh projection matrix
*		program: program in wich to load the projection matrix
*	
*	return: 0 on error or 1 on success
*/
int loadCamera(vec3 cameraPos, GLuint program)
{
	glUseProgram(program);
	GLint loc = glGetUniformLocation(program, "cameraPos");
	if (loc == -1)
	{
		cerr << "Uniform: \"cameraPos\" not found." << endl;
		return 0;
	}
	glUniform3f(loc, cameraPos[0], cameraPos[1], cameraPos[2]);

	return 1;
}
//########################################################################################

//========================================================================================
/*
*	Shader Functions:
*/
//========================================================================================

/*
*	Initialize the fields of a shader object using a glsl shader file
*
*	Params:
*		s: the shader struct into which the info will be loaded
*		file: the file path (relative or absolute) where the shader program is defined 
*		type: the type of shader (e.g vertex,fragment, tesselation...)
*/
void createShader(Shader &s, string file, GLenum type)
{
	s.fileName = file;
	compileShader(s.shaderID, file, type);
	s.type = GL_VERTEX_SHADER;
	s.program = 0;
}

/*
* Delete a shader struct
*	
*	Params: 
*		s: the shader struct to delete
*/
void deleteShader(Shader &s)
{
	glUseProgram(0);
	glDeleteShader(s.shaderID);
	s.program = 0;
}

/*
* Compile a glsl file and generate an OpenGL shading program on teh GPU
*
*	Params:	
*		shader: where the shader ID will be returned
*		filename: file path to the glsl program definition
*		type: the type of shader (e.g vertex,fragment, tesselation...)
*/
void compileShader(GLuint &shader, string &filename, GLenum shaderType)
{
	string source = loadSourceFile(filename);
	const GLchar* s_ptr = source.c_str();//get raw c string (char array)

	shader = glCreateShader(shaderType);//create shader on GPU
	glShaderSource(shader, 1, &s_ptr, NULL);//set shader program source

	glCompileShader(shader);


	//verify compilation
	GLint status;
	glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
	if(status!=GL_TRUE)
	{
		cout << "Shader compilation error. Could not compile: "
		<< filename << "\nShader type: "
		<< shaderType
		<<endl;

		GLint length;
		glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);

		string log(length, ' ');
		glGetShaderInfoLog(shader, log.length(), &length, &log[0]);

		cerr<< endl << source <<endl;
		cerr << endl << log <<endl;
	}
}

/*
* Copy a string from a file into a a string
*	
*	Params:
*		filepath: path to the file
*
*	Return: a string that is the copy of the source file
*/
string loadSourceFile(string &filepath)
{
	string source;

	ifstream input(filepath.c_str());
	if (input) {
		copy(istreambuf_iterator<char>(input),
			istreambuf_iterator<char>(),
			back_inserter(source));
		input.close();
	}

	else {
		cerr << "ERROR: Could not load shader source from file: "
			<< filepath << endl;
	}

	return source;
}
//########################################################################################


//========================================================================================
/*
*	Geometry Functions:
*/
//========================================================================================

//TODO: comment this
void createGeometry(Geometry &g, vector<vec3> vertices,  vector<vec3> normals, 
	vector<vec2> uvs, vector<uint> indices)
{
	//set vertex info
	glEnableVertexAttribArray(0);
	glGenBuffers(1, &(g.vertexBuffer));
	glBindBuffer(GL_ARRAY_BUFFER, g.vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(vec3),
		vertices.data(), GL_DYNAMIC_DRAW);

	//set normals info
	glEnableVertexAttribArray(1);
	glGenBuffers(1, &g.normalsBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, g.normalsBuffer);
	glBufferData(GL_ARRAY_BUFFER, normals.size()*sizeof(vec3),
		normals.data(), GL_DYNAMIC_DRAW);

	//set texture coordinates info
	glEnableVertexAttribArray(2);
	glGenBuffers(1, &g.uvBuffer);
	glBindBuffer(GL_ARRAY_BUFFER, g.uvBuffer);
	glBufferData(GL_ARRAY_BUFFER, uvs.size()*sizeof(vec2),
		uvs.data(), GL_DYNAMIC_DRAW);

	//set element info
	glGenBuffers(1, &(g.elmentBuffer));
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g.elmentBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertices.size()*sizeof(uint),
		indices.data(), GL_DYNAMIC_DRAW);

	//Init VAO
	glGenVertexArrays(1, &(g.vertexArray));
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	g.vertices=vertices;
	g.normals=normals;
	g.uvs=uvs;
	g.indices=indices;
}

/*
*	Initialize the fields of a geometry object using arrays
*
*	Params:
*		g: the geometry struct into which the info will be loaded
*		vertices: the vertex info of the geometry
*		indices: the index information (non-sequencial association of vertices) 
*/
void createGeometry(Geometry &g, vector<vec3> vertices, vector<uint> indices)
{
	//set vertex info
	glEnableVertexAttribArray(0);
	glGenBuffers(1, &(g.vertexBuffer));
	glBindBuffer(GL_ARRAY_BUFFER, g.vertexBuffer);
	glBufferData(GL_ARRAY_BUFFER, vertices.size()*sizeof(vec3),
		vertices.data(), GL_DYNAMIC_DRAW);

	//set element info
	glGenBuffers(1, &(g.elmentBuffer));
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, g.elmentBuffer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, vertices.size()*sizeof(uint),
		indices.data(), GL_DYNAMIC_DRAW);

	//set normals info
	glEnableVertexAttribArray(1);
	glGenBuffers(1, &g.normalsBuffer);

	//set normals info
	glEnableVertexAttribArray(2);
	glGenBuffers(1, &g.uvBuffer);

	//Init VAO
	glGenVertexArrays(1, &(g.vertexArray));
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	g.vertices=vertices;
	g.indices=indices;
}

/*
*	Initialize the fields of a default geometry struct
*
*	Params:
*		g: the geometry struct into which the info will be loaded
*/
void createGeometry(Geometry &g)
{
	//set vertices
	glEnableVertexAttribArray(0);
	glGenBuffers(1, &(g.vertexBuffer));

	//set normals
	glEnableVertexAttribArray(1);
	glGenBuffers(1, &g.normalsBuffer);

	//set uvs
	glEnableVertexAttribArray(2);
	glGenBuffers(1, &g.uvBuffer);

	//set indices
	glGenBuffers(1, &(g.elmentBuffer));

	//init VAO
	glGenVertexArrays(1, &(g.vertexArray));
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}

/*
* Delete a geometry struct
*	
*	Params: 
*		g: the geometry struct to delete
*/
void deleteGeometry(Geometry &g)
{
	glBindVertexArray(0);
	glDeleteVertexArrays(1, &(g.vertexArray));

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glDeleteBuffers(1, &(g.vertexBuffer));

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glDeleteBuffers(1, &(g.elmentBuffer));
}
//########################################################################################

//========================================================================================
/*
*	Texture Functions:
*/
//========================================================================================

/*
*	Initialize the fields of a textureobject using arrays
*
*	Params:
*		texture: a pointer to a texture struct into which the info will be loaded
*		filename: the filepath to the texture file
*		target: the OpenGL texture target (e.g 2D, rectangle...)
*
*	Return: a boolean value indicating whether an error ocurred (true means no error)
*/
bool createTexture(Texture &texture, const char* filename, GLuint target)
{
	int numComponents;
	stbi_set_flip_vertically_on_load(true);
	void *data = stbi_load(filename, &texture.width, &texture.height, &numComponents, 0);
	if (data != nullptr)
	{
		texture.target = target;
		glGenTextures(1, &texture.textureID);
		glBindTexture(texture.target, texture.textureID);
		GLuint format = numComponents == 3 ? GL_RGB : GL_RGBA;
		//cout << numComponents << endl;
		glTexImage2D(texture.target, 0, format, texture.width, texture.height, 0, format, GL_UNSIGNED_BYTE, data);

		// Note: Only wrapping modes supported for GL_TEXTURE_RECTANGLE when defining
		// GL_TEXTURE_WRAP are GL_CLAMP_TO_EDGE or GL_CLAMP_TO_BORDER
		glTexParameteri(texture.target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(texture.target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(texture.target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(texture.target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

		// Clean up
		glBindTexture(texture.target, 0);
		stbi_image_free(data);

		return true;
	}
	else
	{
		cerr << "Problem when loading texture" << endl;
	}
	return false; //error
}

/*
* Delete a texture struct
*	
*	Params: 
*		texture: the texture struct to delete
*/
void DestroyTexture(Texture &texture)
{
	glBindTexture(texture.target, 0);
	glDeleteTextures(1, &texture.textureID);
}

//########################################################################################

//========================================================================================
/*
*	Core functions:
*/
//========================================================================================

//TODO: verify that the following are up to date and well refactored
/*
* The following functions are not final at all, if modifications can be done, do them
*/

void ellipse(vector<vec3> &vertices, vector<uint> &indices, vector<vec3> &normals, float a, float b, float c)
{
	for(uint i=0; i<100; i++)
	{
		float v=(i/99.f)*M_PI;
		for(uint j=0; j<100; j++)
		{
			float u=(j/99.f)*2*M_PI;
			vec3 normal = vec3(a*cos(u)*sin(v), b*sin(u)*sin(v), c*cos(v));
			vertices.push_back(normal);
			indices.push_back(i*100+j);
			indices.push_back((i+1)*100+j);
			indices.push_back((i)*100+j+1);
			indices.push_back((i+1)*100+j+1);
			normals.push_back(normalize(normal));
		}
	}
}

dvec3 lerp(dvec3 p1, dvec3 p2, double t)
{
	return (1-t)*p1+t*p2;
}

dvec3 slerp(dvec3 p1, dvec3 p2, double t)
{
	dvec3 u1 = normalize(p1);
	dvec3 u2 = normalize(p2);

	double omega = acos(dot(p1,p2));

	double c1=sin((1-t)*omega)/(sin(omega));
	double c2=sin(t*omega)/sin(omega);

	return c1*p1+c2*p2;
}

vec3 rectangle_to_sphere(vec2 p, float r)
{
	float u=p[0], v=p[1]; 
	return vec3(r*cos(u)*sin(v), r*sin(u)*sin(v), r*cos(v));
}

vector<dvec3> subdivision(vector<dvec3> points, dvec3(*interp)(dvec3,dvec3,double))
{
	vector<dvec3> new_shape;

	//point duplication
	for(dvec3 point : points)
	{
		new_shape.push_back(point);
		new_shape.push_back(point);
	}

	int n=new_shape.size();

	//G_0
	for(uint i=0; i<n/2; i++)
	{
		new_shape[2*i]=interp(new_shape[2*i], interp(new_shape[((2*i-1)%n+n)%n], new_shape[(2*i+1)%n], 0.5), 0.5);
	}

	//G_1
	for(uint i=0; i<n/2; i++)
	{
		new_shape[2*i+1]=interp(new_shape[2*i+1], interp(new_shape[2*i], new_shape[(2*i+2)%n], 0.5), 0.5);
	}

	return new_shape;
}

void elliptical_P_Decomposition(vector<dvec3> fine, vector<double> w, 
	vector<dvec3> *coarse, vector<dvec4> *details, 
	dvec3(*interp)(dvec3,dvec3,double))
{
	int m=fine.size();
	for(int j=w.size()-1; j>=0; j--)
	{
		if((j&1)==0)
			for(int i=0; i<=m-2; i+=2)
			{
				dvec3 mid = interp(fine[((i-1)%m+m)%m], fine[(i+1)%m], 0.5f);
				fine[i] = interp(fine[i],mid,w[j]/(w[j]-1.f));
			}

		else
			for(int i=1; i<=m-2; i+=2)
			{
				dvec3 mid = interp(fine[((i-1)%m+m)%m], fine[(i+1)%m], 0.5f);
				fine[i] = interp(fine[i],mid,w[j]/(w[j]-1.f));
			}
	}
	for(int i=0; i<=m-2;i+=2)
	{
		dvec3 mid = interp(fine[i], fine[(i+2)%m], 0.5);
		/*(*coarse)[i/2]=*/(*coarse).push_back(fine[i]);
		/*(*details)[i/2]=*/(*details).push_back(dvec4(cross(fine[i],mid),acos(dot((fine[i]), (fine[(i+1)%m])))));
	}
}

void elliptical_P_reconstruction(vector<dvec3> *fine, vector<double> w, 
	vector<dvec3> coarse, vector<dvec4> details, 
	dvec3(*interp)(dvec3,dvec3,double))
{
	int n = coarse.size();

	for(int i=0; i<=n-1; i++)
	{
		dmat4 temp = dmat4(1);
		(*fine)[2*i]=coarse[i];
		(*fine)[2*i+1]= dvec3(rotate(temp, details[i][3], dvec3(details[i][0], details[i][1],details[i][2]))
			* dvec4(interp(coarse[i], coarse[(i+1)%n],1/2),1));
	}
	int l = w.size();
	int f = (*fine).size();
	for(int j=0; j<l-1; j++)
	{
		if(j & 1 ==0 )
		{
			for(int i=0; i<= 2*n-2; i+=2)
			{
				dvec3 mid = interp((*fine)[((i-1)%f+f)%f], 
					(*fine)[(i+1)%f],0.5);
				(*fine)[i]=interp((*fine)[i], mid, w[j]);
			}
		}
		else
		{
			for(int i=1; i<=2*n-1; i+=2)
			{
				dvec3 mid = interp((*fine)[((i-1)%f+f)%f], 
					(*fine)[(i+1)%f],0.5);
				(*fine)[i]=interp((*fine)[i], mid, w[j]);
			}
		}
	}
}

vector<double> weights = {0.5};
vector<dvec3> holder;
//main render loop
void dtof(vector<dvec3> ds, vector<vec3> &fs)
{
	fs = vector<vec3>(ds.size());
	for(uint i=0; i<ds.size(); i++)
	{
		fs[i]=vec3(ds[i]);
	}
}
void ftod(vector<vec3> fs, vector<dvec3> &ds)
{
	ds = vector<dvec3>(fs.size());
	for(uint i=0; i<fs.size(); i++)
	{
		ds[i]=dvec3(fs[i]);
	}
}
void render_loop(GLFWwindow* window)
{
	ellipse(shapes[0].vertices, shapes[0].indices, shapes[0].normals, 1.f,1.f,1.f);

	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.5, 1), 1));
	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.1, 1), 1));
	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.1, 1.5), 1));
	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.5, 1.5), 1));

	//shapes[1].vertices = fine_points;

	loadGeometryArrays(programs[0], shapes[0]);
	loadGeometryArrays(programs[0], shapes[1]);

	ftod(shapes[1].vertices, holder);
    while (!glfwWindowShouldClose(window))
	{
		/*if(temp)
		{
			shapes[1].vertices = subdivision(shapes[1].vertices, slerp);
			loadGeometryArrays(programs[0], shapes[1]);
		}*/

		if(!loadViewProjMatrix(cam, programs[0]))
		{
			cerr << "Error when loading projection matrix!" << endl;
			return;
		}
		glClearColor(0.7f, 0.7f, 0.7f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		loadColor(vec4(0,0.5,0.9,1), programs[0]);
		//loadTexture(programs[0], textures[0]);
		//render(programs[0], shapes[0], GL_TRIANGLE_STRIP);

		glDisable(GL_DEPTH_TEST);
		loadColor(vec4(1.f,0,0,1), programs[0]);
		render(programs[0], shapes[1], GL_POINTS);
		glEnable(GL_DEPTH_TEST);
		
		glfwPollEvents();
		glfwSwapBuffers(window);
	}
}

//cleanup
void end_rendering(GLFWwindow* window)
{
    	//Cleanup
	for(Shader s: shaders)
        deleteShader(s);
    for(GLuint p: programs)
        glDeleteProgram(p);
    for(Geometry g: shapes)
        deleteGeometry(g);

    glfwDestroyWindow(window);
    glfwTerminate();
}
//########################################################################################