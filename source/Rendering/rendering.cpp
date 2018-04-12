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
#include <random>

#define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stb/stb_image_write.h>

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

dvec3 (*INTERP)(dvec3, dvec3, double) = splerp;
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


/*
* The following functions are not final at all, if modifications can be done, do them
*/
Graph::Graph(vector<vec3> *vertices, vector<uint> *indices)
{
	graph = vector<vector<uint>>(vertices->size());
	nodes = (*vertices);

	uint n=indices->size();
	for(uint i=0; i<n; i+=2)
	{
		uint n1 = (*indices)[(i)];
		uint n2 = (*indices)[(i+1)%n];
		bool seen = false;
		for(uint j=0; j<graph[n2].size(); j++)
			if(n1 == graph[n2][j] || n1==n2)
				seen = true;
		if(!seen)
			graph[n2].push_back(n1);

		seen = false;
		for(uint j=0; j<graph[n1].size(); j++)
			if(n2 == graph[n1][j] || n1==n2)
				seen = true;
		if(!seen)
			graph[n1].push_back(n2);
	}
}

uint getMin(vector<double> *distances, vector<bool> *visited)
{
	uint min = -1;
	for(uint i=0; i<distances->size(); i++)
		if(!(*visited)[i])
		{
			min = i;
			break;
		}

	if(min == -1)
		return -1;

	for(uint i=0; i<distances->size(); i++)
	{
		if((*distances)[i] <= (*distances)[min] && !(*visited)[i])
			min = i;
	}

	return min;
}

void Graph::djikstra(uint start)
{
	uint toVisit = graph.size();
	lengths=vector<double>(graph.size());
	vector<bool> visited = vector<bool>(graph.size());
	
	for(uint i=0; i<lengths.size(); i++)
		lengths[i]=1.f/0.f;

	for(uint i=0; i<lengths.size(); i++)
		visited[i]=false;
	
	lengths[start] = 0;

	uint current = start;

	while(toVisit>0)
	{
		for(uint loc_n=0; loc_n<graph[current].size(); loc_n++)
		{
			uint neighbour = (graph[current][loc_n]);
			vec3 edge = (nodes[current])-nodes[neighbour];
			double edge_length = (double)length(edge);

			double n_length = lengths[current]+edge_length;
			//cout << "current: " << current << endl;
			//cout << "Length: " << n_length << " Neighbour: " << lengths[neighbour] <<endl;
			if(n_length <= lengths[neighbour])
			{
				//cout << "check" << endl;
				lengths[neighbour] = n_length;
			}
		}
		visited[current]=true;
		current = getMin(&lengths, &visited);
		toVisit--;
		if(current == -1)
		{
			toVisit = 0;
		}

	}
}

void Graph::toString()
{
	for(uint i=0; i<graph.size(); i++)
	{
		cout << i << ": ";
		for(uint j=0; j<graph[i].size(); j++)
		{
			//cout << j << " ";
			if(graph[i][j] >= 10000)
				cout << "too big" << i << " " << graph[i][j] << endl;
			
			cout << graph[i][j] << ", ";

		}
		cout << endl;
	}

	for(uint i=0; i<nodes.size(); i++)
		cout << "index " << i << ": (" << nodes[i].x << ", "<< nodes[i].y << ", "<< nodes[i].z <<")" << endl;
}

double Graph::node_length(uint node)
{
	return lengths[node];
}

uint inline cap(uint num, uint cap)
{
	if(num > cap-1)
		return cap-1;

	return num;
}

int resolution = 100;
void ellipse(vector<vec3> &vertices, vector<uint> &indices, vector<vec3> &normals, float a, float b, float c)
{
	for(uint i=0; i<resolution; i++)
	{
		float v=(i/((float)resolution-1))*M_PI;
		for(uint j=0; j<resolution; j++)
		{
			float u=(j/((float)resolution-1))*2*M_PI;
			vec3 normal = vec3(a*cos(u)*sin(v), b*sin(u)*sin(v), c*cos(v));
			vertices.push_back(normal);

			indices.push_back(cap(i*resolution+j, pow(resolution, 2)));
			indices.push_back(cap((i+1)*resolution+j, pow(resolution, 2)));

			indices.push_back(cap((i+1)*resolution+j, pow(resolution, 2)));
			indices.push_back(cap((i+1)*resolution+j+1, pow(resolution, 2)));

			indices.push_back(cap((i+1)*resolution+j+1, pow(resolution, 2)));
			indices.push_back(cap(i*resolution+j, pow(resolution, 2)));

			indices.push_back(cap(i*resolution+j, pow(resolution, 2)));
			indices.push_back(cap(i*resolution+j+1, pow(resolution, 2)));

			indices.push_back(cap(i*resolution+j+1, pow(resolution, 2)));
			indices.push_back(cap((i+1)*resolution+j+1, pow(resolution, 2)));

			indices.push_back(cap(i*resolution+j+1, pow(resolution, 2)));
			indices.push_back(cap((i+1)*resolution+j, pow(resolution, 2)));

			/*indices.push_back((i*100+j)%(100*100));
			indices.push_back(((i+1)*100+j)%(100*100));

			indices.push_back(((i)*100+j+1)%(100*100));
			indices.push_back((i*100+j)%(100*100));

			indices.push_back(((i+1)*100+j)%(100*100));
			indices.push_back(((i)*100+j+1)%(100*100));*/

			normals.push_back(normalize(normal));
		}
	}
}

void project_line(vector<vec3> &vertices, vec2 p1, vec2 p2, float a, float b, float c)
{
	float step = 0.01;
	float t = 0;
	for(uint i=0; i<100; i++)
	{
		float f0 = mix(p1[0],p2[0], t);
		float f1 = mix(p1[1],p2[1], t);

		float theta = abs(f0) < 0.0001? 0 : f1/atan(f1/f0);
		float r = f0/cos(theta);

		cout << "theta: " << theta << endl;
		cout << "r: " << r << endl;

		vertices.push_back(vec3(a*cos(theta)*sin(r), b*sin(theta)*sin(r), c*cos(r)));

		t+=step;
	}
}

dvec3 ellipse_project(dvec3 p, float a, float b, float c)
{
	double vx=p[0], vy=p[1], vz=p[2];

	vx*=vx, vy*=vy, vz*=vz;

	a*=a, b*=b, c*=c;

	vx=vx/a, vy=vy/b, vz=vz/c;

	double t = 1.d/sqrt(vx+vy+vz);

	return t*p;
}

dvec3 sphere_project(dvec3 p, float r)
{
	double vx=p[0], vy=p[1], vz=p[2];
	
	vx*=vx, vy*=vy, vz*=vz;

	double t = r/(vx+vy+vz);

	return t*p;
}

dvec3 mlerp(dvec3 p1, dvec3 p2, double t)
{
	return (1-t)*p1+t*p2;
}

dvec3 plerp(dvec3 p1, dvec3 p2, double u)
{
	dvec3 n = normalize(p1);
	dvec3 v1 = cross(p1,p2);
	dvec3 v2 = cross(n,v1);

	v1 = normalize(v1), v2=normalize(v2);

	double t = dot(n,p1)/dot(n,p2);

	dvec3 projected = t*p2;

	dvec3 new_p = mlerp(p1,p2, u);

	return ellipse_project(new_p, a_axis, b_axis, c_axis);
}

dvec3 slerp(dvec3 p1, dvec3 p2, double t)
{
	dvec3 u1 = normalize(p1);
	dvec3 u2 = normalize(p2);

	double omega = acos(dot(u1,u2));

	if(abs(omega) < 0.00000001 || isnan(omega))
		return p1;

	double c1=sin((1-t)*omega)/(sin(omega));
	double c2=sin(t*omega)/sin(omega);

	if(isnan(c1) || isnan(c2))
	{
		cout << "nan error" << endl;
	}

	return c1*p1+c2*p2;
}

dvec3 splerp(dvec3 p1, dvec3 p2, double t)
{
	double r = glm::max(a_axis, b_axis);
	r = glm::max(r, c_axis);

	p1 = sphere_project(p1, r);
	p2 = sphere_project(p2, r);

	dvec3 pi = slerp(p1,p2,t);

	pi = ellipse_project(pi, a_axis, b_axis, c_axis);

	return pi;
}

dvec3 ilerp(dvec3 p1, dvec3 p2, double t)
{
	double big = std::max(a_axis,b_axis);
	big = std::max(big, c_axis);

	p1.x = big/a_axis*p1.x;
	p1.y = big/b_axis*p1.y;
	p1.z = big/c_axis*p1.z;

	p2.x = big/a_axis*p2.x;
	p2.y = big/b_axis*p2.y;
	p2.z = big/c_axis*p2.z;

	dvec3 result = slerp(p1, p2, t);

	result.x = result.x * a_axis/big;
	result.y = result.y * b_axis/big;
	result.z = result.z * c_axis/big;

	return result;
}

dvec3 seedSeeking(dvec3 p, double limit, double *dret,
	double(*implicit)(double x, double y, double z),
	dvec3(*gradient)(double x, double y, double z))
{
	dvec3 dP;
	double distance = 0;
	do 
	{
		dvec3 grad = gradient(p.x, p.y, p.z);
		dP = (implicit(p.x,p.y,p.z)*grad);
		if (dot(grad, grad) > 0)
			dP = dP / dot(grad, grad);
		else
			return p;

		p -= dP;
		distance += length(dP);
		*dret = distance;

		cout << length(dP) << endl;
		cout << dP.x << ", " << dP.y << ", "<< dP.z << endl;
		cout << p.x << ", " << p.y << ", "<< p.z << endl << endl;
	} while (length(dP) > 0.0001);

	return p;
	
}

dvec3 goal = dvec3(0);
double surfaceDistance(double x, double y, double z)
{
	double value = x*x/(a_axis*a_axis) + y*y/(b_axis*b_axis) + z*z/(c_axis+c_axis) -1.0; 

	value + length(dvec3(x,y,z)-goal)*length(dvec3(x,y,z)-goal);

	return value;
}

dvec3 surfaceGrad(double x, double y, double z)
{
	double dx = 2*x/(a_axis*a_axis) + 2*x;
	double dy = 2*y/(b_axis*b_axis) + 2*y;
	double dz = 2*z/(c_axis*c_axis) + 2*z;

	return dvec3(dx,dy,dz);
}

dvec3 glerp(dvec3 p1, dvec3 p2, double t)
{
	if(t<0 || t>1)
		cout << "Here: " << t <<endl;
	goal = p2;

	double distance;

	seedSeeking(p1, 1.f/0.f, &distance, surfaceDistance, surfaceGrad);

	dvec3 ret = seedSeeking(p1, distance*t, &distance, surfaceDistance, surfaceGrad);

	return ret;
}

vec3 rectangle_to_sphere(vec2 p, float a, float b, float c)
{
	float u=p[0], v=p[1]; 
	return vec3(a*cos(u)*sin(v), b*sin(u)*sin(v), c*cos(v));
}

vector<dvec3> subdivision(vector<dvec3> points, dvec3(*interp)(dvec3,dvec3,double))
{
	vector<dvec3> new_shape;

	//point duplication
	for(uint i=0; i<points.size() ;i++)
	{
		new_shape.push_back(points[i]);
		new_shape.push_back(interp(points[i], points[(i+1)%points.size()], 0.5));
	}

	int n=new_shape.size();

	//G_0
	//if('SCHEME' == 'D')
	{
		for(uint i=0; i<n-1; i+=2)
		{
			new_shape[i]=interp(new_shape[i], interp(new_shape[(i-1+n)%n], new_shape[(i+1)%n], 0.5), 0.5);
		}
	}

	for(uint i=0; i<new_shape.size(); i++)
	{
		new_shape[i] = ellipse_project(new_shape[i], a_axis, b_axis, c_axis);
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
				dvec3 mid = interp(fine[(i-1+m)%m], fine[(i+1)%m], 0.5f);
				fine[i] = interp(fine[i],mid,w[j]/(w[j]-1.f));
			}

		else
			for(int i=1; i<=/*m-2*/m-1; i+=2)
			{
				dvec3 mid = interp(fine[(i-1+m)%m], fine[(i+1)%m], 0.5f);
				fine[i] = interp(fine[i],mid,w[j]/(w[j]-1.f));
			}
	}
	for(int i=0; i<=m-2;i+=2)
	{
		dvec3 mid = interp(fine[i], fine[(i+2)%m], 0.5);
		(*coarse).push_back(fine[i]);
		dvec4 dets = dvec4(cross(mid,fine[(i+1)%m]),acos(dot((fine[(i+1)%m]), (mid))));
		if(length(vec3(dets)) <= 0.1)
			dets=dvec4(1,1,1,0);
		(*details).push_back(dets);
	}
}

void elliptical_P_reconstruction(vector<dvec3> *fine, vector<double> w, 
	vector<dvec3> coarse, vector<dvec4> &details, 
	dvec3(*interp)(dvec3,dvec3,double))
{
	int n = coarse.size();

	int offset = details.size() - n;

	for(int i=0; i<=n-1; i++)
	{
		dmat4 temp = dmat4(1);
		(*fine)[2*i]=coarse[i];
		dvec4 det_vec = details[i+offset]; 

		(*fine)[2*i+1]= dvec3(rotate(temp, det_vec[3], dvec3(det_vec[0], det_vec[1], det_vec[2]))
			* dvec4(interp(coarse[i], coarse[(i+1)%n],0.5),1));
	}
	for(int i=0; i<=n-1; i++)
		details.pop_back();
	int l = w.size();
	int f = (*fine).size();
	for(int j=0; j<=l-1; j++)
	{
		if((j & 1) ==0)
		{
			for(int i=0; i<= 2*n-2; i+=2)
			{
				dvec3 mid = interp((*fine)[(i-1+f)%f], (*fine)[(i+1)%f],0.5);
				(*fine)[i]=interp((*fine)[i], mid, w[j]);
			}
		}
		else
		{
			for(int i=1; i<=f-1; i+=2)
			{
				dvec3 mid = interp((*fine)[(i-1+f)%f], (*fine)[(i+1)%f],0.5);
				(*fine)[i]= interp((*fine)[i], mid, w[j]);
			}
		}
		
	}
}

void elliptical_D_Decomposition(vector<dvec3> fine, vector<double> w, 
	vector<dvec3> *coarse, vector<dvec4> *details, 
	dvec3(*interp)(dvec3,dvec3,double))
{
	int m=fine.size();
	for(int j=w.size()-1; j>=0; j--)
	{
		if((j&1)==0)
			for(int i=1; i<=m-1; i+=2)
			{
				dvec3 p = fine[i];
				fine[i] = interp(p, fine[(i+1)%m], w[j]/(2.d*w[j]-2.d));
				fine[(i+1)%m] = interp(fine[(i+1)%m],p,w[j]/(2.d*w[j]-2.d));
			}

		else
			for(int i=0; i<=m-2; i+=2)
			{
				dvec3 p = fine[i];
				fine[i] = interp(p, fine[(i+1)%m], w[j]/2.d);
				fine[i] = interp(fine[(i+1)%m],p,w[j]/(2*w[j]-2.d));
			}
	}
	for(int i=0; i<=m-2;i+=2)
	{
		(*coarse).push_back(interp(fine[i],fine[(i+1)%m], 0.5));
		dvec4 dets = dvec4(cross(fine[i],fine[(i+1)%m]),acos(dot((fine[(i+1)%m]), fine[i]))/0.5d);
		if(length(vec3(dets)) <= 0.1)
			dets=dvec4(1,1,1,0);
		(*details).push_back(dets);
	}
}

void elliptical_D_reconstruction(vector<dvec3> *fine, vector<double> w, 
	vector<dvec3> coarse, vector<dvec4> &details, 
	dvec3(*interp)(dvec3,dvec3,double))
{
	int n = coarse.size();

	int offset = details.size() - n;

	for(int i=0; i<=n-1; i++)
	{
		dmat4 temp = dmat4(1);
		dvec4 det_vec = details[i+offset]; 

		(*fine)[2*i]= dvec3(rotate(temp, -det_vec[3], dvec3(det_vec[0], det_vec[1], det_vec[2]))
		* dvec4(coarse[i],1));

		(*fine)[2*i+1]= dvec3(rotate(temp, det_vec[3], dvec3(det_vec[0], det_vec[1], det_vec[2]))
			* dvec4(coarse[i],1));
	}
	for(int i=0; i<=n-1; i++)
		details.pop_back();
	int l = w.size();
	int f = (*fine).size();
	for(int j=0; j<=l-1; j++)
	{
		if((j & 1) ==0)
		{
			for(int i=1; i<= 2*n-1; i+=2)
			{
				dvec3 p = (*fine)[i];
				(*fine)[i] = interp(p, (*fine)[(i+1)%f],w[j]/2.d);
				(*fine)[(i+1)%f] = interp((*fine)[(i+1)%f], p, w[j]/2.d);
			}
		}
		else
		{
			for(int i=0; i<=2*n-2; i+=2)
			{
				dvec3 p = (*fine)[i];
				(*fine)[i] = interp(p, (*fine)[(i+1)%f],w[j]/2.d);
				(*fine)[(i+1)%f]= interp((*fine)[i+1], p, w[j]/2.d);
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

uint inline getIndex(vector<vec3> *vertices, vec3 v)
{
	double cDistance = 1.f/0.f;
	uint minIndex = 0;

	for(uint i=0; i<vertices->size(); i++)
	{
		double lDistance = length(v-(*vertices)[i]);
		if(lDistance < cDistance)
		{
			minIndex = i;
			cDistance = lDistance;
		}
	}

	return minIndex;
}

double calculateAverageLength(vector<vec3> vertices, Graph g)
{
	double sum = 0;
	uint n = vertices.size();

	vector<double> values;

	for(uint i=0; i<n; i++)
	{
		uint n1 = getIndex(&g.nodes, vertices[i]);
		uint n2 = getIndex(&g.nodes, vertices[(i+1)%n]);
		g.djikstra(n1);
		double v = g.node_length(n2);
		sum+=v;
		values.push_back(v);

		if(i%10==0)
			cout << "Evaluating distance average:" << ((float)i/(float)n)*100 << "\%" <<endl;
	}

	sort(values.begin(), values.end());

	double median;
	if(values.size()%2 == 0)
		median = (values[values.size()/2]+values[values.size()/2+1])/2.f;
	else
		median = values[values.size()/2 + 1];

	sum = sum/vertices.size();

	double stdv = 0;
	for(uint i=0; i<values.size(); i++)
	{
		stdv+=pow(values[i]-sum, 2);
	}
	stdv = sqrt(stdv);
	cout << "Average geodesic distance: " << sum << endl;
	cout << "Median geodesic distance: " << median << endl; 
	cout << "Standard deviation: " << stdv << endl;


	return sum;
}

vector<vec3> uniform_points(uint size)
{
	vector<vec3> return_val;
	std::default_random_engine generator (glfwGetTime());
	std::normal_distribution<double> distributionx(0.0,a_axis);
	std::normal_distribution<double> distributiony(0.0,b_axis);
	std::normal_distribution<double> distributionz(0.0,c_axis);

	for(uint i=0; i<size; i++)
	{
		float x,y,z;
		x = distributionx(generator);
		y = distributiony(generator);
		z = distributionz(generator);

		return_val.push_back(vec3(x,y,z));
	}

	for(uint i=0; i<size; i++)
	{
		return_val[i] = ellipse_project(return_val[i], a_axis, b_axis, c_axis);
	}

	return return_val;
}

int subdivisions = 10;
vector<double> connect_all(vector<vec3> o_set)
{
	vector<double> ret = vector<double>(o_set.size());
	for(uint i=0; i<o_set.size(); i++)
	{
		for(uint j=i+1; j<o_set.size(); j++)
		{
			vector<vec3> temp;
			for(uint k=0; k<=subdivisions; k++)
			{
				temp.push_back(INTERP(o_set[i], o_set[j], (float)k/(float)subdivisions));
			}

			double distance = 0;
			for(uint k=1; k<subdivisions; k++)
			{
				distance += length(temp[k]-temp[k-1]);
			}

			ret.push_back(distance);
		}
	}

	return ret;
}

vector<double> merge(vector<vector<double>> set)
{
	vector<double> ret;
	for(uint i=0; i<set.size(); i++)
	{
		for(uint j=i; j<set[i].size(); j++)
		{
			ret.push_back(set[i][j]);
		}
	}

	return ret;
}

Graph g;
std::ofstream outfile ("tests.csv");

string inline erpToString(dvec3 (*f)(dvec3, dvec3, double))
{
	if (f==mlerp)
	{
		return "Lerp average";
	}

	if (f==plerp)
	{
		return "Plane Projection average";
	}

	if (f==splerp)
	{
		return "Spherical Projection average";
	}

	if (f==ilerp)
	{
		return "Implicit Projection average";
	}
}

void stats(vector<double> data, double &average, double &median, double &std)
{
	sort(data.begin(), data.end());

	if(data.size()%2 ==0)
	{
		median = (data[data.size()/2]+data[(data.size()/2) +1])/2.d;
	}

	else
		median = (data[(data.size()/2) +1]);

	average = 0;
	for(int i=0; i<data.size(); i++)
		 average+=data[i];
	
	average = average/(double)data.size();

	std = 0;
	for(int i=0; i<data.size(); i++)
	{
		std += pow(average-data[i], 2);
	}
	std = sqrt(std/data.size());
}

vector<double> specialMeasure(vector<double> original, vector<double> current)
{
	vector<double> ret;
	for(int i=0; i<original.size(); i++)
		ret.push_back(current[i]/original[i]);

	return ret;
}

vector<double> calcGeodesics(vector<vec3> data)
{
	vector<double> distancesD;
	cout << "Calcualting geodesics" << endl;
	for(uint i=0; i<data.size(); i++)
	{
		uint n1 = getIndex(&g.nodes, data[i]);
		g.djikstra(n1);
		cout << "Evaluation: " << (float)i/(float)data.size() * 100 << "\%" << endl;
		for(uint j=i+1; j<data.size(); j++)
		{
			uint n2 = getIndex(&g.nodes, data[j]);
			distancesD.push_back(g.node_length(n2));
		}	
	}

	return distancesD;
}

vector<double> evaluate(vector<vec3> data)
{
	vector<double> ret;
	for(uint i=0; i<100; i++)
	{
		subdivisions = (i+1)*10;
		double average, median, std;
		vector<double> temp = connect_all(data);
		stats(temp, average, median, std);
		ret.push_back(average);
	}

	return ret;
}

vector<dvec3 (*)(dvec3, dvec3, double)> methods = {mlerp, plerp, splerp, ilerp};
void render_loop(GLFWwindow* window)
{
	ellipse(shapes[0].vertices, shapes[0].indices, shapes[0].normals, a_axis,b_axis,c_axis);

	g = Graph(&(shapes[0].vertices), &(shapes[0].indices));

	/*shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.5, 1), a_axis,b_axis,c_axis));
	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.1, 1), a_axis,b_axis,c_axis));
	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.1, 1.5), a_axis,b_axis,c_axis));
	shapes[1].vertices.push_back(rectangle_to_sphere(vec2(0.5, 1.5), a_axis,b_axis,c_axis));*/

	vector<vec3> data = (uniform_points(10));
	shapes[1].vertices = (data);

	INTERP = mlerp;
	vector<double> geodesics = calcGeodesics(data);
	double average, median, std;
	stats(geodesics, average, median, std);

	vector<vector<double>> tests = vector<vector<double>>(4);
	for(uint i=0; i<4; i++)
	{
		INTERP = methods[i];
		tests[i] = evaluate(data);
	}

	outfile << "Subdivisions, " << "Geodesic average," << erpToString(methods[0]) << ", " << erpToString(methods[1]) <<", "  
		<< erpToString(methods[2]) << ", "<< erpToString(methods[3]) <<endl;
	for(int i=0;  i<100; i++)
	{
		outfile << to_string(i*10)<< ", " << average << ", " << tests[0][i] << ", "<< tests[1][i]<< ", " << tests[2][i] << ", "
			<< tests[3][i] << endl;

	}
	//shapes[1].vertices = uniform_points(1000);

	//evaluate(data_lines);

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
		render(programs[0], shapes[0], GL_LINES);

		//glDisable(GL_DEPTH_TEST);
		/*vector<vec3> temp = shapes[1].vertices;
		for(uint i=0; i<shapes[1].vertices.size(); i++)
		{
			shapes[1].vertices[i] = vec3(ellipse_project(shapes[1].vertices[i],a,b,c));
		}*/
		loadGeometryArrays(programs[0], shapes[1]);
		loadColor(vec4(1.f,0,0,1), programs[0]);
		render(programs[0], shapes[1], GL_POINTS);
		glEnable(GL_DEPTH_TEST);

		//shapes[1].vertices = temp;
		
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