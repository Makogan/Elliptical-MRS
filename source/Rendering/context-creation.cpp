//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
/*
*	Author:	Camilo Talero
*
*
*	Version: 0.0.1
*
*	File defining all relevant OpenGL and GLFW related functions needed to create an 
*	OpenGL context and GLFW window.
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Includes and macros
*/
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "context-creation.hpp"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
*	Global Values
*/
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Camera cam;//Global camera, main camera used for rendering and perspective projection

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//========================================================================================
/*
*	Opengl Context Creation:
*/
//========================================================================================

/*
* Function to create the OpenGL context.
* 
* return GLFWwindow*: the pointer to the GLFW window containing the current context.
*/
GLFWwindow* create_context()
{
    GLFWwindow* window = createWindow();
	if(window==NULL)
	{
		cerr << "Error when creating GLFW window" << endl;
	}
    callBackInit(window);

    // glew initilization, this is so that the program is crossplatform,
    // also things won't work without it
    glewExperimental = GL_TRUE;
    glewInit(); glGetError();
    // An error will always be thrown when initializing glew.
    // It can be safely discarded so we call glGetError() to delete it and move on.
	
	//add 2 new shaders to the global list of shaders
	shaders.push_back(Shader());
	shaders.push_back(Shader());

	//Crete a vertex and a fragment shader. bare minimum for rendering
	createShader(shaders[0], "./Shaders/VertexShader.glsl", GL_VERTEX_SHADER);
	createShader(shaders[1], "./Shaders/FragmentShader.glsl", GL_FRAGMENT_SHADER);

	//Add a new program to the global shading programs list
	programs.push_back(glCreateProgram());
	//Attach both shaders to the program
	glAttachShader(programs[0], shaders[0].shaderID);
	glAttachShader(programs[0], shaders[1].shaderID);
	//Link the program to the current context
	glLinkProgram(programs[0]);
	//Set the program as the current activated shading program
	glUseProgram(programs[0]);
	//Update both shader structures to keep track of the program they have been attached to
	shaders[0].program=programs[0];
	shaders[1].program=programs[0];
	
	//TODO: this should not be here, delete once it's not needed
	createGeometry(shapes[0]);//Create a geometry object 
	createGeometry(shapes[1]);//Create a geometry object 
	createGeometry(shapes[2]);//Create a geometry object 
	createGeometry(shapes[3]);//Create a geometry object 
	
	//Create a new camera object with defined orientation, position, and dimensions
    int width, height;
    glfwGetWindowSize(window, &width, &height);
	cam = *(new Camera(mat3(-1), vec3(3,0,0), width, height));
	cam.setLookDirection(vec3(-1,0,0));

	//Set default OpenGL values for rendering
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
	glPointSize(10.f);
	glLineWidth(10.f);

	GLfloat LineRange[2];

	return window;
	return NULL;
}
//########################################################################################

//========================================================================================
/*
*	Error checking and debugging functions:
*/
//========================================================================================

/*
* Check for OpenGL errors and print the appropriate error message if needed.
*
* return int: the number of the generated error.
*/
int openGLerror()
{
	GLenum errorNum = glGetError();
	switch(errorNum)
	{
		case GL_NO_ERROR:
			break;
		case GL_INVALID_ENUM:
			cerr << "An unacceptable value is specified for an enumerated argument."
				<< endl;
			break;
		case GL_INVALID_VALUE:
			cerr << "A numeric argument is out of range."
				<< endl;
			break;
		case GL_INVALID_OPERATION:
			cerr << "The specified operation is not allowed in the current state."
				<< endl;
			break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			cerr << "The framebuffer object is not complete."
				<< endl;
			break;
		case GL_OUT_OF_MEMORY:
			cerr << "There is not enough memory left to execute the command."
				<< endl;
			break;
		default:
			cerr<< "Undefined error, you messed up big time now."
				<< endl;
	}
	return errorNum;
}
//########################################################################################

//========================================================================================
/*
*	Initialization:
*/
//========================================================================================
//Initialize GLFW callBack Functions
void callBackInit(GLFWwindow* window)
{
	//Set GLFW callback functions
	glfwSetErrorCallback(error_callback);
	glfwSetKeyCallback(window, key_callback);
	glfwSetCursorPosCallback(window, cursor_pos_callback);
	glfwSetMouseButtonCallback(window, mouse_button_callback);

}

/*
* Method to create a GLFW window, window will be maximized and decorated.
*
* return GLFWwindow*: a pointer to teh created window.
*/
GLFWwindow* createWindow()
{
	//Initialize GLFW
	if (!glfwInit())
	{
		cerr<< "Failed to initialize GLFW.\nTerminating program." << endl;
		return NULL;
	}

	//Get the primiray monitor of the current system's info
	const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);//OpenGL major version
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);//OpenGL minor version
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);//Set Forward compatibility
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);//Use GLFW defaults
	glfwWindowHint(GLFW_DECORATED, GLFW_TRUE);//Make the window decorated
	//Create a GLFW window with the main monitors width, reduced height, 
	//name, on windowed mode, not sharing resources with any context
	GLFWwindow* window = glfwCreateWindow(mode->width, mode->height-40, 
		"Voxel World", NULL, NULL);
	if (!window)//Check for errors
	{
		cerr<< "Failed to glfwCreateWindow.\nTerminating program." << endl;
		return NULL;
	}
	//TODO: delete or uncomment at one point
	//glfwMaximizeWindow(window);//Make the window maximized
	//Set the current window to be the current OpenGL context
	glfwMakeContextCurrent(window);

	return window;
}
//########################################################################################

//========================================================================================
/*
*	GLFW helper functions:
*/
//========================================================================================

float max(float a, float b)
{
	return a>=b? a : b;
}

int cursorSelectNode(GLFWwindow *window)
{
	/*double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);

	mat4 view= cam.getViewMatrix();
	mat4 proj= cam.getPerspectiveMatrix();

	uint count = 0;
	int width, height;
	glfwGetWindowSize(window, &width, &height);

	vec3 screenPos;
	vec3 pos;
	float depth;
	vec3 projCursor;*/
	return 0;
}
//########################################################################################

//========================================================================================
/*
*	GLFW callback functions:
*/
//========================================================================================

/*
* GLFW callback functions for Interrupt based event handling.
* More information at http://www.glfw.org/docs/latest/input_guide.html
*/

void error_callback(int error, const char* description)
{
    cout << "Error: " << description << endl;
}

void cursor_pos_callback(GLFWwindow* window, double xpos, double ypos)
{
	
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{

}

vector<dvec4> dets;
char scheme = 'D';
extern Graph g;
#define CAM_SPEED 0.1f
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GLFW_TRUE);

    else if(key == GLFW_KEY_F11 && action == GLFW_PRESS)
    {
    	//Get the primary monitor and the monitor attached to the current window
    	const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    	GLFWmonitor* monitor = glfwGetWindowMonitor(window);

    	//If no monitor is attached to the window attach the primary monitor
    	//This makes the window full screen
    	if(monitor==NULL)
    		glfwSetWindowMonitor(window, glfwGetPrimaryMonitor(), 0, 0,
    			mode->width, mode->height, mode->refreshRate);

    	//If a monitor is attached, detach it and set the dimesnions to that of
    	//The primary's monitor resolution. This makes the window windowed
    	else
    		glfwSetWindowMonitor(window, NULL, 0, 0,
    			mode->width, mode->height-1, mode->refreshRate);
    			//subtract 1 to prevent it from going into full screen mode

    	glfwMaximizeWindow(window);
    }

    else if(key == GLFW_KEY_R && action == GLFW_PRESS)
   	{
		   
   	}

    else if(key == GLFW_KEY_F12 && action == GLFW_PRESS)
    	cout << glfwGetVersionString() << endl;

    else if(key == GLFW_KEY_W)
    	cam.move(normalize(cam.getForward())*CAM_SPEED);

    else if(key == GLFW_KEY_S)
    	cam.move(-normalize(cam.getForward())*CAM_SPEED);

    else if(key == GLFW_KEY_A)
    	cam.move(normalize(cam.getSide())*CAM_SPEED);

    else if(key == GLFW_KEY_D)
		cam.move(-normalize(cam.getSide())*CAM_SPEED);

    else if(key == GLFW_KEY_Q)
		cam.move(normalize(cam.getUp())*CAM_SPEED);

    else if(key == GLFW_KEY_E)
		cam.move(-normalize(cam.getUp())*CAM_SPEED);

    else if(key == GLFW_KEY_KP_6)
    	cam.turnH(radians(-1.f));

    else if(key == GLFW_KEY_KP_4)
    	cam.turnH(radians(1.f));

    else if(key == GLFW_KEY_KP_8)
    	cam.turnV(radians(1.f));

    else if(key == GLFW_KEY_KP_2)
		cam.turnV(radians(-1.f));
	
	else if(key == GLFW_KEY_M && action == GLFW_PRESS)
	{
		for(uint i=0; i<holder.size(); i++)
		{
			float radius = std::max(a_axis, b_axis);
			radius = std::max((double)radius, c_axis);
			holder[i] = sphere_project(holder[i], radius);
		}
		holder = subdivision(holder, slerp);
		for(uint i=0; i<holder.size(); i++)
		{
			holder[i] = ellipse_project(holder[i], a_axis, b_axis, c_axis);
		}
		dtof(holder, shapes[1].vertices);
		loadGeometryArrays(programs[0], shapes[1]);
	}

	else if(key == GLFW_KEY_KP_ADD && action == GLFW_PRESS)
	{
		if(dets.size()>0)
		{
			vector<dvec3> temp=vector<dvec3>(holder.size()*2);
			for(uint i=0; i<holder.size(); i++)
			{
				float radius = std::max(a_axis, b_axis);
				radius = std::max((double)radius, c_axis);
				holder[i] = sphere_project(holder[i], radius);
			}
			if(scheme == 'D')
				elliptical_D_reconstruction(&temp, weights, 
					holder, dets, 
					INTERP);
			else
				elliptical_P_reconstruction(&temp, weights, 
					holder, dets, 
					INTERP);

			holder = temp;
			for(uint i=0; i<holder.size(); i++)
			{
				holder[i] = ellipse_project(holder[i], a_axis, b_axis, c_axis);
			}
			dtof(holder, shapes[1].vertices);
			loadGeometryArrays(programs[0], shapes[1]);
		}
	}

	else if(key == GLFW_KEY_KP_SUBTRACT && action == GLFW_PRESS)
	{
		//glDisable(GL_DEPTH_TEST);
		vector<dvec3> temp;
		for(uint i=0; i<holder.size(); i++)
		{
			float radius = std::max(a_axis, b_axis);
			radius = std::max((double)radius, c_axis);
			holder[i] = sphere_project(holder[i], radius);
		}
		if(scheme == 'D')
			elliptical_D_Decomposition(holder, weights, 
				&temp, &dets, 
				INTERP);
		else
			elliptical_P_Decomposition(holder, weights, 
				&temp, &dets, 
				INTERP);

		holder = temp;
		//loc_dets.insert(loc_dets.end(), dets.begin(), dets.end());
		for(uint i=0; i<holder.size(); i++)
		{
			holder[i] = ellipse_project(holder[i], a_axis, b_axis, c_axis);
		}
		dtof(holder, shapes[1].vertices);
		loadGeometryArrays(programs[0], shapes[1]);
	}

    else if(key == GLFW_KEY_KP_MULTIPLY)
    	cam.resetView();

    else if(key == GLFW_KEY_KP_DIVIDE)
		cam.resetCamera();
	
	else if(key == GLFW_KEY_PRINT_SCREEN && action == GLFW_PRESS)
	{
		vec3 p = cam.getPosition();
		cout << "("<< p[0] << ", " << p[1] << ", " << p[2] << ")" << endl;
	}

	else if(key == GLFW_KEY_C && action == GLFW_PRESS)
	{
		double a = calculateAverageLength(shapes[1].vertices, g);
		cout << endl;
	}
}

//########################################################################################