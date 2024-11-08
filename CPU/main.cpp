///////////////////////////////////////////////////////////////////////
//
// P3D Course
// (c) 2021 by Jo o Madeiras Pereira
//Ray Tracing P3F scenes and drawing points with Modern OpenGL
//
///////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <stdio.h>
#include <chrono>
#include <conio.h>

#include <GL/glew.h>
#include <GL/freeglut.h>
#include <IL/il.h>

#include "scene.h"
#include "rayAccelerator.h"
#include "maths.h"
#include "macros.h"

//Enable OpenGL drawing.  
bool drawModeEnabled = true;

bool P3F_scene = false; //choose between P3F scene or a built-in random scene
bool DEPTH_OF_FIELD = false;
bool ANTI_ALIASING = true;
bool hardS = false;

#define MAX_DEPTH 4  //number of bounces

#define CAPTION "Whitted Ray-Tracer"
#define VERTEX_COORD_ATTRIB 0
#define COLOR_ATTRIB 1
#define FUZZY_REFLECTIONS true

Grid* grid_ptr = NULL;
BVH* bvh_ptr = NULL;
accelerator Accel_Struct = BVH_ACC;

float off_x, off_y = 0;

unsigned int FrameCount = 0;

// Current Camera Position
float camX, camY, camZ;

//Original Camera position;
Vector Eye;

// Mouse Tracking Variables
int startX, startY, tracking = 0;

// Camera Spherical Coordinates
float alpha = 0.0f, beta = 0.0f;
float r = 4.0f;

// Frame counting and FPS computation
long myTime, timebase = 0, frame = 0;
char s[32];


// Points defined by 2 attributes: positions which are stored in vertices array and colors which are stored in colors array
float* colors;
float* vertices;
int size_vertices;
int size_colors;

//Array of Pixels to be stored in a file by using DevIL library
uint8_t* img_Data;

GLfloat m[16];  //projection matrix initialized by ortho function

GLuint VaoId;
GLuint VboId[2];

GLuint VertexShaderId, FragmentShaderId, ProgramId;
GLint UniformId;

Scene* scene = NULL;

int RES_X, RES_Y;

int WindowHandle = 0;

bool SCHLICK_APPROX = false;


/////////////////////////////////////////////////////////////////////// ERRORS

bool isOpenGLError() {
	bool isError = false;
	GLenum errCode;
	const GLubyte* errString;
	while ((errCode = glGetError()) != GL_NO_ERROR) {
		isError = true;
		errString = gluErrorString(errCode);
		std::cerr << "OpenGL ERROR [" << errString << "]." << std::endl;
	}
	return isError;
}

void checkOpenGLError(std::string error)
{
	if (isOpenGLError()) {
		std::cerr << error << std::endl;
		exit(EXIT_FAILURE);
	}
}

/////////////////////////////////////////////////////////////////////// SHADERs

const GLchar* VertexShader =
{
	"#version 430 core\n"

	"in vec2 in_Position;\n"
	"in vec3 in_Color;\n"
	"uniform mat4 Matrix;\n"
	"out vec4 color;\n"

	"void main(void)\n"
	"{\n"
	"	vec4 position = vec4(in_Position, 0.0, 1.0);\n"
	"	color = vec4(in_Color, 1.0);\n"
	"	gl_Position = Matrix * position;\n"

	"}\n"
};

const GLchar* FragmentShader =
{
	"#version 430 core\n"

	"in vec4 color;\n"
	"out vec4 out_Color;\n"

	"void main(void)\n"
	"{\n"
	"	out_Color = color;\n"
	"}\n"
};

void createShaderProgram()
{
	VertexShaderId = glCreateShader(GL_VERTEX_SHADER);
	glShaderSource(VertexShaderId, 1, &VertexShader, 0);
	glCompileShader(VertexShaderId);

	FragmentShaderId = glCreateShader(GL_FRAGMENT_SHADER);
	glShaderSource(FragmentShaderId, 1, &FragmentShader, 0);
	glCompileShader(FragmentShaderId);

	ProgramId = glCreateProgram();
	glAttachShader(ProgramId, VertexShaderId);
	glAttachShader(ProgramId, FragmentShaderId);

	glBindAttribLocation(ProgramId, VERTEX_COORD_ATTRIB, "in_Position");
	glBindAttribLocation(ProgramId, COLOR_ATTRIB, "in_Color");

	glLinkProgram(ProgramId);
	UniformId = glGetUniformLocation(ProgramId, "Matrix");

	checkOpenGLError("ERROR: Could not create shaders.");
}

void destroyShaderProgram()
{
	glUseProgram(0);
	glDetachShader(ProgramId, VertexShaderId);
	glDetachShader(ProgramId, FragmentShaderId);

	glDeleteShader(FragmentShaderId);
	glDeleteShader(VertexShaderId);
	glDeleteProgram(ProgramId);

	checkOpenGLError("ERROR: Could not destroy shaders.");
}

/////////////////////////////////////////////////////////////////////// VAOs & VBOs

void createBufferObjects()
{
	glGenVertexArrays(1, &VaoId);
	glBindVertexArray(VaoId);
	glGenBuffers(2, VboId);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);

	/* S  se faz a aloca  o dos arrays glBufferData (NULL), e o envio dos pontos para a placa gr fica
	  feito na drawPoints com GlBufferSubData em tempo de execu  o pois os arrays s o GL_DYNAMIC_DRAW */
	glBufferData(GL_ARRAY_BUFFER, size_vertices, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glVertexAttribPointer(VERTEX_COORD_ATTRIB, 2, GL_FLOAT, 0, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferData(GL_ARRAY_BUFFER, size_colors, NULL, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(COLOR_ATTRIB);
	glVertexAttribPointer(COLOR_ATTRIB, 3, GL_FLOAT, 0, 0, 0);

	// unbind the VAO
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	//	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB); 
	//	glDisableVertexAttribArray(COLOR_ATTRIB);
	checkOpenGLError("ERROR: Could not create VAOs and VBOs.");
}

void destroyBufferObjects()
{
	glDisableVertexAttribArray(VERTEX_COORD_ATTRIB);
	glDisableVertexAttribArray(COLOR_ATTRIB);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	glDeleteBuffers(1, VboId);
	glDeleteVertexArrays(1, &VaoId);
	checkOpenGLError("ERROR: Could not destroy VAOs and VBOs.");
}

void drawPoints()
{
	FrameCount++;
	glClear(GL_COLOR_BUFFER_BIT);

	glBindVertexArray(VaoId);
	glUseProgram(ProgramId);

	glBindBuffer(GL_ARRAY_BUFFER, VboId[0]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_vertices, vertices);
	glBindBuffer(GL_ARRAY_BUFFER, VboId[1]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, size_colors, colors);

	glUniformMatrix4fv(UniformId, 1, GL_FALSE, m);
	glDrawArrays(GL_POINTS, 0, RES_X * RES_Y);
	glFinish();

	glUseProgram(0);
	glBindVertexArray(0);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	checkOpenGLError("ERROR: Could not draw scene.");
}

ILuint saveImgFile(const char* filename) {
	ILuint ImageId;

	ilEnable(IL_FILE_OVERWRITE);
	ilGenImages(1, &ImageId);
	ilBindImage(ImageId);

	ilTexImage(RES_X, RES_Y, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, img_Data /*Texture*/);
	ilSaveImage(filename);

	ilDisable(IL_FILE_OVERWRITE);
	ilDeleteImages(1, &ImageId);
	if (ilGetError() != IL_NO_ERROR)return ilGetError();

	return IL_NO_ERROR;
}

/////////////////////////////////////////////////////////////////////// CALLBACKS

void timer(int value)
{
	std::ostringstream oss;
	oss << CAPTION << ": " << FrameCount << " FPS @ (" << RES_X << "x" << RES_Y << ")";
	std::string s = oss.str();
	glutSetWindow(WindowHandle);
	glutSetWindowTitle(s.c_str());
	FrameCount = 0;
	glutTimerFunc(1000, timer, 0);
}


// Callback function for glutCloseFunc
void cleanup()
{
	destroyShaderProgram();
	destroyBufferObjects();
}

void ortho(float left, float right, float bottom, float top,
	float nearp, float farp)
{
	m[0 * 4 + 0] = 2 / (right - left);
	m[0 * 4 + 1] = 0.0;
	m[0 * 4 + 2] = 0.0;
	m[0 * 4 + 3] = 0.0;
	m[1 * 4 + 0] = 0.0;
	m[1 * 4 + 1] = 2 / (top - bottom);
	m[1 * 4 + 2] = 0.0;
	m[1 * 4 + 3] = 0.0;
	m[2 * 4 + 0] = 0.0;
	m[2 * 4 + 1] = 0.0;
	m[2 * 4 + 2] = -2 / (farp - nearp);
	m[2 * 4 + 3] = 0.0;
	m[3 * 4 + 0] = -(right + left) / (right - left);
	m[3 * 4 + 1] = -(top + bottom) / (top - bottom);
	m[3 * 4 + 2] = -(farp + nearp) / (farp - nearp);
	m[3 * 4 + 3] = 1.0;
}

void reshape(int w, int h)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glViewport(0, 0, w, h);
	ortho(0, (float)RES_X, 0, (float)RES_Y, -1.0, 1.0);
}

void processKeys(unsigned char key, int xx, int yy)
{
	switch (key) {

	case 27:
		glutLeaveMainLoop();
		break;

	case 'r':
		camX = Eye.x;
		camY = Eye.y;
		camZ = Eye.z;
		r = Eye.length();
		beta = asinf(camY / r) * 180.0f / 3.14f;
		alpha = atanf(camX / camZ) * 180.0f / 3.14f;
		break;

	case 'c':
		printf("Camera Spherical Coordinates (%f, %f, %f)\n", r, beta, alpha);
		printf("Camera Cartesian Coordinates (%f, %f, %f)\n", camX, camY, camZ);
		break;
	}
}


// ------------------------------------------------------------
//
// Mouse Events
//

void processMouseButtons(int button, int state, int xx, int yy)
{
	// start tracking the mouse
	if (state == GLUT_DOWN) {
		startX = xx;
		startY = yy;
		if (button == GLUT_LEFT_BUTTON)
			tracking = 1;
		else if (button == GLUT_RIGHT_BUTTON)
			tracking = 2;
	}

	//stop tracking the mouse
	else if (state == GLUT_UP) {
		if (tracking == 1) {
			alpha -= (xx - startX);
			beta += (yy - startY);
		}
		else if (tracking == 2) {
			r += (yy - startY) * 0.01f;
			if (r < 0.1f)
				r = 0.1f;
		}
		tracking = 0;
	}
}

// Track mouse motion while buttons are pressed

void processMouseMotion(int xx, int yy)
{

	int deltaX, deltaY;
	float alphaAux, betaAux;
	float rAux;

	deltaX = -xx + startX;
	deltaY = yy - startY;

	// left mouse button: move camera
	if (tracking == 1) {


		alphaAux = alpha + deltaX;
		betaAux = beta + deltaY;

		if (betaAux > 85.0f)
			betaAux = 85.0f;
		else if (betaAux < -85.0f)
			betaAux = -85.0f;
		rAux = r;
	}
	// right mouse button: zoom
	else if (tracking == 2) {

		alphaAux = alpha;
		betaAux = beta;
		rAux = r + (deltaY * 0.01f);
		if (rAux < 0.1f)
			rAux = 0.1f;
	}

	camX = rAux * sin(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
	camZ = rAux * cos(alphaAux * 3.14f / 180.0f) * cos(betaAux * 3.14f / 180.0f);
	camY = rAux * sin(betaAux * 3.14f / 180.0f);
}

void mouseWheel(int wheel, int direction, int x, int y) {

	r += direction * 0.1f;
	if (r < 0.1f)
		r = 0.1f;

	camX = r * sin(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	camZ = r * cos(alpha * 3.14f / 180.0f) * cos(beta * 3.14f / 180.0f);
	camY = r * sin(beta * 3.14f / 180.0f);
}


void setupGLEW() {
	glewExperimental = GL_TRUE;
	GLenum result = glewInit();
	if (result != GLEW_OK) {
		std::cerr << "ERROR glewInit: " << glewGetString(result) << std::endl;
		exit(EXIT_FAILURE);
	}
	GLenum err_code = glGetError();
	printf("Vendor: %s\n", glGetString(GL_VENDOR));
	printf("Renderer: %s\n", glGetString(GL_RENDERER));
	printf("Version: %s\n", glGetString(GL_VERSION));
	printf("GLSL: %s\n", glGetString(GL_SHADING_LANGUAGE_VERSION));
}

void setupGLUT(int argc, char* argv[])
{
	glutInit(&argc, argv);

	glutInitContextVersion(4, 3);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	glutInitContextProfile(GLUT_CORE_PROFILE);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	glutInitWindowPosition(100, 250);
	glutInitWindowSize(RES_X, RES_Y);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glDisable(GL_DEPTH_TEST);
	WindowHandle = glutCreateWindow(CAPTION);
	if (WindowHandle < 1) {
		std::cerr << "ERROR: Could not create a new rendering window." << std::endl;
		exit(EXIT_FAILURE);
	}
}


/////////////////////////////////////////////////////YOUR CODE HERE///////////////////////////////////////////////////////////////////////////////////////


float fresnel(Vector I, Vector N, const float& ior, float& kr)
{
	float cosi = clamp(-1, 1, I * N);
	float etai = 1, etat = ior;
	if (cosi > 0) { std::swap(etai, etat); }
	// Compute sini using Snell's law
	float sint = etai / etat * sqrtf(std::max(0.f, 1 - cosi * cosi));
	// Total internal reflection
	if (sint >= 1) {
		kr = 1;
	}
	else {
		float cost = sqrtf(std::max(0.f, 1 - sint * sint));
		cosi = fabsf(cosi);
		float Rs = ((etat * cosi) - (etai * cost)) / ((etat * cosi) + (etai * cost));
		float Rp = ((etai * cosi) - (etat * cost)) / ((etai * cosi) + (etat * cost));
		kr = (Rs * Rs + Rp * Rp) / 2;
	}
	return kr;
	// As a consequence of the conservation of energy, the transmittance is given by:
	// kt = 1 - kr;
}



Color rayTracing(Ray ray, int depth, float ior_1)  //index of refraction of medium 1 where the ray is travelling
{

	int num_objects = scene->getNumObjects();

	bool one_interception = false, interception = false;

	Object* final_obj = NULL;

	Color color, rColor, rfColor;

	float check_dist;
	float min_dist = FLT_MAX;

	Vector hit_point, hit;

	// NO ACCELERATION
	if (Accel_Struct == NONE) {
		for (int i = 0; i < num_objects; i++) {

			Object* obj = scene->getObject(i);
			interception = obj->intercepts(ray, check_dist);

			if (interception) {
				one_interception = true;
				if (check_dist < min_dist) {
					final_obj = obj;
					min_dist = check_dist;
				}
			}
		}

		hit = ray.origin + ray.direction * min_dist;

	}
	else if (Accel_Struct == BVH_ACC) { // BVH ACCELERATION
		if (bvh_ptr->Traverse(ray, &final_obj, hit)) {
			one_interception = true;

		}

	}
	else if (Accel_Struct == GRID_ACC) {
		//num_objects = grid_ptr->getNumObjects();
		if (grid_ptr->Traverse(ray, &final_obj, hit)) {
			one_interception = true;
		}
	}

	if (!one_interception) {
		return scene->GetBackgroundColor();
	}

	hit_point = hit + final_obj->getNormal(hit) * EPSILON;

	int num_lights = scene->getNumLights();

	color = Color(0.0f, 0.0f, 0.0f);

	for (int j = 0; j < num_lights; j++) {
		Light* light = scene->getLight(j);

		Vector Lpos;

		if (!ANTI_ALIASING) {
			Lpos = light->position;
		}
		else {
			// Pick random position within area of light
			Lpos = Vector(
				light->position.x + 0.5f * ((off_x + rand_float()) / 4),
				light->position.y + 0.5f * ((off_y + rand_float()) / 4),
				light->position.z);

		}

		Vector L = (Lpos - hit).normalize();
		Vector h = (L - ray.direction).normalize();

		if (L * final_obj->getNormal(hit_point) > 0) {

			one_interception = false;
			float d;

			float ShadowRaysX = 2.0;
			float ShadowRaysY = 2.0;
			Vector ShadowCast = Lpos;

			if (!hardS) { 	//SOFT SHADOWS: assumes area light source; light as many point lights

				for (int s = 0; s <= ShadowRaysX * ShadowRaysY; s++)
				{
					float x_off = s % int(ShadowRaysX);
					float z_off = round((s / ShadowRaysX));

					ShadowCast.x += (x_off * 0.03);
					ShadowCast.z += (z_off * 0.03);


					Vector SVector = (ShadowCast - hit_point);
					SVector.normalize();

					Ray rShadow = Ray(hit_point, SVector);

					if (Accel_Struct == BVH_ACC) { // BVH
						if (bvh_ptr->Traverse(rShadow)) {
							one_interception = true;
						}
					}
					else if (Accel_Struct == GRID_ACC) { // GRID
						if (grid_ptr->Traverse(rShadow)) {
							one_interception = true;
						}
					}
					else { // NO ACCELERATION
						for (size_t k = 0; k < num_objects; k++)
						{
							Object* obj = scene->getObject(k);
							if (obj->intercepts(rShadow, d)) {
								one_interception = true;
								break;
							}
						}
					}

					if (!one_interception) {
						Color LC = light->color * (1 / (ShadowRaysX * ShadowRaysY));

						// PHONG ILLUMINATION
						color += LC * final_obj->GetMaterial()->GetDiffColor() * final_obj->GetMaterial()->GetDiffuse() * max(final_obj->getNormal(hit_point) * L, 0.0f) + LC * final_obj->GetMaterial()->GetSpecColor() * final_obj->GetMaterial()->GetSpecular() * pow(max(h * final_obj->getNormal(hit_point), 0.0f), final_obj->GetMaterial()->GetShine());
					}

				}
			}
			else { // HARD SHADOWS: assume an infinitely small (point) light source
				Ray rShadow = Ray(hit_point, L);

				if (Accel_Struct == BVH_ACC) { //BVH
					if (bvh_ptr->Traverse(rShadow)) {
						one_interception = true;
					}
				}
				else if (Accel_Struct == GRID_ACC) {
					if (grid_ptr->Traverse(rShadow)) {
						one_interception = true;
					}
				}
				else { // NO ACCELERATION
					for (size_t k = 0; k < num_objects; k++)
					{
						Object* obj = scene->getObject(k);
						if (obj->intercepts(rShadow, d)) {
							one_interception = true;
							break;
						}
					}
				}

				if (!one_interception) {
					color += final_obj->GetMaterial()->GetDiffColor() * final_obj->GetMaterial()->GetDiffuse() * max(final_obj->getNormal(hit_point) * L, 0.0f) + final_obj->GetMaterial()->GetSpecColor() * final_obj->GetMaterial()->GetSpecular() * pow(max(h * final_obj->getNormal(hit_point), 0.0f), final_obj->GetMaterial()->GetShine());
				}
			}
		}
	}

	if (depth > MAX_DEPTH) {
		return color.clamp();
	}

	Vector V = ray.direction * -1; //VIEW
	Vector n = final_obj->getNormal(hit_point);
	Vector bias = n * EPSILON;
	bool Inside = (V * n) < 0;

	if (final_obj->GetMaterial()->GetReflection() > 0) {

		if (Inside) {
			hit_point = hit - bias;
			n = n * -1;
		}
		else {
			hit_point = hit + bias;
		}
		Vector reflected = n * (2 * (V * n)) - V;
		if (FUZZY_REFLECTIONS) {
			// generate a random vector within a unit sphere
			Vector p;
			p = Vector(rand_float(), rand_float(), rand_float()) * 2 - Vector(1.0, 1.0, 1.0);
			while (p * p >= 1.0) {
				p = Vector(rand_float(), rand_float(), rand_float()) * 2 - Vector(1.0, 1.0, 1.0);
			}

			Ray rRay = Ray(hit_point, (reflected + (p * 0.1f)).normalize());
			rColor = rayTracing(rRay, depth + 1, ior_1);

		}
		else {
			Ray rRay = Ray(hit_point, reflected);
			rColor = rayTracing(rRay, depth + 1, ior_1);
		}
	}

	// REFRACTION //
	float Kr = 1; // Reflection mix value 
	//Vector normal = n;
	Vector Vt = n * (V * n) - V;

	if (final_obj->GetMaterial()->GetTransmittance() > 0) {

		float index;
		if (!Inside) {  // SNELL LAW
			index = ior_1 / final_obj->GetMaterial()->GetRefrIndex();
		}
		else {
			index = ior_1;
		}

		float sin = (index)*Vt.length();
		float temp = 1 - pow(sin, 2);

		if (temp >= 0) {
			float cos = sqrt(temp);

			// refracted direction 
			Vector tDir = (Vt.normalize() * sin + n * (-cos)).normalize();

			// Fresnel Equations 
			Kr = fresnel(V, n, final_obj->GetMaterial()->GetRefrIndex(), Kr);

			Vector refracted = hit + tDir * EPSILON;//avoid self-intersection

			Ray tRay = Ray(refracted, tDir);

			//if inside new index is set to 1
			float newIoR = !Inside ? final_obj->GetMaterial()->GetRefrIndex() : 1; // New Indice of Refraction

			// compute refracted color using recursion (tColor = rayTracing(refracted ray direction, depth+1)
			rfColor = rayTracing(tRay, depth + 1, newIoR);
		}
	}
	else { // is opaque
		Kr = final_obj->GetMaterial()->GetSpecular();
	}

	// join colors
	color += rColor * Kr * final_obj->GetMaterial()->GetSpecColor() + rfColor * (1 - Kr);

	return color.clamp();

}

// Render function by primary ray casting from the eye towards the scene's objects
void renderScene()
{
	int index_pos = 0;
	int index_col = 0;
	unsigned int counter = 0;

	if (drawModeEnabled) {
		glClear(GL_COLOR_BUFFER_BIT);
		scene->GetCamera()->SetEye(Vector(camX, camY, camZ));  //Camera motion
	}


	for (int y = 0; y < RES_Y; y++)
	{
		for (int x = 0; x < RES_X; x++)
		{
			Color color, jitter;

			if (ANTI_ALIASING) {
				int n = 4;

				for (int p = 0; p < n; p++)
				{
					for (int q = 0; q < n; q++) {

						off_x = p;
						off_y = q;

						Vector pixel;  //viewport coordinates
						pixel.x = x + (p + EPSILON) / n;
						pixel.y = y + (q + EPSILON) / n;

						if (DEPTH_OF_FIELD) {
							Vector lens;
							float aperture = scene->GetCamera()->GetAperture();
							//float aperture = 0.8;

							Vector p;
							p = Vector(rand_float(), rand_float(), 0.0) * 2 - Vector(1.0, 1.0, 0.0);
							while (p * p >= 1.0) {
								p = Vector(rand_float(), rand_float(), 0.0) * 2 - Vector(1.0, 1.0, 0.0);
							}

							// Compute the sample point on the "thin lens"
							lens = p * aperture;

							Ray ray = scene->GetCamera()->PrimaryRay(lens, pixel);
							jitter += rayTracing(ray, 1, 1.0).clamp();

						}
						else {
							Ray ray = scene->GetCamera()->PrimaryRay(pixel);   //function from camera.h
							jitter += rayTracing(ray, 1, 1.0).clamp();

						}

					}
				}

				color = jitter * (1 / pow(n, 2));

			}
			else {
				Vector pixel;  //viewport coordinates
				pixel.x = x + 0.5f;
				pixel.y = y + 0.5f;

				if (DEPTH_OF_FIELD) {
					Vector lens;
					float aperture = scene->GetCamera()->GetAperture();

					Vector p;
					p = Vector(rand_float(), rand_float(), 0.0) * 2 - Vector(1.0, 1.0, 0.0);
					while (p * p >= 1.0) {
						p = Vector(rand_float(), rand_float(), 0.0) * 2 - Vector(1.0, 1.0, 0.0);
					}

					// Compute the sample point on the "thin lens"
					lens = p * aperture;

					Ray ray = scene->GetCamera()->PrimaryRay(lens, pixel);
					color = rayTracing(ray, 1, 1.0).clamp();
				}
				else {
					Ray ray = scene->GetCamera()->PrimaryRay(pixel);   //function from camera.h
					color = rayTracing(ray, 1, 1.0).clamp();
				}
			}

			img_Data[counter++] = u8fromfloat((float)color.r());
			img_Data[counter++] = u8fromfloat((float)color.g());
			img_Data[counter++] = u8fromfloat((float)color.b());

			if (drawModeEnabled) {
				vertices[index_pos++] = (float)x;
				vertices[index_pos++] = (float)y;
				colors[index_col++] = (float)color.r();

				colors[index_col++] = (float)color.g();

				colors[index_col++] = (float)color.b();
			}
		}

	}
	if (drawModeEnabled) {
		drawPoints();
		glutSwapBuffers();
	}
	else {
		printf("Terminou o desenho!\n");
		if (saveImgFile("RT_Output.png") != IL_NO_ERROR) {
			printf("Error saving Image file\n");
			exit(0);
		}
		printf("Image file created\n");
	}
}


///////////////////////////////////////////////////////////////////////  SETUP     ///////////////////////////////////////////////////////

void setupCallbacks()
{
	glutKeyboardFunc(processKeys);
	glutCloseFunc(cleanup);
	glutDisplayFunc(renderScene);
	glutReshapeFunc(reshape);
	glutMouseFunc(processMouseButtons);
	glutMotionFunc(processMouseMotion);
	glutMouseWheelFunc(mouseWheel);

	glutIdleFunc(renderScene);
	glutTimerFunc(0, timer, 0);
}
void init(int argc, char* argv[])
{
	// set the initial camera position on its spherical coordinates
	Eye = scene->GetCamera()->GetEye();
	camX = Eye.x;
	camY = Eye.y;
	camZ = Eye.z;
	r = Eye.length();
	beta = asinf(camY / r) * 180.0f / 3.14f;
	alpha = atanf(camX / camZ) * 180.0f / 3.14f;

	setupGLUT(argc, argv);
	setupGLEW();
	std::cerr << "CONTEXT: OpenGL v" << glGetString(GL_VERSION) << std::endl;
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	createShaderProgram();
	createBufferObjects();
	setupCallbacks();
}


void init_scene(void)
{
	char scenes_dir[70] = "P3D_Scenes/";
	char input_user[50];
	char scene_name[70];

	scene = new Scene();

	if (P3F_scene) {  //Loading a P3F scene

		while (true) {
			cout << "Input the Scene Name: ";
			cin >> input_user;
			strcpy_s(scene_name, sizeof(scene_name), scenes_dir);
			strcat_s(scene_name, sizeof(scene_name), input_user);

			ifstream file(scene_name, ios::in);
			if (file.fail()) {
				printf("\nError opening P3F file.\n");
			}
			else
				break;
		}

		scene->load_p3f(scene_name);
		printf("Scene loaded.\n\n");
	}
	else {
		printf("Creating a Random Scene.\n\n");
		scene->create_random_scene();
	}


	RES_X = scene->GetCamera()->GetResX();
	RES_Y = scene->GetCamera()->GetResY();
	printf("\nResolutionX = %d  ResolutionY= %d.\n", RES_X, RES_Y);

	// Pixel buffer to be used in the Save Image function
	img_Data = (uint8_t*)malloc(3 * RES_X * RES_Y * sizeof(uint8_t));
	if (img_Data == NULL) exit(1);

	//Accel_Struct = scene->GetAccelStruct();   //Type of acceleration data structure

	if (Accel_Struct == GRID_ACC) {
		grid_ptr = new Grid();
		vector<Object*> objs;
		int num_objects = scene->getNumObjects();

		for (int o = 0; o < num_objects; o++) {
			objs.push_back(scene->getObject(o));
		}
		grid_ptr->Build(objs);
		printf("Grid built.\n\n");
	}
	else if (Accel_Struct == BVH_ACC) {
		printf("enters here");
		vector<Object*> objs;
		int num_objects = scene->getNumObjects();
		bvh_ptr = new BVH();

		for (int o = 0; o < num_objects; o++) {
			objs.push_back(scene->getObject(o));
		}
		bvh_ptr->Build(objs);
		printf("BVH built.\n\n");
	}
	else
		printf("No acceleration data structure.\n\n");

	unsigned int spp = scene->GetSamplesPerPixel();
	if (spp == 0)
		printf("Whitted Ray-Tracing\n");
	else
		printf("Distribution Ray-Tracing\n");

}

int main(int argc, char* argv[])
{
	//Initialization of DevIL 
	if (ilGetInteger(IL_VERSION_NUM) < IL_VERSION)
	{
		printf("wrong DevIL version \n");
		exit(0);
	}
	ilInit();

	int
		ch;
	if (!drawModeEnabled) {

		do {
			init_scene();

			auto timeStart = std::chrono::high_resolution_clock::now();
			renderScene();  //Just creating an image file
			auto timeEnd = std::chrono::high_resolution_clock::now();
			auto passedTime = std::chrono::duration<double, std::milli>(timeEnd - timeStart).count();
			printf("\nDone: %.2f (sec)\n", passedTime / 1000);
			if (!P3F_scene) break;
			cout << "\nPress 'y' to render another image or another key to terminate!\n";
			delete(scene);
			free(img_Data);
			ch = _getch();
		} while ((toupper(ch) == 'Y'));
	}

	else {   //Use OpenGL to draw image in the screen
		printf("OPENGL DRAWING MODE\n\n");
		init_scene();
		size_vertices = 2 * RES_X * RES_Y * sizeof(float);
		size_colors = 3 * RES_X * RES_Y * sizeof(float);
		vertices = (float*)malloc(size_vertices);
		if (vertices == NULL) exit(1);
		colors = (float*)malloc(size_colors);
		if (colors == NULL) exit(1);
		memset(colors, 0, size_colors);

		/* Setup GLUT and GLEW */
		init(argc, argv);
		glutMainLoop();
	}

	free(colors);
	free(vertices);
	printf("Program ended normally\n");
	exit(EXIT_SUCCESS);
}
///////////////////////////////////////////////////////////////////////