
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "SDL.h" 
#include "SDL2_gfxPrimitives.h" 
#include "SDL_mixer.h"
#include <vector>
#include <algorithm>
#include "SDL_image.h" 
#include "SDL_ttf.h"
#include <iostream>
#include "SDLgeneral.h"
#include "WorldMap.h"


using namespace std;

int initSDL();
void closeSDL();
void mouseHandleEvent(SDL_Event* e);

struct vec3d
{
	double x, y, z, w = 1;
};
struct triangle
{
	vec3d p[3];
	double shadow;
};
struct mat4x4
{
	double m[4][4] = { 0 };
};
struct mesh
{
	vector<triangle> tris;
};
struct boundaryData {
	float d;
	float dot;
};
struct enemy {
	float x;
	float y;
	int img;
};
struct {
	int x, y;
	int phase;
} port[2];

struct antiInfo {
	int c;
	int r;
	float dis;
	int move_count;
	int wall_count;
};

//jesus perspective variables
triangle triRotatedZ, triRotatedZX;
double cx = -0.5, cy = -0.5, cz = -0.5;
//store projected triangle coordinate
const int n = 12;//2*6 faces
Sint16 vx[n][3];
Sint16 vy[n][3];
double fThetax = 0.0, fThetay = 0.0, fThetaz = 0.0f, fThetaA = 0.0;;
vector<double> shadow[n];

//first perspective variables
triangle triProjected, triTranslated, triViewed;
mat4x4 matProj;
mesh meshCube, land, originalCube;
vec3d vCamera = { 0, 0, -20 };
//vec3d firstview = { 1.0,-1.0 ,-1.0 };
vec3d vLookDir;
vec3d vRight;
vector<triangle> listTriangles;
double fYaw;

vec3d normal, line1, line2;
double Distance;

bool LoadFromObjectFile(const char* sFilename)
{
	FILE* fp;
	errno_t err;
	err = fopen_s(&fp, sFilename, "r");
	if (err)
		return false;
	// Local cache of verts
	vector<vec3d> verts;
	char junk;
	while (!feof(fp) && fscanf_s(fp, " %c", &junk, 1))
	{
		//printf("junk:%c\n", junk);
		if (junk == 'v')
		{
			vec3d v;
			fscanf_s(fp, "%lf%lf%lf", &v.x, &v.y, &v.z);
			verts.push_back(v);
			//printf("%lf %lf %lf %lf\n", v.x, v.y, v.z, v.w);
		}
		else if (junk == 'f')
		{
			int f[3];
			fscanf_s(fp, "%d%d%d", &f[0], &f[1], &f[2]);
			//printf("%d %d %d\n", f[0], f[1], f[2]);
			meshCube.tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
			originalCube.tris.push_back({ verts[f[0] - 1], verts[f[1] - 1], verts[f[2] - 1] });
		}
		else {
			char line[128];
			fgets(line, 128, fp);
			//printf("%s\n", line);
		}
	}
	return true;
}

vec3d Vector_Add(vec3d& v1, vec3d& v2)
{
	return { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
}
vec3d Vector_Sub(vec3d& v1, vec3d& v2)
{
	return { v1.x - v2.x, v1.y - v2.y, v1.z - v2.z };
}
vec3d Vector_Mul(vec3d& v1, double k)
{
	return { v1.x * k, v1.y * k, v1.z * k };
}
vec3d Vector_Div(vec3d& v1, double k)
{
	return { v1.x / k, v1.y / k, v1.z / k };
}
double Vector_DotProduct(vec3d& v1, vec3d& v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}
double Vector_Length(vec3d& v)
{
	return sqrtf(Vector_DotProduct(v, v));
}
vec3d Vector_Normalise(vec3d& v)
{
	double l = Vector_Length(v);
	return { v.x / l, v.y / l, v.z / l };
}
vec3d Vector_CrossProduct(vec3d& v1, vec3d& v2)
{
	vec3d v;
	v.x = v1.y * v2.z - v1.z * v2.y;
	v.y = v1.z * v2.x - v1.x * v2.z;
	v.z = v1.x * v2.y - v1.y * v2.x;
	return v;
}
vec3d Vector_IntersectPlane(vec3d& plane_p, vec3d& plane_n, vec3d& lineStart, vec3d& lineEnd)
{
	plane_n = Vector_Normalise(plane_n);
	float plane_d = -Vector_DotProduct(plane_n, plane_p);//the distance of plane to the origin
	float ad = Vector_DotProduct(lineStart, plane_n);//distance of point a at vector n
	float bd = Vector_DotProduct(lineEnd, plane_n);//distance of point b at vecto n
	float t = (-plane_d - ad) / (bd - ad);
	vec3d lineStartToEnd = Vector_Sub(lineEnd, lineStart);
	vec3d lineToIntersect = Vector_Mul(lineStartToEnd, t);
	return Vector_Add(lineStart, lineToIntersect);
}

// Return signed shortest distance from point to plane, plane normal must be normalised
//https://stackoverflow.com/questions/12262019/c-operator

double dist(vec3d plane_n, vec3d plane_p, vec3d& p) {
	//vec3d n = Vector_Normalise(p);
	return (plane_n.x * p.x + plane_n.y * p.y + plane_n.z * p.z - Vector_DotProduct(plane_n, plane_p));
}

int Triangle_ClipAgainstPlane(vec3d plane_p, vec3d plane_n, triangle& in_tri, triangle& out_tri1, triangle& out_tri2)
{
	//clear triangle components to avoid error

	for (int i = 0; i < 3; i++) {
		out_tri1.p[i].x = 0;
		out_tri1.p[i].y = 0;
		out_tri1.p[i].z = 0;
		out_tri1.p[i].w = 0;
		out_tri2.p[i].x = 0;
		out_tri2.p[i].y = 0;
		out_tri2.p[i].z = 0;
		out_tri2.p[i].w = 0;
	}


	//return how many triangles crated in the end
	// Make sure plane normal is indeed normal
	plane_n = Vector_Normalise(plane_n);

	// Create two temporary storage arrays to classify points either side of plane
	// If distance sign is positive, point lies on "inside" of plane
	vec3d* inside_points[3];  int nInsidePointCount = 0;
	vec3d* outside_points[3]; int nOutsidePointCount = 0;

	// Get signed distance of each point in triangle to plane
	double d0 = dist(plane_n, plane_p, in_tri.p[0]);
	double d1 = dist(plane_n, plane_p, in_tri.p[1]);
	double d2 = dist(plane_n, plane_p, in_tri.p[2]);

	//Why &in_tri.p[0] because inside_points is pointer
	//store corresponding point of the triangle
	if (d0 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[0]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[0]; }

	if (d1 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[1]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[1]; }

	if (d2 >= 0) { inside_points[nInsidePointCount++] = &in_tri.p[2]; }
	else { outside_points[nOutsidePointCount++] = &in_tri.p[2]; }

	// Now classify triangle points, and break the input triangle into 
	// smaller output triangles if required. There are four possible
	// outcomes...
	if (nInsidePointCount == 0)
	{
		// All points lie on the outside of plane, so clip whole triangle
		// It ceases to exist
		return 0; // No returned triangles are valid
	}

	if (nInsidePointCount == 3)
	{
		// All points lie on the inside of plane, so do nothing
		// and allow the triangle to simply pass through
		out_tri1 = in_tri;

		return 1; // Just the one returned original triangle is valid
	}

	if (nInsidePointCount == 1 && nOutsidePointCount == 2)
	{
		// Triangle should be clipped. As two points lie outside
		// the plane, the triangle simply becomes a smaller triangle

		// Copy appearance info to new triangle
		out_tri1.shadow = in_tri.shadow;

		// The inside point is valid, so keep that...
		out_tri1.p[0] = *inside_points[0];

		// but the two new points are at the locations where the 
		// original sides of the triangle (lines) intersect with the plane
		out_tri1.p[1] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);
		out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[1]);

		return 1; // Return the newly formed single triangle
	}

	if (nInsidePointCount == 2 && nOutsidePointCount == 1)
	{
		// Triangle should be clipped. As two points lie inside the plane,
		// the clipped triangle becomes a "quad". Fortunately, we can
		// represent a quad with two new triangles

		// Copy appearance info to new triangles
		out_tri1.shadow = in_tri.shadow;
		out_tri2.shadow = in_tri.shadow;

		// The first triangle consists of the two inside points and a new
		// point determined by the location where one side of the triangle
		// intersects with the plane
		out_tri1.p[0] = *inside_points[0];
		out_tri1.p[1] = *inside_points[1];
		out_tri1.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[0], *outside_points[0]);

		// The second triangle is composed of one of he inside points, a
		// new point determined by the intersection of the other side of the 
		// triangle and the plane, and the newly created point above
		out_tri2.p[0] = *inside_points[1];
		out_tri2.p[1] = out_tri1.p[2];
		out_tri2.p[2] = Vector_IntersectPlane(plane_p, plane_n, *inside_points[1], *outside_points[0]);

		return 2; // Return two newly formed triangles which form a quad
	}
}

vec3d MultiplyMatrixVector(mat4x4& m, vec3d& i)
{
	vec3d v;
	v.x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
	v.y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
	v.z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
	v.w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];
	return v;
}
mat4x4 Matrix_MakeIdentity()
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0;
	matrix.m[1][1] = 1.0;
	matrix.m[2][2] = 1.0;
	matrix.m[3][3] = 1.0;
	return matrix;
}
mat4x4 Matrix_MakeRotationX(double fAngleRad)
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0;
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[1][2] = sinf(fAngleRad);
	matrix.m[2][1] = -sinf(fAngleRad);
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0;
	return matrix;
}
mat4x4 Matrix_MakeRotationY(double fAngleRad)
{
	mat4x4 matrix;
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][2] = sinf(fAngleRad);
	matrix.m[2][0] = -sinf(fAngleRad);
	matrix.m[1][1] = 1.0;
	matrix.m[2][2] = cosf(fAngleRad);
	matrix.m[3][3] = 1.0;
	return matrix;
}
mat4x4 Matrix_MakeRotationZ(double fAngleRad)
{
	mat4x4 matrix;
	matrix.m[0][0] = cosf(fAngleRad);
	matrix.m[0][1] = sinf(fAngleRad);
	matrix.m[1][0] = -sinf(fAngleRad);
	matrix.m[1][1] = cosf(fAngleRad);
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_MakeRotationA(double fAngleRad, vec3d a)
{
	mat4x4 matrix;

	a = { a.x - WIDTH / 2, (a.y - HEIGHT / 2), 0, 0 };
	double d = sqrtf(powf(a.x, 2.0) + powf(a.y, 2.0));
	//if (d < 100) return Matrix_MakeIdentity();
	vec3d axis = { a.y, -(a.x), 0, 0 };
	axis = Vector_Normalise(axis);
	//printf("asix = {%lf, %lf, %lf}\n", axis.x, axis.y, axis.z);
	//system("pause");
	double u = axis.x, v = axis.y, w = axis.z;
	matrix.m[0][0] = u * u + (v * v + w * w) * cosf(fAngleRad);
	matrix.m[0][1] = u * v * (1 - cosf(fAngleRad)) + w * sinf(fAngleRad);
	matrix.m[0][2] = u * w * (1 - cosf(fAngleRad)) - v * sinf(fAngleRad);
	matrix.m[0][3] = 0;
	matrix.m[1][0] = u * v * (1 - cosf(fAngleRad)) - w * sinf(fAngleRad);
	matrix.m[1][1] = v * v + (u * u + w * w) * cosf(fAngleRad);
	matrix.m[1][2] = v * w * (1 - cosf(fAngleRad)) + u * sinf(fAngleRad);
	matrix.m[1][3] = 0;
	matrix.m[2][0] = u * w * (1 - cosf(fAngleRad)) + v * sinf(fAngleRad);
	matrix.m[2][1] = v * w * (1 - cosf(fAngleRad)) - u * sinf(fAngleRad);
	matrix.m[2][2] = w * w + (u * u + v * v) * cosf(fAngleRad);
	matrix.m[2][3] = 0;
	/*matrix.m[3][0] = (a.x * (v * v + w * w) - u * (a.y * v + a.z * w)) * (1 - cosf(fAngleRad)) + (a.y * w - a.z * v) * sinf(fAngleRad);
	matrix.m[3][1] = (a.y * (u * u + w * w) - v * (a.x * u + a.z * w)) * (1 - cosf(fAngleRad)) + (a.z * u - a.x * w) * sinf(fAngleRad);
	matrix.m[3][2] = (a.z * (v * v + u * u) - w * (a.x * u + a.y * v)) * (1 - cosf(fAngleRad)) + (a.x * v - a.y * u) * sinf(fAngleRad);
	matrix.m[3][3] = 1.0f;*/
	matrix.m[3][0] = 0;
	matrix.m[3][1] = 0;
	matrix.m[3][2] = 0;
	matrix.m[3][3] = 1.0f;
	return matrix;
}
mat4x4 Matrix_MakeTranslation(double x, double y, double z)
{
	mat4x4 matrix;
	matrix.m[0][0] = 1.0f;
	matrix.m[1][1] = 1.0f;
	matrix.m[2][2] = 1.0f;
	matrix.m[3][3] = 1.0f;
	matrix.m[3][0] = x;
	matrix.m[3][1] = y;
	matrix.m[3][2] = z;
	return matrix;
}
mat4x4 Matrix_MakeProjection(double fFovDegrees, double fAspectRatio, double fNear, double fFar)
{
	float fFovRad = 1.0 / tanf(fFovDegrees * 0.5 / 180.0 * M_PI);
	mat4x4 matrix;
	matrix.m[0][0] = fAspectRatio * fFovRad;
	matrix.m[1][1] = fFovRad;
	matrix.m[2][2] = fFar / (fFar - fNear);
	matrix.m[3][2] = (-fFar * fNear) / (fFar - fNear);
	matrix.m[2][3] = 1.0;
	matrix.m[3][3] = 0.0;
	return matrix;
}
mat4x4 Matrix_MultiplyMatrix(mat4x4& m1, mat4x4& m2)
{
	mat4x4 matrix;
	for (int c = 0; c < 4; c++)
		for (int r = 0; r < 4; r++)
			matrix.m[r][c] = m1.m[r][0] * m2.m[0][c] + m1.m[r][1] * m2.m[1][c] + m1.m[r][2] * m2.m[2][c] + m1.m[r][3] * m2.m[3][c];
	return matrix;
}
mat4x4 Matrix_QuickInverse(mat4x4& m) // Only for Rotation/Translation Matrices
{
	mat4x4 matrix;
	matrix.m[0][0] = m.m[0][0]; matrix.m[0][1] = m.m[1][0]; matrix.m[0][2] = m.m[2][0]; matrix.m[0][3] = 0.0;
	matrix.m[1][0] = m.m[0][1]; matrix.m[1][1] = m.m[1][1]; matrix.m[1][2] = m.m[2][1]; matrix.m[1][3] = 0.0;
	matrix.m[2][0] = m.m[0][2]; matrix.m[2][1] = m.m[1][2]; matrix.m[2][2] = m.m[2][2]; matrix.m[2][3] = 0.0;
	matrix.m[3][0] = -(m.m[3][0] * matrix.m[0][0] + m.m[3][1] * matrix.m[1][0] + m.m[3][2] * matrix.m[2][0]);
	matrix.m[3][1] = -(m.m[3][0] * matrix.m[0][1] + m.m[3][1] * matrix.m[1][1] + m.m[3][2] * matrix.m[2][1]);
	matrix.m[3][2] = -(m.m[3][0] * matrix.m[0][2] + m.m[3][1] * matrix.m[1][2] + m.m[3][2] * matrix.m[2][2]);
	matrix.m[3][3] = 1.0;
	return matrix;
}
mat4x4 Matrix_PointAt(vec3d& pos, vec3d& target, vec3d& up)
{
	// Calculate new forward direction
	vec3d newForward = Vector_Sub(target, pos);
	newForward = Vector_Normalise(newForward);

	// Calculate new Up direction
	vec3d a = Vector_Mul(newForward, Vector_DotProduct(up, newForward));
	vec3d newUp = Vector_Sub(up, a);
	newUp = Vector_Normalise(newUp);

	// New Right direction is easy, its just cross product
	vec3d newRight = Vector_CrossProduct(newUp, newForward);

	// Construct Dimensioning and Translation Matrix	
	mat4x4 matrix;
	matrix.m[0][0] = newRight.x;	matrix.m[0][1] = newRight.y;	matrix.m[0][2] = newRight.z;	matrix.m[0][3] = 0.0f;
	matrix.m[1][0] = newUp.x;		matrix.m[1][1] = newUp.y;		matrix.m[1][2] = newUp.z;		matrix.m[1][3] = 0.0f;
	matrix.m[2][0] = newForward.x;	matrix.m[2][1] = newForward.y;	matrix.m[2][2] = newForward.z;	matrix.m[2][3] = 0.0f;
	matrix.m[3][0] = pos.x;			matrix.m[3][1] = pos.y;			matrix.m[3][2] = pos.z;			matrix.m[3][3] = 1.0f;
	return matrix;

}

void CubeGenerate() {
	// Load object file
	if (c_simple.click) {
		if (LoadFromObjectFile("Simple(3D).obj") == false) {
			printf("cube loading fail!\n");
			exit(1);
		}
	}
	else {
		if (LoadFromObjectFile("Classic(3D).obj") == false) {
			printf("cube loading fail!\n");
			exit(1);
		}
	}

	matProj = Matrix_MakeProjection(fFov, (double)HEIGHT / WIDTH, 0.1, 1000.0);
}

vector<triangle> queue;

void drawCube()
{

	matProj = Matrix_MakeProjection(fFov, (double)HEIGHT / WIDTH, 0.1, 1000.0);
	vector<triangle> vecTrianglesToRaster;
	// Set up rotation matrices
	// Set up "World Tranmsform" though not updating theta 
	// makes this a bit redundant

	int x, y;
	SDL_GetMouseState(&x, &y);
	//printf("x = %d, y = %d\n", x, y);

	mat4x4 matRotZ, matRotY, matRotX, matRotA, matRot;

	// Rotation Arbitrary axis
	matRotA = Matrix_MakeRotationA(fThetaA * 0.5, { (double)x, (double)y, 0.0, 0.0 });
	// Rotation Z
	matRotZ = Matrix_MakeRotationZ(fThetaz);
	// Rotation Y
	matRotY = Matrix_MakeRotationY(fThetay);
	// Rotation X
	matRotX = Matrix_MakeRotationX(fThetax);

	//printf("%lf %lf\n", fThetax, fThetay);
	matRot = Matrix_MakeIdentity();
	matRot = Matrix_MultiplyMatrix(matRotX, matRotY);
	matRot = Matrix_MultiplyMatrix(matRot, matRotZ);
	matRot = Matrix_MultiplyMatrix(matRot, matRotA);


	mat4x4 matTrans;
	matTrans = Matrix_MakeTranslation(0.0, 0.0, 5.0);

	mat4x4 matWorld;

	// Form World Matrix
	matWorld = Matrix_MakeIdentity();

	// Transform by rotation
	//matWorld = matRotA; 
	//matWorld = Matrix_MultiplyMatrix(matRotX, matRotY);

	// Transform by translation
	//matWorld = Matrix_MultiplyMatrix(matWorld, matTrans); 

	// Create "Point At" Matrix for camera
	vec3d vUp = { 0,1,0 };
	vec3d vTarget = { 0,0,1 };
	mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
	vLookDir = MultiplyMatrixVector(matCameraRot, vTarget);
	vRight = Vector_CrossProduct(vLookDir, vUp);
	vRight = Vector_Normalise(vRight);
	vTarget = Vector_Add(vCamera, vLookDir);
	mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);

	// Make view matrix from camera
	mat4x4 matView = Matrix_QuickInverse(matCamera);

	//listTriangles.clear();

	int c = 0;
	for (int i = 0; i < meshCube.tris.size(); i++) {

		meshCube.tris[i].p[0] = MultiplyMatrixVector(matRot, meshCube.tris[i].p[0]);
		meshCube.tris[i].p[1] = MultiplyMatrixVector(matRot, meshCube.tris[i].p[1]);
		meshCube.tris[i].p[2] = MultiplyMatrixVector(matRot, meshCube.tris[i].p[2]);


		triTranslated.p[0] = MultiplyMatrixVector(matTrans, meshCube.tris[i].p[0]);
		triTranslated.p[1] = MultiplyMatrixVector(matTrans, meshCube.tris[i].p[1]);
		triTranslated.p[2] = MultiplyMatrixVector(matTrans, meshCube.tris[i].p[2]);


		// Use Cross-Product to get surface normal
		vec3d normal, line1, line2;

		//get two sides of the triangle
		line1 = Vector_Sub(triTranslated.p[1], triTranslated.p[0]);
		line2 = Vector_Sub(triTranslated.p[2], triTranslated.p[0]);

		//cross-product to find the normal vector
		normal = Vector_CrossProduct(line1, line2);
		normal = Vector_Normalise(normal);

		// Get Ray from triangle to camera
		vec3d vCameraRay = Vector_Sub(triTranslated.p[0], vCamera);

		if (Vector_DotProduct(normal, vCameraRay) < 0.0)
		{
			//direction of light: object to eye
			vec3d light_direction = { 0.0f, 0.0f, -1.0f };
			//vec3d light_direction = Vector_Add(vCamera, vLookDir);
			light_direction = Vector_Normalise(light_direction);
			light_direction = Vector_Mul(light_direction, -1);

			// Convert World Space --> View Space
			triViewed.p[0] = MultiplyMatrixVector(matView, triTranslated.p[0]);
			triViewed.p[1] = MultiplyMatrixVector(matView, triTranslated.p[1]);
			triViewed.p[2] = MultiplyMatrixVector(matView, triTranslated.p[2]);

			// Clip Viewed Triangle against near plane, this could form two additional
			//additional triangles. 
			int nClippedTriangles = 0;//only 1 or 2 or 0  
			triangle clipped[2] = { 0 };
			nClippedTriangles = Triangle_ClipAgainstPlane({ 0.0, 0.0, 0.1 }, { 0.0, 0.0, 1.0 }, triViewed, clipped[0], clipped[1]);

			//printf("ININ\n\n");
			// We may end up with multiple triangles form the clip, so project as required
			//0 <= n <= 2
			for (int n = 0; n < nClippedTriangles; n++)
			{
				// Project triangles from 3D --> 2D
				triProjected.p[0] = MultiplyMatrixVector(matProj, clipped[n].p[0]);
				triProjected.p[1] = MultiplyMatrixVector(matProj, clipped[n].p[1]);
				triProjected.p[2] = MultiplyMatrixVector(matProj, clipped[n].p[2]);

				// Scale into view, we moved the normalising into cartesian space
				// out of the matrix.vector function from the previous videos, so
				// do this manually
				triProjected.p[0] = Vector_Div(triProjected.p[0], triProjected.p[0].w);
				triProjected.p[1] = Vector_Div(triProjected.p[1], triProjected.p[1].w);
				triProjected.p[2] = Vector_Div(triProjected.p[2], triProjected.p[2].w);

				// X/Y are inverted so put them back
				triProjected.p[0].x *= -1.0;
				triProjected.p[1].x *= -1.0;
				triProjected.p[2].x *= -1.0;
				triProjected.p[0].y *= -1.0;
				triProjected.p[1].y *= -1.0;
				triProjected.p[2].y *= -1.0;

				// Offset verts into visible normalised space
				vec3d vOffsetView = { 1,1,0 };
				triProjected.p[0] = Vector_Add(triProjected.p[0], vOffsetView);
				triProjected.p[1] = Vector_Add(triProjected.p[1], vOffsetView);
				triProjected.p[2] = Vector_Add(triProjected.p[2], vOffsetView);
				triProjected.p[0].x *= 0.5 * WIDTH;
				triProjected.p[0].y *= 0.5 * HEIGHT;
				triProjected.p[1].x *= 0.5 * WIDTH;
				triProjected.p[1].y *= 0.5 * HEIGHT;
				triProjected.p[2].x *= 0.5 * WIDTH;
				triProjected.p[2].y *= 0.5 * HEIGHT;

				triProjected.shadow = 0.8 * (normal.x * light_direction.x + normal.y * light_direction.y + normal.z * light_direction.z) * 0.5;

				// Store triangle for sorting
				vecTrianglesToRaster.push_back(triProjected);
			}
		}

	}


	sort(vecTrianglesToRaster.begin(), vecTrianglesToRaster.end(), [](triangle& t1, triangle& t2)
		{
			double z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0f;
			double z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0f;
			return z1 > z2;
		});



	for (int tri = 0; tri < vecTrianglesToRaster.size(); tri++) {
		queue.push_back(vecTrianglesToRaster[tri]);
		//filledTrigonRGBA(renderer, listTriangles[tri].p[0].x, listTriangles[tri].p[0].y, listTriangles[tri].p[1].x, listTriangles[tri].p[1].y, listTriangles[tri].p[2].x, listTriangles[tri].p[2].y, 255, 255, 255, 255);
		//trigonRGBA(renderer, listTriangles[tri].p[0].x, listTriangles[tri].p[0].y, listTriangles[tri].p[1].x, listTriangles[tri].p[1].y, listTriangles[tri].p[2].x, listTriangles[tri].p[2].y, 255, 255, 255, 255);
		//filledTrigonRGBA(renderer, listTriangles[i].p[0].x, listTriangles[i].p[0].y, listTriangles[i].p[1].x, listTriangles[i].p[1].y, listTriangles[i].p[2].x, listTriangles[i].p[2].y, cR, cG, cB, 255);
		//FilledPolygonRGBA(renderer, vx[i], vy[i], 3, 255, 255, 255, 255);
	}
	//printf("clip out\n");
	listTriangles.clear();
	vecTrianglesToRaster.clear();

}

//up down left right

int dr[4] = { 0,1,0,-1 };
int dc[4] = { 1,0,-1,0 };
vector<int> rQ, cQ;
int nodes_left_in_layer = 1;
int nodes_in_next_layer = 0;
bool visited[N][N] = { 0 };
bool reachEnd = false;
int pathR[N][N] = { 0 };
int pathC[N][N] = { 0 };

void explore_neighbors(int r, int c) {
	int rr, cc;
	for (int i = 0; i < 4; i++) {
		rr = r + dr[i];
		cc = c + dc[i];

		if (rr < 0 || rr >= N || cc < 0 || cc >= N) continue;
		if (visited[rr][cc] || map[rr][cc] == '#') continue;

		rQ.push_back(rr);
		cQ.push_back(cc);
		visited[rr][cc] = true;
		pathR[rr][cc] = r;
		pathC[rr][cc] = c;
		nodes_in_next_layer++;
	}
}

int solve(int& move_count, int sR, int sC, int eR, int eC) {
	int r, c;
	rQ.push_back(sR);
	cQ.push_back(sC);
	visited[sR][sC] = true;
	pathR[sR][sC] = -1;
	pathC[sR][sC] = -1;

	while (!(rQ.empty())) {
		r = rQ[0];
		c = cQ[0];
		rQ.erase(rQ.begin());
		cQ.erase(cQ.begin());

		if (r == eR && c == eC) {
			reachEnd = true;
			break;
		}
		explore_neighbors(r, c);
		nodes_left_in_layer--;
		if (nodes_left_in_layer == 0) {
			nodes_left_in_layer = nodes_in_next_layer;
			nodes_in_next_layer = 0;
			move_count++;
		}
	}
	if (reachEnd) {
		return move_count;
	}
	return -1;
}

float lastX = enemyX, lastY = enemyY;
int enemy_move_count;
int count_steps = 0;
int routeR[100] = { 0 };
int routeC[100] = { 0 };

//Uint32 enemy_move(Uint32 interval, void* param)
void enemy_move(void)
{

	reachEnd = false;
	nodes_left_in_layer = 1;
	nodes_in_next_layer = 0;
	for (int i = 0; i < mapSize; i++) {
		for (int j = 0; j < mapSize; j++) {
			visited[i][j] = false;
			pathR[i][j] = 0;
			pathC[i][j] = 0;
		}
	}
	int eR = (int)userY, eC = (int)userX;
	int sR = (int)enemyY, sC = (int)enemyX;

	int move_Count = 0;

	if (solve(move_Count, sR, sC, eR, eC) != -1) {

		//printf("total move: %d\n\n", move_Count);
		/*
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				printf("%d ", pathR[j][k]);
			}
			printf("\n");
		}
		printf("\n");
		for (int j = 0; j < N; j++) {
			for (int k = 0; k < N; k++) {
				printf("%d ", pathC[j][k]);
			}
			printf("\n");
		}

		printf("\n");
		*/

		//if (move_Count > 5) {
		count_steps = 0;
		enemy_move_count = move_Count;
		int ansR = eR, ansC = eC;
		int nextR = eR, nextC = eC;
		int i = 1;

		while (pathR[ansR][ansC] != -1) {
			ansR = nextR;
			ansC = nextC;
			routeR[move_Count - 1] = nextR;
			routeC[move_Count - 1] = nextC;
			move_Count--;
			//printf("%d: %d %d\n", i++, nextR, nextC);
			nextR = pathR[ansR][ansC];
			nextC = pathC[ansR][ansC];
		}
		//}

		/*
		printf("enemyX = %lf enemyY = %lf\n", enemyX, enemyY);
		for (int i = 0; i < enemy_move_count; i++) {
			printf("%d %d\n", routeR[i], routeC[i]);
		}
		printf("\n\n");
		*/
	}
	rQ.clear();
	cQ.clear();
	//if ((int)userX == (int)enemyX && (int)userY == (int)enemyY) return 0;
	//return interval;
}
int lastR, lastC;

void enemy_trace(void)
{
	//printf("\n\n");
	//map[(int)enemyY][(int)enemyX] = '0';
	if (count_steps < enemy_move_count) {
		//printf("countsteps: %d\n", count_steps);
		//printf("%lf %lf\n", enemyX, enemyY);

		/*
		enemyR = atan2f(enemyY - (float)routeR[count_steps], (float)routeC[count_steps] - enemyX)-M_PI/2.f;
		float repX, repY;
		repX = (float)routeC[count_steps] - enemyX;
		repY = (float)routeR[count_steps] - enemyY;
		if ((repX > 0 && repY > 0) || (repX > 0 && repY < 0)) enemyR = atan2f(repY, repX) + M_PI / 2.f;
		else if ((repX < 0 && repY > 0) || (repX < 0 || repY < 0)) enemyR = atan2f(repY, repX) - M_PI / 2.f;
		else enemyR = 0.f;
		EshiftX = -sinf(enemyR) * 4.f * elapsedTime;
		EshiftY = cosf(-enemyR) * 4.f * elapsedTime;
		enemyX += EshiftX;
		enemyY += EshiftY;
		*/
		if ((int)enemyX == routeC[count_steps] && enemyY < (float)routeR[count_steps]) {
			//printf("1\n");
			enemyY += Espeed * elapsedTime;
		}
		else if ((int)enemyX == routeC[count_steps] && (float)enemyY > routeR[count_steps]) {
			//printf("2\n");
			enemyY -= Espeed * elapsedTime;
		}
		else if (enemyX < (float)routeC[count_steps] && (int)enemyY == routeR[count_steps]) {
			//printf("3\n");
			enemyX += Espeed * elapsedTime;
		}
		else if (enemyX > (float)routeC[count_steps] && (int)enemyY == routeR[count_steps]) {
			//printf("4\n");
			enemyX -= Espeed * elapsedTime;
		}
		/*
		enemyX = (float)routeC[count_steps];
		enemyY = (float)routeR[count_steps];
		count_steps++;


		enemyX +=  4.f * elapsedTime;
		enemyY +=  4.f * elapsedTime;
		*/
		if (map[(int)enemyY][(int)enemyX] == '#' || enemyX < 1.f || enemyX >= (float)mapSize - 1.f || enemyY < 1.f || enemyY >= (float)mapSize - 1.f)
		{
			enemyX -= Espeed * EshiftX;
			enemyY -= Espeed * EshiftY;
		}


		//map[(int)enemyY][(int)enemyX] = 'E';
		if ((int)enemyX == routeC[count_steps] && (int)enemyY == routeR[count_steps]) {
			count_steps++;
			//printf("in countsteps: %d\n", count_steps);
		}
	}
	//if ((int)userX == (int)enemyX && (int)userY == (int)enemyY) return 0;
	//return interval;
}

float test_enemyR = enemyR;
//Uint32 enemy_check(Uint32 interval, void* param) {
void enemy_check(void) {
	do {
		hitUser = false;
		test_enemyR += 5.f;
		//ErayAngle = (enemyR - fov / 2);
		EeyeX = sinf(test_enemyR);
		EeyeY = cosf(test_enemyR);
		EhitWall = false;
		EdToWall = 0.f;
		dToUser = 0.f;
		while (!EhitWall && !hitUser && EdToWall < 10.f) {
			EdToWall += rayStepSize;
			dToUser += rayStepSize;
			EtestX = (int)(enemyX + EeyeX * EdToWall);
			EtestY = (int)(enemyY + EeyeY * EdToWall);
			if (EtestX < 2 || EtestX > mapSize - 2 || EtestY < 2 || EtestY > mapSize - 2 || map[EtestY][EtestX] == '#') {
				EhitWall = true;
			}
			else if ((EtestX == (int)userX) && (EtestY == (int)(userY))) {
				hitUser = true;
			}
		}

	} while (!hitUser && test_enemyR - enemyR < M_PI * 2.f);
	if (hitUser) {
		enemyR = fmod(enemyR, M_PI * 2.f);
		enemyR = test_enemyR - M_PI;
	}
	//else enemyR = fmod((float)rand(), M_PI * 2);
	//return interval;
}
void enemy_run(void) {
	//printf("in\n\n");
	enemyR = fmod(enemyR, M_PI * 2.f);
	test_enemyR = enemyR;
	do {
		//enemyR = fmod(enemyR, M_PI * 2.f);
		//printf("fmod R = %lf\n\n", fmod(enemyR, M_PI * 2));
		if (hitUser) test_enemyR = M_PI / 2.f + enemyR - fmod((float)rand(), M_PI);
		else test_enemyR += 1.f;
		//printf("test enemyR = %lf\n\n", test_enemyR);
		EeyeX = sinf(test_enemyR);
		EeyeY = cosf(test_enemyR);
		EhitWall = false;
		EdToWall = 0.f;
		while (!EhitWall && EdToWall < 0.5f) {
			EdToWall += rayStepSize;
			EtestX = (int)(enemyX + EeyeX * EdToWall);
			EtestY = (int)(enemyY + EeyeY * EdToWall);
			if (EtestX < 2 || EtestX > mapSize - 2 || EtestY < 2 || EtestY > mapSize - 2 || map[EtestY][EtestX] == '#') {
				EhitWall = true;
			}
		}
		if (hitUser && test_enemyR - enemyR < M_PI * 2.f) break;
	} while (EhitWall);
	//printf("out\n\n");
	enemyX += sinf(test_enemyR) * Espeed * elapsedTime;
	enemyY += cosf(test_enemyR) * Espeed * elapsedTime;
	if (map[(int)enemyY][(int)enemyX] == '#' || enemyX<1.f || enemyX>mapSize - 1.f || enemyY<1.f || enemyY>mapSize - 1.f) {
		enemyX -= sinf(test_enemyR) * Espeed * elapsedTime;
		enemyY -= cosf(test_enemyR) * Espeed * elapsedTime;
	}
	//printf("enemyX = %lf, enemyY = %lf\n\n", enemyX, enemyY);
}


int last_move_Count = mapSize * mapSize;
float last_dis = 0.0f;

bool GPS[6] = { 1, 0, 0, 0, 0, 0 };

int rotate(void) {
	//bool rotateFlag = false;
	float junk, size = (float)mapSize;

	if (c_classic.click) {
		//up
		if (userY < 0.5f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = m3;
				GPS[0] = 0;
				GPS[2] = 1; //change to m5
				userR = userR;
				userY = (int)(size - 1.f - 0.8f);
				return 13;
			}
			else if (GPS[1]) {
				map = m3;
				GPS[1] = 0;
				GPS[2] = 1; //change to m5
				userR = userR + (M_PI / 2.f);
				userY = (size - 1.f - userX);
				userX = (size - 1.f - 0.8f);
				return 23;
			}
			else if (GPS[4]) {
				map = m1;
				GPS[4] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userY = size - 1.f - 0.8f;
				return 51;
			}
			else if (GPS[3]) {
				map = m3;
				GPS[3] = 0;
				GPS[2] = 1; //change to m5
				userR = userR - M_PI / 2.f;
				userY = userX;
				userX = 0.8f;
				return 43;
			}
			else if (GPS[2]) {
				map = m6;
				GPS[2] = 0;
				GPS[5] = 1; //change to m6
				userR = userR + M_PI;
				userX = size - 1.f - userX;
				userY = 0.8f;
				return 36;
			}
			else if (GPS[5]) {
				map = m3;
				GPS[5] = 0;
				GPS[2] = 1; //change to m5
				userR = userR + M_PI;
				userY = 0.8f;
				userX = size - 1.f - userX;
				return 63;
			}
		}
		//down
		else if (userY > size - 0.5f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = m5;
				GPS[0] = 0;
				GPS[4] = 1; //change to m5
				userR = userR;
				userY = 0.8f;
				return 15;
			}
			else if (GPS[1]) {
				map = m5;
				GPS[1] = 0;
				GPS[4] = 1; //change to m5
				userR = userR - (M_PI / 2.f);
				userY = userX;
				userX = size - 1.f - 0.8f;
				return 25;
			}
			else if (GPS[4]) {
				map = m6;
				GPS[4] = 0;
				GPS[5] = 1; //change to m6
				userR = userR + M_PI;
				userY = size - 1.f - 0.8f;
				userX = size - 1.f - userX;
				return 56;
			}
			else if (GPS[3]) {
				map = m5;
				GPS[3] = 0;
				GPS[4] = 1; //change to m5
				userR = userR + M_PI / 2.f;
				userY = size - 1.f - userX;
				userX = 0.8f;
				return 45;
			}
			else if (GPS[2]) {
				map = m1;
				GPS[2] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userY = 0.8f;
				return 31;
			}
			else if (GPS[5]) {
				map = m5;
				GPS[5] = 0;
				GPS[4] = 1; //change to m3
				userR = userR - M_PI;
				userY = size - 1.f - 0.8f;;
				userX = size - 1.f - userX;
				return 65;
			}
		}
		//left
		else if (userX < 0.8f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = m4;
				GPS[0] = 0;
				GPS[3] = 1; //change to m4
				userR = userR;
				userX = size - 1.f - 0.8f;
				userY = userY;
				return 14;
			}
			else if (GPS[1]) {
				map = m1;
				GPS[1] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userX = size - 1.f - 0.8f;
				return 21;
			}
			else if (GPS[4]) {
				map = m4;
				GPS[4] = 0;
				GPS[3] = 1; //change to m4
				userR = userR - (M_PI / 2.f);
				userY = size - 1.f - 0.8f;
				userX = size - 1.f - userY;
				return 54;
			}
			else if (GPS[3]) {
				map = m6;
				GPS[3] = 0;
				GPS[5] = 1; //change to m6
				userR = userR;
				userX = size - 1.f - 0.8f;
				return 46;
			}
			else if (GPS[2]) {
				map = m4;
				GPS[2] = 0;
				GPS[3] = 1; //change to m4
				userR = userR + (M_PI / 2.f);
				userY = 0.8f;
				return 34;
			}
			else if (GPS[5]) {
				map = m2;
				GPS[5] = 0;
				GPS[1] = 1; //change to m2
				userR = userR;
				userX = size - 1.f - 0.8f;
				return 62;
			}
		}
		//right
		else if (userX > size - 0.8f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = m2;
				GPS[0] = 0;
				GPS[1] = 1; //change to m2
				userR = userR;
				userX = 0.8f;
				return 12;
			}
			else if (GPS[1]) {
				map = m6;
				GPS[1] = 0;
				GPS[5] = 1; //change to m6
				userR = userR;
				userX = 0.8f;
				return 26;
			}
			else if (GPS[4]) {
				map = m2;
				GPS[4] = 0;
				GPS[1] = 1; //change to m2
				userR = userR + (M_PI / 2.f);
				userX = userY;
				userY = size - 1.f - 0.8f;
				return 52;
			}
			else if (GPS[3]) {
				map = m1;
				GPS[3] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userX = 0.8f;
				userY = userY;
				return 41;
			}
			else if (GPS[2]) {
				map = m2;
				GPS[2] = 0;
				GPS[1] = 1; //change to m2
				userR = userR - M_PI / 2.f;
				userX = size - 1.f - userY;
				userY = 0.8f;
				return 32;
			}
			else if (GPS[5]) {
				map = m4;
				GPS[5] = 0;
				GPS[3] = 1; //change to m4
				userR = userR;
				userX = 0.8f;
				return 64;
			}
		}
	}
	else {
		//up
		if (userY < 0.5f) {
			//printf("user x = %lf, user y = %lf\n\n", userX, userY);
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = s_m3;
				GPS[0] = 0;
				GPS[2] = 1; //change to m5
				userR = userR;
				userY = (int)(size - 1.f - 0.8f);
				return 13;
			}
			else if (GPS[1]) {
				map = s_m3;
				GPS[1] = 0;
				GPS[2] = 1; //change to m5
				userR = userR + (M_PI / 2.f);
				userY = (size - 1.f - userX);
				userX = (size - 1.f - 0.8f);
				return 23;
			}
			else if (GPS[4]) {
				map = s_m1;
				GPS[4] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userY = size - 1.f - 0.8f;
				return 51;
			}
			else if (GPS[3]) {
				map = s_m3;
				GPS[3] = 0;
				GPS[2] = 1; //change to m5
				userR = userR - M_PI / 2.f;
				userY = userX;
				userX = 0.8f;
				return 43;
			}
			else if (GPS[2]) {
				map = s_m6;
				GPS[2] = 0;
				GPS[5] = 1; //change to m6
				userR = userR + M_PI;
				shiftX *= -1;
				shiftY *= -1;
				userX = size - 1.f - userX;
				userY = 0.8f;
				return 36;
			}
			else if (GPS[5]) {
				map = s_m3;
				GPS[5] = 0;
				GPS[2] = 1; //change to m5
				userR = userR + M_PI;
				shiftX *= -1;
				shiftY *= -1;
				userY = 0.8f;
				userX = size - 1.f - userX;
				return 63;
			}
		}
		//down
		else if (userY > size - 0.5f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = s_m5;
				GPS[0] = 0;
				GPS[4] = 1; //change to m5
				userR = userR;
				userY = 0.8f;
				return 15;
			}
			else if (GPS[1]) {
				map = s_m5;
				GPS[1] = 0;
				GPS[4] = 1; //change to m5
				userR = userR - (M_PI / 2.f);
				userY = userX;
				userX = size - 1.f - 0.8f;
				return 25;
			}
			else if (GPS[4]) {
				map = s_m6;
				GPS[4] = 0;
				GPS[5] = 1; //change to m6
				userR = userR + M_PI;
				shiftX *= -1;
				shiftY *= -1;
				userY = size - 1.f - 0.8f;
				userX = size - 1.f - userX;
				return 56;
			}
			else if (GPS[3]) {
				map = s_m5;
				GPS[3] = 0;
				GPS[4] = 1; //change to m5
				userR = userR + M_PI / 2.f;
				userY = size - 1.f - userX;
				userX = 0.8f;
				return 45;
			}
			else if (GPS[2]) {
				map = s_m1;
				GPS[2] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userY = 0.8f;
				return 31;
			}
			else if (GPS[5]) {
				map = s_m5;
				GPS[5] = 0;
				GPS[4] = 1; //change to m3
				userR = userR - M_PI;
				shiftX *= -1;
				shiftY *= -1;
				userY = size - 1.f - 0.8f;;
				userX = size - 1.f - userX;
				return 65;
			}
		}
		//left
		else if (userX < 0.8f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = s_m4;
				GPS[0] = 0;
				GPS[3] = 1; //change to m4
				userR = userR;
				userX = size - 1.f - 0.8f;
				userY = userY;
				return 14;
			}
			else if (GPS[1]) {
				map = s_m1;
				GPS[1] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userX = size - 1.f - 0.8f;
				return 21;
			}
			else if (GPS[4]) {
				map = s_m4;
				GPS[4] = 0;
				GPS[3] = 1; //change to m4
				userR = userR - (M_PI / 2.f);
				userY = size - 1.f - 0.8f;
				userX = size - 1.f - userY;
				return 54;
			}
			else if (GPS[3]) {
				map = s_m6;
				GPS[3] = 0;
				GPS[5] = 1; //change to m6
				userR = userR;
				userX = size - 1.f - 0.8f;
				return 46;
			}
			else if (GPS[2]) {
				map = s_m4;
				GPS[2] = 0;
				GPS[3] = 1; //change to m4
				userR = userR + (M_PI / 2.f);
				userY = 0.8f;
				return 34;
			}
			else if (GPS[5]) {
				map = s_m2;
				GPS[5] = 0;
				GPS[1] = 1; //change to m2
				userR = userR;
				userX = size - 1.f - 0.8f;
				return 62;
			}
		}
		//right
		else if (userX > size - 0.8f) {
			//rotateFlag = true;
			//userR has no need to change
			if (GPS[0]) {
				map = s_m2;
				GPS[0] = 0;
				GPS[1] = 1; //change to m2
				userR = userR;
				userX = 0.8f;
				return 12;
			}
			else if (GPS[1]) {
				map = s_m6;
				GPS[1] = 0;
				GPS[5] = 1; //change to m6
				userR = userR;
				userX = 0.8f;
				return 26;
			}
			else if (GPS[4]) {
				map = s_m2;
				GPS[4] = 0;
				GPS[1] = 1; //change to m2
				userR = userR + (M_PI / 2.f);
				userX = userY;
				userY = size - 1.f - 0.8f;
				return 52;
			}
			else if (GPS[3]) {
				map = s_m1;
				GPS[3] = 0;
				GPS[0] = 1; //change to m1
				userR = userR;
				userX = 0.8f;
				userY = userY;
				return 41;
			}
			else if (GPS[2]) {
				map = s_m2;
				GPS[2] = 0;
				GPS[1] = 1; //change to m2
				userR = userR - M_PI / 2.f;
				userX = size - 1.f - userY;
				userY = 0.8f;
				return 32;
			}
			else if (GPS[5]) {
				map = s_m4;
				GPS[5] = 0;
				GPS[3] = 1; //change to m4
				userR = userR;
				userX = 0.8f;
				return 64;
			}
		}
	}
	return 0;
}

enemy enemyINFO[6];
void move() {

	userX += shiftX;
	userY += shiftY;
	//up

	userR += shiftR;
	userR = fmod(userR, M_PI * 2.f);
	int rotateINFO = rotate();
	if (rotateINFO) {
		enemyINFO[rotateINFO / 10 - 1].x = enemyX;
		enemyINFO[rotateINFO / 10 - 1].y = enemyY;

		enemyX = enemyINFO[rotateINFO % 10 - 1].x;
		enemyY = enemyINFO[rotateINFO % 10 - 1].y;

		enemypaper = loadImgTexture(enemyimg[enemyINFO[rotateINFO % 10 - 1].img], 1, 1, 1, true, 0xFF, 0xFF, 0xFF);
		enemytexture = IMG_LoadTexture(renderer, enemyimg[enemyINFO[rotateINFO % 10 - 1].img]);
		if (SDL_SetTextureBlendMode(enemytexture, SDL_BLENDMODE_BLEND) == -1)
		{
			printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
			exit(1);
		}
	}
	while (map[(int)userY][(int)userX] == '#')
	{
		userX -= 1.1 * shiftX;
		userY -= 1.1 * shiftY;
	}
}

void raycast(void)
{
	for (int x = 0; x < sWIDTH; x++) {

		rayAngle = (userR - fov / 2) + ((float)x / (float)sWIDTH) * fov;
		dToWall = 0.0f;
		dToEnemy = 0.0f;
		dToPort = 0.0f;

		hitWall = false;
		hitEnemy = false;
		hitPort = false;
		EhitWall = false;
		boundary = false;
		eyeX = sinf(rayAngle);
		eyeY = cosf(rayAngle);

		WsampleX = 0.0f;

		//background
		//while (!hitWall && !hitEnemy && dToWall < depth) {
		while (!hitWall && dToWall < depth) {
			dToWall += rayStepSize;
			if (!hitEnemy)
				dToEnemy += rayStepSize;
			//if(!hitPort)
				//dToPort += rayStepSize;
			testX = (int)(userX + eyeX * dToWall);
			testY = (int)(userY + eyeY * dToWall);

			if (testX < 0 || testX > mapSize - 1 || testY < 0 || testY > mapSize - 1) {
				hitWall = true;
				//hitEnemy = false;
				dToWall = depth;
			}
			else {
				if (map[testY][testX] == '#') {
					hitWall = true;
					/*vector<boundaryData> p;
					float vx, vy, d, dot;
					float bound = 0.005f;

					for (int tx = 0; tx < 2; tx++) {
						for (int ty = 0; ty < 2; ty++) {
							vx = (float)testX + tx - userX;
							vy = (float)testY + ty - userY;
							d = sqrt(vx * vx + vy * vy);
							dot = eyeX * vx / d + eyeY * vy / d;
							p.push_back({ d, dot });
						}
					}
					sort_boundary(p);

					if (acos(p[0].dot) < bound) boundary = true;
					if (acos(p[1].dot) < bound) boundary = true;
					if (acos(p[2].dot) < bound) boundary = true;
					p.clear();*/


					float blockMidX = (float)testX + 0.5f;
					float blockMidY = (float)testY + 0.5f;

					float testPointX = userX + eyeX * dToWall;
					float testPointY = userY + eyeY * dToWall;

					float testAngle = atan2f((testPointY - blockMidY), (testPointX - blockMidX));

					if (testAngle >= -3.14159f * 0.25f && testAngle < 3.14159f * 0.25f)
						WsampleX = testPointY - (float)testY;
					if (testAngle >= 3.14159f * 0.25f && testAngle < 3.14159f * 0.75f)
						WsampleX = testPointX - (float)testX;
					if (testAngle < -3.14159f * 0.25f && testAngle >= -3.14159f * 0.75f)
						WsampleX = testPointX - (float)testX;
					if (testAngle >= 3.14159f * 0.75f || testAngle < -3.14159f * 0.75f)
						WsampleX = testPointY - (float)testY;

				}
				/*if (testX == (int)enemyX && testY == (int)enemyY) {
					hitEnemy = true;

					float blockMidX = (float)testX + 0.5f;
					float blockMidY = (float)testY + 0.5f;

					float testPointX = userX + eyeX * dToWall;
					float testPointY = userY + eyeY * dToWall;

					float testAngle = atan2f(-(testPointY - blockMidY), (testPointX - blockMidX));
					/*if (testAngle >= -3.14159f * 0.0f && testAngle < 3.14159f * 0.5f)
						EsampleX = 1-(testPointY - (float)testY);
					if (testAngle >= 3.14159f * 0.5f && testAngle < 3.14159f * 1.0f)
						EsampleX = testPointX - (float)testX;
					if (testAngle >= -3.14159f * 1.0f || testAngle < -3.14159f * 0.5f)
						EsampleX = testPointY - (float)testY;
					if (testAngle >= -3.14159f * 0.5f && testAngle < -3.14159f * 0.0f)
						EsampleX = 1-(testPointX - (float)testX);

					if (testAngle >= -3.14159f * 0.25f && testAngle < 3.14159f * 0.25f)
						EsampleX = 1-(testPointY - (float)testY);
					else if (testAngle >= 3.14159f * 0.25f && testAngle < 3.14159f * 0.75f)
						EsampleX = 1-(testPointX - (float)testX);
					else if (testAngle < -3.14159f * 0.25f && testAngle >= -3.14159f * 0.75f)
						EsampleX = testPointX - (float)testX;
					else if (testAngle >= 3.14159f * 0.75f || testAngle < -3.14159f * 0.75f)
						EsampleX = testPointY - (float)testY;
				}*/
				/*if (testX == portX && testY == portY) {
					hitPort = true;

					float blockMidX = (float)testX + 0.5f;
					float blockMidY = (float)testY + 0.5f;

					float testPointX = userX + eyeX * dToWall;
					float testPointY = userY + eyeY * dToWall;

					float testAngle = atan2f(-(testPointY - blockMidY), (testPointX - blockMidX));

					if (testAngle >= -3.14159f * 0.25f && testAngle < 3.14159f * 0.25f)
						PsampleX = 1-(testPointY - (float)testY);
					else if (testAngle >= 3.14159f * 0.25f && testAngle < 3.14159f * 0.75f)
						PsampleX = 1-(testPointX - (float)testX);
					else if (testAngle < -3.14159f * 0.25f && testAngle >= -3.14159f * 0.75f)
						PsampleX = testPointX - (float)testX;
					else if (testAngle >= 3.14159f * 0.75f || testAngle < -3.14159f * 0.75f)
						PsampleX = testPointY - (float)testY;
				}*/


			}
		}


		float shade = 0.0f, Eshade = 0.0f;
		/*	if (hitEnemy) {

			if (dToEnemy <= depth / 10.0f) Eshade = 1.0f;
			else if (dToEnemy < depth / 9.0f) Eshade = 0.9f;
			else if (dToEnemy < depth / 8.0f) Eshade = 0.8f;
			else if (dToEnemy < depth / 7.0f) Eshade = 0.7f;
			else if (dToEnemy < depth / 6.0f) Eshade = 0.6f;
			else if (dToEnemy < depth / 5.0f) Eshade = 0.5f;
			else if (dToEnemy < depth / 4.0f) Eshade = 0.4f;
			else if (dToEnemy < depth / 3.0f) Eshade = 0.3f;
			else if (dToEnemy < depth / 2.0f) Eshade = 0.2f;
			else if (dToEnemy < depth) Eshade = 0.1f;
			else Eshade = 0.0f;

			if (dToEnemy > 1.f) Eshade = 1.f / dToEnemy;
			else Eshade = 1.f;
		}*/


		if (dToWall <= depth / 10.0f) shade = 1.0f;
		else if (dToWall < depth / 9.0f) shade = 0.9f;
		else if (dToWall < depth / 8.0f) shade = 0.8f;
		else if (dToWall < depth / 7.0f) shade = 0.7f;
		else if (dToWall < depth / 6.0f) shade = 0.6f;
		else if (dToWall < depth / 5.0f) shade = 0.5f;
		else if (dToWall < depth / 4.0f) shade = 0.4f;
		else if (dToWall < depth / 3.0f) shade = 0.3f;
		else if (dToWall < depth / 2.0f) shade = 0.2f;
		else if (dToWall < depth) shade = 0.1f;
		else shade = 0.0f;


		shade = 1;

		/*if (dToWall > 1.f) shade = 1.f / dToWall;
		else shade = 1.f;
		if (boundary) shade = 0.0f;*/


		int ceiling = (float)sHEIGHT / 2.0f - sHEIGHT / (float)dToWall;
		int floor = sHEIGHT - ceiling;
		int Eceiling = (float)sHEIGHT / 2.0f - (sHEIGHT / (float)dToEnemy) * 0.5;
		int Efloor = sHEIGHT - Eceiling;
		int Pceiling = (float)sHEIGHT / 2.0f - (sHEIGHT / (float)dToPort);
		int Pfloor = sHEIGHT - Pceiling;

		// Update Depth Buffer
		depthBuffer[x] = dToWall;


		for (int y = 0; y < sHEIGHT; y++) {

			//ceiling
			if (y <= ceiling) shadowC[y][x] = 1.0f;
			//wall
			if (y > ceiling && y < floor) {
				if (dToWall < depth) {
					float WsampleY = ((float)y - (float)ceiling) / ((float)floor - (float)ceiling);
					int w = int((float)wallpaper.height / ((float)floor - (float)ceiling));
					if (w == 0) w = 1;
					wallPos[y][x].src = { int(WsampleX * (float)wallpaper.width), int(WsampleY * (float)wallpaper.height),FirstW / sWIDTH, w };//int((float)wallpaper.height / ((float)floor - (float)ceiling))
					wallPos[y][x].dst = { 200 + x * (FirstW / sWIDTH), y * (FirstH / sHEIGHT), FirstW / sWIDTH, FirstH / sHEIGHT };
					//printf("%d %d", wallPos[y][x].src.x);
					//shadowW[y][x] = shade
					shadowW[y][x] = ((float)y - sHEIGHT / 2.f) / sHEIGHT;
				}
				else
					shadowW[y][x] = 0;
			}

			//ground
			else {
				float b = 1.0f - (((float)y - sHEIGHT / 2.0f) / ((float)sHEIGHT / 2.0f));
				/*if (b < 1. / 10) shadowG[y][x] = 1;
				else if (b < 1. / 9) shadowG[y][x] = 0.9f;
				else if (b < 1. / 8) shadowG[y][x] = 0.8f;
				else if (b < 1. / 7) shadowG[y][x] = 0.7f;
				else if (b < 1. / 6) shadowG[y][x] = 0.6f;
				else if (b < 1. / 5) shadowG[y][x] = 0.5f;
				else if (b < 1. / 4) shadowG[y][x] = 0.4f;
				else if (b < 1. / 3) shadowG[y][x] = 0.3f;
				else if (b < 1. / 2) shadowG[y][x] = 0.2f;
				else if (b < 1. / 1) shadowG[y][x] = 0.1f;
				else shadowG[y][x] = 0.0f;
				*/

				shadowG[y][x] = ((float)y - sHEIGHT / 2.f) / sHEIGHT;
			}

			/*if (y > Eceiling && y < Efloor && hitEnemy) {
				if (dToEnemy < depth) {
					float EsampleY = ((float)y - (float)Eceiling) / ((float)Efloor - (float)Eceiling);
					EImgPos[y][x].src = { int(EsampleX * (float)enemypaper.width), int(EsampleY * (float)enemypaper.height), int((float)wallpaper.width / ((float)Efloor - (float)Eceiling)), HEIGHT / sHEIGHT };
					EImgPos[y][x].dst = { x * (WIDTH / sWIDTH), y * (HEIGHT / sHEIGHT), WIDTH / sWIDTH, HEIGHT / sHEIGHT };
					shadowE[y][x] = 1;
				}
				else
					shadowE[y][x] = 0;
				//shadowE[y][x] = Eshade;
			}*/
			//port
			/*if (y > Pceiling && y < Pfloor && hitPort) {
				if (dToPort < depth) {
					float PsampleY = ((float)y - (float)Pceiling) / ((float)Pfloor - (float)Pceiling);
					PImgPos[y][x].src = { int(PsampleX * (float)portpaper.width), int(PsampleY * (float)portpaper.height), int((float)portpaper.width / ((float)Pfloor - (float)Pceiling)), HEIGHT / sHEIGHT };
					PImgPos[y][x].dst = { x * (WIDTH / sWIDTH), y * (HEIGHT / sHEIGHT), WIDTH / sWIDTH, HEIGHT / sHEIGHT };
					shadowP[y][x] = 1;
				}
				else
					shadowP[y][x] = 0;
				//shadowE[y][x] = Eshade;
			}*/

		}
		/*
		//enemy object
		eyeX = sinf(rayAngle);
		eyeY = -cosf(rayAngle);
		hitWall = false;
		hitEnemy = false;
		dToEnemy = 0.0f;
		//while (!hitEnemy && dToEnemy < depth) {
		while (!hitEnemy && !hitWall && dToEnemy < depth) {
			dToEnemy += rayStepSize;
			testX = (int)(userX + eyeX * dToEnemy);
			testY = (int)(userY + eyeY * dToEnemy);
			if (testX < 0 || testX >= mapSize || testY < 0 || testY >= mapSize) {
				hitEnemy = false;
				hitWall = true;
				dToEnemy = 0.0f;
			}
			else if (map[testY][testX] == 'E') {
				hitEnemy = true;

				vector<boundaryData> p;
				float vx, vy, d, dot;
				float bound = 0.005f;

				for (int tx = 0; tx < 2; tx++) {
					for (int ty = 0; ty < 2; ty++) {
						vx = (float)testX + tx - userX;
						vy = (float)testY + ty - userY;
						d = sqrt(vx * vx + vy * vy);
						dot = eyeX * vx / d + eyeY * vy / d;
						p.push_back({ d, dot });
					}
				}
				sort_boundary(p);

				if (acos(p[0].dot) < bound) boundary = true;
				if (acos(p[1].dot) < bound) boundary = true;
				if (acos(p[2].dot) < bound) boundary = true;
				p.clear();
			}
			else if (map[testY][testX] == '#') {
				hitEnemy = false;
				hitWall = true;
				//dToEnemy = 0.0f;
			}
		}
		float Eshade = 0.0f;
		if (hitEnemy) {
			Eshade = 1.0f;
			if (boundary) Eshade = 0.0f;

			if (dToEnemy <= depth / 10.0f) Eshade = 1.0f;
			else if (dToEnemy < depth / 9.0f) Eshade = 0.9f;
			else if (dToEnemy < depth / 8.0f) Eshade = 0.8f;
			else if (dToEnemy < depth / 7.0f) Eshade = 0.7f;
			else if (dToEnemy < depth / 6.0f) Eshade = 0.6f;
			else if (dToEnemy < depth / 5.0f) Eshade = 0.5f;
			else if (dToEnemy < depth / 4.0f) Eshade = 0.4f;
			else if (dToEnemy < depth / 3.0f) Eshade = 0.3f;
			else if (dToEnemy < depth / 2.0f) Eshade = 0.2f;
			else if (dToEnemy < depth) Eshade = 0.1f;
			else Eshade = 0.0f;

		}
		else {
			Eshade = 0.f;
		}

		for (int y = 0; y < sHEIGHT; y++) {
			//&& (y > (ceiling + 50) && y < (floor - 50))
			if ((y > (ceiling -20)) && (y < (floor +20))) {
				shadowE[y][x] = Eshade;
			}
			else shadowE[y][x] = 0.f;
		}
		*/
	}
}

void enemycast(void) {

	// Can object be seen?
	float vecX;
	float vecY;
	float a = 0.4;
	if (map[(int)enemyY][int(enemyX + a)] == '#') {
		if (map[(int)(enemyY + a)][int(enemyX + a)] == '#' || map[(int)(enemyY - a)][int(enemyX + a)] == '#') {
			//vecX = ((int)enemyX - userX)+0.5 ;
			//vecY = ((int)enemyY - userY) + 0.5;
			enemyX = (int)enemyX + 0.5;
			enemyY = (int)enemyY + 0.5;
		}
		else {
			//vecX = ((int)enemyX - userX) + 0.5;
			//vecY = enemyY - userY;
			enemyX = (int)enemyX + 0.5;
		}
	}
	else if (map[(int)enemyY][int(enemyX - a)] == '#') {
		if (map[(int)(enemyY + a)][int(enemyX - a)] == '#' || map[(int)(enemyY - a)][int(enemyX - a)] == '#') {
			//vecX = ((int)enemyX - userX) + 0.5;
			//vecY = ((int)enemyY - userY) + 0.5;
			enemyX = (int)enemyX + 0.5;
			enemyY = (int)enemyY + 0.5;
		}
		else {
			//vecX = ((int)enemyX - userX) + 0.5;
			//vecY = enemyY - userY;
			enemyX = (int)enemyX + 0.5;
		}
	}
	else if (map[(int)(enemyY + a)][int(enemyX)] == '#') {
		if (map[(int)(enemyY + a)][int(enemyX + a)] == '#' || map[(int)(enemyY + a)][int(enemyX - a)] == '#') {
			//vecX = ((int)enemyX - userX) + 0.5;
			//vecY = ((int)enemyY - userY) + 0.5;
			enemyX = (int)enemyX + 0.5;
			enemyY = (int)enemyY + 0.5;
		}
		else {
			//vecX = enemyX - userX;
			//vecY = ((int)enemyY - userY) + 0.5;
			enemyY = (int)enemyY + 0.5;
		}
	}
	else if (map[(int)(enemyY - a)][int(enemyX)] == '#') {
		if (map[(int)(enemyY - a)][int(enemyX + a)] == '#' || map[(int)(enemyY - a)][int(enemyX - a)] == '#') {
			//vecX = ((int)enemyX - userX) + 0.5;
			//vecY = ((int)enemyY - userY) + 0.5;
			enemyX = (int)enemyX + 0.5;
			enemyY = (int)enemyY + 0.5;
		}
		else {
			//vecX = enemyX - userX;
			//vecY = ((int)enemyY - userY) + 0.5;
			enemyY = (int)enemyY + 0.5;
		}
	}
	/*else {
		vecX = enemyX - userX;
		vecY = enemyY - userY;
	}*/

	vecX = enemyX - userX;
	vecY = enemyY - userY;

	//float vecX = ((int)enemyX - userX) + 0.5 + (enemyX - (int)enemyX - 0.5f) * 0.5f;
	//float vecY = ((int)enemyY - userY) + 0.5 + (enemyY - (int)enemyY - 0.5f) * 0.5f;

	/*if (map[(int)(enemyX+0.5)][(int)(enemyX + 0.5)] == '#') {
		vecX -= 0.5;
		vecY -= 0.5;
	}*/

	float dFromUser = sqrtf(vecX * vecX + vecY * vecY);

	float eyeX = sinf(userR);
	float eyeY = cosf(userR);

	// Calculate angle between enemy and players feet, and players looking direction
	// to determine if the enemy is in the players field of view
	float objectAngle = atan2f(eyeY, eyeX) - atan2f(vecY, vecX);
	//float objectAngle = M_PI / 2.f - userR - fov / 2.f - atan2f(vecY, vecX);
	if (objectAngle < -3.14159f)
		objectAngle += 2.0f * 3.14159f;
	if (objectAngle > 3.14159f)
		objectAngle -= 2.0f * 3.14159f;


	bool bInUserFOV = fabs(objectAngle) < fov / 2.0f;

	if (bInUserFOV && dFromUser >= 0.5f && dFromUser < depth)
	{
		/*float ceiling = (float)(sHEIGHT / 2.0) - (float)sHEIGHT / (dFromUser);
		//printf("c=%f\n", ceiling);
		float floor = sHEIGHT - ceiling;
		float height = floor - ceiling;
		float objectAspectRatio = (float)enemypaper.height / (float)enemypaper.width;
		float width = (float)height / objectAspectRatio;
		float middleOfObject = (0.5f * (objectAngle / (fov / 2.0f)) + 0.5f) * (float)sWIDTH;
		*/
		int ceiling = int(((float)sHEIGHT / 2.0) - ((float)sHEIGHT / (dFromUser)));
		int Wceiling = int(((float)sHEIGHT / 2.0) - ((float)sHEIGHT / (dFromUser)));
		int floor = sHEIGHT - ceiling;
		int Wfloor = sHEIGHT - Wceiling;
		int height = floor - ceiling;
		int Wheight = Wfloor - Wceiling;
		float objectAspectRatio = (float)enemypaper.height / (float)enemypaper.width;
		int width = int((float)height / objectAspectRatio);
		int Wwidth = int((float)Wheight / objectAspectRatio);
		int middleOfObject = int((0.5f * (objectAngle / (fov / 2.0f)) + 0.5f) * (float)sWIDTH);// + Wwidth/2.f/(float)sWIDTH

		//printf("%d %d\n", width, height);
		/*for (float lx = 0; lx < width; lx++)
		{
			for (float ly = 0; ly < height; ly++)
			{
				float EsampleX = lx / width;
				float EsampleY = ly / height;
				int objectColumn = (int)(middleOfObject + lx - (width / 2.0f));

				if (objectColumn >= 0 && objectColumn < sWIDTH)
					if (depthBuffer[objectColumn] >= dFromUser)
					{
						if (int(ceiling + ly) >= 0 && int(ceiling + ly) < sHEIGHT) {
							printf("in\n");
							shadowE[int(ceiling + ly)][objectColumn] = 1.0;
							EImgPos[int(ceiling + ly)][objectColumn].src = { int(EsampleX * (float)enemypaper.width), int(EsampleY * (float)enemypaper.height), WIDTH / sWIDTH, int((float)wallpaper.height / height) };
							EImgPos[int(ceiling + ly)][objectColumn].dst = { objectColumn * (WIDTH / sWIDTH), (int)(ceiling + ly) * (HEIGHT / sHEIGHT), WIDTH / sWIDTH, HEIGHT / sHEIGHT };
							depthBuffer[objectColumn] = dFromUser;
						}
						//else shadowE[int(ceiling + ly)][objectColumn] = 0;
					}
			}

		}
		*/
		for (int lx = 0; lx < width; lx++)
		{
			for (int ly = 0; ly < height; ly++)
			{
				float EsampleX = (float)lx / (float)width;
				float EsampleY = (float)ly / (float)height;
				int objectColumn = (int)(middleOfObject + lx - (width / 2.0f));

				if (objectColumn >= 0 && objectColumn < sWIDTH)
					if (depthBuffer[objectColumn] >= dFromUser)
					{
						if (ceiling + ly >= 0 && ceiling + ly < sHEIGHT) {
							shadowE[ceiling + ly][objectColumn] = 1.0;
							int h = int((float)wallpaper.height / (float)height);
							if (h == 0)h = 1;
							int w = int((float)wallpaper.width / (float)width);
							if (w == 0)w = 1;
							EImgPos[ceiling + ly][objectColumn].src = { int(EsampleX * (float)enemypaper.width), int(EsampleY * (float)enemypaper.height), w, h };// FirstW / sWIDTHint((float)wallpaper.height / (float)height)
							EImgPos[ceiling + ly][objectColumn].dst = { objectColumn * (FirstW / sWIDTH) + 200, (int)(ceiling + ly) * (FirstH / sHEIGHT), FirstW / sWIDTH, FirstH / sHEIGHT };
							depthBuffer[objectColumn] = dFromUser;
						}
						//else shadowE[int(ceiling + ly)][objectColumn] = 0;
					}
			}
		}


	}
}

void portcast(int num) {

	// Can object be seen?
	float vecX = (port[num].x - userX) + 0.5;
	float vecY = (port[num].y - userY) + 0.5;
	/*if (map[(int)(enemyX+0.5)][(int)(enemyX + 0.5)] == '#') {
		vecX -= 0.5;
		vecY -= 0.5;
	}*/

	float dFromUser = sqrtf(vecX * vecX + vecY * vecY);

	float eyeX = sinf(userR);
	float eyeY = cosf(userR);

	float objectAngle = atan2f(eyeY, eyeX) - atan2f(vecY, vecX);
	//float objectAngle = M_PI / 2.f - userR - fov / 2.f - atan2f(vecY, vecX);
	if (objectAngle < -3.14159f)
		objectAngle += 2.0f * 3.14159f;
	if (objectAngle > 3.14159f)
		objectAngle -= 2.0f * 3.14159f;


	bool bInUserFOV = fabs(objectAngle) < fov / 2.0f;

	if (bInUserFOV && dFromUser >= 0.5f && dFromUser < depth)
	{
		/*float ceiling = (float)(sHEIGHT / 2.0) - (float)sHEIGHT / (dFromUser);
		//printf("c=%f\n", ceiling);
		float floor = sHEIGHT - ceiling;
		float height = floor - ceiling;
		float objectAspectRatio = (float)enemypaper.height / (float)enemypaper.width;
		float width = (float)height / objectAspectRatio;
		float middleOfObject = (0.5f * (objectAngle / (fov / 2.0f)) + 0.5f) * (float)sWIDTH;
		*/
		int ceiling = int(((float)sHEIGHT / 2.0) - ((float)sHEIGHT / (dFromUser)) * 0.5);
		int Wceiling = int(((float)sHEIGHT / 2.0) - ((float)sHEIGHT / (dFromUser)));
		int floor = sHEIGHT - ceiling;
		int Wfloor = sHEIGHT - Wceiling;
		int height = floor - ceiling;
		int Wheight = Wfloor - Wceiling;
		float objectAspectRatio = (float)portpaper.height / (float)portpaper.width;
		int width = int((float)height / objectAspectRatio);
		int Wwidth = int((float)Wheight / objectAspectRatio);
		int middleOfObject = int((0.5f * (objectAngle / (fov / 2.0f)) + 0.5f) * (float)sWIDTH);// + Wwidth/2.f/(float)sWIDTH

		//printf("%d %d\n", width, height);
		/*for (float lx = 0; lx < width; lx++)
		{
			for (float ly = 0; ly < height; ly++)
			{
				float EsampleX = lx / width;
				float EsampleY = ly / height;
				int objectColumn = (int)(middleOfObject + lx - (width / 2.0f));

				if (objectColumn >= 0 && objectColumn < sWIDTH)
					if (depthBuffer[objectColumn] >= dFromUser)
					{
						if (int(ceiling + ly) >= 0 && int(ceiling + ly) < sHEIGHT) {
							printf("in\n");
							shadowE[int(ceiling + ly)][objectColumn] = 1.0;
							EImgPos[int(ceiling + ly)][objectColumn].src = { int(EsampleX * (float)enemypaper.width), int(EsampleY * (float)enemypaper.height), WIDTH / sWIDTH, int((float)wallpaper.height / height) };
							EImgPos[int(ceiling + ly)][objectColumn].dst = { objectColumn * (WIDTH / sWIDTH), (int)(ceiling + ly) * (HEIGHT / sHEIGHT), WIDTH / sWIDTH, HEIGHT / sHEIGHT };
							depthBuffer[objectColumn] = dFromUser;
						}
						//else shadowE[int(ceiling + ly)][objectColumn] = 0;
					}
			}

		}
		*/
		for (int lx = 0; lx < width; lx++)
		{
			for (int ly = 0; ly < height; ly++)
			{
				float PsampleX = (float)lx / (float)width;
				float PsampleY = (float)ly / (float)height;
				int objectColumn = (int)(middleOfObject + lx - (width / 2.0f));

				if (objectColumn >= 0 && objectColumn < sWIDTH)
					if (depthBuffer[objectColumn] >= dFromUser)
					{
						if (ceiling + ly >= 0 && ceiling + ly < sHEIGHT) {
							shadowP[ceiling + ly][objectColumn] = 1.0;
							int h = int((float)portpaper.height / (float)height);
							if (h == 0)h = 1;
							int w = int((float)portpaper.width / (float)width);
							if (w == 0)w = 1;
							PImgPos[ceiling + ly][objectColumn].src = { int(PsampleX * (float)portpaper.width), int(PsampleY * (float)portpaper.height), w, h };//FirstW / sWIDTH
							PImgPos[ceiling + ly][objectColumn].dst = { objectColumn * (FirstW / sWIDTH) + 200, (int)(ceiling + ly) * (FirstH / sHEIGHT), FirstW / sWIDTH, FirstH / sHEIGHT };
							depthBuffer[objectColumn] = dFromUser;
						}
						//else shadowP[int(ceiling + ly)][objectColumn] = 0;
					}
			}
		}


	}
}


double camshiftX = 0, camshiftY = 0, yawshift = 0;
vec3d vHor = { 0,0,0 }, vFor = { 0,0,0 };

void viewChange() {
	//vHor = Vector_Mul(vHor, 0.01);
	//vFor = Vector_Mul(vFor, 0.01);
	vCamera = Vector_Add(vCamera, vHor);
	vCamera = Vector_Add(vCamera, vFor);
	//printf("%lf\n", vCamera.x);
	//vCamera = Vector_Mul(vCamera, 0.1);
	//vCamera.x += camshiftX;
	vCamera.y += camshiftY;
	//printf("%d %d\n", vCamera.x, vCamera.y);
	fYaw += yawshift;
}


int main(int argc, char* args[])
{
	if (initSDL())
	{
		printf("Failed to initialize SDL!\n");
		return -1;
	}


	bool quit = false;
	SDL_Event e;
	MouseState mouseState;
	int mouseX, mouseY;

	int midMouseX = WIDTH / 2, midMouseY = HEIGHT / 2;

	SDL_RendererFlip no = SDL_FLIP_NONE;
	SDL_RendererFlip ho = SDL_FLIP_HORIZONTAL;
	SDL_RendererFlip ve = SDL_FLIP_VERTICAL;
	SDL_RendererFlip hove = (SDL_RendererFlip)(SDL_FLIP_HORIZONTAL | SDL_FLIP_VERTICAL);

	//****************** ENTRANCE *******************

	TextData title1 = loadTextTexture("BATTLE", "fonts/Melted Monster.ttf", 80, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData title2 = loadTextTexture("ROYALE", "fonts/Melted Monster.ttf", 80, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData title3 = loadTextTexture("3D", "fonts/Melted Monster.ttf", 80, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData press = loadTextTexture("--- PRESS SPACE ---", "fonts/PressStart2P.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	int titleINFO[4] = { 0,6, WIDTH / 2 - title1.width / 2 - 100, -title1.height };
	int dropINFO1[3] = { WIDTH / 2 - title1.width / 2 - title2.width / 2 - title3.width / 2 - 50,-title1.height,title1.height };//y
	int dropINFO2[3] = { WIDTH / 2 - title2.width / 2 + 60,-title1.height - title2.height,title1.height };//y
	int dropINFO3[3] = { WIDTH / 2 + title2.width / 2 + title3.width / 2 + 50,-title1.height - title2.height - title3.height,title1.height };//y
	bool initial = true;



	//****************** MENU *******************

	TextData HOME = loadTextTexture("HOME", "fonts/Impacted2.0.ttf", 60, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_HOME = loadTextTexture("HOME", "fonts/Impacted2.0.ttf", 60, 255, 255, 0, BLENDED, NULL, NULL, NULL);
	TextData menuTitle = loadTextTexture("BATTLE ROYALE 3D", "fonts/Impacted2.0.ttf", 70, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData MODE = loadTextTexture("MODE", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData WORLD = loadTextTexture("WORLD", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData MUSIC = loadTextTexture("MUSIC", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData survival = loadTextTexture("SURVIVAL", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_survival = loadTextTexture("SURVIVAL", "fonts/Impacted2.0.ttf", 30, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData TIME = loadTextTexture("TIMER", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_TIME = loadTextTexture("TIMER", "fonts/Impacted2.0.ttf", 30, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData simple = loadTextTexture("SIMPLE", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_simple = loadTextTexture("SIMPLE", "fonts/Impacted2.0.ttf", 30, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData classic = loadTextTexture("CLASSIC", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_classic = loadTextTexture("CLASSIC", "fonts/Impacted2.0.ttf", 30, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData areUready = loadTextTexture("ARE U READY ... ?", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData on = loadTextTexture("ON", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_on = loadTextTexture("ON", "fonts/Impacted2.0.ttf", 30, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData off = loadTextTexture("OFF", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_off = loadTextTexture("OFF", "fonts/Impacted2.0.ttf", 30, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData play = loadTextTexture("P L A Y", "fonts/Impacted2.0.ttf", 50, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_play = loadTextTexture("P L A Y", "fonts/Impacted2.0.ttf", 50, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData QUIT = loadTextTexture("EXIT", "fonts/Impacted2.0.ttf", 50, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData H_QUIT = loadTextTexture("EXIT", "fonts/Impacted2.0.ttf", 50, 255, 255, 0, BLENDED, NULL, NULL, NULL);

	t_HOME = { 64,58, HOME.width, HOME.height };
	t_menuTitle = { WIDTH / 2 - menuTitle.width / 2,  182 - 8 - menuTitle.height / 2 };
	t_MODE = { WIDTH / 2 - 182 * 3 / 2 - 80 + 182 / 2 - MODE.width / 2 ,438 + 182 - areUready.height - 40 - 290 - MODE.height / 2 };
	t_WORLD = { WIDTH / 2 - WORLD.width / 2,438 + 182 - areUready.height - 40 - 290 - MODE.height / 2 };
	t_MUSIC = { WIDTH / 2 + 182 + 80 - MUSIC.width / 2,438 + 182 - areUready.height - 40 - 290 - MODE.height / 2 };
	t_QUIT = { WIDTH - 60 - QUIT.width,HEIGHT - 58 - QUIT.height + 20, QUIT.width, QUIT.height };
	t_PLAY = { WIDTH / 2 - play.width / 2, 182 + 438 - play.height / 2 - 3, play.width, play.height };

	c_TIME;
	c_TIME.w = survival.width + 20;
	c_TIME.h = survival.height + 20;
	c_TIME.x = WIDTH / 2 - 182 * 3 / 2 - 80 + (182) / 2 - survival.width / 2 - 10;
	c_TIME.y = 6 + 438 + 182 - areUready.height - 40 - 290 + (438 + 182 - areUready.height - 40 - (438 + 182 - areUready.height - 40 - 290)) / 2 - 30 - survival.height - 20;

	c_survival = c_TIME;
	c_survival.y = c_TIME.y + c_TIME.h + 60;

	c_simple = c_TIME;
	c_simple.x = WIDTH / 2 - 182 * 1 / 2 + 182 / 2 - survival.width / 2 - 10;
	c_simple.y = c_TIME.y;
	c_classic = c_simple;
	c_classic.y = c_survival.y;

	c_on = c_TIME;
	c_on.x = WIDTH / 2 + 182 * 1 / 2 + 80 + (182) / 2 - survival.width / 2 - 10;
	c_off = c_on;
	c_off.y = c_survival.y;

	//****************** JESUS *******************

	TextData jesusTitle_1 = loadTextTexture("BATTLE", "fonts/Impacted2.0.ttf", 70, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData jesusTitle_2 = loadTextTexture("ROYALE", "fonts/Impacted2.0.ttf", 70, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData jesusTitle_3 = loadTextTexture("3D", "fonts/Impacted2.0.ttf", 70, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	way_img = loadImgTexture(way, 1, 1, 1, true, 0xFF, 0xFF, 0xFF);
	waytexture = IMG_LoadTexture(renderer, way);
	if (SDL_SetTextureBlendMode(waytexture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}

	ponpon_1 = loadImgTexture(ponpon_Path1, 1, 1, 1, false, 12, 12, 12);
	ponpon_2 = loadImgTexture(ponpon_Path2, 1, 1, 1, false, 12, 12, 12);
	ponpon_1_texture = IMG_LoadTexture(renderer, ponpon_Path1);
	ponpon_2_texture = IMG_LoadTexture(renderer, ponpon_Path2);
	if (SDL_SetTextureBlendMode(ponpon_1_texture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}
	if (SDL_SetTextureBlendMode(ponpon_2_texture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}
	int frame = -1;

	//****************** FIRST *******************

	wallpaper = loadImgTexture(wallimg, 1, 1, 1, false, 0xFF, 0xFF, 0xFF);
	walltexture = IMG_LoadTexture(renderer, wallimg);

	//enemypaper = loadImgTexture(enemyimg, 1, 1, 1, true, 0xFF, 0xFF, 0xFF);
	//enemytexture = IMG_LoadTexture(renderer, enemyimg);
	/*if (SDL_SetTextureBlendMode(enemytexture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}*/
	if (walltexture == NULL)
	{
		printf("SDL_CreateTextureFromSurface failed: %s\n", SDL_GetError());
	}

	skeletonpaper = loadImgTexture(skeletonimg, 1, 1, 1, true, 0xFF, 0xFF, 0xFF);
	skeletontexture = IMG_LoadTexture(renderer, skeletonimg);
	if (SDL_SetTextureBlendMode(skeletontexture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}
	skeleton_redpaper = loadImgTexture(skeleton_redimg, 1, 1, 1, true, 0xFF, 0xFF, 0xFF);
	skeleton_redtexture = IMG_LoadTexture(renderer, skeleton_redimg);
	if (SDL_SetTextureBlendMode(skeleton_redtexture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}

	portpaper = loadImgTexture(portimg, 1, 1, 1, false, 0xFF, 0xFF, 0xFF);
	porttexture = IMG_LoadTexture(renderer, portimg);
	if (SDL_SetTextureBlendMode(porttexture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}

	SDL_Rect leftBar;
	leftBar.x = 0;
	leftBar.y = 0;
	leftBar.w = 200;
	leftBar.h = HEIGHT;

	SDL_Rect firstView;
	firstView.x = 200;
	firstView.y = 0;
	firstView.w = FirstW;
	leftBar.h = FirstH;

	TextData f_title = loadTextTexture("BATTLE ROYALE 3D", "fonts/Impacted2.0.ttf", 23, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData f_position = loadTextTexture("Position", "fonts/Impacted2.0.ttf", 46, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData f_time = loadTextTexture("Time", "fonts/Impacted2.0.ttf", 50, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData f_portkey = loadTextTexture("Port Key", "fonts/Impacted2.0.ttf", 37, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData f_press = loadTextTexture("[ or press ' p ' ]", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData showtime = loadTextTexture("[ or press ' p ' ]", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);


	//enemy enemyINFO[6];

	//****************** TRANSFORM *******************

	TextData loading = loadTextTexture("LOADING", "fonts/PressStart2P.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData dot = loadTextTexture(".", "fonts/PressStart2P.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	clock_t timesUP;

	//****************** SCRIPT *******************

	text script_background = { WIDTH / 2 - 700 / 2, HEIGHT - 80 - 600 - loading.height, 700, 600 };
	TextData instruction = loadTextTexture("INSTRUCTION", "fonts/Impacted2.0.ttf", 50, 58, 72, 50, BLENDED, NULL, NULL, NULL);
	TextData god_pers = loadTextTexture("GOD'S PERSPECTIVE", "fonts/Impacted2.0.ttf", 30, 127, 96, 0, BLENDED, NULL, NULL, NULL);
	TextData god_pers_1 = loadTextTexture("WHEEL UP: ZOOM IN", "fonts/Impacted2.0.ttf", 25, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData god_pers_2 = loadTextTexture("WHEEL DOWN: ZOOM OUT", "fonts/Impacted2.0.ttf", 25, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	TextData man_pers = loadTextTexture("MAN'S PERSPECTIVE", "fonts/Impacted2.0.ttf", 30, 127, 96, 0, BLENDED, NULL, NULL, NULL);
	TextData man_pers_1 = loadTextTexture("CLICK MOUSE[RIGHT]: TURN BACK", "fonts/Impacted2.0.ttf", 25, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	TextData switch_pers_1 = loadTextTexture("SWITCH", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData switch_pers_2 = loadTextTexture("PERSPECTIVE", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	TextData turn_1 = loadTextTexture("TURN TO", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData turn_2 = loadTextTexture("CURRENT", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData turn_3 = loadTextTexture("PHASE", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	TextData pause = loadTextTexture("PAUSE", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	TextData forward = loadTextTexture("FORWARD", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData backward = loadTextTexture("BACKWARD", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData left = loadTextTexture("LEFT", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData right = loadTextTexture("RIGHT", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	TextData  tab = loadTextTexture("TAB", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  O = loadTextTexture("O", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  P = loadTextTexture("P", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  W = loadTextTexture("W", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  A = loadTextTexture("A", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  S = loadTextTexture("S", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  D = loadTextTexture("D", "fonts/Impacted2.0.ttf", 40, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	//CIRCLE PART
	TextData  operation = loadTextTexture("MOVE YOUR MOUSE TO ROTATE", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  operation_1 = loadTextTexture("1. ROTATE TOWARD UPPER LEFT", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  operation_2 = loadTextTexture("2. NO ROTATION", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  operation_3 = loadTextTexture("3. ROTATE TOWARD RIGHT", "fonts/Impacted2.0.ttf", 20, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  center = loadTextTexture("CENTER", "fonts/Impacted2.0.ttf", 25, 0, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData  B_center = loadTextTexture("CENTER", "fonts/Impacted2.0.ttf", 27, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  center_1 = loadTextTexture("1", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  center_2 = loadTextTexture("2", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData  center_3 = loadTextTexture("3", "fonts/Impacted2.0.ttf", 30, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	/*
	char cur[100] = "images/cursor.png";
	ImageData cursor;
	cursor = loadImgTexture(cur, 1, 1, 1, 0, 0xFF, 0xFF, 0xFF);
	cursortexture = IMG_LoadTexture(renderer, cur);
	if (SDL_SetTextureBlendMode(cursortexture, SDL_BLENDMODE_BLEND) == -1)
	{
		printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
		return -1;
	}
	if (cursortexture == NULL)
	{
		printf("SDL_CreateTextureFromSurface failed: %s\n", SDL_GetError());
	}
	*/

	//****************** PAUSE  SURFACE *******************

	TextData  CONTINUE = loadTextTexture("CONTINUE", "fonts/Impacted2.0.ttf", 40, 58, 72, 50, BLENDED, NULL, NULL, NULL);
	TextData  H_CONTINUE = loadTextTexture("CONTINUE", "fonts/Impacted2.0.ttf", 40, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData  BACK = loadTextTexture("Q U I T", "fonts/Impacted2.0.ttf", 40, 58, 72, 50, BLENDED, NULL, NULL, NULL);
	TextData  H_BACK = loadTextTexture("Q U I T", "fonts/Impacted2.0.ttf", 40, 255, 0, 0, BLENDED, NULL, NULL, NULL);

	t_continue = { (WIDTH + 200) / 2 - 700 / 2 , HEIGHT - 80 - loading.height - 30, 700 / 2, loading.height + 40 + 30 };
	t_back = { (WIDTH + 200) / 2  , HEIGHT - 80 - loading.height - 30, 700 / 2, loading.height + 40 + 30 };

	//****************** FIRST & END *******************

	TextData hash = loadTextTexture("#", "fonts/PressStart2P.ttf", 12, 230, 230, 230, BLENDED, NULL, NULL, NULL);
	TextData H_hash = loadTextTexture("#", "fonts/PressStart2P.ttf", 12, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData AND = loadTextTexture("&", "fonts/PressStart2P.ttf", 12, 150, 150, 150, BLENDED, NULL, NULL, NULL);
	TextData zero = loadTextTexture("0", "fonts/PressStart2P.ttf", 12, 40, 40, 40, BLENDED, NULL, NULL, NULL);

	t_hash.x = (WIDTH - (game_over_SizeW * hash.width)) / 2;
	t_hash.y = 403;

	t_AND.x = (WIDTH - (game_over_SizeW * hash.width)) / 2;
	t_AND.y = 403;

	t_zero.x = (WIDTH - (game_over_SizeW * hash.width)) / 2;
	t_zero.y = 403;

	TextData BADGE_1 = loadTextTexture("BATTLE", "fonts/Impacted2.0.ttf", 45, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData BADGE_2 = loadTextTexture("ROYALE", "fonts/Impacted2.0.ttf", 45, 255, 255, 255, BLENDED, NULL, NULL, NULL);
	TextData BADGE_3 = loadTextTexture("3D", "fonts/Impacted2.0.ttf", 45, 255, 255, 255, BLENDED, NULL, NULL, NULL);

	Sint16 BADGE_X[5] = { t_hash.x + 30, t_hash.x + 30, t_hash.x + (184) / 2 + 30, t_hash.x + 184 + 30,t_hash.x + 184 + 30 };
	Sint16 BADGE_Y[5] = { 58 + 40, 58 + 200 + 40, 58 + 267 + 40 , 58 + 200 + 40,58 + 40 };
	TextData SUCCESS_1 = loadTextTexture("SUCCESSFULLY", "fonts/Impacted2.0.ttf", 70, 255, 255, 111, BLENDED, NULL, NULL, NULL);
	TextData SUCCESS_2 = loadTextTexture("ESCAPED !", "fonts/Impacted2.0.ttf", 70, 255, 255, 111, BLENDED, NULL, NULL, NULL);
	TextData FAIL_1 = loadTextTexture("F A I L   T O", "fonts/Impacted2.0.ttf", 70, 255, 0, 0, BLENDED, NULL, NULL, NULL);
	TextData FAIL_2 = loadTextTexture("E  S  C  A  P  E", "fonts/Impacted2.0.ttf", 70, 255, 0, 0, BLENDED, NULL, NULL, NULL);

	c_TryAgain.x = (t_AND.x + 12);
	c_TryAgain.y = (t_AND.y + 12 * (game_over_SizeH / 2 + 1));
	c_TryAgain.w = (t_AND.x + (game_over_SizeW - 1) * 12) - c_TryAgain.x;
	c_TryAgain.h = (t_AND.y + 12 * (game_over_SizeH - 1)) - c_TryAgain.y;

	c_BackToMenu.x = (t_AND.x + 12);
	c_BackToMenu.y = (t_AND.y + 12 * (game_over_SizeH / 2 + 1));
	c_BackToMenu.w = (t_AND.x + (game_over_SizeW - 1) * 12) - c_BackToMenu.x;
	c_BackToMenu.h = (t_AND.y + 12 * (game_over_SizeH - 1)) - c_BackToMenu.y;

	SDL_Rect junksrc, junkdst;
	mouseState = NONE;
	loadAudio();
	srand(time(NULL));

	//****************** BETA *******************
	TextData INVERSION = loadTextTexture("INVERSION", "fonts/Impacted2.0.ttf", 25, 150, 150, 150, BLENDED, NULL, NULL, NULL);
	TextData H_INVERSION = loadTextTexture("INVERSION", "fonts/Impacted2.0.ttf", 25, 150, 0, 0, BLENDED, NULL, NULL, NULL);
	//c_INVERSION = {WIDTH/2-310, 550, INVERSION.width, INVERSION.height};
	c_INVERSION = { WIDTH - 200, 100, INVERSION.width, INVERSION.height };



	//**************** INITIALISZE ALL OF THE FLAGS ***************
	//*** ENTRANCE ***
	entFlag = true;
	entbgm = true;
	c_on.click = true;
	c_off.click = false;


	//*** MENU ***
	menuFlag = false;
	menubgm = false;

	//*** TRANSFORM ***
	transFlag = false;

	//*** FIRST ***
	firstFlag = false;
	firstDisplay = false;
	pauseSurface = false;
	H_pause = false;
	pauseSurface = false;
	approach = false;
	quickFlag = false;
	reveal = false;

	//*** JESUS ***
	jesusFlag = false;
	jesusDisplay = false;

	//for testing
	/*
	firstFlag = true;
	entFlag = false;
	firstDisplay = true;
	*/

	while (!quit)
	{
		if (firstFlag) t = clock();
		while (SDL_PollEvent(&e) != 0)
		{
			switch (e.type)
			{
			case SDL_QUIT: // User requests quit
				quit = true;
				break;
			}
			mouseHandleEvent(&e, &mouseState, &mouseX, &mouseY);
			handleEvent(e);

			//*********** BGM ***********
			if (c_on.click) {
				if (entFlag && entbgm) {
					Mix_PlayChannel(1, entMusic, -1);
					entbgm = false;
				}
				else if (menuFlag && menubgm) {
					Mix_PlayChannel(2, menuMusic, -1);
					menubgm = false;
				}
				else if (failbgm || successbgm) {
					/*if (c_survival.click) {
						Mix_HaltChannel(5);
					}
					else Mix_HaltChannel(11);*/
					if (failbgm && fail) {
						Mix_HaltChannel(-1);
						Mix_PlayChannel(4, failMusic, -1);
						failbgm = false;
					}
					else if (successbgm && succeed) {
						Mix_HaltChannel(-1);
						Mix_PlayChannel(7, successMusic, -1);
						successbgm = false;
					}
				}
				else if (firstFlag) {
					//printf("in\n\n\n\n");
					if (survivalbgm || timerbgm) {
						if (survivalbgm) {
							Mix_PlayChannel(5, survivalMusic, -1);
							//survivalflag = true;
							//c_TIME.click = false;
							survivalbgm = false;
						}
						else if (timerbgm) {
							Mix_PlayChannel(11, timerMusic, -1);
							//timerflag = true, survivalflag = false;
							timerbgm = false;
						}

					}
					if (alertbgm) {
						//printf("alert\n");
						/*if (c_survival.click) {
							printf("volume change\n");
							//Mix_Volume(5, (float)MIX_MAX_VOLUME * 0.5);
							//Mix_PlayChannel(5, survivalMusic, -1);
						}
						else if (c_TIME.click) {
							printf("volume change\n");
							//Mix_Volume(11, (float)MIX_MAX_VOLUME * 0.5);
							//Mix_PlayChannel(11, timerMusic, -1);
						}*/
						Mix_Volume(8, MIX_MAX_VOLUME);
						Mix_PlayChannel(8, alertMusic, -1);
						alertbgm = false;
					}
					else if (!approach) {
						if (c_survival.click)
							Mix_Volume(5, MIX_MAX_VOLUME);
						else
							Mix_Volume(11, MIX_MAX_VOLUME);
						Mix_FadeOutChannel(8, 1000);
					}
					/*if (fail)printf("fail\n");
					if (failbgm)printf("failbgm\n");
					if (successbgm && succeed) {
						printf("succed!!!!!!!!!\n");
						Mix_HaltChannel(-1);
						Mix_PlayChannel(7, successMusic, -1);
						successbgm = false;
					}
					if (failbgm && fail) {
						printf("fail!!!!!!!!!\n");
						Mix_HaltChannel(-1);
						Mix_PlayChannel(4, failMusic, -1);
						failbgm = false;
					}*/


				}

			}
			else if (c_off.click) {
				Mix_HaltChannel(-1);
				menubgm = true;
			}
		}


		//entFlag = true;
		//firstFlag = false;
		//jesusFlag = false;

		//************ ENTRANCE ************

		if (entFlag) {
			if (initial) {
				SDL_TimerID timerID_Title1 = SDL_AddTimer(10, TitleDrop, dropINFO1);
				SDL_TimerID timerID_Title2 = SDL_AddTimer(10, TitleDrop, dropINFO2);
				SDL_TimerID timerID_Title3 = SDL_AddTimer(10, TitleDrop, dropINFO3);
				initial = false;
			}
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xFF);
			SDL_RenderClear(renderer);
			textRender(renderer, title1, dropINFO1[0], dropINFO1[1], title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, title2, dropINFO2[0], dropINFO2[1], title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, title3, dropINFO3[0], dropINFO3[1], title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, press, WIDTH / 2 - press.width / 2, HEIGHT / 2 + 150, title1.width / 2, title1.height / 2, NULL, no, 255);
			SDL_RenderPresent(renderer);
		}

		//************ MENU ************
		else if (!entFlag && menuFlag) {
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xFF);
			SDL_RenderClear(renderer);

			if (t_HOME.hover) textRender(renderer, H_HOME, t_HOME.x, t_HOME.y, 0, 0, NULL, no, 255);
			else textRender(renderer, HOME, t_HOME.x, t_HOME.y, 0, 0, NULL, no, 255);

			for (int i = 0; i < 6; i++) {
				roundedRectangleRGBA(renderer, WIDTH / 2 - 872 / 2 - i, 182 - i, WIDTH / 2 + 872 / 2 + i, 182 + 438 + i, 10, 255, 255, 255, 255);
			}
			junksrc = { WIDTH / 2 - menuTitle.width / 2 - 10,182 - 6 - menuTitle.height / 2 - 10, menuTitle.width + 20,menuTitle.height + 20 };
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, menuTitle, WIDTH / 2 - menuTitle.width / 2, 182 - 8 - menuTitle.height / 2, title1.width / 2, title1.height / 2, NULL, no, 255);

			roundedRectangleRGBA(renderer, WIDTH / 2 - 182 * 3 / 2 - 80, 438 + 182 - areUready.height - 40 - 290, WIDTH / 2 - 182 * 3 / 2 - 80 + 182, 438 + 182 - areUready.height - 40, 10, 255, 255, 255, 255);
			roundedRectangleRGBA(renderer, WIDTH / 2 - 182 * 3 / 2 - 80 + 6, 438 + 182 - areUready.height - 40 - 290 + 6, WIDTH / 2 - 182 * 3 / 2 - 80 + 182 - 6, 438 + 182 - areUready.height - 40 - 6, 10, 255, 255, 255, 255);

			roundedRectangleRGBA(renderer, WIDTH / 2 - 182 * 1 / 2, 438 + 182 - areUready.height - 40 - 290, WIDTH / 2 + 182 * 1 / 2, 438 + 182 - areUready.height - 40, 10, 255, 255, 255, 255);
			roundedRectangleRGBA(renderer, WIDTH / 2 - 182 * 1 / 2 + 6, 438 + 182 - areUready.height - 40 - 290 + 6, WIDTH / 2 + 182 * 1 / 2 - 6, 438 + 182 - areUready.height - 40 - 6, 10, 255, 255, 255, 255);

			roundedRectangleRGBA(renderer, WIDTH / 2 + 182 * 1 / 2 + 80, 438 + 182 - areUready.height - 40 - 290, WIDTH / 2 + 182 * 1 / 2 + 80 + 182, 438 + 182 - areUready.height - 40, 10, 255, 255, 255, 255);
			roundedRectangleRGBA(renderer, WIDTH / 2 + 182 * 1 / 2 + 80 + 6, 438 + 182 - areUready.height - 40 - 290 + 6, WIDTH / 2 + 182 * 1 / 2 + 80 + 182 - 6, 438 + 182 - areUready.height - 40 - 6, 10, 255, 255, 255, 255);

			junksrc = { WIDTH / 2 - 182 * 3 / 2 - 80 + 182 / 2 - MODE.width / 2 - 10, 438 + 182 - areUready.height - 40 - 290 - WORLD.height / 2 - 10 , MODE.width + 20, WORLD.height + 20 };
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			junksrc = { WIDTH / 2 - WORLD.width / 2 - 10, 438 + 182 - areUready.height - 40 - 290 - WORLD.height / 2 - 10 , WORLD.width + 20, WORLD.height + 20 };
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			junksrc = { WIDTH / 2 + 182 + 80 - MUSIC.width / 2 - 10, 438 + 182 - areUready.height - 40 - 290 - WORLD.height / 2 - 10 , MUSIC.width + 20, WORLD.height + 20 };
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderFillRect(renderer, &junksrc);

			textRender(renderer, MODE, WIDTH / 2 - 182 * 3 / 2 - 80 + 182 / 2 - MODE.width / 2, 438 + 182 - areUready.height - 40 - 290 - MODE.height / 2, title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, WORLD, WIDTH / 2 - WORLD.width / 2, 438 + 182 - areUready.height - 40 - 290 - MODE.height / 2, title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, MUSIC, WIDTH / 2 + 182 + 80 - MUSIC.width / 2, 438 + 182 - areUready.height - 40 - 290 - MODE.height / 2, title1.width / 2, title1.height / 2, NULL, no, 255);

			if (c_TIME.hover || c_TIME.click) {
				roundedRectangleRGBA(renderer, c_TIME.x, c_TIME.y, c_TIME.x + c_TIME.w, c_TIME.y + c_TIME.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_TIME, c_TIME.x + c_TIME.w / 2 - TIME.width / 2, c_TIME.y + c_TIME.h / 2 - TIME.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			else {
				roundedRectangleRGBA(renderer, c_TIME.x, c_TIME.y, c_TIME.x + c_TIME.w, c_TIME.y + c_TIME.h, 8, 255, 255, 255, 255);
				textRender(renderer, TIME, c_TIME.x + c_TIME.w / 2 - TIME.width / 2, c_TIME.y + c_TIME.h / 2 - TIME.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}

			if (c_survival.hover || c_survival.click) {
				roundedRectangleRGBA(renderer, c_survival.x, c_survival.y, c_survival.x + c_survival.w, c_survival.y + c_survival.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_survival, c_survival.x + c_survival.w / 2 - survival.width / 2, c_survival.y + c_survival.h / 2 - survival.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			else {
				roundedRectangleRGBA(renderer, c_survival.x, c_survival.y, c_survival.x + c_survival.w, c_survival.y + c_survival.h, 8, 255, 255, 255, 255);
				textRender(renderer, survival, c_survival.x + c_survival.w / 2 - survival.width / 2, c_survival.y + c_survival.h / 2 - survival.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}

			if (c_simple.hover || c_simple.click) {
				roundedRectangleRGBA(renderer, c_simple.x, c_simple.y, c_simple.x + c_simple.w, c_simple.y + c_simple.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_simple, c_simple.x + c_simple.w / 2 - simple.width / 2, c_simple.y + c_simple.h / 2 - survival.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			else {
				roundedRectangleRGBA(renderer, c_simple.x, c_simple.y, c_simple.x + c_simple.w, c_simple.y + c_simple.h, 8, 255, 255, 255, 255);
				textRender(renderer, simple, c_simple.x + c_simple.w / 2 - simple.width / 2, c_simple.y + c_simple.h / 2 - survival.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			if (c_classic.hover || c_classic.click) {
				roundedRectangleRGBA(renderer, c_classic.x, c_classic.y, c_classic.x + c_classic.w, c_classic.y + c_classic.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_classic, c_classic.x + c_classic.w / 2 - classic.width / 2, c_classic.y + c_classic.h / 2 - classic.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			else {
				roundedRectangleRGBA(renderer, c_classic.x, c_classic.y, c_classic.x + c_classic.w, c_classic.y + c_classic.h, 8, 255, 255, 255, 255);
				textRender(renderer, classic, c_classic.x + c_classic.w / 2 - classic.width / 2, c_classic.y + c_classic.h / 2 - classic.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}

			if (c_on.hover || c_on.click) {
				roundedRectangleRGBA(renderer, c_on.x, c_on.y, c_on.x + c_on.w, c_on.y + c_on.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_on, c_on.x + c_on.w / 2 - on.width / 2, c_on.y + c_on.h / 2 - on.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			else {
				roundedRectangleRGBA(renderer, c_on.x, c_on.y, c_on.x + c_on.w, c_on.y + c_on.h, 8, 255, 255, 255, 255);
				textRender(renderer, on, c_on.x + c_on.w / 2 - on.width / 2, c_on.y + c_on.h / 2 - on.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}

			if (c_off.hover || c_off.click) {
				roundedRectangleRGBA(renderer, c_off.x, c_off.y, c_off.x + c_off.w, c_off.y + c_off.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_off, c_off.x + c_off.w / 2 - off.width / 2, c_off.y + c_off.h / 2 - off.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			else {
				roundedRectangleRGBA(renderer, c_off.x, c_off.y, c_off.x + c_off.w, c_off.y + c_off.h, 8, 255, 255, 255, 255);
				textRender(renderer, off, c_off.x + c_off.w / 2 - off.width / 2, c_off.y + c_off.h / 2 - off.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			}
			if (c_INVERSION.hover || c_INVERSION.click) {
				//roundedRectangleRGBA(renderer, c_INVERSION.x, c_INVERSION.y, c_INVERSION.x + c_INVERSION.w, c_INVERSION.y + c_INVERSION.h, 8, 255, 0, 0, 255);
				textRender(renderer, H_INVERSION, c_INVERSION.x, c_INVERSION.y, NULL, NULL, NULL, no, 255);
			}
			else {
				//roundedRectangleRGBA(renderer, c_INVERSION.x, c_INVERSION.y, c_INVERSION.x + c_INVERSION.w, c_INVERSION.y + c_INVERSION.h, 8, 255, 255, 255, 255);
				textRender(renderer, INVERSION, c_INVERSION.x, c_INVERSION.y, NULL, NULL, NULL, no, 255);
			}


			junksrc.x = WIDTH / 2 - play.width / 2 - 10;
			junksrc.y = 182 + 438 - play.height / 2 - 10;
			junksrc.w = play.width + 20;
			junksrc.h = play.height + 20;
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			if (t_PLAY.hover) textRender(renderer, H_play, WIDTH / 2 - play.width / 2, 182 + 438 - play.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);
			else textRender(renderer, play, WIDTH / 2 - play.width / 2, 182 + 438 - play.height / 2 - 3, title1.width / 2, title1.height / 2, NULL, no, 255);

			if (t_QUIT.hover) textRender(renderer, H_QUIT, t_QUIT.x, t_QUIT.y, title1.width / 2, title1.height / 2, NULL, no, 255);
			else textRender(renderer, QUIT, t_QUIT.x, t_QUIT.y, title1.width / 2, title1.height / 2, NULL, no, 255);

			textRender(renderer, areUready, WIDTH / 2 + 872 / 2 - areUready.width - 15, 438 + 182 - areUready.height - 20, title1.width / 2, title1.height / 2, NULL, no, 255);
			SDL_RenderPresent(renderer);
		}

		//************ TRANSFORM ************
		else if (transFlag) {
			//record the start time
			if (trans_start == 0) trans_start = clock();

			//loading for 10 secs
			if (clock() - trans_start <= duration) {
				SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xFF);
				SDL_RenderClear(renderer);

				textRender(renderer, loading, WIDTH / 2 - loading.width / 2 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
				if (((clock() - trans_start) / 500) % 3 == 0) textRender(renderer, dot, WIDTH / 2 + loading.width / 2 + 10 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
				else if (((clock() - trans_start) / 500) % 3 == 1) {
					textRender(renderer, dot, WIDTH / 2 + loading.width / 2 + 10 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
					textRender(renderer, dot, WIDTH / 2 + loading.width / 2 + 10 + 20 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
				}
				else if (((clock() - trans_start) / 500) % 3 == 2) {
					textRender(renderer, dot, WIDTH / 2 + loading.width / 2 + 10 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
					textRender(renderer, dot, WIDTH / 2 + loading.width / 2 + 10 + 20 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
					textRender(renderer, dot, WIDTH / 2 + loading.width / 2 + 10 + 20 + 20 - 35, HEIGHT - 40 - loading.height, 0, 0, NULL, no, 255);
				}

				textRender(renderer, instruction, WIDTH / 2 - instruction.width / 2, script_background.y + 10, 0, 0, NULL, no, 255);


				//GOD'S PERSPECTIVE: total block height = script_background.y + 10 + 6 + god_pers.height + god_pers_1.height * 2
				int god_block_h = 6 + god_pers.height + god_pers_1.height * 2;
				textRender(renderer, god_pers, script_background.x + 20, instruction.height + 6 + script_background.y + 10, 0, 0, NULL, no, 255);
				textRender(renderer, god_pers_1, script_background.x + 20 + 20, instruction.height + 6 + script_background.y + 10 + god_pers.height + 3, 0, 0, NULL, no, 255);
				textRender(renderer, god_pers_2, script_background.x + 20 + 20, instruction.height + 6 + script_background.y + 10 + god_pers.height + god_pers_1.height + 6, 0, 0, NULL, no, 255);

				//MAN'S PERSPECTIVE: total block height = script_background.y + 10 + 6 + god_pers.height + god_pers_1.height * 2
				int man_block_h = 6 + man_pers.height + man_pers_1.height;
				textRender(renderer, man_pers, script_background.x + 20, instruction.height + 6 + script_background.y + 10 + god_block_h + 30, 0, 0, NULL, no, 255);
				textRender(renderer, man_pers_1, script_background.x + 20 + 20, instruction.height + 6 + script_background.y + 10 + god_block_h + 30 + man_pers.height + 3, 0, 0, NULL, no, 255);

				text keyboard;
				keyboard.x = script_background.x + 20;
				keyboard.y = 5 + script_background.y + 10 + instruction.height + 6 + god_block_h + man_block_h + 100;
				keyboard.w = 20 + man_pers_1.width;
				keyboard.h = script_background.y + 600 - 10 - keyboard.y;

				junksrc.x = keyboard.x + keyboard.w / 2 - 50 / 2;
				junksrc.y = keyboard.y + keyboard.h / 2 - 50;
				junksrc.w = 50;
				junksrc.h = junksrc.w;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, W, keyboard.x + keyboard.w / 2 - W.width / 2, keyboard.y + keyboard.h / 2 - 50 / 2 - W.height / 2 - 3, 0, 0, NULL, no, 255);
				textRender(renderer, forward, keyboard.x + keyboard.w / 2 - forward.width / 2, keyboard.y + keyboard.h / 2 - 50 - 3 - forward.height, 0, 0, NULL, no, 255);

				junksrc.x = keyboard.x + keyboard.w / 2 - 50 / 2;
				junksrc.y = keyboard.y + keyboard.h / 2 + 5;
				junksrc.w = 50;
				junksrc.h = junksrc.w;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, S, keyboard.x + keyboard.w / 2 - S.width / 2, keyboard.y + keyboard.h / 2 + 50 / 2 - S.height / 2 - 3 + 5, 0, 0, NULL, no, 255);
				textRender(renderer, backward, keyboard.x + keyboard.w / 2 - backward.width / 2, keyboard.y + keyboard.h / 2 + 5 + 50 + 3, 0, 0, NULL, no, 255);

				junksrc.x = keyboard.x + keyboard.w / 2 - backward.width / 2 - 10 - 50;
				junksrc.y = keyboard.y + keyboard.h / 2 + 5;
				junksrc.w = 50;
				junksrc.h = junksrc.w;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, A, keyboard.x + keyboard.w / 2 - backward.width / 2 - 10 - 50 / 2 - A.width / 2, keyboard.y + keyboard.h / 2 + 50 / 2 - S.height / 2 - 3 + 5, 0, 0, NULL, no, 255);
				textRender(renderer, left, keyboard.x + keyboard.w / 2 - 10 - backward.width / 2 - 50 / 2 - left.width / 2, keyboard.y + keyboard.h / 2 + 5 + 50 + 3, 0, 0, NULL, no, 255);

				junksrc.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10;
				junksrc.y = keyboard.y + keyboard.h / 2 + 5;
				junksrc.w = 50;
				junksrc.h = junksrc.w;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, D, keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 - D.width / 2, keyboard.y + keyboard.h / 2 + 50 / 2 - S.height / 2 - 3 + 5, 0, 0, NULL, no, 255);
				textRender(renderer, right, keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 - right.width / 2, keyboard.y + keyboard.h / 2 + 5 + 50 + 3, 0, 0, NULL, no, 255);

				junksrc.x = keyboard.x + keyboard.w / 2 + forward.width / 2 + 50 / 2;
				junksrc.y = keyboard.y + keyboard.h / 2 - 50 - 50 - 3;
				junksrc.w = 50;
				junksrc.h = junksrc.w;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, O, junksrc.x + 50 / 2 - O.width / 2, junksrc.y + 50 / 2 - O.height / 2 - 3, 0, 0, NULL, no, 255);
				textRender(renderer, turn_3, junksrc.x + 50 / 2 - turn_3.width / 2, junksrc.y - turn_3.height, 0, 0, NULL, no, 255);
				textRender(renderer, turn_2, junksrc.x + 50 / 2 - turn_2.width / 2, junksrc.y - turn_3.height - turn_2.height, 0, 0, NULL, no, 255);
				textRender(renderer, turn_1, junksrc.x + 50 / 2 - turn_1.width / 2, junksrc.y - turn_3.height - turn_2.height - turn_1.height, 0, 0, NULL, no, 255);

				junksrc.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50;
				junksrc.y = keyboard.y + keyboard.h / 2 - 50 - 50 - 3;
				junksrc.w = 50;
				junksrc.h = junksrc.w;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, P, junksrc.x + 50 / 2 - P.width / 2, junksrc.y + 50 / 2 - P.height / 2 - 3, 0, 0, NULL, no, 255);
				textRender(renderer, pause, junksrc.x + 50 / 2 - pause.width / 2, junksrc.y - pause.height, 0, 0, NULL, no, 255);

				//TAB
				junksrc.x = keyboard.x + keyboard.w / 2 - forward.width / 2 - 50 / 2 - 50 * 2;
				junksrc.w = 100;
				junksrc.h = 50;

				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, tab, junksrc.x + junksrc.w / 2 - tab.width / 2, junksrc.y + 50 / 2 - tab.height / 2 - 3, 0, 0, NULL, no, 255);
				textRender(renderer, switch_pers_2, junksrc.x + junksrc.w / 2 - switch_pers_2.width / 2, junksrc.y - switch_pers_2.height, 0, 0, NULL, no, 255);
				textRender(renderer, switch_pers_1, junksrc.x + junksrc.w / 2 - switch_pers_1.width / 2, junksrc.y - switch_pers_2.height - switch_pers_1.height, 0, 0, NULL, no, 255);

				//circle part
				text circle;
				circle.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50 + 50 / 2 + pause.width / 2;
				circle.y = -30 + keyboard.y + keyboard.h / 2 - 100 - 3 - pause.height / 2; //MIDDLE
				circle.w = script_background.x + script_background.w - 20 - circle.x;
				circle.h = circle.y * 2;
				/*
				junksrc.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50 + 50 / 2 + pause.width / 2;
				junksrc.y = -30 + keyboard.y + keyboard.h / 2 - 100 - 3 - pause.height / 2; //MIDDLE
				junksrc.w = cursor.width;
				junksrc.h = cursor.height;


				SDL_Rect junkdst;
				junkdst.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50 + 50 / 2 + pause.width / 2 + 30;
				junkdst.y = -30 + keyboard.y + keyboard.h / 2 - 100 - 3 - pause.height / 2 + 30; //MIDDLE
				junkdst.w = cursor.width;
				junkdst.h = cursor.height;
				*/

				filledCircleRGBA(renderer, circle.x + circle.w / 2, circle.y, 130, 58, 72, 50, 255);
				filledCircleRGBA(renderer, circle.x + circle.w / 2, circle.y, 20, 255, 255, 255, 255);
				textRender(renderer, B_center, circle.x + circle.w / 2 - B_center.width / 2 + 1, circle.y - B_center.height / 2, 0, 0, NULL, no, 255);
				textRender(renderer, center, circle.x + circle.w / 2 - center.width / 2, circle.y - center.height / 2, 0, 0, NULL, no, 255);
				//SDL_RenderCopy(renderer, cursortexture, &junksrc, &junkdst);
				//imgRender(renderer, cursor, junksrc.x, junksrc.y);
				//imgRender(renderer, sp, WIDTH / 2 - sp.width / 2, 100);

				textRender(renderer, center_1, circle.x + circle.w / 2 - center.width / 2 - 20, circle.y - center.height / 2 - 40, 0, 0, NULL, no, 255);
				textRender(renderer, center_2, circle.x + circle.w / 2, circle.y - center.height / 2 + 30, 0, 0, NULL, no, 255);
				textRender(renderer, center_3, circle.x + circle.w / 2 + center.width / 2 + 30, circle.y - center_3.height / 2, 0, 0, NULL, no, 255);

				junksrc.x = circle.x + circle.w / 2 - operation.width / 2 - 10;
				junksrc.y = circle.y - 130 - 20 - operation.height - 20;
				junksrc.w = operation.width + 20;
				junksrc.h = operation.height + 20;
				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
				SDL_RenderFillRect(renderer, &junksrc);
				textRender(renderer, operation, junksrc.x + junksrc.w / 2 - operation.width / 2, junksrc.y + junksrc.h / 2 - operation.height / 2 - 3, 0, 0, NULL, no, 255);

				junksrc.x = circle.x + circle.w / 2 - operation_3.width / 2;
				junksrc.y = circle.y + 130 + 10;
				junksrc.w = operation_3.width;
				junksrc.h = operation_3.height;
				textRender(renderer, operation_1, junksrc.x, junksrc.y, 0, 0, NULL, no, 255);
				textRender(renderer, operation_2, junksrc.x, junksrc.y + operation_3.height, 0, 0, NULL, no, 255);
				textRender(renderer, operation_3, junksrc.x, junksrc.y + operation_3.height * 2, 0, 0, NULL, no, 255);


				junksrc = { WIDTH / 2 - 700 / 2, HEIGHT - 80 - 600 - loading.height, 700, 600 };
				SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
				SDL_SetRenderDrawColor(renderer, 255, 255, 255, 20);
				SDL_RenderFillRect(renderer, &junksrc);

			}

			// show first perspective
			else {
				Mix_HaltChannel(9);
				if (c_survival.click) survivalbgm = true;
				else timerbgm = true;
				transFlag = false;
				firstFlag = true;
				trans_start = 0;
				firstDisplay = true;
				//jesusDisplay = true;
				duration = 3000; //change to 3 sec
			}
			SDL_RenderPresent(renderer);
		}

		//************ JESUS ************
		else if ((!menuFlag) && (!firstFlag) && (!transFlag) && !(fail || succeed)) {
			//mouseState = NONE;

			//reset
			fThetax = 0.f;
			fThetay = 0.f;
			fThetaz = 0.f;

			if (jesusDisplay) {
				meshCube.tris.clear();
				CubeGenerate();
				if (GPS[0]) {
					fThetay = -M_PI / 2;
				}
				else if (GPS[2]) {
					fThetax = M_PI / 2;
					fThetaz = M_PI / 2;
				}
				else if (GPS[3]) {
					fThetay = M_PI;
				}
				else if (GPS[4]) {
					fThetax = -M_PI / 2;
					fThetaz = -M_PI / 2;
				}
				else if (GPS[5]) {
					fThetay = M_PI / 2;
				}
				jesusDisplay = false;
			}


			//rotate to recent face
			if (recentFace) {
				//printf("in\n\n\n");
				meshCube = originalCube;
				if (GPS[0]) {
					fThetay = -M_PI / 2;
				}
				else if (GPS[2]) {
					fThetax = M_PI / 2;
					fThetaz = M_PI / 2;
				}
				else if (GPS[3]) {
					fThetay = M_PI;
				}
				else if (GPS[4]) {
					fThetax = -M_PI / 2;
					fThetaz = -M_PI / 2;
				}
				else if (GPS[5]) {
					fThetay = M_PI / 2;
				}
				recentFace = false;
			}

			if (mouseState == HOVER) {
				if (!pauseSurface) {
					double d = sqrtf(powf(mouseX - midMouseX, 2.0) + powf(mouseY - midMouseY, 2.0));
					if ((d) > 50) {
						//fThetaA = fmod(fThetaA - 0.0005 * (d - 300), 360);
						fThetaA = -0.0005 * (d - 50);
					}
					else fThetaA = 0.0;
				}
				else fThetaA = 0.0;
			}

			drawCube();
			vCamera = { 0, 0, -20 };
			if (backgroundDisplay) {
				SDL_TimerID timerID_background = SDL_AddTimer(100, background, NULL);
				backgroundDisplay = false;
			}

			SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderClear(renderer);


			int cR, cG, cB;
			cR = rand() % 256;
			cG = rand() % 256;
			cB = rand() % 256;

			for (int i = 0; i < 100; i++) {
				filledCircleRGBA(renderer, starX[i], starY[i], 1, rand() % 256, rand() % 256, rand() % 256, 255);
			}
			for (int tri = 0; tri < queue.size(); tri++) {
				filledTrigonRGBA(renderer, queue[tri].p[0].x, queue[tri].p[0].y, queue[tri].p[1].x, queue[tri].p[1].y, queue[tri].p[2].x, queue[tri].p[2].y, 255 * queue[tri].shadow, 255 * queue[tri].shadow, 255 * queue[tri].shadow, 255);
				//trigonRGBA(renderer, queue[tri].p[0].x, queue[tri].p[0].y, queue[tri].p[1].x, queue[tri].p[1].y, queue[tri].p[2].x, queue[tri].p[2].y, 0 , 0 , 0 , 255);
			}
			//SDL_RenderPresent(renderer);
			queue.clear();

			//set title
			textRender(renderer, jesusTitle_1, 50, t_HOME.y, 0, 0, NULL, no, 255);
			textRender(renderer, jesusTitle_2, 50, t_HOME.y + jesusTitle_1.height, 0, 0, NULL, no, 255);
			textRender(renderer, jesusTitle_3, 50, t_HOME.y + jesusTitle_1.height * 2, 0, 0, NULL, no, 255);

			//set press
			filledCircleRGBA(renderer, 100, 620, 50, 0, 0, 0, 255);
			if (H_pause) {
				filledTrigonRGBA(renderer, 70, 600, 95, 620, 70, 640, 255, 0, 0, 100);
				thickLineRGBA(renderer, 110, 600, 110, 640, 10, 255, 0, 0, 100);
				thickLineRGBA(renderer, 130, 600, 130, 640, 10, 255, 0, 0, 100);
			}
			else {
				filledTrigonRGBA(renderer, 70, 600, 95, 620, 70, 640, 255, 255, 255, 255);
				thickLineRGBA(renderer, 110, 600, 110, 640, 10, 255, 255, 255, 255);
				thickLineRGBA(renderer, 130, 600, 130, 640, 10, 255, 255, 255, 255);
			}
			t_f_press = { 200 / 2 - f_press.width / 2, 680 };
			if (t_f_press.hover)  textRender(renderer, f_press, t_f_press.x, t_f_press.y, 0, 0, NULL, no, 255);
			else  textRender(renderer, f_press, t_f_press.x, t_f_press.y, 0, 0, NULL, no, 255);



			//set ponpon box
			junksrc.x = 980;
			junksrc.y = 545;
			junksrc.w = 170;
			junksrc.h = 170;
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			//junksrc = { 0, 0, way_img.width, way_img.height };
			//junkdst = { 1010-15,575-15, 140, 140 };
			//SDL_RenderCopy(renderer, waytexture, &junksrc, &junkdst);


			//set ponpon
			frame = (++frame) % 2;
			if (frame) {
				junksrc = { 0, 0, ponpon_1.width, ponpon_1.height };
				junkdst = { 1010 - 15,575 - 15, 140, 140 };
				SDL_RenderCopy(renderer, ponpon_1_texture, &junksrc, &junkdst);
			}
			else {
				junksrc = { 0, 0, ponpon_2.width, ponpon_2.height };
				junkdst = { 1010 - 15,575 - 15, 140, 140 };
				SDL_RenderCopy(renderer, ponpon_2_texture, &junksrc, &junkdst);
			}

		}

		//************ FIRST ************
		else if (firstFlag && !(succeed || fail)) {
			if (firstDisplay) {
				reveal = false;
				quickFlag = false;
				userR = fmod((float)rand(), M_PI * 2.f);
				shiftX = shiftY = 0;
				for (int i = 0; i < 6; i++) GPS[i] = 0;
				//get random world phase 
				int phase = rand() % 6;
				//test
				//phase = 1;
				if (!c_simple.click) {
					map = classicWorld[phase];
				}
				else {
					map = simpleWorld[phase];
				}
				GPS[phase] = 1;
				//initialize enemy position on 6 surface
				//assign user position

				for (int i = 0; i < 6; i++) {
					if (c_classic.click) {
						do {
							enemyX = (float)(rand() % mapSize) + 0.5;
							enemyY = (float)(rand() % mapSize) + 0.5;
						} while (classicWorld[i][(int)enemyY][(int)enemyX] == '#' || enemyX<2.f || enemyX>mapSize - 2.f || enemyY < 2.f || enemyX > mapSize - 2.f);//|| (int)enemyX == (int)userX || (int)enemyY == (int)userY
					}
					else {
						do {
							enemyX = (float)(rand() % mapSize) + 0.5;
							enemyY = (float)(rand() % mapSize) + 0.5;
						} while (simpleWorld[i][(int)enemyY][(int)enemyX] == '#' || enemyX<2.f || enemyX>mapSize - 2.f || enemyY < 2.f || enemyX > mapSize - 2.f);//|| (int)enemyX == (int)userX || (int)enemyY == (int)userY
					}
					enemyINFO[i].x = enemyX;
					enemyINFO[i].y = enemyY;
					enemyINFO[i].img = rand() % 4;
				}
				do {
					userX = (float)(rand() % mapSize);
					userY = (float)(rand() % mapSize);
				} while (map[(int)userY][(int)userX] == '#' || userX < 2.f || userX > mapSize - 2.f || userY < 2.f || userY > mapSize - 2.f || (int)enemyINFO[phase].x == (int)userX || (int)enemyINFO[phase].y == (int)userY);
				
				enemyX = enemyINFO[phase].x;
				enemyY = enemyINFO[phase].y;

				enemypaper = loadImgTexture(enemyimg[enemyINFO[phase].img], 1, 1, 1, true, 0xFF, 0xFF, 0xFF);
				enemytexture = IMG_LoadTexture(renderer, enemyimg[enemyINFO[phase].img]);
				if (SDL_SetTextureBlendMode(enemytexture, SDL_BLENDMODE_BLEND) == -1)
				{
					printf("SDL_SetTextureBlendMode failed: %s\n", SDL_GetError());
					return -1;
				}
				speed = 30.f;
				startTime = clock();
				jesusDisplay = true;
				firstDisplay = false;
				alertflag = true;
			}
			else {
				//record latest time
				if (!pauseSurface) {
					now = clock();
				}

				//the time passed
				int timePass = (now - startTime) / CLOCKS_PER_SEC;//total second

				//fail or success
				if (c_survival.click) {
					float d = sqrt(((int)enemyX - (int)userX) * ((int)enemyX - (int)userX) + ((int)enemyY - (int)userY) * ((int)enemyY - (int)userY));
					//printf("%f\n", d);
					if (d <= 0.5) {
						fail = true;
						//firstFlag = false;
						//firstDisplay = false;
						jesusFlag = false;
						jesusDisplay = false;
						failbgm = true;
						continue;
					}
					//enter key port
					for (int i = 0; i < 2; i++) {
						if ((int)userX == port[i].x && (int)userY == port[i].y && GPS[port[i].phase] == 1 && reveal) {
							succeed = true;
							successbgm = true;
						}
					}
					if (succeed) continue;
					//reset port key position
					if (timePass % 75 == 15) {
						char(*checkmap)[mapSize];
						do {
							port[0].phase = rand() % 6;
							port[1].phase = rand() % 6;
						} while (port[0].phase == port[1].phase);
						//test
						//port[0].phase = 1;
						//port[1].phase = 2;
						for (int i = 0; i < 2; i++) {
							do {
								port[i].x = (rand() % (mapSize - 2)) + 1;
								port[i].y = (rand() % (mapSize - 2)) + 1;
								if (c_simple.click)
									checkmap = simpleWorld[port[i].phase];
								else
									checkmap = classicWorld[port[i].phase];
							} while (checkmap[port[i].y][port[i].x] == '#' || port[i].x == (int)userX || port[i].y == (int)userY);//|| !(port[i].x == 1 || port[i].y == mapSize - 2 || port[i].y == 1 || port[i].y == mapSize - 2)
						}
						reveal = true;
					}
					else if (timePass % 75 == 0) reveal = false;

				}
				else if (c_TIME.click) {
					float d = sqrt(((int)enemyX - (int)userX) * ((int)enemyX - (int)userX) + ((int)enemyY - (int)userY) * ((int)enemyY - (int)userY));
					//printf("%f\n", d);
					if (d <= 0.5) {
						fail = true;
						firstFlag = false;
						firstDisplay = false;
						jesusFlag = false;
						jesusDisplay = false;
						failbgm = true;
						continue;
					}
					if (timePass == timeLimit) {
						succeed = true;
						firstFlag = false;
						firstDisplay = false;
						jesusFlag = false;
						jesusDisplay = false;
						successbgm = true;
						continue;
					}
				}
				else if (c_INVERSION.click) {
					float d = sqrt(((int)enemyX - (int)userX) * ((int)enemyX - (int)userX) + ((int)enemyY - (int)userY) * ((int)enemyY - (int)userY));
					//printf("%f\n", d);
					if (d <= 0.5) {
						succeed = true;
						firstFlag = false;
						firstDisplay = false;
						jesusFlag = false;
						jesusDisplay = false;
						successbgm = true;
						continue;
					}
					if (timePass == timeLimit) {
						fail = true;
						firstFlag = false;
						firstDisplay = false;
						jesusFlag = false;
						jesusDisplay = false;
						failbgm = true;
						continue;
					}
				}

				//reset
				shiftR = 0;
				for (int i = 0; i < sWIDTH; i++) {
					depthBuffer[i] = depth;
					for (int j = 0; j < sHEIGHT; j++) {
						shadowC[j][i] = 0.f;
						shadowG[j][i] = 0.f;
						shadowW[j][i] = 0.f;
						shadowE[j][i] = 0.f;
						shadowP[j][i] = 0.f;
					}
				}

				//if ((userX >= 1 && userX <= mapSize - 2) && (userY >= 1 && userY <= mapSize - 2)) {
				if (!pauseSurface && mouseX > 0 && midMouseX < WIDTH) {
					if ((mouseX - (200 + FirstW / 2)) > 80) {
						shiftR = (speed / 2.f * abs(mouseX - (200 + FirstW / 2)) / (FirstW / 2) * 2) * elapsedTime;
					}
					else if ((mouseX - (200 + FirstW / 2)) < -80) {
						shiftR = -(speed/ 2.f * abs(mouseX - (200 + FirstW / 2)) / (FirstW / 2) * 2) * elapsedTime;
					}
				}

				//enemy BFS
				if (!pauseSurface) {
					if (c_survival.click || c_TIME.click) {
						enemy_move();
						enemy_trace();
					}
					else if (c_INVERSION.click) {
						enemy_check();
						enemy_run();
					}
				}

				t = clock() - t;
				//elapsedTime = ((float)t) / CLOCKS_PER_SEC;
				elapsedTime = 0.001;

				raycast();
				enemycast();
				move();
				int showport = 0;
				for (int i = 0; i < 2; i++)
					if (GPS[port[i].phase] && reveal) {
						portcast(i);
						showport = i + 1;
						break;
					}

				SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
				SDL_RenderClear(renderer);
				SDL_Rect srcrect;

				for (int y = 0; y < sHEIGHT; y++) {
					for (int x = 0; x < sWIDTH; x++) {
						srcrect.x = x * (FirstW / sWIDTH) + firstView.x;
						srcrect.y = y * (HEIGHT / sHEIGHT) + firstView.y;
						srcrect.h = (HEIGHT / sHEIGHT);
						srcrect.w = (FirstW / sWIDTH);

						//ceiling
						if (shadowC[y][x] != 0 && shadowG[y][x] == 0 && shadowW[y][x] == 0) {
							//printf("c\n");
							SDL_SetRenderDrawColor(renderer, 0 * shadowC[y][x], 0 * shadowC[y][x], 0 * shadowC[y][x], 255);
							SDL_RenderFillRect(renderer, &srcrect);
							//SDL_RenderDrawRects(renderer, &srcrect, 1);
							//fillRectangleRGBA(renderer, x * (WIDTH / sWIDTH), y * (HEIGHT / sHEIGHT), (x + 1) * (WIDTH / sWIDTH), (y + 1) * (HEIGHT / sHEIGHT), 30 * shadowW[y][x], 30 * shadowW[y][x], 180 * shadowW[y][x], 180);
							//pixelRGBA(renderer, x, y, 30, 30, 180, 200);
						}

						//wall
						else if (shadowC[y][x] == 0 && shadowW[y][x] != 0 && shadowG[y][x] == 0) {
							//printf("w\n");
							//SDL_RenderCopy(renderer, texture, &imgsrcrect, &imgdstrect);
							//SDL_SetRenderDrawColor(renderer, 44 * shadowW[y][x], 106 * shadowW[y][x], 110 * shadowW[y][x], 255);
							//SDL_RenderFillRect(renderer, &srcrect);
							//SDL_SetRenderDrawColor(renderer, 0, 0, 255, int(255.f*(1-shadowW[y][x])));
							//SDL_RenderFillRect(renderer, &srcrect);
							//printf("%f\n", 255.f * (1 - shadowW[y][x]));
							/*if (SDL_SetTextureAlphaMod(walltexture, shadowW[y][x]) == -1)//255.f*
							{
								printf("SDL_SetTextureAlphaMod failed: %s\n", SDL_GetError());
								return -1;
							}*/
							SDL_RenderCopy(renderer, walltexture, &(wallPos[y][x].src), &(wallPos[y][x].dst));
							//SDL_RenderDrawRects(renderer, &srcrect, 1);
							//rectangleRGBA(renderer, x * (WIDTH / sWIDTH), y * (HEIGHT / sHEIGHT), (x + 1) * (WIDTH / sWIDTH), (y + 1) * (HEIGHT / sHEIGHT), 100 * shadowW[y][x], 100 * shadowW[y][x], 100 * shadowW[y][x], 255);
							//pixelRGBA(renderer, x, y, 100*shadowW[y][x], 100 * shadowW[y][x], 100 * shadowW[y][x], 255);
						}

						//ground
						else if (shadowC[y][x] == 0 && shadowW[y][x] == 0 && shadowG[y][x] != 0) {
							//printf("g\n");
							SDL_SetRenderDrawColor(renderer, 201 * shadowG[y][x], 191 * shadowG[y][x], 145 * shadowG[y][x], 255);
							SDL_RenderFillRect(renderer, &srcrect);
							//rectangleRGBA(renderer, firstView.x + x * (FirstW / sWIDTH), y * (FirstH / sHEIGHT), firstView.x + (x + 1) * (FirstW / sWIDTH), (y + 1) * (FirstH / sHEIGHT), 0 * shadowG[y][x], 0 * shadowG[y][x], 0 * shadowG[y][x], 255);
							//SDL_RenderFillRect(renderer, &srcrect);
						}
						//else if (shadowC[y][x] == 0 && shadowW[y][x] == 0 && shadowG[y][x] == 0 && shadowE[y][x] != 0) {
						if (shadowE[y][x] != 0) {
							//printf("e\n");
							/*
							SDL_SetRenderDrawColor(renderer, 201 * shadowE[y][x], 0 * shadowE[y][x], 0 * shadowE[y][x], 180);
							SDL_RenderFillRect(renderer, &srcrect);
							*/
							//SDL_SetRenderDrawColor(renderer, 201 * shadowE[y][x], 0 * shadowE[y][x], 0 * shadowE[y][x], 180);
							//SDL_RenderFillRect(renderer, &srcrect);
							SDL_RenderCopy(renderer, enemytexture, &(EImgPos[y][x].src), &(EImgPos[y][x].dst));
						}
						//port
						if (shadowP[y][x] != 0 && showport) {
							/*if (SDL_SetTextureAlphaMod(walltexture, shadowW[y][x]) == -1)//255.f*
							{
								printf("SDL_SetTextureAlphaMod failed: %s\n", SDL_GetError());
								return -1;
							}*/
							SDL_RenderCopy(renderer, porttexture, &(PImgPos[y][x].src), &(PImgPos[y][x].dst));
						}
					}
				}

				//*************set bar*****************
				SDL_SetRenderDrawColor(renderer, 58, 72, 50, 0xFF);
				SDL_RenderFillRect(renderer, &leftBar);

				//set title
				SDL_Rect titlePort = { 10, 30, 180, 50 };
				SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
				SDL_RenderFillRect(renderer, &titlePort);

				t_f_title = { 15, 40 };
				textRender(renderer, f_title, t_f_title.x, t_f_title.y, 0, 0, NULL, no, 255);

				//set map
				SDL_Rect posPort = { 10, 100, 180, 240 };
				SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
				SDL_RenderFillRect(renderer, &posPort);

				t_f_pos = { 15, 100 };
				textRender(renderer, f_position, t_f_pos.x, t_f_pos.y, 0, 0, NULL, no, 255);

				SDL_Rect mapViewport = { 20, 170, 160, 160 };
				SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
				SDL_RenderFillRect(renderer, &mapViewport);

				SDL_Rect mapBlock;
				mapBlock.w = (float)mapViewport.w / mapSize;
				mapBlock.h = (float)mapViewport.h / mapSize;
				for (int y = 1; y < mapSize - 1; y++) {
					mapBlock.y = (mapSize - 1 - y) * ((float)mapViewport.w / mapSize) + mapViewport.y;
					for (int x = 1; x < mapSize - 1; x++) {
						mapBlock.x = x * ((float)mapViewport.w / mapSize) + mapViewport.x;
						if (map[y][x] == '#') {
							SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
							SDL_RenderFillRect(renderer, &mapBlock);
						}
						if (showport) {
							if (x == port[showport - 1].x && y == port[showport - 1].y) {
								SDL_SetRenderDrawColor(renderer, 0xff, 0xd3, 0x06, 255);
								SDL_RenderFillRect(renderer, &mapBlock);
							}
						}
						if (x == (int)userX && y == (int)userY) {
							SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
							SDL_RenderFillRect(renderer, &mapBlock);
							//printf("user=%f\n", userR);
							//right
							if ((M_PI * 0.25 <= userR && userR < M_PI * 0.75) || (M_PI * -1.75 <= userR && userR < M_PI * -1.25))
								thickLineRGBA(renderer, mapBlock.x + mapBlock.w, mapBlock.y, mapBlock.x + mapBlock.w, mapBlock.y + mapBlock.h, 5, 0, 0, 255, 255);
							//down
							else if ((M_PI * 0.75 <= userR && userR < M_PI * 1.25) || (M_PI * -1.25 <= userR && userR < M_PI * -0.75))
								thickLineRGBA(renderer, mapBlock.x, mapBlock.y + mapBlock.h, mapBlock.x + mapBlock.w, mapBlock.y + mapBlock.h, 5, 0, 0, 255, 255);
							//left
							else if ((M_PI * 1.25 <= userR && userR < M_PI * 1.75) || (M_PI * -0.75 <= userR && userR < M_PI * -0.25))
								thickLineRGBA(renderer, mapBlock.x, mapBlock.y, mapBlock.x, mapBlock.y + mapBlock.h, 5, 0, 0, 255, 255);
							//up
							else if ((M_PI * 1.75 <= userR || userR < M_PI * 0.25) || (M_PI * -0.25 <= userR || userR < M_PI * -1.75))
								thickLineRGBA(renderer, mapBlock.x, mapBlock.y, mapBlock.x + mapBlock.w, mapBlock.y, 5, 0, 0, 255, 255);
						}
						/*if (x == (int)enemyX && y == (int)enemyY) {
							SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
							SDL_RenderFillRect(renderer, &mapBlock);
						}*/
					}
				}

				//timer
				if (c_TIME.click || c_INVERSION.click) {
					//set timeport
					SDL_Rect timePort = { 10, 380, 180, 160 };
					SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
					SDL_RenderFillRect(renderer, &timePort);

					t_f_time = { 50, 380 };
					textRender(renderer, f_time, 50, 380, 0, 0, NULL, no, 255);

					thickLineRGBA(renderer, 15, 445, 185, 445, 5, 191, 191, 191, 255);

					SDL_Rect countPort[4] = { { 20, 460, 30, 60 }, {55, 460, 30, 60 }, {115, 460, 30, 60}, {150, 460, 30, 60} };
					SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
					for (int i = 0; i < 4; i++)
						SDL_RenderFillRect(renderer, &countPort[i]);
					SDL_Rect clockdot[2] = { {95, 470, 10, 10}, {95, 495, 10, 10} };
					for (int i = 0; i < 2; i++)
						SDL_RenderFillRect(renderer, &clockdot[i]);

					//show time
					int show = 0;
					show = timeLimit - timePass;

					if ((c_TIME.click || c_INVERSION.click) && show < 15) {
						//minute
						char t1[2] = { '0' + show / 600, '\0' };
						//printf("%s\n", t1);
						showtime = loadTextTexture(t1, "fonts/Impacted2.0.ttf", 50, 255, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 21, 455, 0, 0, NULL, no, 255);

						char t2[2] = { '0' + show / 60, '\0' };
						showtime = loadTextTexture(t2, "fonts/Impacted2.0.ttf", 50, 255, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 56, 455, 0, 0, NULL, no, 255);

						//second
						char t3[2] = { '0' + (show % 60) / 10, '\0' };
						//printf("%s\n", t1);
						showtime = loadTextTexture(t3, "fonts/Impacted2.0.ttf", 50, 255, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 116, 455, 0, 0, NULL, no, 255);

						char t4[2] = { '0' + (show % 60) % 10, '\0' };
						showtime = loadTextTexture(t4, "fonts/Impacted2.0.ttf", 50, 255, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 151, 455, 0, 0, NULL, no, 255);
					}
					else {
						//minute
						char t1[2] = { '0' + show / 600, '\0' };
						//printf("%s\n", t1);
						showtime = loadTextTexture(t1, "fonts/Impacted2.0.ttf", 50, 0, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 21, 455, 0, 0, NULL, no, 255);

						char t2[2] = { '0' + show / 60, '\0' };
						showtime = loadTextTexture(t2, "fonts/Impacted2.0.ttf", 50, 0, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 56, 455, 0, 0, NULL, no, 255);

						//second
						char t3[2] = { '0' + (show % 60) / 10, '\0' };
						//printf("%s\n", t1);
						showtime = loadTextTexture(t3, "fonts/Impacted2.0.ttf", 50, 0, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 116, 455, 0, 0, NULL, no, 255);

						char t4[2] = { 48 + (show % 60) % 10, '\0' };
						showtime = loadTextTexture(t4, "fonts/Impacted2.0.ttf", 50, 0, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, showtime, 151, 455, 0, 0, NULL, no, 255);
					}
				}

				//survival
				else {
					//set port key port
					SDL_Rect PKPort = { 10, 355, 180, 205 };
					SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
					SDL_RenderFillRect(renderer, &PKPort);

					textRender(renderer, f_portkey, 28, 350, 0, 0, NULL, no, 255);

					SDL_Rect block = { 20, 400, 160, 150 };
					SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
					SDL_RenderFillRect(renderer, &block);

					block = { 30, 410, 140, 130 };
					SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
					SDL_RenderFillRect(renderer, &block);
					if (reveal) {
						TextData f_conceal = loadTextTexture("Conceal", "fonts/Impacted2.0.ttf", 38, 255, 255, 255, BLENDED, NULL, NULL, NULL);
						TextData f_reveal = loadTextTexture("Reveal", "fonts/Impacted2.0.ttf", 37, 255, 0, 0, BLENDED, NULL, NULL, NULL);
						textRender(renderer, f_reveal, 49, 415, 0, 0, NULL, no, 255);
						textRender(renderer, f_conceal, 36, 480, 0, 0, NULL, no, 255);
					}
					else {
						TextData f_conceal = loadTextTexture("Conceal", "fonts/Impacted2.0.ttf", 38, 255, 0, 0, BLENDED, NULL, NULL, NULL);
						TextData f_reveal = loadTextTexture("Reveal", "fonts/Impacted2.0.ttf", 37, 255, 255, 255, BLENDED, NULL, NULL, NULL);
						textRender(renderer, f_reveal, 49, 415, 0, 0, NULL, no, 255);
						textRender(renderer, f_conceal, 36, 480, 0, 0, NULL, no, 255);

					}
					thickLineRGBA(renderer, 35, 475, 160, 475, 5, 255, 255, 255, 255);

				}

				//set press bottum
				filledCircleRGBA(renderer, 100, 620, 50, 0, 0, 0, 255);
				if (H_pause) {
					filledTrigonRGBA(renderer, 70, 600, 95, 620, 70, 640, 255, 0, 0, 100);
					thickLineRGBA(renderer, 110, 600, 110, 640, 10, 255, 0, 0, 100);
					thickLineRGBA(renderer, 130, 600, 130, 640, 10, 255, 0, 0, 100);
				}
				else {
					filledTrigonRGBA(renderer, 70, 600, 95, 620, 70, 640, 255, 255, 255, 255);
					thickLineRGBA(renderer, 110, 600, 110, 640, 10, 255, 255, 255, 255);
					thickLineRGBA(renderer, 130, 600, 130, 640, 10, 255, 255, 255, 255);
				}
				t_f_press = { 200 / 2 - f_press.width / 2, 680 };
				textRender(renderer, f_press, t_f_press.x, t_f_press.y, 0, 0, NULL, no, 255);

				//*************skeleton**************
				float d = sqrtf((enemyX - userX) * (enemyX - userX) + (enemyY - userY) * (enemyY - userY));

				if (d < 5.f) {
					if (c_survival.click)
						Mix_Volume(5, (float)MIX_MAX_VOLUME * (d / 5));
					else Mix_Volume(11, (float)MIX_MAX_VOLUME * (d / 5));
					junksrc = { 0, 0, skeleton_redpaper.width, skeleton_redpaper.height };
					junkdst = { WIDTH - 120, 10, 120, 120 };
					SDL_RenderCopy(renderer, skeleton_redtexture, &junksrc, &junkdst);
					approach = true;
					if (alertflag) {
						alertbgm = true;
						alertflag = false;
					}
				}
				else {
					junksrc = { 0, 0, skeletonpaper.width, skeletonpaper.height };
					junkdst = { WIDTH - 120, 10, 120, 120 };
					SDL_RenderCopy(renderer, skeletontexture, &junksrc, &junkdst);
					approach = false;
					alertflag = true;
				}
			}
		}

		//************ SUCCESS OR FAIL ************
		else if (succeed || fail) {
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 0xFF);
			SDL_RenderClear(renderer);
			for (int i = 0; i < 4; i++) {
				roundedRectangleRGBA(renderer, 64 - i, 58 - i, WIDTH - 64 + i, HEIGHT - 58 + i, 8, 255, 255, 255, 255);
			}

			//bottom
			for (int i = 0; i < game_over_SizeH; i++) {
				for (int j = 0; j < game_over_SizeW; j++) {
					if (game_over[i][j] == '&') {
						junksrc.x = t_AND.x + j * AND.width;
						junksrc.y = t_AND.y + i * AND.height;
						junksrc.w = AND.width;
						junksrc.h = AND.height;
						SDL_SetRenderDrawColor(renderer, 192, 192, 192, 200);
						SDL_RenderFillRect(renderer, &junksrc);

						textRender(renderer, AND, t_AND.x + j * AND.width, t_AND.y + i * AND.height, 0, 0, 0, no, 255);
					}
					else if (game_over[i][j] == '#') {
						junksrc.x = t_AND.x + j * AND.width;
						junksrc.y = t_AND.y + i * AND.height;
						junksrc.w = AND.width;
						junksrc.h = AND.height;

						//if(c_TryAgain.hover&&i>1&&i<game_over_SizeH/2-1) SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
						//else if(c_BackToMenu.hover&&i>game_over_SizeH/2+1&&i<game_over_SizeH-1) SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
						SDL_SetRenderDrawColor(renderer, 80, 80, 80, 200);
						SDL_RenderFillRect(renderer, &junksrc);

						if ((c_TryAgain.hover && i > 1 && i < game_over_SizeH / 2 - 1) || (c_BackToMenu.hover && i > game_over_SizeH / 2 + 1 && i < game_over_SizeH - 1)) textRender(renderer, H_hash, t_hash.x + j * hash.width, t_hash.y + i * hash.height, 0, 0, 0, no, 255);
						else textRender(renderer, hash, t_hash.x + j * hash.width, t_hash.y + i * hash.height, 0, 0, 0, no, 255);

					}
					else if (game_over[i][j] == '0') {
						textRender(renderer, zero, t_zero.x + j * zero.width, t_zero.y + i * zero.height, 0, 0, 0, no, 255);
					}
				}
			}

			filledPolygonRGBA(renderer, BADGE_X, BADGE_Y, 5, 58, 72, 50, 155);
			textRender(renderer, BADGE_1, BADGE_X[0] + 184 / 2 - BADGE_1.width / 2, BADGE_Y[0] + 100 - BADGE_2.height / 2 - 5 - BADGE_2.height, title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, BADGE_2, BADGE_X[0] + 184 / 2 - BADGE_2.width / 2, BADGE_Y[0] + 100 - BADGE_2.height / 2, title1.width / 2, title1.height / 2, NULL, no, 255);
			textRender(renderer, BADGE_3, BADGE_X[0] + 184 / 2 - BADGE_3.width / 2, BADGE_Y[0] + 100 + BADGE_2.height / 2 + 5, title1.width / 2, title1.height / 2, NULL, no, 255);

			//success
			if (succeed) {
				textRender(renderer, SUCCESS_1, ((WIDTH - t_hash.x - BADGE_X[3] - 60 * 2) / 2 + BADGE_X[3] + 60 - SUCCESS_1.width / 2), BADGE_Y[3] - 70 - 10 - SUCCESS_1.height, title1.width / 2, title1.height / 2, NULL, no, 255);
				textRender(renderer, SUCCESS_2, ((WIDTH - t_hash.x - BADGE_X[3] - 60 * 2) / 2 + BADGE_X[3] + 60 - SUCCESS_2.width / 2), BADGE_Y[3] - 70 + 10, title1.width / 2, title1.height / 2, NULL, no, 255);

			}

			//fail
			else {
				textRender(renderer, FAIL_1, ((WIDTH - t_hash.x - BADGE_X[3] - 60 * 2) / 2 + BADGE_X[3] + 60 - FAIL_1.width / 2), BADGE_Y[3] - 70 - 10 - SUCCESS_1.height, title1.width / 2, title1.height / 2, NULL, no, 255);
				textRender(renderer, FAIL_2, ((WIDTH - t_hash.x - BADGE_X[3] - 60 * 2) / 2 + BADGE_X[3] + 60 - FAIL_2.width / 2), BADGE_Y[3] - 70 + 10, title1.width / 2, title1.height / 2, NULL, no, 255);

			}
		}

		//************ PAUSE SURFACE ************
		if (pauseSurface && (firstFlag || jesusFlag)) {
			junksrc = { (WIDTH + 200) / 2 - 700 / 2, HEIGHT - 80 - 600 - loading.height, 700, 600 };
			SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
			SDL_SetRenderDrawColor(renderer, 20, 20, 20, 200);
			SDL_RenderFillRect(renderer, &junksrc);

			textRender(renderer, instruction, WIDTH / 2 + 100 - instruction.width / 2, script_background.y + 10, 0, 0, NULL, no, 255);

			//GOD'S PERSPECTIVE: total block height = script_background.y + 10 + 6 + god_pers.height + god_pers_1.height * 2
			int god_block_h = 6 + god_pers.height + god_pers_1.height * 2;
			textRender(renderer, god_pers, 100 + script_background.x + 20, instruction.height + 6 + script_background.y + 10, 0, 0, NULL, no, 255);
			textRender(renderer, god_pers_1, 100 + script_background.x + 20 + 20, instruction.height + 6 + script_background.y + 10 + god_pers.height + 3, 0, 0, NULL, no, 255);
			textRender(renderer, god_pers_2, 100 + script_background.x + 20 + 20, instruction.height + 6 + script_background.y + 10 + god_pers.height + god_pers_1.height + 6, 0, 0, NULL, no, 255);

			//MAN'S PERSPECTIVE: total block height = script_background.y + 10 + 6 + god_pers.height + god_pers_1.height * 2
			int man_block_h = 6 + man_pers.height + man_pers_1.height;
			textRender(renderer, man_pers, 100 + script_background.x + 20, instruction.height + 6 + script_background.y + 10 + god_block_h + 30, 0, 0, NULL, no, 255);
			textRender(renderer, man_pers_1, 100 + script_background.x + 20 + 20, instruction.height + 6 + script_background.y + 10 + god_block_h + 30 + man_pers.height + 3, 0, 0, NULL, no, 255);

			text keyboard;
			keyboard.x = 100 + script_background.x + 20;
			keyboard.y = 5 + script_background.y + 10 + instruction.height + 6 + god_block_h + man_block_h + 100;
			keyboard.w = 20 + man_pers_1.width;
			keyboard.h = script_background.y + 600 - 10 - keyboard.y;

			junksrc.x = keyboard.x + keyboard.w / 2 - 50 / 2;
			junksrc.y = keyboard.y + keyboard.h / 2 - 50;
			junksrc.w = 50;
			junksrc.h = junksrc.w;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, W, keyboard.x + keyboard.w / 2 - W.width / 2, keyboard.y + keyboard.h / 2 - 50 / 2 - W.height / 2 - 3, 0, 0, NULL, no, 255);
			textRender(renderer, forward, keyboard.x + keyboard.w / 2 - forward.width / 2, keyboard.y + keyboard.h / 2 - 50 - 3 - forward.height, 0, 0, NULL, no, 255);

			junksrc.x = keyboard.x + keyboard.w / 2 - 50 / 2;
			junksrc.y = keyboard.y + keyboard.h / 2 + 5;
			junksrc.w = 50;
			junksrc.h = junksrc.w;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, S, keyboard.x + keyboard.w / 2 - S.width / 2, keyboard.y + keyboard.h / 2 + 50 / 2 - S.height / 2 - 3 + 5, 0, 0, NULL, no, 255);
			textRender(renderer, backward, keyboard.x + keyboard.w / 2 - backward.width / 2, keyboard.y + keyboard.h / 2 + 5 + 50 + 3, 0, 0, NULL, no, 255);

			junksrc.x = keyboard.x + keyboard.w / 2 - backward.width / 2 - 10 - 50;
			junksrc.y = keyboard.y + keyboard.h / 2 + 5;
			junksrc.w = 50;
			junksrc.h = junksrc.w;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, A, keyboard.x + keyboard.w / 2 - backward.width / 2 - 10 - 50 / 2 - A.width / 2, keyboard.y + keyboard.h / 2 + 50 / 2 - S.height / 2 - 3 + 5, 0, 0, NULL, no, 255);
			textRender(renderer, left, keyboard.x + keyboard.w / 2 - 10 - backward.width / 2 - 50 / 2 - left.width / 2, keyboard.y + keyboard.h / 2 + 5 + 50 + 3, 0, 0, NULL, no, 255);

			junksrc.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10;
			junksrc.y = keyboard.y + keyboard.h / 2 + 5;
			junksrc.w = 50;
			junksrc.h = junksrc.w;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, D, keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 - D.width / 2, keyboard.y + keyboard.h / 2 + 50 / 2 - S.height / 2 - 3 + 5, 0, 0, NULL, no, 255);
			textRender(renderer, right, keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 - right.width / 2, keyboard.y + keyboard.h / 2 + 5 + 50 + 3, 0, 0, NULL, no, 255);

			junksrc.x = keyboard.x + keyboard.w / 2 + forward.width / 2 + 50 / 2;
			junksrc.y = keyboard.y + keyboard.h / 2 - 50 - 50 - 3;
			junksrc.w = 50;
			junksrc.h = junksrc.w;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, O, junksrc.x + 50 / 2 - O.width / 2, junksrc.y + 50 / 2 - O.height / 2 - 3, 0, 0, NULL, no, 255);
			textRender(renderer, turn_3, junksrc.x + 50 / 2 - turn_3.width / 2, junksrc.y - turn_3.height, 0, 0, NULL, no, 255);
			textRender(renderer, turn_2, junksrc.x + 50 / 2 - turn_2.width / 2, junksrc.y - turn_3.height - turn_2.height, 0, 0, NULL, no, 255);
			textRender(renderer, turn_1, junksrc.x + 50 / 2 - turn_1.width / 2, junksrc.y - turn_3.height - turn_2.height - turn_1.height, 0, 0, NULL, no, 255);

			junksrc.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50;
			junksrc.y = keyboard.y + keyboard.h / 2 - 50 - 50 - 3;
			junksrc.w = 50;
			junksrc.h = junksrc.w;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, P, junksrc.x + 50 / 2 - P.width / 2, junksrc.y + 50 / 2 - P.height / 2 - 3, 0, 0, NULL, no, 255);
			textRender(renderer, pause, junksrc.x + 50 / 2 - pause.width / 2, junksrc.y - pause.height, 0, 0, NULL, no, 255);

			//TAB
			junksrc.x = keyboard.x + keyboard.w / 2 - forward.width / 2 - 50 / 2 - 50 * 2;
			junksrc.w = 100;
			junksrc.h = 50;

			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, tab, junksrc.x + junksrc.w / 2 - tab.width / 2, junksrc.y + 50 / 2 - tab.height / 2 - 3, 0, 0, NULL, no, 255);
			textRender(renderer, switch_pers_2, junksrc.x + junksrc.w / 2 - switch_pers_2.width / 2, junksrc.y - switch_pers_2.height, 0, 0, NULL, no, 255);
			textRender(renderer, switch_pers_1, junksrc.x + junksrc.w / 2 - switch_pers_1.width / 2, junksrc.y - switch_pers_2.height - switch_pers_1.height, 0, 0, NULL, no, 255);

			//circle part
			text circle;
			circle.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50 + 50 / 2 + pause.width / 2;
			circle.y = -30 + keyboard.y + keyboard.h / 2 - 100 - 3 - pause.height / 2; //MIDDLE
			circle.w = 100 + script_background.x + script_background.w - 20 - circle.x;
			circle.h = circle.y * 2;
			/*
			junksrc.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50 + 50 / 2 + pause.width / 2;
			junksrc.y = -30 + keyboard.y + keyboard.h / 2 - 100 - 3 - pause.height / 2; //MIDDLE
			junksrc.w = cursor.width;
			junksrc.h = cursor.height;


			SDL_Rect junkdst;
			junkdst.x = keyboard.x + keyboard.w / 2 + backward.width / 2 + 10 + 50 / 2 + 3 + 50 + 50 / 2 + pause.width / 2 + 30;
			junkdst.y = -30 + keyboard.y + keyboard.h / 2 - 100 - 3 - pause.height / 2 + 30; //MIDDLE
			junkdst.w = cursor.width;
			junkdst.h = cursor.height;
			*/

			filledCircleRGBA(renderer, circle.x + circle.w / 2, circle.y, 130, 58, 72, 50, 255);
			filledCircleRGBA(renderer, circle.x + circle.w / 2, circle.y, 20, 255, 255, 255, 255);
			textRender(renderer, B_center, circle.x + circle.w / 2 - B_center.width / 2 + 1, circle.y - B_center.height / 2, 0, 0, NULL, no, 255);
			textRender(renderer, center, circle.x + circle.w / 2 - center.width / 2, circle.y - center.height / 2, 0, 0, NULL, no, 255);
			//SDL_RenderCopy(renderer, cursortexture, &junksrc, &junkdst);
			//imgRender(renderer, cursor, junksrc.x, junksrc.y);
			//imgRender(renderer, sp, WIDTH / 2 - sp.width / 2, 100);

			textRender(renderer, center_1, circle.x + circle.w / 2 - center.width / 2 - 20, circle.y - center.height / 2 - 40, 0, 0, NULL, no, 255);
			textRender(renderer, center_2, circle.x + circle.w / 2, circle.y - center.height / 2 + 30, 0, 0, NULL, no, 255);
			textRender(renderer, center_3, circle.x + circle.w / 2 + center.width / 2 + 30, circle.y - center_3.height / 2, 0, 0, NULL, no, 255);

			junksrc.x = circle.x + circle.w / 2 - operation.width / 2 - 10;
			junksrc.y = circle.y - 130 - 20 - operation.height - 20;
			junksrc.w = operation.width + 20;
			junksrc.h = operation.height + 20;
			SDL_SetRenderDrawColor(renderer, 58, 72, 50, 255);
			SDL_RenderFillRect(renderer, &junksrc);
			textRender(renderer, operation, junksrc.x + junksrc.w / 2 - operation.width / 2, junksrc.y + junksrc.h / 2 - operation.height / 2 - 3, 0, 0, NULL, no, 255);

			junksrc.x = circle.x + circle.w / 2 - operation_3.width / 2;
			junksrc.y = circle.y + 130 + 10;
			junksrc.w = operation_3.width;
			junksrc.h = operation_3.height;
			textRender(renderer, operation_1, junksrc.x, junksrc.y, 0, 0, NULL, no, 255);
			textRender(renderer, operation_2, junksrc.x, junksrc.y + operation_3.height, 0, 0, NULL, no, 255);
			textRender(renderer, operation_3, junksrc.x, junksrc.y + operation_3.height * 2, 0, 0, NULL, no, 255);

			junksrc = { (WIDTH + 200) / 2 - 700 / 2, HEIGHT - 80 - loading.height, 700, loading.height + 40 };
			SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);
			SDL_SetRenderDrawColor(renderer, 20, 20, 20, 200);
			SDL_RenderFillRect(renderer, &junksrc);

			thickLineRGBA(renderer, junksrc.x, junksrc.y - 30, junksrc.x + junksrc.w, junksrc.y - 30, 4, 58, 72, 50, 255);
			thickLineRGBA(renderer, junksrc.x + junksrc.w / 2, junksrc.y - 30, junksrc.x + junksrc.w / 2, junksrc.y + junksrc.h, 4, 58, 72, 50, 255);
			if (t_continue.hover) textRender(renderer, H_CONTINUE, junksrc.x + 700 / 4 - CONTINUE.width / 2, junksrc.y - CONTINUE.height / 2 + 20 - 3, 0, 0, NULL, no, 255);
			else textRender(renderer, CONTINUE, junksrc.x + 700 / 4 - CONTINUE.width / 2, junksrc.y - CONTINUE.height / 2 + 20 - 3, 0, 0, NULL, no, 255);
			if (t_back.hover) textRender(renderer, H_BACK, junksrc.x + 700 * 3 / 4 - BACK.width / 2, junksrc.y - BACK.height / 2 + 20 - 3, 0, 0, NULL, no, 255);
			else textRender(renderer, BACK, junksrc.x + 700 * 3 / 4 - BACK.width / 2, junksrc.y - BACK.height / 2 + 20 - 3, 0, 0, NULL, no, 255);
		}

		SDL_RenderPresent(renderer);
	}

	//Free resources and close SDL
	//************ ENTRANCE ************
	SDL_DestroyTexture(title1.texture);
	SDL_DestroyTexture(title2.texture);
	SDL_DestroyTexture(title3.texture);
	SDL_DestroyTexture(press.texture);

	//************ MENU ************
	SDL_DestroyTexture(HOME.texture);
	SDL_DestroyTexture(H_HOME.texture);
	SDL_DestroyTexture(menuTitle.texture);
	SDL_DestroyTexture(MODE.texture);
	SDL_DestroyTexture(WORLD.texture);
	SDL_DestroyTexture(MUSIC.texture);
	SDL_DestroyTexture(survival.texture);
	SDL_DestroyTexture(H_survival.texture);
	SDL_DestroyTexture(TIME.texture);
	SDL_DestroyTexture(H_TIME.texture);
	SDL_DestroyTexture(simple.texture);
	SDL_DestroyTexture(H_simple.texture);
	SDL_DestroyTexture(classic.texture);
	SDL_DestroyTexture(H_classic.texture);
	SDL_DestroyTexture(areUready.texture);
	SDL_DestroyTexture(on.texture);
	SDL_DestroyTexture(H_on.texture);
	SDL_DestroyTexture(off.texture);
	SDL_DestroyTexture(H_off.texture);
	SDL_DestroyTexture(play.texture);
	SDL_DestroyTexture(H_play.texture);
	SDL_DestroyTexture(QUIT.texture);
	SDL_DestroyTexture(H_QUIT.texture);

	//************ SUCCESS & FAIL ************
	SDL_DestroyTexture(hash.texture);
	SDL_DestroyTexture(AND.texture);
	SDL_DestroyTexture(zero.texture);
	SDL_DestroyTexture(BADGE_1.texture);
	SDL_DestroyTexture(BADGE_2.texture);
	SDL_DestroyTexture(BADGE_3.texture);
	SDL_DestroyTexture(SUCCESS_1.texture);
	SDL_DestroyTexture(SUCCESS_2.texture);
	SDL_DestroyTexture(FAIL_1.texture);
	SDL_DestroyTexture(FAIL_2.texture);

	//************ PAUSE SURFACE ************
	SDL_DestroyTexture(instruction.texture);
	SDL_DestroyTexture(god_pers.texture);
	SDL_DestroyTexture(god_pers_1.texture);
	SDL_DestroyTexture(god_pers_2.texture);

	SDL_DestroyTexture(man_pers.texture);
	SDL_DestroyTexture(man_pers_1.texture);
	SDL_DestroyTexture(switch_pers_1.texture);
	SDL_DestroyTexture(switch_pers_2.texture);

	SDL_DestroyTexture(turn_1.texture);
	SDL_DestroyTexture(turn_2.texture);
	SDL_DestroyTexture(turn_3.texture);

	SDL_DestroyTexture(pause.texture);

	SDL_DestroyTexture(forward.texture);
	SDL_DestroyTexture(backward.texture);
	SDL_DestroyTexture(left.texture);
	SDL_DestroyTexture(right.texture);

	SDL_DestroyTexture(tab.texture);
	SDL_DestroyTexture(O.texture);
	SDL_DestroyTexture(P.texture);
	SDL_DestroyTexture(W.texture);
	SDL_DestroyTexture(A.texture);
	SDL_DestroyTexture(S.texture);
	SDL_DestroyTexture(D.texture);

	SDL_DestroyTexture(operation.texture);
	SDL_DestroyTexture(operation_1.texture);
	SDL_DestroyTexture(operation_2.texture);
	SDL_DestroyTexture(operation_3.texture);

	SDL_DestroyTexture(center.texture);
	SDL_DestroyTexture(B_center.texture);
	SDL_DestroyTexture(center_1.texture);
	SDL_DestroyTexture(center_2.texture);
	SDL_DestroyTexture(center_3.texture);

	//************ JESUS ************
	SDL_DestroyTexture(jesusTitle_1.texture);
	SDL_DestroyTexture(jesusTitle_2.texture);
	SDL_DestroyTexture(jesusTitle_3.texture);

	//************ FIRST ************
	SDL_DestroyTexture(f_title.texture);
	SDL_DestroyTexture(f_position.texture);
	SDL_DestroyTexture(f_time.texture);
	SDL_DestroyTexture(f_portkey.texture);
	SDL_DestroyTexture(f_press.texture);
	SDL_DestroyTexture(showtime.texture);

	closeSDL();
	//delete[]shadow;
	return 0;
}
