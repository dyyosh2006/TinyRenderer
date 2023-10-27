#include <vector>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red = TGAColor(255, 0, 0, 255);

Model* model = NULL;
const int width = 800;
const int height = 800;

Vec3f		 L(1, 2, 3);
Vec3f cameraPos(1, 0, 3);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);
float pi = 3.1415926;

float* zbuffer;
float* shadowbuffer;

struct Shader : public IShader 
{
	mat<2, 3, float> varying_uv;  // triangle uv coordinates, written by the vertex shader, read by the fragment shader
	mat<4, 3, float> varying_tri; // triangle coordinates (clip coordinates), written by VS, read by FS
	mat<3, 3, float> varying_nrm; // normal per vertex to be interpolated by FS
	mat<3, 3, float> ndc_tri;     // triangle in normalized device coordinates

	virtual Vec4f vertex(int iface, int nthvert) 
	{
		varying_uv.set_col(nthvert, model->uv(iface, nthvert));
		varying_nrm.set_col(nthvert, proj<3>(/*(Projection * ModelView).invert_transpose() **/ embed<4>(model->normal(iface, nthvert), 0.f)));
		Vec4f gl_Vertex = Projection * ModelView * embed<4>(model->vert(iface, nthvert));
		varying_tri.set_col(nthvert, gl_Vertex);
		ndc_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor& color, Vec3f currentScreenPixel)
	{
		Vec3f bn = (varying_nrm * bar).normalize();
		Vec2f uv = varying_uv * bar;

		/*mat<3, 3, float> A;
		A[0] = ndc_tri.col(1) - ndc_tri.col(0);
		A[1] = ndc_tri.col(2) - ndc_tri.col(0);
		A[2] = bn;

		mat<3, 3, float> AI = A.invert();
		Vec3f i = AI * Vec3f(varying_uv[0][1] - varying_uv[0][0], varying_uv[0][2] - varying_uv[0][0], 0);
		Vec3f j = AI * Vec3f(varying_uv[1][1] - varying_uv[1][0], varying_uv[1][2] - varying_uv[1][0], 0);

		mat<3, 3, float> B;
		B.set_col(0, i.normalize());
		B.set_col(1, j.normalize());
		B.set_col(2, bn);
		Vec3f n = (B * model->normal(uv)).normalize();*/

		//tbn
		float u10 = varying_uv[0][1] - varying_uv[0][0];
		float u20 = varying_uv[0][2] - varying_uv[0][0];
		float v10 = varying_uv[1][1] - varying_uv[1][0];
		float v20 = varying_uv[1][2] - varying_uv[1][0];
		Vec3f p10 = ndc_tri.col(1) - ndc_tri.col(0);
		Vec3f p20 = ndc_tri.col(2) - ndc_tri.col(0);

		Vec3f T = (p20 * v10 - p10 * v20) / (v10 * u20 - v20 * u10);
		Vec3f N = bn;
		T = T - N * (T * N);
		Vec3f Bi = cross(N, T);

		mat<3, 3, float> TBN;
		TBN.set_col(0, T.normalize());
		TBN.set_col(1, Bi.normalize());
		TBN.set_col(2, N.normalize());

		Vec3f normalTexWorld = (TBN * model->normal(uv)).normalize();
		float ndotL = std::max(0.f, normalTexWorld * L.normalize());

		//shadow
		Vec4f modelP = (Screen * Projection * ModelView).invert() * embed<4>(currentScreenPixel);	
		Vec4f lightP = (lightScreen * lightProjection * lightModelView) * modelP;
		Vec3f lightPNDC = Vec3f((lightP[0] / lightP[3]), (lightP[1] / lightP[3]), (lightP[2] / lightP[3]));
		int index = (int)lightPNDC[0] + (int)lightPNDC[1] * width;
		float currentPixelInShadowmap = shadowbuffer[index];
		float depthBias = 2.f;
		float inShadow = lightPNDC[2] + depthBias <  currentPixelInShadowmap;
		float shadowCol = 0.3 + 0.7 *(1 - inShadow);
		color = model->diffuse(uv) * shadowCol * ndotL;

		//color = TGAColor(shadowCol * 255, shadowCol * 255, shadowCol * 255, 255);
		return false;
	}
};

//void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) 
//{
//	for (float t = 0.0; t < 1.; t += 0.01f) 
//	{
//		int x = x0 + (x1 - x0) * t;
//		int y = y0 + (y1 - y0) * t;
//		image.set(x, y, color);
//	}
//}

//void line02(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color)
//{
//	for (int xD = x0; xD <= x1; xD += 1) 
//	{
//		float t = (float)(xD - x0) / (float)(x1 - x0);
//		int yD = y0 * (1 - t) + y1 * t;
//		image.set(xD, yD, color);
//	}
//}

//取变化率大的轴作为循环，避免空点
//x0,x1大小不应该受限制
//void line03(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color)
//{
//	//判断x0 x1 y0 y1大小
//	if (x0 > x1)
//	{
//		std::swap(x0, x1);
//	}
//
//	if (y0 > y1)
//	{
//		std::swap(y0, y1);
//	}
//	//判断变化率
//	if (std::abs(x1 - x0) > std::abs(y1 - y0)) 
//	{
//		for (int xD = x0; xD <= x1; xD++)
//		{
//			float t = (xD - x0) / (float)(x1 - x0);
//			int yD = y0 * (1. - t) + y1 * t;
//			image.set(xD, yD, color);
//		}
//	}
//	else
//	{
//		for (int yD = y0; yD <= y1; yD++)
//		{
//			float t = (yD - y0) / (float)(y1 - y0);
//			int xD = x0 * (1. - t) + x1 * t;
//			image.set(xD, yD, color);
//		}
//	}
//}

//void line04(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color) {
//	bool steep = false;
//	if (std::abs(x0 - x1) < std::abs(y0 - y1)) {
//		std::swap(x0, y0);
//		std::swap(x1, y1);
//		steep = true;
//	}
//	if (x0 > x1) {
//		std::swap(x0, x1);
//		std::swap(y0, y1);
//	}
//
//	for (int x = x0; x <= x1; x++) {
//		float t = (x - x0) / (float)(x1 - x0);
//		int y = y0 * (1. - t) + y1 * t;
//		if (steep) {
//			image.set(y, x, color);
//		}
//		else {
//			image.set(x, y, color);
//		}
//	}
//}


////matrix
//Matrix Convert2Matrix4x1(Vec3f v) 
//{
//	Matrix m(4, 1);
//	m[0][0] = v.x;
//	m[1][0] = v.y;
//	m[2][0] = v.z;
//	m[3][0] = 1.f;
//	return m;
//}
//
//Vec3f Convert2Vec3f(Matrix m)
//{
//	return Vec3f(m[0][0] / m[3][0], m[1][0] / m[3][0], m[2][0] / m[3][0]);
//}
//
//
//Matrix Convert2Screen(int xleft, int xright, int ybottom, int ytop) 
//{
//	Matrix m = Matrix::identity(4);
//	m[0][3] = xleft + xright / 2.f;
//	m[1][3] = ybottom + ytop / 2.f;
//	m[2][3] = 255.f / 2.f;
//
//	m[0][0] = xright / 2.f;
//	m[1][1] = ytop / 2.f;
//	m[2][2] = 255.f / 2.f;
//	return m;
//}
//
//Matrix lookat(Vec3f eye, Vec3f target, Vec3f up) 
//{
//	Vec3f z = (eye - target).normalize();
//	Vec3f x = (up ^ z).normalize();
//	Vec3f y = (z ^ x).normalize();
//	Matrix res = Matrix::identity(4);
//	for (int i = 0; i < 3; i++) 
//	{
//		res[0][i] = x[i];
//		res[1][i] = y[i];
//		res[2][i] = z[i];
//		res[i][3] = -target[i];
//	}
//	return res;
//}
//
////重心坐标
//bool barycentricInside(Vec2i* pts, Vec2i P)
//{ 
//	float x = P.x;
//	float y = P.y;
//	float xa = pts[0].x;
//	float xb = pts[1].x;
//	float xc = pts[2].x;
//	float ya = pts[0].y;
//	float yb = pts[1].y;
//	float yc = pts[2].y;
//	float a = (-(x - xb) * (yc - yb) + (y - yb) * (xc - xb)) / (-(xa - xb) * (yc - yb) + (ya - yb) * (xc - xb));
//	float b = (-(x - xc) * (ya - yc) + (y - yc) * (xa - xc)) / (-(xb - xc) * (ya - yc) + (yb - yc) * (xa - xc));
//	float c = 1 - a - b;
//	return a >= 0 && b >= 0 && c >= 0;
//}
//
////重心坐标
//Vec3f barycentricCoordinate(Vec2i* pts, Vec2i P)
//{
//	float x = P.x;
//	float y = P.y;
//	float xa = pts[0].x;
//	float xb = pts[1].x;
//	float xc = pts[2].x;
//	float ya = pts[0].y;
//	float yb = pts[1].y;
//	float yc = pts[2].y;
//	float a = (-(x - xb) * (yc - yb) + (y - yb) * (xc - xb)) / (-(xa - xb) * (yc - yb) + (ya - yb) * (xc - xb));
//	float b = (-(x - xc) * (ya - yc) + (y - yc) * (xa - xc)) / (-(xb - xc) * (ya - yc) + (yb - yc) * (xa - xc));
//	float c = 1 - a - b;
//	return Vec3f(a, b, c);
//}
//
//Vec3f barycentricCoordinateVec3(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
//{
//	float x = P.x;
//	float y = P.y;
//	float xa = A.x;
//	float xb = B.x;
//	float xc = C.x;
//	float ya = A.y;
//	float yb = B.y;
//	float yc = C.y;
//	float a = (-(x - xb) * (yc - yb) + (y - yb) * (xc - xb)) / (-(xa - xb) * (yc - yb) + (ya - yb) * (xc - xb));
//	float b = (-(x - xc) * (ya - yc) + (y - yc) * (xa - xc)) / (-(xb - xc) * (ya - yc) + (yb - yc) * (xa - xc));
//	float c = 1 - a - b;
//	return Vec3f(a, b, c);
//}


//void triangle(Vec2i* pts, TGAImage& image,TGAColor col)
//{
//	//AABB
//	Vec2i AABBX(0, image.get_width() - 1);
//	Vec2i AABBY(0, image.get_height() - 1);
//
//	int tempXMin = image.get_width() - 1;
//	int tempYMin = image.get_height() - 1;
//	int tempXMax = 0;
//	int tempYMax = 0;
//
//	for (int i = 0; i < 3; i++)
//	{
//		AABBX.x = std::min(pts[i].x, tempXMin);
//		tempXMin = AABBX.x;
//		AABBY.x = std::min(pts[i].y, tempYMin);
//		tempYMin = AABBY.x;
//
//		AABBX.y = std::max(pts[i].x, tempXMax);
//		tempXMax = AABBX.y;
//		AABBY.y = std::max(pts[i].y, tempYMax);
//		tempYMax = AABBY.y;
//	}
//
//	Vec2i P;
//	for (P.x = AABBX.x; P.x <= AABBX.y; P.x++)
//	{
//		for (P.y = AABBY.x; P.y <= AABBY.y; P.y++)
//		{
//			bool inside = barycentricInside(pts, P);
//			if (inside)
//			{
//				image.set(P.x, P.y, col);
//			}
//
//			/*Vec3f s = barycentric(pts, P);
//			if (s.x >= 0 && s.y >= 0 && s.z >= 0)
//			{
//				image.set(P.x, P.y, col);
//			}*/
//		}
//	}
//}
//
//void triangleZTest(Vec3f* pts, float* zbuffer, TGAImage& image, TGAColor col)
//{
//	//AABB
//	Vec2f AABBX(0, image.get_width() - 1);
//	Vec2f AABBY(0, image.get_height() - 1);
//
//	float tempXMin = image.get_width() - 1;
//	float tempYMin = image.get_height() - 1;
//	float tempXMax = 0;
//	float tempYMax = 0;
//
//	for (int i = 0; i < 3; i++)
//	{
//		AABBX.x = std::min(pts[i].x, tempXMin);
//		tempXMin = AABBX.x;
//		AABBY.x = std::min(pts[i].y, tempYMin);
//		tempYMin = AABBY.x;
//
//		AABBX.y = std::max(pts[i].x, tempXMax);
//		tempXMax = AABBX.y;
//		AABBY.y = std::max(pts[i].y, tempYMax);
//		tempYMax = AABBY.y;
//
//		//限制越界
//		AABBX.x = std::max(AABBX.x, 0.0f);
//		AABBX.y = std::min(AABBX.y, image.get_width() - 1.0f);
//		AABBY.x = std::max(AABBY.x, 0.0f);
//		AABBY.y = std::min(AABBY.y, image.get_height() - 1.0f);
//	}
//
//	Vec3f P;
//	for (P.x = AABBX.x; P.x <= AABBX.y; P.x++)
//	{
//		for (P.y = AABBY.x; P.y <= AABBY.y; P.y++)
//		{
//			Vec3f barycentricCoordinate = barycentricCoordinateVec3(pts[0], pts[1], pts[2], P);
//			if (barycentricCoordinate.x >= 0 && barycentricCoordinate.y >= 0 && barycentricCoordinate.z >= 0)
//			{
//				P.z = pts[0].z * barycentricCoordinate.x + pts[1].z * barycentricCoordinate.y + pts[2].z * barycentricCoordinate.z;
//				if (zbuffer[int(P.x + P.y * width)] < P.z)
//				{
//					zbuffer[int(P.x + P.y * width)] = P.z;
//					image.set(P.x, P.y, col);
//				}
//			}
//		}	
//	}
//}
//
//void triangleZTestDiffuse(Vec3f* pts,TGAColor* diffuse,float ndotl, float* zbuffer, TGAImage& image)
//{
//	//AABB
//	Vec2f AABBX(0, image.get_width() - 1);
//	Vec2f AABBY(0, image.get_height() - 1);
//
//	float tempXMin = image.get_width() - 1;
//	float tempYMin = image.get_height() - 1;
//	float tempXMax = 0;
//	float tempYMax = 0;
//
//	for (int i = 0; i < 3; i++)
//	{
//		AABBX.x = std::min(pts[i].x, tempXMin);
//		tempXMin = AABBX.x;
//		AABBY.x = std::min(pts[i].y, tempYMin);
//		tempYMin = AABBY.x;
//
//		AABBX.y = std::max(pts[i].x, tempXMax);
//		tempXMax = AABBX.y;
//		AABBY.y = std::max(pts[i].y, tempYMax);
//		tempYMax = AABBY.y;
//
//		//限制越界
//		AABBX.x = std::max(AABBX.x, 0.0f);
//		AABBX.y = std::min(AABBX.y, image.get_width() - 1.0f);
//		AABBY.x = std::max(AABBY.x, 0.0f);
//		AABBY.y = std::min(AABBY.y, image.get_height() - 1.0f);
//	}
//
//	Vec3f P;
//	for (P.x = AABBX.x; P.x <= AABBX.y; P.x++)
//	{
//		for (P.y = AABBY.x; P.y <= AABBY.y; P.y++)
//		{
//			Vec3f barycentricCoordinate = barycentricCoordinateVec3(pts[0], pts[1], pts[2], P);
//			if (barycentricCoordinate.x >= 0 && barycentricCoordinate.y >= 0 && barycentricCoordinate.z >= 0)
//			{
//				P.z = pts[0].z * barycentricCoordinate.x + pts[1].z * barycentricCoordinate.y + pts[2].z * barycentricCoordinate.z;
//				TGAColor dif = diffuse[0] * barycentricCoordinate.x + diffuse[1] * barycentricCoordinate.y + diffuse[2] * barycentricCoordinate.z;
//				if (zbuffer[int(P.x + P.y * width)] < P.z)
//				{
//					zbuffer[int(P.x + P.y * width)] = P.z;
//					image.set(P.x, P.y, dif * ndotl);
//				}
//			}
//		}
//	}
//}

//Vec3f world2screen(Vec3f v)
//{
//	return Vec3f(int((v.x + 1.) * width / 2. ), int((v.y + 1.) * height / 2.), v.z);
//}

//lesson8
float max_elevation_angle(float* zbuffer, Vec2f p, Vec2f dir) {
	float maxangle = 0;
	for (float t = 0.; t < 20.; t += 1.) {
		Vec2f cur = p + dir * t;
		if (cur.x >= width || cur.y >= height || cur.x < 0 || cur.y < 0) return maxangle;

		float distance = sqrt((p.x - cur.x) * (p.x - cur.x) + (p.y - cur.y) * (p.y - cur.y));
		if (distance < 1.f) continue;
		float elevation = zbuffer[int(cur.x) + int(cur.y) * width] - zbuffer[int(p.x) + int(p.y) * width];
		maxangle = std::max(maxangle, atanf(elevation / distance));
	}
	return maxangle;
}

int main(int argc, char** argv) {
	//TGAImage image(100, 100, TGAImage::RGB);
	//lesson 0 绘制像素点
	//image.set(52, 41, red);

	//lesson 1 绘制线段
	//line(13, 20, 80, 40, image, white);

	/*line02(13, 20, 80, 40, image, white);
	line02(13, 20, 80, 80, image, white);
	line02(13, 20, 20, 80, image, white);*/
	/*line03(13, 20, 80, 40, image, white);
	line03(13, 20, 80, 80, image, white);
	line03(13, 20, 20, 80, image, white);
	line03(50, 80, 10, 10, image, red);*/

	//lesson 1 绘制模型线框
	//if (2 == argc) {
	//	model = new Model(argv[1]);
	//}
	//else {
	//	model = new Model("obj/african_head.obj");
	//}
	//
	//TGAImage image(width, height, TGAImage::RGB);
	//for (int i = 0; i < model->nfaces(); i++) {
	//	std::vector<int> face = model->face(i);
	//	for (int j = 0; j < 3; j++) {
	//		Vec3f v0 = model->vert(face[j]);
	//		Vec3f v1 = model->vert(face[(j + 1) % 3]);
	//		int x0 = (v0.x + 1.) * width / 2.;
	//		int y0 = (v0.y + 1.) * height / 2.;
	//		int x1 = (v1.x + 1.) * width / 2.;
	//		int y1 = (v1.y + 1.) * height / 2.;
	//		line04(x0, y0, x1, y1, image, white);
	//	}
	//}
	//image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//image.write_tga_file("output.tga");
	//delete model;

	////lesson 2 绘制三角形
	//TGAImage img(200, 200, TGAImage::RGB);
	//Vec2i pts[3] = { Vec2i(10,10),Vec2i(100,30),Vec2i(190,160) };
	//triangle(pts, img, white);
	//img.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//img.write_tga_file("output.tga");

	////lesson 2 绘制兰伯特着色
	//Vec3f L(0, 0, -1);
	//TGAImage img(800, 800, TGAImage::RGB);
	//model = new Model("obj/african_head.obj");
	//for (int i = 0; i < model->nfaces(); i++)
	//{
	//	std::vector<int> face = model->face(i);
	//	Vec2i screen_coords[3];
	//	Vec3f world_coords[3];
	//	for (int j = 0; j < 3; j++) {
	//		Vec3f v = model->vert(face[j]);
	//		screen_coords[j] = Vec2i((v.x + 1.) * width / 2., (v.y + 1.) * height / 2.);
	//		world_coords[j] = v;
	//	}
	//	Vec3f n = (world_coords[2] - world_coords[0]) ^ (world_coords[1] - world_coords[0]);
	//	n.normalize();
	//	float NdotL = n * L;
	//	if (NdotL > 0)
	//	{
	//		triangle(screen_coords, img, TGAColor(255 * NdotL, 255 * NdotL, 255 * NdotL, 255));
	//	}
	//}
	//img.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//img.write_tga_file("output.tga");
	//delete model;


	/*//lesson4 vp matrix
	//if (2 == argc) {
	//		model = new Model(argv[1]);
	//	}
	//	else {
	//		model = new Model("obj/african_head.obj");
	//	}
	//	float* zbuffer = new float[width * height];
	//	for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
	//	TGAImage image(width, height, TGAImage::RGB);
	//	

	//	Matrix Projection = Matrix::identity();
	//	Projection[3][2] = -1.f / cameraPos.z;
	//	Projection[3][2] = -1.f / (cameraPos).norm();
	//	
	//	Matrix M_ViewP = lookat(cameraPos, Vec3f(0, 0, 0), Vec3f(0, 1, 0));
	//	Matrix M_ProjP = Projection;
	//	Matrix M_ScreenP = Convert2Screen(100, 600, 100, 600);

	//	for (int i = 0; i < model->nfaces(); i++)
	//	{
	//		std::vector<int> face = model->face(i);
	//		Vec3f screenXYworldZ[3];
	//		Vec3f vertexnormal[3];

	//		Vec2i uv[3];
	//		TGAColor dif[3];
	//		for (int j = 0; j < 3; j++) 
	//		{
	//			///////
	//			Vec3f worldP = model->vert(face[j]);
	//			Matrix M_worldP = Convert2Matrix4x1(worldP);

	//			Matrix finalM = M_ScreenP * M_ProjP * M_ViewP * M_worldP;
	//			screenXYworldZ[j] = Vec3f((int)Convert2Vec3f(finalM).x, (int)Convert2Vec3f(finalM).y, Convert2Vec3f(finalM).z);
	//			vertexnormal[j] = model->norm(i, j);

	//			uv[j] = model->uv(i, j);
	//			dif[j] = model->diffuse(uv[j]);
	//		}
	//		Vec3f n = (model->vert(face[2]) - model->vert(face[0])) ^ (model->vert(face[1]) - model->vert(face[0]));
	//		n.normalize();

	//		
	//		float NdotL = n * L;

	//		
	//		
	//		if(NdotL >0)
	//		{
	//			//triangleZTest(screenXYworldZ, zbuffer, image, TGAColor(255 * NdotL, 255 * NdotL, 255 * NdotL, 255));
	//			//triangleZTest(screenXYworldZ, zbuffer, image, dif * NdotL);
	//			triangleZTestDiffuse(screenXYworldZ, dif, NdotL, zbuffer, image);
	//		}
	//	}
	//	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	//	image.write_tga_file("output.tga");
	//	delete model;
	//return 0;*/

	//lesson6 
	if (2 == argc) 
	{
		model = new Model(argv[1]);
	}
	else 
	{
		model = new Model("obj/african_head.obj");
	}
	zbuffer = new float[width * height];
	for (int i = width * height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
	shadowbuffer = new float[width * height];
	for (int i = width * height; i--; shadowbuffer[i] = -std::numeric_limits<float>::max());
	TGAImage frame(width, height, TGAImage::RGB);
	Shader shader;

	TGAImage shadowMap(width, height, TGAImage::RGB);
	Shader shadowmapShader;


	//shadowmapmvps
	lookat(L, center, up);
	projection(0);
	Convert2Screen(100, 600, 100, 600);

	lightModelView = ModelView;
	lightProjection = Projection;
	lightScreen = Screen;   

	for (int i = 0; i < model->nfaces(); i++)
	{
		Vec3f shadowScreenV3f[3];
		for (int j = 0; j < 3; j++)
		{
			Vec4f gl_Screen = lightScreen * lightProjection * lightModelView * embed<4>(model->vert(i, j));
			shadowScreenV3f[j] = Vec3f((int)(gl_Screen[0] / gl_Screen[3]), (int)(gl_Screen[1] / gl_Screen[3]), (gl_Screen[2] / gl_Screen[3]));
		}
		triangleShadowmap(shadowScreenV3f, shadowMap, shadowbuffer);
	}


	//mvps
	lookat(cameraPos, center, up);
	projection(-1.f / (cameraPos - center).norm());
	Convert2Screen(100, 600, 100, 600);

	for (int i = 0; i < model->nfaces(); i++)
	{
		Vec3f screenV3f[3];
		for (int j = 0; j < 3; j++) 
		{
			shader.vertex(i, j);

			Vec4f gl_Screen = Screen * Projection * ModelView * embed<4>(model->vert(i, j));
			Vec3f w = model->vert(i, j);
		    screenV3f[j] = Vec3f((int)(gl_Screen[0] / gl_Screen[3]), (int)(gl_Screen[1] / gl_Screen[3]), (gl_Screen[2] / gl_Screen[3]));
		}
		triangle(screenV3f, shader, frame, zbuffer);
	}

	//lesson8
	//Post SSAO
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < height; j++)
		{
			if (zbuffer[i + j * width] < -10000) continue;
			float ssao = 0;
			for (float k = 0; k < pi * 2; k += pi / 4)
			{
				ssao += pi/2 - max_elevation_angle(zbuffer,Vec2f(i,j),Vec2f(cos(k),sin(k)));
			}
			ssao /= (pi / 2) * 8;
			ssao = pow(ssao, 4.f);
			//frame.set(i, j, white * ssao);
		}
	}

	//output
	frame.flip_vertically(); // to place the origin in the bottom left corner of the image
	frame.write_tga_file("framebuffer.tga");
	shadowMap.flip_vertically();
	shadowMap.write_tga_file("shadowmap.tga");
	delete[] shadowbuffer;
	delete model;
	delete[] zbuffer;
	return 0;
}
