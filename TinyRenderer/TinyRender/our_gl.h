#ifndef __OUR_GL_H__
#define __OUR_GL_H__
#include "tgaimage.h"
#include "geometry.h"

extern Matrix ModelView;
extern Matrix Projection;
extern Matrix Screen;

extern Matrix lightModelView;
extern Matrix lightProjection;
extern Matrix lightScreen;

void Convert2Screen(int x, int y, int w, int h);
void projection(float coeff = 0.f); // coeff = -1/c
void lookat(Vec3f eye, Vec3f center, Vec3f up);

struct IShader {
    virtual ~IShader();
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    virtual bool fragment(Vec3f bar, TGAColor& color,Vec3f P) = 0;
};

//void triangle(Vec4f *pts, IShader &shader, TGAImage &image, float *zbuffer);
void triangle(Vec3f* pts, IShader& shader, TGAImage& image, float* zbuffer);
void triangleShadowmap(Vec3f* pts, TGAImage& image, float* zbuffer);
#endif //__OUR_GL_H__