#include <cmath>
#include <limits>
#include <cstdlib>
#include "our_gl.h"

Matrix ModelView;
Matrix Screen;
Matrix Projection;

Matrix lightModelView;
Matrix lightProjection;
Matrix lightScreen;

IShader::~IShader() {}

void Convert2Screen(int xleft, int xright, int ybottom, int ytop) 
{
    Screen = Matrix::identity();
    Screen[0][3] = xleft + xright / 2.f;
    Screen[1][3] = ybottom + ytop / 2.f;
    Screen[2][3] = 255.f / 2.f;

    Screen[0][0] = xright / 2.f;
    Screen[1][1] = ytop / 2.f;
    Screen[2][2] = 255.f / 2.f;
}

void projection(float coeff) {
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
}

void lookat(Vec3f eye, Vec3f target, Vec3f up) {
    Vec3f z = (eye - target).normalize();
    Vec3f x = cross(up, z).normalize();
    Vec3f y = cross(z, x).normalize();
    ModelView = Matrix::identity();
    for (int i = 0; i < 3; i++) {
        ModelView[0][i] = x[i];
        ModelView[1][i] = y[i];
        ModelView[2][i] = z[i];
        ModelView[i][3] = -target[i];
    }
}


Vec3f barycentricCoordinateVec3(Vec3f A, Vec3f B, Vec3f C, Vec3f P)
{
    float x = P.x;
    float y = P.y;
    float xa = A.x;
    float xb = B.x;
    float xc = C.x;
    float ya = A.y;
    float yb = B.y;
    float yc = C.y;
    float a = (-(x - xb) * (yc - yb) + (y - yb) * (xc - xb)) / (-(xa - xb) * (yc - yb) + (ya - yb) * (xc - xb));
    float b = (-(x - xc) * (ya - yc) + (y - yc) * (xa - xc)) / (-(xb - xc) * (ya - yc) + (yb - yc) * (xa - xc));
    float c = 1 - a - b;
    return Vec3f(a, b, c);
}

void triangle(Vec3f* pts, IShader& shader, TGAImage& image, float* buffer) 
{
    //AABB
    Vec2f AABBX(0, image.get_width() - 1);
    Vec2f AABBY(0, image.get_height() - 1);

    float tempXMin = image.get_width() - 1;
    float tempYMin = image.get_height() - 1;
    float tempXMax = 0;
    float tempYMax = 0;

    for (int i = 0; i < 3; i++)
    {
        AABBX.x = std::min(pts[i].x, tempXMin);
        tempXMin = AABBX.x;
        AABBY.x = std::min(pts[i].y, tempYMin);
        tempYMin = AABBY.x;

        AABBX.y = std::max(pts[i].x, tempXMax);
        tempXMax = AABBX.y;
        AABBY.y = std::max(pts[i].y, tempYMax);
        tempYMax = AABBY.y;

        //限制越界
        AABBX.x = std::max(AABBX.x, 0.0f);
        AABBX.y = std::min(AABBX.y, image.get_width() - 1.0f);
        AABBY.x = std::max(AABBY.x, 0.0f);
        AABBY.y = std::min(AABBY.y, image.get_height() - 1.0f);
    }
    Vec3f P;
    TGAColor color;
    for (P.x = AABBX.x; P.x <= AABBX.y; P.x++) 
    {
        for (P.y = AABBY.x; P.y <= AABBY.y; P.y++)
        {
            Vec3f a = pts[0];
            Vec3f b = pts[1];
            Vec3f c = pts[2];
            Vec3f barycentricCoordinate = barycentricCoordinateVec3(a, b, c, P);
            Vec2i f;
            if (barycentricCoordinate.x >= 0 && barycentricCoordinate.y >= 0 && barycentricCoordinate.z >= 0)
            {
                P.z = pts[0].z * barycentricCoordinate.x + pts[1].z * barycentricCoordinate.y + pts[2].z * barycentricCoordinate.z;
                if (buffer[int(P.x + P.y * image.get_width())] < P.z)
                {
                    buffer[(int)(P.x + P.y * image.get_width())] = P.z;
                    bool discard = shader.fragment(barycentricCoordinate, color, P);
                    if (!discard)
                    {
                        image.set(P.x, P.y, color);
                    }
                }
            }
        }
    }
}

void triangleShadowmap(Vec3f* pts, TGAImage& image, float* buffer)
{
    //AABB
    Vec2f AABBX(0, image.get_width() - 1);
    Vec2f AABBY(0, image.get_height() - 1);

    float tempXMin = image.get_width() - 1;
    float tempYMin = image.get_height() - 1;
    float tempXMax = 0;
    float tempYMax = 0;

    for (int i = 0; i < 3; i++)
    {
        AABBX.x = std::min(pts[i].x, tempXMin);
        tempXMin = AABBX.x;
        AABBY.x = std::min(pts[i].y, tempYMin);
        tempYMin = AABBY.x;

        AABBX.y = std::max(pts[i].x, tempXMax);
        tempXMax = AABBX.y;
        AABBY.y = std::max(pts[i].y, tempYMax);
        tempYMax = AABBY.y;

        //限制越界
        AABBX.x = std::max(AABBX.x, 0.0f);
        AABBX.y = std::min(AABBX.y, image.get_width() - 1.0f);
        AABBY.x = std::max(AABBY.x, 0.0f);
        AABBY.y = std::min(AABBY.y, image.get_height() - 1.0f);
    }
    Vec3f P;
    TGAColor color;
    for (P.x = AABBX.x; P.x <= AABBX.y; P.x++)
    {
        for (P.y = AABBY.x; P.y <= AABBY.y; P.y++)
        {
            Vec3f a = pts[0];
            Vec3f b = pts[1];
            Vec3f c = pts[2];
            Vec3f barycentricCoordinate = barycentricCoordinateVec3(a, b, c, P);
            Vec2i f;
            if (barycentricCoordinate.x >= 0 && barycentricCoordinate.y >= 0 && barycentricCoordinate.z >= 0)
            {
                P.z = pts[0].z * barycentricCoordinate.x + pts[1].z * barycentricCoordinate.y + pts[2].z * barycentricCoordinate.z;
                if (buffer[int(P.x + P.y * image.get_width())] < P.z)
                {
                    buffer[(int)(P.x + P.y * image.get_width())] = P.z;
                    image.set(P.x, P.y, TGAColor(255, 255, 255, 255) * (P.z / 255));
                }
            }
        }
    }
}

