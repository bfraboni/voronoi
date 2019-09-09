
#include <cassert>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "mat.h"


float radians( const float deg )
{
    return ((float) M_PI  / 180.f) * deg;
}

float degrees( const float rad )
{
    return (180.f / (float) M_PI) * rad;
}

Transform::Transform (
    const float t00, const float t01, const float t02, const float t03,
    const float t10, const float t11, const float t12, const float t13,
    const float t20, const float t21, const float t22, const float t23,
    const float t30, const float t31, const float t32, const float t33 )
{
    m[0][0]= t00; m[0][1]= t01; m[0][2]= t02; m[0][3]= t03;
    m[1][0]= t10; m[1][1]= t11; m[1][2]= t12; m[1][3]= t13;
    m[2][0]= t20; m[2][1]= t21; m[2][2]= t22; m[2][3]= t23;
    m[3][0]= t30; m[3][1]= t31; m[3][2]= t32; m[3][3]= t33;
}

Transform::Transform( const Vector& x, const Vector& y, const Vector& z, const Vector& w )
{
	m[0][0] = x.x;	m[0][1] = y.x;	m[0][2] = z.x;	m[0][3] = w.x;
	m[1][0] = x.y;	m[1][1] = y.y;	m[1][2] = z.y;	m[1][3] = w.y;
	m[2][0] = x.z;	m[2][1] = y.z;	m[2][2] = z.z;	m[2][3] = w.z;
	m[3][0] = 0;	m[3][1] = 0;	m[3][2] = 0;	m[3][3] = 1;
}

Vector Transform::operator[]( int c ) const
{
	assert(c >= 0 && c <= 3);
	return Vector(m[0][c], m[1][c], m[2][c]);
}


//! renvoie le point transforme.
Point Transform::operator() ( const Point& p ) const
{
    float x= p.x;
    float y= p.y;
    float z= p.z;

    float xt= m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];        // dot(vec4(m[0]), vec4(p, 1))
    float yt= m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];        // dot(vec4(m[1]), vec4(p, 1))
    float zt= m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];        // dot(vec4(m[2]), vec4(p, 1))
    float wt= m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];        // dot(vec4(m[3]), vec4(p, 1))

    assert(wt != 0);
    float w= 1.f / wt;
    if(wt == 1.f)
        return Point(xt, yt, zt);
    else
        return Point(xt*w, yt*w, zt*w);
}

//! renvoie le vecteur transforme.
Vector Transform::operator() ( const Vector& v ) const
{
    float x= v.x;
    float y= v.y;
    float z= v.z;

    float xt= m[0][0] * x + m[0][1] * y + m[0][2] * z;                  // dot(vec4(m[0]), vec4(v, 0))
    float yt= m[1][0] * x + m[1][1] * y + m[1][2] * z;                  // dot(vec4(m[1]), vec4(v, 0))
    float zt= m[2][0] * x + m[2][1] * y + m[2][2] * z;                  // dot(vec4(m[2]), vec4(v, 0))
    // dot(vec4(m[3]), vec4(v, 0)) == dot(vec4(0, 0, 0, 1), vec4(v, 0)) == 0 par definition

    return Vector(xt, yt, zt);
}

//! renvoie le point/vecteur homogene transforme.
vec4 Transform::operator() ( const vec4& v ) const
{
    float x= v.x;
    float y= v.y;
    float z= v.z;
    float w= v.w;

    float xt= m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3] * w;    // dot(vec4(m[0]), v)
    float yt= m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3] * w;    // dot(vec4(m[1]), v)
    float zt= m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3] * w;    // dot(vec4(m[2]), v)
    float wt= m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3] * w;    // dot(vec4(m[3]), v)

    return vec4(xt, yt, zt, wt);
}

//! renvoie la transposee de la matrice.
Transform Transform::transpose( ) const
{
    return Transform(
        m[0][0], m[1][0], m[2][0], m[3][0],
        m[0][1], m[1][1], m[2][1], m[3][1],
        m[0][2], m[1][2], m[2][2], m[3][2],
        m[0][3], m[1][3], m[2][3], m[3][3]);
}


Transform Transform::operator() ( const Transform& b ) const
{
    return compose_transform(*this, b);
}

//! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
Transform Transform::normal( ) const
{
    return inverse().transpose();
}


Transform Identity( )
{
    return Transform();
}

Transform Transpose( const Transform& m )
{
    return m.transpose();
}

Transform Inverse( const Transform& m )
{
    return m.inverse();
}

Transform Normal( const Transform& m )
{
    return m.normal();
}

Transform Scale( const float x, const float y, const float z )
{
    return Transform(
        x, 0, 0, 0,
        0, y, 0, 0,
        0, 0, z, 0,
        0, 0, 0, 1);
}

Transform Translation( const Vector& v )
{
    return Transform(
        1, 0, 0, v.x,
        0, 1, 0, v.y,
        0, 0, 1, v.z,
        0, 0, 0, 1);
}

Transform Translation( const float x, const float y, const float z )
{
    return Translation( Vector(x, y, z) );
}

Transform RotationX( const float a )
{
    float sin_t= sinf(radians(a));
    float cos_t= cosf(radians(a));

    return Transform(
        1,     0,      0, 0,
        0, cos_t, -sin_t, 0,
        0, sin_t,  cos_t, 0,
        0,     0,      0, 1 );
}

Transform RotationY( const float a )
{
    float sin_t= sinf(radians(a));
    float cos_t= cosf(radians(a));

    return Transform(
         cos_t, 0, sin_t, 0,
             0, 1,     0, 0,
        -sin_t, 0, cos_t, 0,
             0, 0,     0, 1 );
}

Transform RotationZ( const float a )
{
    float sin_t= sinf(radians(a));
    float cos_t= cosf(radians(a));

    return Transform(
        cos_t, -sin_t, 0, 0,
        sin_t,  cos_t, 0, 0,
            0,      0, 1, 0,
            0,      0, 0, 1 );
}

Transform Rotation( const Vector& axis, const float angle )
{
    Vector a= normalize(axis);
    float s= sinf(radians(angle));
    float c= cosf(radians(angle));

    return Transform(
        a.x * a.x + (1 - a.x * a.x ) * c,
        a.x * a.y * (1 - c ) - a.z * s,
        a.x * a.z * (1 - c ) + a.y * s,
        0,

        a.x * a.y * (1 - c ) + a.z * s,
        a.y * a.y + (1 - a.y * a.y ) * c,
        a.y * a.z * (1 - c ) - a.x * s,
        0,

        a.x * a.z * (1 - c ) - a.y * s,
        a.y * a.z * (1 - c ) + a.x * s,
        a.z * a.z + (1 - a.z * a.z ) * c,
        0,

        0, 0, 0, 1);
}


Transform Perspective( const float fov, const float aspect, const float znear, const float zfar )
{
    // perspective, openGL version
    float itan= 1 / tanf(radians(fov) * 0.5f);
    float id= 1 / (znear - zfar);

    return Transform(
        itan/aspect,    0,               0,                 0,
                  0, itan,               0,                 0,
                  0,    0, (zfar+znear)*id, 2.f*zfar*znear*id,
                  0,    0,              -1,                 0);
}

Transform Viewport( const float width, const float height )
{
    float w= width / 2.f;
    float h= height / 2.f;

    return Transform(
        w, 0,   0,   w,
        0, h,   0,   h,
        0, 0, .5f, .5f,
        0, 0,   0,   1);
}

Transform Lookat( const Point& from, const Point& to, const Vector& up )
{
    Vector dir= normalize( Vector(from, to) );
    Vector right= normalize( cross(dir, normalize(up)) );
    Vector newUp= normalize( cross(right, dir) );

    Transform m(
        right.x, newUp.x, -dir.x, from.x,
        right.y, newUp.y, -dir.y, from.y,
        right.z, newUp.z, -dir.z, from.z,
        0,       0,        0,     1);

    return m.inverse();
}

Transform compose_transform( const Transform& a, const Transform& b )
{
    Transform m;
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            m.m[i][j]= a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j] + a.m[i][3] * b.m[3][j];

    return m;
}

Transform operator* ( const Transform& a, const Transform& b )
{
    return compose_transform(a, b);
}

Transform Transform::inverse( ) const
{
    Transform minv= *this;

    int indxc[4], indxr[4];
    int ipiv[4] = { 0, 0, 0, 0 };

    for (int i = 0; i < 4; i++) {
        int irow = -1, icol = -1;
        float big = 0.f;

        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (fabsf(minv.m[j][k]) >= big) {
                            big = std::abs(minv.m[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        printf("singular matrix in make_inverse()\n");
                }
            }
        }

        assert(irow >= 0 && irow < 4);
        assert(icol >= 0 && icol < 4);

        ++ipiv[icol];
        // Swap rows _irow_ and _icol_ for pivot
        if (irow != icol) {
            for (int k = 0; k < 4; ++k)
                std::swap(minv.m[irow][k], minv.m[icol][k]);
        }

        indxr[i] = irow;
        indxc[i] = icol;
        if (minv.m[icol][icol] == 0.)
            printf("singular matrix in make_inverse()\n");

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        float pivinv = 1.f / minv.m[icol][icol];
        minv.m[icol][icol] = 1.f;
        for (int j = 0; j < 4; j++)
            minv.m[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                float save = minv.m[j][icol];
                minv.m[j][icol] = 0;
                for (int k = 0; k < 4; k++)
                    minv.m[j][k] -= minv.m[icol][k]*save;
            }
        }
    }

    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(minv.m[k][indxr[j]], minv.m[k][indxc[j]]);
        }
    }

    return minv;
}
