
#include <algorithm>
#include <cmath>

#include "vec.h"


Point Origin( )
{
    return Point(0, 0, 0);
}

float distance( const Point& a, const Point& b )
{
    return length(a - b);
}

float distance2( const Point& a, const Point& b )
{
    return length2(a - b);
}

Point center( const Point& a, const Point& b )
{
    return Point((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
}

Point min( const Point& a, const Point& b )
{ 
    return Point( std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z) ); 
}

Point max( const Point& a, const Point& b ) 
{ 
    return Point( std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z) ); 
}

 
Vector operator- ( const Point& a, const Point& b )
{
    return Vector(a.x - b.x, a.y - b.y, a.z - b.z);
}

Point operator* ( const float k, const Point& a )
{
    return Point(k * a.x, k * a.y, k * a.z);
}

Point operator* ( const Point& a, const float k )
{
    return k * a;
}

Point operator/ ( const Point& a, const float k )
{ 
    float kk= 1.f / k; 
    return kk * a; 
}

Point operator+ ( const Point& a, const Point& b )
{
    return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

Vector operator- ( const Vector& v )
{
    return Vector(-v.x, -v.y, -v.z);
}

Point operator+ ( const Point& a, const Vector& v )
{
    return Point(a.x + v.x, a.y + v.y, a.z + v.z);
}

Point operator+ ( const Vector& v, const Point& a )
{
    return a + v;
}

Point operator- ( const Vector& v, const Point& a )
{
    return a + (-v);
}

Point operator- ( const Point& a, const Vector& v )
{
    return a + (-v);
}

Vector operator+ ( const Vector& u, const Vector& v )
{
    return Vector(u.x + v.x, u.y + v.y, u.z + v.z);
}

Vector operator- ( const Vector& u, const Vector& v )
{
    return Vector(u.x - v.x, u.y - v.y, u.z - v.z);
}

Vector operator* ( const float k, const Vector& v )
{
    return Vector(k * v.x, k * v.y, k * v.z);
}

Vector operator* ( const Vector& v, const float k )
{
    return k * v;
}

Vector operator* ( const Vector& a, const Vector& b )
{
    return Vector(a.x * b.x, a.y * b.y, a.z * b.z);
}

Vector operator/ ( const Vector& v, const float k )
{
    float kk= 1 / k;
    return kk * v;
}

Vector normalize( const Vector& v )
{
    float kk= 1 / length(v);
    return kk * v;
}

Vector cross( const Vector& u, const Vector& v )
{
    return Vector(
        (u.y * v.z) - (u.z * v.y),
        (u.z * v.x) - (u.x * v.z),
        (u.x * v.y) - (u.y * v.x));
}

float dot( const Vector& u, const Vector& v )
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

float length( const Vector& v )
{
    return std::sqrt(length2(v));
}

float length2( const Vector& v )
{
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

float dot(const vec2& a, const vec2& b) 
{
    return a.x * b.x + a.y * b.y;
}

vec2 lerp(const vec2& a, const vec2& b, const float t) 
{   
    if( t == 0 ) return a;
    if( t == 1 ) return b;

    return vec2(a.x - t * (b.x - a.x), a.y - t * (b.y - a.y));
}

vec2 operator-(const vec2& a, const vec2& b) 
{
    return vec2(a.x - b.x, a.y - b.y);
}

vec2 operator-(const vec2& a) 
{
    return vec2(-a.x, -a.y);
}

vec2 operator+(const vec2& a, const vec2& b) 
{
    return vec2(a.x + b.x, a.y + b.y);
}

vec2 operator+(const float v, const vec2& b) 
{
    return vec2(v + b.x, v + b.y);
}

vec2 operator+(const vec2& b, const float v) 
{
    return vec2(v + b.x, v + b.y);
}

vec2 operator/(const vec2& a, float v)
{
    return v != 0 ? vec2(a.x / v, a.y / v) : vec2(0,0);
}

vec2 operator*(const vec2& a, float v) 
{
    return vec2(a.x * v, a.y * v);
}

vec2 operator*(float v, const vec2& a) 
{
    return vec2(a.x * v, a.y * v);
}

float length2( const vec2& v )
{
    return v.x * v.x + v.y * v.y;
}

float length( const vec2& v )
{
    return std::sqrt(length2(v));
}

vec2 normalize(const vec2& n) 
{
    float v = std::sqrt(dot(n,n));
    return v > 0 ? n / v : n;
}