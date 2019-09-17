#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "vec.h"
#include "dunavant.hpp"

struct Dunavant
{
    int degree, order;
    double * points;
    double * weights;

    Dunavant( int degree = 2 ) : degree(degree)
    {
        order = dunavant_order_num( degree );
        points = new double[2*order];
        weights = new double[order];
        dunavant_rule ( degree, order, points, weights );

        for ( int i = 0; i < order_num; i++ )
        {
            printf("p%d %f %f %f\n", i, points[2*i], points[2*i+1], weights[i] );
        }
    }

    const int size() const {return order;}

    vec2 point(int i) const {return vec2(points[2*i],points[2*i+1]);}
    float weight(int i) const {return weights[i];}
};

struct Triangle
{
    vec2 a, b, c;
    
    Triangle( ) = default;
    Triangle( const vec2& a, const vec2& b, const vec2& c ) : a(a), b(b), c(c) {}
    Triangle& operator=(const Triangle&) = default;
    
    float area() const
    {
        return std::abs( 0.5f * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)));
    }

    vec2 point( const float u, const float v ) const
    {
        float w = 1.f - u - v;
        return vec2(a * w + b * u + c * v);
    }
};
#endif