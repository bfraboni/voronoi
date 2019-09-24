
#include "color.h"
#include <cmath>
#include <math.h>

float Color::power( ) const
{
    return (r+g+b) / 3.f;
}

float Color::sum( ) const
{
    return r+g+b;
}

float Color::length( ) const
{
    return std::sqrt(this->length2());
}

float Color::length2( ) const
{
    return r * r + g * g + b * b;
}

Color lerp(const Color& a, const Color& b, const float t) 
{   
    if( t == 0 ) return a;
    if( t == 1 ) return b;

    return Color(   a.r - t * (b.r - a.r), 
                    a.g - t * (b.g - a.g),
                    a.b - t * (b.b - a.b)  );
}

Color clamp(const Color& c, const float cmin, const float cmax) 
{   
    return Color(   std::fmax(cmin, std::fmin(cmax, c.r)), 
                    std::fmax(cmin, std::fmin(cmax, c.g)),
                    std::fmax(cmin, std::fmin(cmax, c.b)),
                    std::fmax(cmin, std::fmin(cmax, c.a))  );
}

Color Black( )
{
    return Color(0, 0, 0);
}

Color White( )
{
    return Color(1, 1, 1);
}

Color Red( )
{
    return Color(1, 0, 0);
}

Color Green( )
{
    return Color(0, 1, 0);
}

Color Blue( )
{
    return Color(0, 0, 1);
}

Color Yellow( )
{
    return Color(1, 1, 0);
}

Color operator+ ( const Color& a, const Color& b )
{
    return Color(a.r + b.r, a.g + b.g, a.b + b.b, a.a + b.a);
}

Color operator- ( const Color& c )
{
    return Color(-c.r, -c.g, -c.b, -c.a);
}

Color operator- ( const Color& a, const Color& b )
{
    return a + (-b);
}

Color operator* ( const Color& a, const Color& b )
{
    return Color(a.r * b.r, a.g * b.g, a.b * b.b, a.a * b.a);
}

Color operator* ( const float k, const Color& c )
{
    return Color(c.r * k, c.g * k, c.b * k, c.a * k);
}

Color operator* ( const Color& c, const float k )
{
    return k * c;
}

Color operator/ ( const Color& a, const Color& b )
{
    return Color(a.r / b.r, a.g / b.g, a.b / b.b, a.a / b.a);
}

Color operator/ ( const float k, const Color& c )
{
    return Color(k / c.r, k / c.g, k / c.b, k / c.a);
}

Color operator/ ( const Color& c, const float k )
{
    float kk= 1 / k;
    return kk * c;
}
