#ifndef COLORSPACE_H
#define COLORSPACE_H

#include "color.h"

Color rgb2hsv( const Color& in )
{
    Color         out;
    double      min, max, delta;

    min = in.r < in.g ? in.r : in.g;
    min = min  < in.b ? min  : in.b;

    max = in.r > in.g ? in.r : in.g;
    max = max  > in.b ? max  : in.b;

    out.b = max;                                // v
    delta = max - min;
    if (delta < 0.00001)
    {
        out.g = 0;
        out.r = 0; // undefined, maybe nan?
        return out;
    }
    if( max > 0.0 ) 
    {   // NOTE: if Max is == 0, this divide would cause a crash
        out.g = (delta / max);                  // s
    } 
    else 
    {
        // if max is 0, then r = g = b = 0              
        // s = 0, h is undefined
        out.g = 0.0;
        out.r = NAN;                            // its now undefined
        return out;
    }
    if( in.r >= max )                           // > is bogus, just keeps compilor happy
        out.r = ( in.g - in.b ) / delta;        // between yellow & magenta
    else if( in.g >= max )
        out.r = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
    else
        out.r = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

    out.r *= 60.0;                              // degrees

    if( out.r < 0.0 )
        out.r += 360.0;

    return out;
}

Color hsv2rgb( const Color& in )
{
    double hh, p, q, t, ff;
    long i;
    Color out;

    if(in.g <= 0.0) 
    {       // < is bogus, just shuts up warnings
        out.r = in.b;
        out.g = in.b;
        out.b = in.b;
        return out;
    }

    hh = in.r;
    if(hh >= 360.0) hh = 0.0;
    hh /= 60.0;
    i = (long)hh;
    ff = hh - i;
    p = in.b * (1.0 - in.g);
    q = in.b * (1.0 - (in.g * ff));
    t = in.b * (1.0 - (in.g * (1.0 - ff)));

    switch(i) {
    case 0:
        out = Color(in.b, t, p);
        break;
    case 1:
        out = Color(q, in.b, p);
        break;
    case 2:
        out = Color(p, in.b, t);
        break;

    case 3:
        out = Color(p, q, in.b);
        break;
    case 4:
        out = Color(t, p, in.b);
        break;
    case 5:
    default:
        out = Color(in.b, p, q);
        break;
    }

    return out;     
}

Color rgb2yuv( const Color& rgb ) 
{
    double y = rgb.r * .299000 + rgb.g * .587000 + rgb.b * .114000;
    double u = rgb.r * -.168736 + rgb.g * -.331264 + rgb.b * .500000 + 128;
    double v = rgb.r * .500000 + rgb.g * -.418688 + rgb.b * -.081312 + 128;
    return Color(y, u, v);
}

Color yuv2rgb( const Color& yuv ) 
{
    float r = yuv.r + 1.4075 * (yuv.b - 128);
    float g = yuv.r - 0.3455 * (yuv.g - 128) - (0.7169 * (yuv.b - 128));
    float b = yuv.r + 1.7790 * (yuv.g - 128);
    return Color(r, g, b);
}

Color lab2rgb( const Color& lab)
{
    float   y = (lab.r + 16) / 116,
            x = lab.g / 500 + y,
            z = y - lab.b / 200,
            r, g, b;

    x = 0.95047 * ((x * x * x > 0.008856) ? x * x * x : (x - 16.0 / 116.0) / 7.787);
    y = 1.00000 * ((y * y * y > 0.008856) ? y * y * y : (y - 16.0 / 116.0) / 7.787);
    z = 1.08883 * ((z * z * z > 0.008856) ? z * z * z : (z - 16.0 / 116.0) / 7.787);

    r = x *  3.2406 + y * -1.5372 + z * -0.4986;
    g = x * -0.9689 + y *  1.8758 + z *  0.0415;
    b = x *  0.0557 + y * -0.2040 + z *  1.0570;

    r = (r > 0.0031308) ? (1.055 * std::pow(r, 1.0/2.4) - 0.055) : 12.92 * r;
    g = (g > 0.0031308) ? (1.055 * std::pow(g, 1.0/2.4) - 0.055) : 12.92 * g;
    b = (b > 0.0031308) ? (1.055 * std::pow(b, 1.0/2.4) - 0.055) : 12.92 * b;

    return clamp(Color(r, g, b, 1));
}

Color rgb2lab( const Color& rgb )
{
    float   r = rgb.r,
            g = rgb.g,
            b = rgb.b,
            x, y, z;

    r = (r > 0.04045) ? std::pow((r + 0.055) / 1.055, 2.4) : r / 12.92;
    g = (g > 0.04045) ? std::pow((g + 0.055) / 1.055, 2.4) : g / 12.92;
    b = (b > 0.04045) ? std::pow((b + 0.055) / 1.055, 2.4) : b / 12.92;

    x = (r * 0.4124 + g * 0.3576 + b * 0.1805) / 0.95047;
    y = (r * 0.2126 + g * 0.7152 + b * 0.0722) / 1.00000;
    z = (r * 0.0193 + g * 0.1192 + b * 0.9505) / 1.08883;

    x = (x > 0.008856) ? std::pow(x, 1.0 / 3.0) : (7.787 * x) + 16.0 / 116.0;
    y = (y > 0.008856) ? std::pow(y, 1.0 / 3.0) : (7.787 * y) + 16.0 / 116.0;
    z = (z > 0.008856) ? std::pow(z, 1.0 / 3.0) : (7.787 * z) + 16.0 / 116.0;

    return Color((116 * y) - 16, 500 * (x - y), 200 * (y - z));
}

#endif