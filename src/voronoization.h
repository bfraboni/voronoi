#ifndef VORONOIZATION_H
#define VORONOIZATION_H

#include <omp.h>

#include "vec.h"
#include "color.h"
#include "image_io.h"

#include "kdtree.h"
#include "line.h"
#include "triangle.h"

#include "voronoi.hpp"
#include "legendre.hpp"

std::vector<vec2> init( const int w, const int h, const int size, std::mt19937& rng, std::uniform_real_distribution<float>& distribution )
{
    // generate sites position randomly
    std::vector<vec2> sites;
    for(int i = 0; i < size; i++)
        sites.push_back(vec2(distribution(rng) * w, distribution(rng) * h));

    return sites;
}

cinekine::voronoi::Graph build_graph( const std::vector<vec2>& sites, const int w, const int h )
{
    // copy sites for nearest neighbour search and voronoi construction
    cinekine::voronoi::Sites vsites( sites.size() );
    // #pragma omp for 
    for(int i = 0; i < (int)sites.size(); ++i) 
        vsites[i] = cinekine::voronoi::Vertex(sites[i].x, sites[i].y);

    // construct geometric voronoi diagram
    return cinekine::voronoi::build(std::move(vsites), w, h);
}

void evaluate_colors(   const cinekine::voronoi::Graph& graph, 
                        const Image& image,
                        std::vector<Color>& colors )
{   
    const float min_area = 1.f; 
    Dunavant dunavant(2);
    // #pragma omp for
    for(int i = 0; i < (int)graph.sites().size(); ++i)
    {   
        // get current site / cell info
        const int site_id = i;
        const auto& site = graph.sites()[site_id];
        const int cell_id = site.cell;
        const auto& cell = graph.cells()[cell_id];

        float wcell = 0;
        Color ccell;
        for(const auto& halfedge: cell.halfEdges)
        {
            // get current edge info
            const int edge_id = halfedge.edge;
            const auto& edge = graph.edges()[edge_id];

            // get triangle informations
            vec2 p0(edge.p0.x, edge.p0.y);
            vec2 p1(edge.p1.x, edge.p1.y);
            vec2 p2(site.x, site.y);
            Triangle t(p0, p1, p2);
            float area = t.area();

            // evaluate the triangle color integral using Dunavant quadrature
            Color ctriangle;
            float wtriangle = 0;
            for( int j = 0; j < dunavant.size(); ++j )
            {
                float w = dunavant.weight(j);
                vec2 uv = dunavant.point(j);
                vec2 p = t.point( uv );
                ctriangle += image.sample(p.x, p.y) * w;
                wtriangle += w;
            }
            if( wtriangle > 0 ) ctriangle /= wtriangle;

            ccell += ctriangle * area;
            wcell += area;
        }
        // normalize color
        if( wcell > 0 ) ccell /= wcell;
        
        colors[site_id] = ccell;
    }
}

void evaluate_areas( const cinekine::voronoi::Graph& graph, std::vector<float>& areas )
{   
    Dunavant dunavant(5);
    // #pragma omp for
    for(int i = 0; i < (int)graph.sites().size(); ++i)
    {   
        // get current site / cell info
        const int site_id = i;
        const auto& site = graph.sites()[site_id];
        const int cell_id = site.cell;
        const auto& cell = graph.cells()[cell_id];

        for(const auto& halfedge: cell.halfEdges)
        {
            // get current edge info
            const int edge_id = halfedge.edge;
            const auto& edge = graph.edges()[edge_id];

            // get triangle informations
            vec2 p0(edge.p0.x, edge.p0.y);
            vec2 p1(edge.p1.x, edge.p1.y);
            vec2 p2(site.x, site.y);
            Triangle t(p0, p1, p2);
            areas[site_id] += t.area();
        }
    }
}

void evaluate_gradient(     const cinekine::voronoi::Graph& graph, 
                            const std::vector<Color>& colors, 
                            const Image& image,
                            std::vector<vec2>& gradients )
{
    // foreach cell of the voronoi diagram compute the gradient of the objective function over edges
    typedef Rosetta::GaussLegendreQuadrature<vec2, 3> Legendre;
    Legendre legendre;
    // std::vector<vec2> gradients((int)graph.sites().size(), vec2::zero());

    // #pragma omp for
    for(int i = 0; i < (int)graph.sites().size(); ++i)
    {   
        // get current site / cell info
        const int site_id = i;
        const auto& site = graph.sites()[site_id];
        const int cell_id = site.cell;
        const auto& cell = graph.cells()[cell_id];

        for(const auto& halfedge: cell.halfEdges)
        {
            // get current edge info
            const int edge_id = halfedge.edge;
            const auto& edge = graph.edges()[edge_id];
            const int neighbor_id = edge.leftSite == site_id ? edge.rightSite : edge.leftSite;
            if( neighbor_id < 0 ) continue; // no neighbor

            vec2 pl( graph.sites()[site_id].x, graph.sites()[site_id].y);
            vec2 pk( graph.sites()[neighbor_id].x, graph.sites()[neighbor_id].y);

            // get segment info
            float dx = std::abs(edge.p1.x - edge.p0.x), dy = std::abs(edge.p1.y - edge.p0.y);
            vec2 n = normalize(vec2(-dy, dx));
            vec2 p0(edge.p0.x, edge.p0.y);
            vec2 p1(edge.p1.x, edge.p1.y);

            // if normal is ill-oriented
            if( dot(n, (pk - p0)) < 0 ) n = -n;

            // lambda function for edge integration 
            auto objective = [&pl, &pk, &n, site_id, neighbor_id, &image, &colors]( const vec2& p ) -> vec2 
            { 
                if( p.x < 0 || p.x >= image.width() || p.y < 0 || p.y >= image.height() ) return vec2(0, 0);

                // image gradient term
                const Color& color_pixel = image.sample(p.x, p.y); 
                 
                const Color& color_current = colors[site_id];  
                float gl = (color_pixel - color_current).length2();
                const Color& color_neighbor = colors[neighbor_id]; 
                float gk = (color_pixel - color_neighbor).length2();

                float g = gl - gk;

                // speed vector term times gradient term
                return g * (pl - p) / dot(n, (pl - pk));
            };

            // integral of F(x) over the edge segment using quadrature
            gradients[cell.site] = gradients[cell.site] + legendre.integrate<>(p0, p1, objective);
        }
        // normalize gradient
        // gradients[cell.site] = normalize(gradients[cell.site]);
    }
}

struct Voronoization
{
    std::vector<vec2> sites;
    Voronoization(){};

    Voronoization( const Image& image, const float size, const int max_iter )
    {
        printf("voronoization %d...\n", (int)omp_get_thread_num());
        const int w = image.width();
        const int h = image.height();

        // gradient descent params
        const float delta0 = 0.02f * std::sqrt(float(w * w + h * h)); 
        const float sigma = 0.5f; 

        // init sites randomly
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_real_distribution<float> distribution(0.f,1.f);
        sites = init(w, h, size, rng, distribution);

        // image voronoisation
        // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams, Nivoliers V, Lévy B, 2013
        cinekine::voronoi::Graph graph;
        std::vector<Color> colors( size, Black() );
        std::vector<Color> old_colors( size, Black() );
        std::vector<vec2> gradients( size, vec2::zero() );
        std::vector<vec2> old_gradients( size, vec2::zero() );
        std::vector<float> factors( size, 1 );
        std::vector<float> areas( size, 0 );

        std::vector<vec2> means( size, vec2::zero() );
        std::vector<float> variances( size, 0 );
        const float beta1 = 0.9;
        const float beta2 = 0.999;
        const float eps = 1e-6;

        for(int iter = 0; iter < max_iter; ++iter)
        {
            // compute geometric Voronoi graph
            cinekine::voronoi::Graph graph = build_graph( sites, w, h );

            // compute sites colors
            evaluate_colors( graph, image, colors );

            // compute sites colors
            evaluate_areas( graph, areas );

            // compute gradients
            evaluate_gradient( graph, colors, image, gradients );

            // kdtree construction
            typedef kdtree::Point<5> Point5;
            typedef kdtree::KDTree<float, 5> KDTree;
            std::vector<Point5> points( graph.sites().size() );
            #pragma omp parallel for 
            for( int i = 0; i < (int)graph.sites().size(); ++i)
                points[i] = {sites[i].x, sites[i].y, colors[i].r, colors[i].g, colors[i].b};

            KDTree kdtree( points, areas );

            // move sites in the gradient direction
            const float exp = iter / (max_iter - iter);
            const float delta = delta0 * std::pow(sigma, exp); 

        #if 0
            // classic gradient descenet
            #pragma omp for
            // for(int i = 0; i < this->size; ++i)
            //     sites[i] = sites[i] - delta * gradients[i];
        #else
            // Nesterov optimistaion
            const float gamma = 0.9;

            // #pragma omp for
            for(int i = 0; i < size; ++i)
            {
                // Nesterov simple + optimisation w.r.t to cell area (small cells moves faster)
                sites[i] = sites[i] - delta * 2.f / std::sqrt(areas[i]) * (gamma * normalize(old_gradients[i]) + normalize(gradients[i]));
                
                // take out of bounds cells back in
                if( sites[i].x < 0 || sites[i].x >= w )
                {
                    float rd = distribution(rng) * delta0;
                    sites[i].x = std::min(w - rd, std::max(rd, sites[i].x));
                }

                // take out of bounds cells back in
                if( sites[i].y < 0 || sites[i].y >= h )
                {
                    float rd = distribution(rng) * delta0;
                    sites[i].y = std::min(w - rd, std::max(rd, sites[i].y));
                }
            }
        #endif

            std::swap(gradients, old_gradients);
            std::swap(colors, old_colors);
            std::fill(gradients.begin(), gradients.end(), vec2::zero());
            std::fill(colors.begin(), colors.end(), Black());
            std::fill(factors.begin(), factors.end(), 1);
            std::fill(areas.begin(), areas.end(), 0);

            // printf("iteration %d\n", iter);
        }
        printf("voronoization %d done\n", (int)omp_get_thread_num());
    }
};

Image draw_graph( const Image& image, const std::vector<vec2>& sites )
{
    const int w = image.width();
    const int h = image.height();

    auto graph = build_graph( sites, w, h );

    // draw graph
    Image gimage(image);    
    for (auto& edge: graph.edges()) 
        WuLine(image, gimage, edge.p0.x, edge.p0.y, edge.p1.x, edge.p1.y, Color(0,0,0));

    return gimage;
}

Image draw_cells_kd( const Image& image, const std::vector<vec2>& sites )
{
    typedef kdtree::Point<2> Point2;
    typedef kdtree::KDTree<Color, 2> KDTree2Color;

    const int w = image.width();
    const int h = image.height();

    auto graph = build_graph( sites, w, h );
    std::vector<Color> colors( sites.size(), Black() );
    evaluate_colors( graph, image, colors );

    std::vector<Point2> points( graph.sites().size() );
    // #pragma omp for 
    for( int i = 0; i < (int)graph.sites().size(); ++i)
        points[i] = {graph.sites()[i].x, graph.sites()[i].y};

    KDTree2Color kdtree( points, colors );

    Image output(w, h);
    // #pragma omp for
    for( int i = 0; i < w; ++i )
    for( int j = 0; j < h; ++j )
    {   
        Point2 pixel = {float(i), float(j)};
        output(i, j) = kdtree.nearest_value( pixel );
    }
    return output;
}

Image draw_cells_kd( const std::vector<vec2>& sites, const std::vector<Color>& colors, const int w, const int h )
{
    typedef kdtree::Point<2> Point2;
    typedef kdtree::KDTree<Color, 2> KDTree2Color;

    std::vector<Point2> points( sites.size() );
    // #pragma omp for 
    for( int i = 0; i < (int)sites.size(); ++i)
        points[i] = { sites[i].x, sites[i].y };

    KDTree2Color kdtree( points, colors );

    Image output(w, h);
    // #pragma omp for
    for( int i = 0; i < w; ++i )
    for( int j = 0; j < h; ++j )
    {   
        Point2 pixel = {float(i), float(j)};
        output(i, j) = kdtree.nearest_value( pixel );
    }
    return output;
}

Color rgb2hsv(const Color& in)
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
    if( max > 0.0 ) { // NOTE: if Max is == 0, this divide would cause a crash
        out.g = (delta / max);                  // s
    } else {
        // if max is 0, then r = g = b = 0              
        // s = 0, h is undefined
        out.g = 0.0;
        out.r = NAN;                            // its now undefined
        return out;
    }
    if( in.r >= max )                           // > is bogus, just keeps compilor happy
        out.r = ( in.g - in.b ) / delta;        // between yellow & magenta
    else
    if( in.g >= max )
        out.r = 2.0 + ( in.b - in.r ) / delta;  // between cyan & yellow
    else
        out.r = 4.0 + ( in.r - in.g ) / delta;  // between magenta & cyan

    out.r *= 60.0;                              // degrees

    if( out.r < 0.0 )
        out.r += 360.0;

    return out;
}


Color hsv2rgb(const Color& in)
{
    double      hh, p, q, t, ff;
    long        i;
    Color         out;

    if(in.g <= 0.0) {       // < is bogus, just shuts up warnings
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
        out.r = in.b;
        out.g = t;
        out.b = p;
        break;
    case 1:
        out.r = q;
        out.g = in.b;
        out.b = p;
        break;
    case 2:
        out.r = p;
        out.g = in.b;
        out.b = t;
        break;

    case 3:
        out.r = p;
        out.g = q;
        out.b = in.b;
        break;
    case 4:
        out.r = t;
        out.g = p;
        out.b = in.b;
        break;
    case 5:
    default:
        out.r = in.b;
        out.g = p;
        out.b = q;
        break;
    }
    return out;     
}

Color rgb2hsl(const Color& rgb) {
    Color hsl;

    // float r = (rgb.r / 255.0f);
    // float g = (rgb.g / 255.0f);
    // float b = (rgb.b / 255.0f);
    float r = rgb.r;
    float g = rgb.g;
    float b = rgb.b;

    float min = std::min(std::min(r, g), b);
    float max = std::max(std::max(r, g), b);
    float delta = max - min;

    hsl.b = (max + min) / 2;

    if (delta == 0)
    {
        hsl.r = 0;
        hsl.g = 0.0f;
    }
    else
    {
        hsl.g = (hsl.b <= 0.5) ? (delta / (max + min)) : (delta / (2 - max - min));

        float hue;

        if (r == max)
        {
            hue = ((g - b) / 6) / delta;
        }
        else if (g == max)
        {
            hue = (1.0f / 3) + ((b - r) / 6) / delta;
        }
        else
        {
            hue = (2.0f / 3) + ((r - g) / 6) / delta;
        }

        if (hue < 0)
            hue += 1;
        if (hue > 1)
            hue -= 1;

        hsl.r = static_cast<int>(hue * 360.f);
    }

    return hsl;
}

float hue2rgb(float v1, float v2, float vH) 
{
    if (vH < 0)
        vH += 1;

    if (vH > 1)
        vH -= 1;

    if ((6 * vH) < 1)
        return (v1 + (v2 - v1) * 6 * vH);

    if ((2 * vH) < 1)
        return v2;

    if ((3 * vH) < 2)
        return (v1 + (v2 - v1) * ((2.0f / 3) - vH) * 6);

    return v1;
}

Color hsl2rgb(const Color& hsl) 
{
    float r = 0;
    float g = 0;
    float b = 0;

    if (hsl.g == 0)
    {
        r = g = b = (unsigned char)(hsl.b * 255);
    }
    else
    {
        float v1, v2;
        float hue = (float)hsl.r / 360.f;

        v2 = (hsl.b < 0.5) ? (hsl.b * (1 + hsl.g)) : ((hsl.b + hsl.g) - (hsl.b * hsl.g));
        v1 = 2 * hsl.b - v2;

        // r = 255.f * hue2rgb(v1, v2, hue + (1.f / 3.f));
        // g = 255.f * hue2rgb(v1, v2, hue);
        // b = 255.f * hue2rgb(v1, v2, hue - (1.f / 3.f));
        r = hue2rgb(v1, v2, hue + (1.f / 3.f));
        g = hue2rgb(v1, v2, hue);
        b = hue2rgb(v1, v2, hue - (1.f / 3.f));
    }

    return Color(r, g, b);
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
