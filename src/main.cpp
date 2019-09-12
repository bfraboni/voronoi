#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>
#include <functional>

#include "vec.h"
#include "color.h"
#include "image.h"
#include "image_io.h"

#include "kdtree.h"
#include "voronoi.hpp"

typedef kdtree::Point<2> Point2;
typedef kdtree::KDTree<2> KDTree2;

void BresenhamLine( Image& image, float x1, float y1, float x2, float y2, const Color& color )
{
    const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
    if(steep)
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
    }

    if(x1 > x2)
    {
        std::swap(x1, x2);
        std::swap(y1, y2);
    }

    const float dx = x2 - x1;
    const float dy = fabs(y2 - y1);

    float error = dx / 2.0f;
    const int ystep = (y1 < y2) ? 1 : -1;
    int y = (int)y1;

    const int maxX = (int)x2;

    for(int x=(int)x1; x<maxX; x++)
    {
        if(steep)
        {
            image(y,x) = color;
        }
        else
        {
            image(x,y) = color;
        }

        error -= dy;
        if(error < 0)
        {
            y += ystep;
            error += dx;
        }
    }
}

void WuLine( const Image& image, Image& output, float x0, float y0, float x1, float y1, const Color& color ) 
{
    auto ipart = [](float x) -> int     {return int(std::floor(x));};
    auto round = [](float x) -> float   {return std::round(x);};
    auto fpart = [](float x) -> float   {return x - std::floor(x);};
    auto rfpart = [=](float x) -> float {return 1 - fpart(x);};
 
    const bool steep = abs(y1 - y0) > abs(x1 - x0);
    if (steep) 
    {
        std::swap(x0,y0);
        std::swap(x1,y1);
    }
    
    if (x0 > x1) 
    {
        std::swap(x0,x1);
        std::swap(y0,y1);
    }
 
    const float dx = x1 - x0;
    const float dy = y1 - y0;
    const float gradient = (dx == 0) ? 1 : dy/dx;
 
    int xpx11;
    float intery;
    {
        const float xend = round(x0);
        const float yend = y0 + gradient * (xend - x0);
        const float xgap = rfpart(x0 + 0.5);
        xpx11 = int(xend);
        const int ypx11 = ipart(yend);
        if (steep) 
        {
            output(ypx11, xpx11) = image(ypx11, xpx11) + rfpart(yend) * xgap * (color - image(ypx11, xpx11));
            output(ypx11 + 1, xpx11) = image(ypx11 + 1, xpx11) + fpart(yend) * xgap * (color - image(ypx11+1, xpx11));
        } 
        else 
        {
            output(xpx11, ypx11) = image(xpx11, ypx11) + rfpart(yend) * xgap * (color - image(xpx11, ypx11));
            output(xpx11, ypx11 + 1) = image(xpx11, ypx11 + 1) + fpart(yend) * xgap * (color - image(xpx11, ypx11 + 1));
        }
        intery = yend + gradient;
    }
 
    int xpx12;
    {
        const float xend = round(x1);
        const float yend = y1 + gradient * (xend - x1);
        const float xgap = rfpart(x1 + 0.5);
        xpx12 = int(xend);
        const int ypx12 = ipart(yend);
        if (steep) 
        {
            output(ypx12, xpx12) = image(ypx12, xpx12) + rfpart(yend) * xgap *(color - image(ypx12, xpx12));
            output(ypx12 + 1, xpx12) = image(ypx12 + 1, xpx12) +  fpart(yend) * xgap *(color - image(ypx12 + 1, xpx12));
        } 
        else 
        {
            output(xpx12, ypx12) = image(xpx12, ypx12) + rfpart(yend) * xgap * (color - image(xpx12, ypx12));
            output(xpx12, ypx12 + 1) = image(xpx12, ypx12 + 1) + fpart(yend) * xgap * (color - image(xpx12, ypx12 + 1));
        }
    }
 
    if (steep) 
    {
        for (int x = xpx11 + 1; x < xpx12; x++) 
        {
            output(ipart(intery), x) = image(ipart(intery), x) + rfpart(intery) * (color - image(ipart(intery), x));
            output(ipart(intery) + 1, x) = image(ipart(intery) + 1, x) +  fpart(intery) * (color - image(ipart(intery) + 1, x));
            intery += gradient;
        }
    } 
    else 
    {
        for (int x = xpx11 + 1; x < xpx12; x++) 
        {
            output(x, ipart(intery)) = image(x, ipart(intery)) + rfpart(intery) * (color - image(x, ipart(intery)));
            output(x, ipart(intery) + 1) = image(x, ipart(intery) + 1) + fpart(intery) * (color - image(x, ipart(intery) + 1));
            intery += gradient;
        }
    }
}

// traverse all pixels under a line
// https://stackoverflow.com/questions/3233522/elegant-clean-special-case-straight-line-grid-traversal-algorithm
// https://playtechs.blogspot.com/2007/03/raytracing-on-grid.html
struct Traverse
{
    float dx, dy, error;
    int x, y, n, x_inc, y_inc;

    explicit Traverse( float x0, float y0, float x1, float y1 )
    {
        // line slope
        dx = std::fabs(x1 - x0), dy = std::fabs(y1 - y0);
        // pixel start
        x = int(std::floor(x0)), y = int(std::floor(y0));
        // number of pixels to traverse
        n = 1;

        if (dx == 0)
        {
            x_inc = 0;
            error = std::numeric_limits<float>::infinity();
        }
        else if (x1 > x0)
        {
            x_inc = 1;
            n += int(floor(x1)) - x;
            error = (floor(x0) + 1 - x0) * dy;
        }
        else
        {
            x_inc = -1;
            n += x - int(floor(x1));
            error = (x0 - floor(x0)) * dy;
        }

        if (dy == 0)
        {
            y_inc = 0;
            error -= std::numeric_limits<float>::infinity();
        }
        else if (y1 > y0)
        {
            y_inc = 1;
            n += int(floor(y1)) - y;
            error -= (floor(y0) + 1 - y0) * dx;
        }
        else
        {
            y_inc = -1;
            n += y - int(floor(y1));
            error -= (y0 - floor(y0)) * dx;
        }
    }

    Traverse& operator++()
    {
        if (error > 0)
        {
            y += y_inc;
            error -= dx;
        }
        else
        {
            x += x_inc;
            error += dy;
        }

        return *this;
    }
};

int main( int argc, char * argv[] )
{
    if(argc < 4 || argc > 5)
        printf("usage : voronoi <size> <input> <output> [metric]");

    int size = atoi(argv[1]);
    Image image = read_image(argv[2]);
    int dtype = argc > 4 ? atoi(argv[4]) : 2;

    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

    // generate sites position randomly
    std::vector<Point2> sites;
    for(int i = 0; i < size; i++)
    {
        float y = distribution(rng) * image.width();
        float x = distribution(rng) * image.height();
        sites.push_back({x, y});
    } 

    // build kdtree for nearest neighbour search
    KDTree2 kdtree( sites );

    // output image
    Image tmp( image.width(), image.height() );

    // compute sites colors
    std::vector<Color> colors( sites.size() );
    #pragma omp for
    for( int i = 0; i < image.width(); i++ )
    for( int j = 0; j < image.height(); j++ )
    {   
        Point2 pixel = {(float)i, (float)j};
        int id = kdtree.nearest( pixel );
        assert(id >= 0);
        
        // maj site mean
        for( int k = 0; k < 4; k++ )
        {
            #pragma omp atomic
            colors[id](k) += image(i,j)(k);
        }

        // maj image
        tmp(i,j).r = id;
    }

    // normalize cells colors
    #pragma omp for
    for( int i = 0; i < (int)colors.size(); i++ )
        colors[i] = colors[i].a > 0 ? colors[i] / colors[i].a : Black();

    // compute voronoi image
    #pragma omp for
    for( int i = 0; i < image.width(); i++ )
    for( int j = 0; j < image.height(); j++ )
    {
        int id = tmp(i,j).r;
        tmp(i,j) = colors[id];
    }

    Image out = tmp;

    // construct geometric voronoi diagram
    cinekine::voronoi::Sites vsites;
    for(int i = 0; i < (int)sites.size(); i++) 
        vsites.push_back(cinekine::voronoi::Vertex(sites[i][0], sites[i][1]));
    cinekine::voronoi::Graph graph = cinekine::voronoi::build(std::move(vsites), image.width(), image.height());
    
    // optimisation process
    // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams
    //      Nivoliers V, Lévy B
    //
    // optimisation process
    // foreach cell compute the gradient of the objective function over edges
    std::vector<vec2> gradients(graph.sites().size(), vec2(0,0));
    #pragma omp for
    for(int i = 0; i < (int)graph.cells().size(); i++)
    {   
        const auto& cell = graph.cells()[i];
        for(const auto& halfedge: cell.halfEdges)
        {
            // get current edge
            const auto& edge = graph.edges()[halfedge.edge];
            
            // integrate gradients over edge
            Traverse tr(edge.p0.x, edge.p0.y, edge.p1.x, edge.p1.y);
            for( int nb = tr.n; nb > 0; --nb, ++tr )
            {
                int current = cell.site;
                int neighbor = edge.leftSite == current ? edge.rightSite : edge.leftSite;

                Color color_pixel = image(tr.x, tr.y);  
                Color color_current = colors[current];  
                Color color_neighbor = colors[neighbor]; 

                // color gradient
                float g_current = ((color_pixel - color_current) * (color_pixel - color_current)).sum() ;
                float g_neighbor = ((color_pixel - color_neighbor) * (color_pixel - color_neighbor)).sum() ;
                
                // angular speed term
                vec2 n = normalize(vec2(tr.dy, -tr.dx));
                vec2 p(tr.x, tr.y);
                vec2 pc(graph.sites()[current].x, graph.sites()[current].y);
                vec2 pn(graph.sites()[neighbor].x, graph.sites()[neighbor].y);
                float dot_prod = dot(n, (pn - p));
                vec2 vg = dot_prod > 0 ? (pc - p) / dot_prod : vec2(0,0);

                // integrate value
                gradients[cell.site] = gradients[cell.site] + vg * (g_current - g_neighbor);
            }
        }
    }

    Image grad(image.width(), image.height());
    // compute new image
    #pragma omp for
    for( int i = 0; i < image.width(); i++ )
    for( int j = 0; j < image.height(); j++ )
    {
        Point2 pixel = {(float)i, (float)j};
        int id = kdtree.nearest( pixel );
        assert(id >= 0);
        grad(i,j).r = gradients[id].x;
        grad(i,j).g = gradients[id].y;
    }

    // draw graph
    for (auto& edge: graph.edges()) 
        WuLine(tmp, out, edge.p0.x, edge.p0.y, edge.p1.x, edge.p1.y, Color(0.75,0.75,0.75));

    write_image(out, argv[3]);
    write_image(grad, "gradients.png");

    return 0;
}