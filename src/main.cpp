#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

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


void WuLine( Image& image, float x0, float y0, float x1, float y1, const Color& color ) 
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
            image(ypx11, xpx11) = image(ypx11, xpx11) + rfpart(yend) * xgap * (color - image(ypx11, xpx11));
            image(ypx11 + 1, xpx11) = image(ypx11 + 1, xpx11) + fpart(yend) * xgap * (color - image(ypx11+1, xpx11));
        } 
        else 
        {
            image(xpx11, ypx11) = image(xpx11, ypx11) + rfpart(yend) * xgap * (color - image(xpx11, ypx11));
            image(xpx11, ypx11 + 1) = image(xpx11, ypx11 + 1) + fpart(yend) * xgap * (color - image(xpx11, ypx11 + 1));
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
            image(ypx12, xpx12) = image(ypx12, xpx12) + rfpart(yend) * xgap *(color - image(ypx12, xpx12));
            image(ypx12 + 1, xpx12) = image(ypx12 + 1, xpx12) +  fpart(yend) * xgap *(color - image(ypx12 + 1, xpx12));
        } 
        else 
        {
            image(xpx12, ypx12) = image(xpx12, ypx12) + rfpart(yend) * xgap * (color - image(xpx12, ypx12));
            image(xpx12, ypx12 + 1) = image(xpx12, ypx12 + 1) + fpart(yend) * xgap * (color - image(xpx12, ypx12 + 1));
        }
    }
 
    if (steep) 
    {
        for (int x = xpx11 + 1; x < xpx12; x++) 
        {
            image(ipart(intery), x) = image(ipart(intery), x) + rfpart(intery) * (color - image(ipart(intery), x));
            image(ipart(intery) + 1, x) = image(ipart(intery) + 1, x) +  fpart(intery) * (color - image(ipart(intery) + 1, x));
            intery += gradient;
        }
    } 
    else 
    {
        for (int x = xpx11 + 1; x < xpx12; x++) 
        {
            image(x, ipart(intery)) = image(x, ipart(intery)) + rfpart(intery) * (color - image(x, ipart(intery)));
            image(x, ipart(intery) + 1) = image(x, ipart(intery) + 1) + fpart(intery) * (color - image(x, ipart(intery) + 1));
            intery += gradient;
        }
    }
}

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
    Image out( image.width(), image.height() );

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
        out(i,j).r = id;
    }

    // compute new image
    #pragma omp for
    for( int i = 0; i < image.width(); i++ )
    for( int j = 0; j < image.height(); j++ )
    {
        int id = out(i,j).r;
        if( colors[id].a > 0 )
            out(i,j) = colors[id] / colors[id].a;
    }

    // construct geometric voronoi diagram
    cinekine::voronoi::Sites vsites;
    for(int i = 0; i < (int)sites.size(); i++) 
        vsites.push_back(cinekine::voronoi::Vertex(sites[i][0], sites[i][1]));
    cinekine::voronoi::Graph graph = cinekine::voronoi::build(std::move(vsites), image.width(), image.height());
    
    // draw edges
    for (auto& cell: graph.cells()) 
    for (auto& halfedge : cell.halfEdges) 
    {
        auto& edge = graph.edges()[halfedge.edge];

        WuLine(out, edge.p0.x, edge.p0.y, edge.p1.x, edge.p1.y, Black());
    }

    write_image(out, argv[3]);

    return 0;
}