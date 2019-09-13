#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>
#include <functional>

#include "vec.h"
#include "color.h"
#include "image_io.h"

#include "kdtree.h"
#include "line.h"

#include "voronoi.hpp"

typedef kdtree::Point<2> Point2;
typedef kdtree::KDTree<2> KDTree2;

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