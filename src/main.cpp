#include <iostream>
#include <limits>
#include <sstream>
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

#include "legendre.hpp"

typedef Rosetta::GaussLegendreQuadrature<vec2, 2> Legendre;

typedef kdtree::Point<2> Point2;
typedef kdtree::KDTree<2> KDTree2;

void init(std::vector<vec2>& sites, const int w, const int h, const int size)
{
    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

    // generate sites position randomly
    for(int i = 0; i < size; i++)
    {
        float y = distribution(rng) * w;
        float x = distribution(rng) * h;
        sites.push_back(vec2(x, y));
    } 
}

Image reconstruct( const Image& image, const std::vector<vec2>& sites )
{
    int w = image.width(); 
    int h = image.height(); 
    int size = sites.size();

    // output image
    Image out( w, h );

    // copy sites for nearest neighbour search
    std::vector<Point2> ksites( size );
    #pragma omp for 
    for(int i = 0; i < size; ++i)
        ksites[i] = {sites[i].x, sites[i].y};
    
    // build kdtree for nearest neighbour search
    KDTree2 kdtree( ksites );

    // compute sites color means
    std::vector<Color> colors( size );
    #pragma omp for
    for( int i = 0; i < image.width(); ++i )
    for( int j = 0; j < image.height(); ++j )
    {   
        Point2 pixel = {(float)i, (float)j};
        int id = kdtree.nearest( pixel );
        // maj site mean
        for( int k = 0; k < 4; k++ )
        {
            #pragma omp atomic
            colors[id](k) += image(i,j)(k);
        }

        out(i, j).r = id;
    }

    // normalize sites color means
    #pragma omp for
    for( int i = 0; i < size; ++i )
        colors[i] = colors[i].a > 0 ? colors[i] / colors[i].a : Black();

    // compute voronoi image
    #pragma omp for
    for( int i = 0; i < w; i++ )
    for( int j = 0; j < h; j++ )
    {
        int id = out(i,j).r;
        out(i,j) = colors[id];
    }

    return out;
}

int main( int argc, char * argv[] )
{
    if(argc < 4 || argc > 5)
        printf("usage : voronoi <site_number> <input> <output> [metric]");

    // init parameters
    const int size = atoi(argv[1]);
    const int dtype = argc > 4 ? atoi(argv[4]) : 2;
    Image image = read_image(argv[2]);
    const int w = image.width();
    const int h = image.height();

    // init sites randomly
    std::vector<vec2> sites;
    init(sites, w, h, size);

    // image voronoisation
    // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams
    //      Nivoliers V, Lévy B
    const int max_iter = 200;
    for(int iter = 0; iter < max_iter; ++iter)
    {
        // copy sites for nearest neighbour search
        std::vector<Point2> ksites( size );
        #pragma omp for 
        for(int i = 0; i < size; ++i)
            ksites[i] = {sites[i].x, sites[i].y};
        
        // build kdtree for nearest neighbour search
        KDTree2 kdtree( ksites );

        // compute sites color means
        std::vector<Color> colors( size );
        #pragma omp for
        for( int i = 0; i < image.width(); ++i )
        for( int j = 0; j < image.height(); ++j )
        {   
            Point2 pixel = {(float)i, (float)j};
            int id = kdtree.nearest( pixel );
            // assert(id >= 0);
            for( int k = 0; k < 4; k++ )
            {
                #pragma omp atomic
                colors[id](k) += image(i, j)(k);
            }
        }

        // normalize sites color means
        #pragma omp for
        for( int i = 0; i < size; ++i )
            colors[i] = colors[i].a > 0 ? colors[i] / colors[i].a : Black();

        // construct geometric voronoi diagram
        cinekine::voronoi::Sites vsites( size );
        #pragma omp for
        for(int i = 0; i < size; ++i) 
            vsites[i] = cinekine::voronoi::Vertex(sites[i].x, sites[i].y);
        cinekine::voronoi::Graph graph = cinekine::voronoi::build(std::move(vsites), w, h);

        // foreach cell of the voronoi diagram compute the gradient of the objective function over edges
        Legendre legendre;
        std::vector<vec2> gradients(size, vec2(0,0));
        #pragma omp for
        for(int i = 0; i < size; ++i)
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
                const auto& neighbor = graph.sites()[neighbor_id];

                // integral of F(x) over the edge segment using quadrature
                {   
                    // get segment info
                    float dx = std::fabs(edge.p1.x - edge.p0.x), dy = std::fabs(edge.p1.y - edge.p0.y);
                    vec2 n = normalize(vec2(-dy, dx));
                    const vec2& pl = sites[site_id];
                    const vec2& pk = sites[neighbor_id];

                    vec2 p0(edge.p0.x, edge.p0.y);
                    // if( p0.x < 0 || p0.x >= w || p0.y < 0 || p0.y >= h ) continue;
                    vec2 p1(edge.p1.x, edge.p1.y);
                    // if( p1.x < 0 || p1.x >= w || p1.y < 0 || p1.y >= h ) continue;

                    // if normal is ill-oriented return backward
                    if( dot(n, (pk - p0)) < 0 ) n = -n;

                    // lambda function for edge integration 
                    auto objective = [pl, pk, n, site_id, neighbor_id, &image, &colors]( vec2 p ) -> vec2 
                    { 
                        // image gradient term
                        Color& color_pixel = image(p.x, p.y);  
                        
                        Color& color_current = colors[site_id];  
                        float gl = (color_pixel - color_current).length2();

                        Color& color_neighbor = colors[neighbor_id]; 
                        float gk = (color_pixel - color_neighbor).length2();

                        float g = gl - gk;

                        // speed vector term times gradient term
                        return g * (pl - p) / dot(n, (pl - pk));
                    };

                    gradients[cell.site] = gradients[cell.site] + legendre.integrate<>(p0, p1, objective);
                }
            }
            // normalize gradient
            gradients[cell.site] = normalize(gradients[cell.site]);
        }

        // move sites in the gradient direction
        const float exp = iter / (max_iter - iter); 
        const float delta0 = 0.02f * std::sqrt((float)w * (float)w + (float)h * (float)h);
        const float sigma = 0.5f; 
        const float delta = delta0 * std::pow(sigma, exp); 

        #pragma omp for
        for(int i = 0; i < size; ++i)
            sites[i] = sites[i] + delta * gradients[i];

        {
            // std::stringstream ss;
            // ss << iter << "-" << argv[3];
            // write_image(reconstruct(image, sites), ss.str().c_str());
        }

        printf("iteration %d\n", iter);
    }

    // draw graph
    // for (auto& edge: graph.edges()) 
    //     WuLine(tmp, out, edge.p0.x, edge.p0.y, edge.p1.x, edge.p1.y, Color(0.75,0.75,0.75));

    write_image(reconstruct(image, sites), argv[3]);
    // write_image(grad, "gradients.png");

    return 0;
}