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

typedef Rosetta::GaussLegendreQuadrature<vec2, 3> Legendre;

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
    const int max_iter = 300;
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

        
        Legendre legendre;
        // foreach cell of the voronoi diagram compute the gradient of the objective function over edges
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

                // line integral using quadrature : order 2
                {
                    float dx = std::fabs(edge.p1.x - edge.p0.x), dy = std::fabs(edge.p1.y - edge.p0.y);
                    vec2 n = normalize(vec2(dy, -dx));
                    vec2 ps = sites[site_id];
                    vec2 pn = sites[neighbor_id];

                    vec2 p(edge.p0.x, edge.p0.y);
                    if( p.x < 0 || p.x >= w || p.y < 0 || p.y >= h ) continue;
                    
                    vec2 q(edge.p1.x, edge.p1.y);
                    if( q.x < 0 || q.x >= w || q.y < 0 || q.y >= h ) continue;

                    auto project = [&dy, &dx](const vec2& p0, const vec2& p1, const vec2& p) -> vec2
                    {
                        float m = dy / dy;
                        float b = p0.y - (m * p0.x);

                        float x = (m * p.y + p.x - m * b) / (m * m + 1);
                        float y = (m * m * p.y + m * p.x + b) / (m * m + 1);

                        return vec2(x, y);
                    };
                    
                    vec2 proj = project(vec2(edge.p0.x, edge.p0.y), vec2(edge.p1.x, edge.p1.y), ps);
                    vec2 n2 = ps - proj;
                    if( dot(n, n2) > 0 ) 
                    {
                        n = -n;
                    }

                    auto objective = [ps, pn, n, site_id, neighbor_id, &image, &colors]( vec2 p ) -> vec2 
                    { 
                        // image gradient term
                        Color& color_pixel = image(p.x, p.y);  
                        
                        Color& color_current = colors[site_id];  
                        float g_current = (color_pixel - color_current).length2();

                        Color& color_neighbor = colors[neighbor_id]; 
                        float g_neighbor = (color_pixel - color_neighbor).length2();

                        float g = g_current - g_neighbor;

                        // speed vector term
                        float dot_prod = dot(n, (pn - p));
                        vec2 vs = dot_prod > 0 ? (ps - p) * g / dot_prod : vec2(0,0);

                        return vs;
                    };

                    gradients[cell.site] = gradients[cell.site] + legendre.integrate<>(p, q, objective);
                    
                    // {
                    //     // image gradient term
                    //     Color& color_pixel = image(p.x, p.y);  
                        
                    //     Color& color_current = colors[site_id];  
                    //     float g_current = (color_pixel - color_current).length2();

                    //     Color& color_neighbor = colors[neighbor_id]; 
                    //     float g_neighbor = (color_pixel - color_neighbor).length2();

                    //     float g = g_current - g_neighbor;

                    //     // speed vector term
                    //     float dot_prod = dot(n, (pn - p));
                    //     vec2 vs = dot_prod > 0 ? (ps - p) * g / dot_prod : vec2(0,0);

                    //     // integrate value
                    //     gradients[cell.site] = gradients[cell.site] + legendre;
                    // }

                    // {
                    //     // image gradient term
                    //     Color& color_pixel = image(q.x, q.y);  
                        
                    //     Color& color_current = colors[site_id];  
                    //     float g_current = (color_pixel - color_current).length2();

                    //     Color& color_neighbor = colors[neighbor_id]; 
                    //     float g_neighbor = (color_pixel - color_neighbor).length2();

                    //     float g = g_current - g_neighbor;

                    //     // speed vector term
                    //     float dot_prod = dot(n, (pn - q));
                    //     vec2 vs = dot_prod > 0 ? (ps - q) * g / dot_prod : vec2(0,0);

                    //     // integrate value
                    //     gradients[cell.site] = gradients[cell.site] + vs;
                    // }
                }
                
                // integrate gradients over edge
                // float x0 = edge.p0.x; 
                // float y0 = edge.p0.y; 
                // float x1 = edge.p1.x; 
                // float y1 = edge.p1.y; 
                // Traverse tr(x0, y0, x1, y1);
                // for( int nb = tr.n; nb > 0; --nb, ++tr )
                // {
                //     int x = tr.x, y = tr.y;
                //     if( x < 0 || x >= w || y < 0 || y >= h ) continue; // out of bound pixel

                //     // image gradient term
                //     Color& color_pixel = image(x, y);  
                    
                //     Color& color_current = colors[site_id];  
                //     float g_current = (color_pixel - color_current).length2();

                //     Color& color_neighbor = colors[neighbor_id]; 
                //     float g_neighbor = (color_pixel - color_neighbor).length2();

                //     float g = g_current - g_neighbor;

                //     // speed vector term
                //     vec2 n = normalize(vec2(tr.dy, -tr.dx));
                //     vec2 p(tr.x, tr.y);
                //     vec2 ps = sites[site_id];
                //     vec2 pn = sites[neighbor_id];

                //     assert( pn.x == sites[neighbor_id].x );
                //     assert( pn.y == sites[neighbor_id].y );
                //     float dot_prod = dot(n, (pn - p));
                //     vec2 vs = dot_prod > 0 ? (ps - p) / dot_prod : vec2(0,0);

                //     // integrate value
                //     gradients[cell.site] = gradients[cell.site] + vs * g;
                // }
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
            std::stringstream ss;
            ss << iter << "-" << argv[3];
            write_image(reconstruct(image, sites), ss.str().c_str());
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