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
#include "triangle.h"

#include "voronoi.hpp"
#include "legendre.hpp"


std::vector<vec2> init( const int w, const int h, const int size )
{
    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

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
    #pragma omp for 
    for(int i = 0; i < (int)sites.size(); ++i) 
        vsites[i] = cinekine::voronoi::Vertex(sites[i].x, sites[i].y);

    // construct geometric voronoi diagram
    return cinekine::voronoi::build(std::move(vsites), w, h);
}

std::vector<Color> evaluate_colors( const cinekine::voronoi::Graph& graph, const Image& image )
{   
    Dunavant dunavant;
    std::vector<Color> colors( (int)graph.sites().size() );
    #pragma omp for
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
            float area = t.area();

            // evaluate the triangle color integral using Dunavant quadrature
            for( int j = 0; j < dunavant.size(); ++j )
            {
                vec2 uv = dunavant.point(j);
                float w = dunavant.weight(j) * area;
                vec2 p = t.point( uv );
                colors[site_id] = colors[site_id] + image.sample(p.x, p.y) * w;
            }
        }
        // normalize color
        if( colors[site_id].a > 0 ) colors[site_id] = colors[site_id] / colors[site_id].a;
    }

    return colors;
}

std::vector<vec2> evaluate_gradient( const cinekine::voronoi::Graph& graph, const std::vector<Color>& colors, const Image& image )
{
    // foreach cell of the voronoi diagram compute the gradient of the objective function over edges
    typedef Rosetta::GaussLegendreQuadrature<vec2, 2> Legendre;
    Legendre legendre;
    std::vector<vec2> gradients((int)graph.sites().size(), vec2(0,0));

    #pragma omp for
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
                // int x = p.x, y = p.y;
                // if( x < 0 || x >= image.width() || y < 0 || y >= image.height() ) return vec2(0, 0);

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
            gradients[cell.site] = gradients[cell.site] - legendre.integrate<>(p0, p1, objective);
        }
        // normalize gradient
        gradients[cell.site] = normalize(gradients[cell.site]);
    }

    return gradients;
}

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
    std::vector<Color> colors = evaluate_colors( graph, image );

    std::vector< Point2 > points( graph.sites().size() );
    #pragma omp for 
    for( int i = 0; i < (int)graph.sites().size(); ++i)
        points[i] = {graph.sites()[i].x, graph.sites()[i].y};

    KDTree2Color kdtree( points, colors );

    Image output(w, h);
    #pragma omp for
    for( int i = 0; i < w; ++i )
    for( int j = 0; j < h; ++j )
    {   
        Point2 pixel = {float(i), float(j)};
        output(i, j) = kdtree.nearest_value( pixel );
    }
    return output;
}

int main( int argc, char * argv[] )
{
    if(argc < 5 || argc > 6)
        printf("usage : voronoi <sites> <iteration> <input> <output> [metric]");

    // parameters
    const int size = atoi(argv[1]);
    const int max_iter = atoi(argv[2]);
    const int dtype = argc > 5 ? atoi(argv[5]) : 2;

    // i/o info 
    Image image = read_image(argv[3]);
    const int w = image.width();
    const int h = image.height();
    std::string filename(argv[4]);

    // gradient descent params
    const float delta0 = 0.02f * std::sqrt(float(w * w + h * h)); 
    const float sigma = 0.5f; 

    // init sites randomly
    std::vector<vec2> sites = init(w, h, size);

    // image voronoisation
    // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams, Nivoliers V, Lévy B, 2013
    for(int iter = 0; iter < max_iter; ++iter)
    {
        // compute geometric Voronoi graph
        auto graph = build_graph( sites, w, h );
        
        // compute sites colors
        std::vector<Color> colors = evaluate_colors( graph, image );

        // compute gradients
        std::vector<vec2> gradients = evaluate_gradient( graph, colors, image );

        // move sites in the gradient direction
        const float exp = iter / (max_iter - iter);
        const float delta = delta0 * std::pow(sigma, exp); 

        #pragma omp for
        for(int i = 0; i < size; ++i)
            sites[i] = sites[i] + delta * gradients[i];

        printf("iteration %d\n", iter);
    }

    Image voronoi = draw_cells_kd( image, sites );
    Image graph = draw_graph( voronoi, sites );

    std::stringstream ss;
    ss << filename;
    write_image(voronoi, ss.str().c_str());
    ss.str("");
    ss << "graph-" << filename;
    write_image(graph, ss.str().c_str());

    return 0;
}