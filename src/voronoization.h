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

// generate sites position randomly
std::vector<vec2> init( const int w, const int h, const int size, std::mt19937& rng, std::uniform_real_distribution<float>& distribution )
{
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

void evaluate_colors_areas( const cinekine::voronoi::Graph& graph, 
                            const Image& image,
                            std::vector<Color>& colors,
                            std::vector<float>& areas )
{   
    Dunavant dunavant(3);
    for(int i = 0; i < (int)graph.sites().size(); ++i)
    {   
        // get current site - cell info
        const int site_id = i;
        const auto& site = graph.sites()[site_id];
        const int cell_id = site.cell;
        const auto& cell = graph.cells()[cell_id];

        Color ccell = Color::zero();
        for(const auto& halfedge: cell.halfEdges)
        {
            // get current edge info
            const int edge_id = halfedge.edge;
            const auto& edge = graph.edges()[edge_id];

            // get triangle informations
            Triangle t( vec2(edge.p0.x, edge.p0.y),
                        vec2(edge.p1.x, edge.p1.y),
                        vec2(site.x, site.y) );
            float area = t.area();

            // evaluate the triangle color using Dunavant triangle quadrature
            // good property of Dunavant quadrature: sum_i w_i = 1
            // no normalization needed
            Color ctriangle = Color::zero();
            for( int j = 0; j < dunavant.size(); ++j )
            {
                float w = dunavant.weight(j);
                vec2 uv = dunavant.point(j);
                vec2 p = t.point( uv );
                ctriangle += image.sample(p.x, p.y) * w;
            }

            if( ctriangle.a <= 0 ) ctriangle = Black();

            ccell += Color(ctriangle * area, area);
        }
        
        // cell area
        areas[site_id] = ccell.a;

        // normalize cell color
        ccell = ccell.a > 0 ? Color(ccell / ccell.a, 1) : Black();

        // cell color
        colors[site_id] = ccell;
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

            // lambda function for edge gradient integration 
            auto objective = [&pl, &pk, &n, site_id, neighbor_id, &image, &colors]( const vec2& p ) -> vec2 
            { 
                if( p.x < 0 || p.x >= image.width() || p.y < 0 || p.y >= image.height() ) return vec2(0, 0);

                // image gradient term
                const Color& color_pixel = image.sample(p.x, p.y); 
                
                // site color gradient
                const Color& color_current = colors[site_id];  
                float gl = (color_pixel - color_current).length2();
                
                // neighbor color gradient
                const Color& color_neighbor = colors[neighbor_id]; 
                float gk = (color_pixel - color_neighbor).length2();

                // speed vector term times gradient term
                return (gl - gk) * (pl - p) / dot(n, (pl - pk));
            };

            // integral of F(x) over the edge segment using quadrature
            gradients[cell.site] = gradients[cell.site] + legendre.integrate<>(p0, p1, objective);
        }
    }
}

struct Voronoization
{
    std::vector<vec2> sites;
    std::vector<Color> colors;
    cinekine::voronoi::Graph graph;

    Voronoization(){};

    // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams, Nivoliers V, Lévy B, 2013
    Voronoization( const Image& image, const float size, const int max_iter, const int seed = 1234 )
    {        
        const int w = image.width();
        const int h = image.height();
        // init sites randomly
        std::mt19937 rng( seed );
        std::uniform_real_distribution<float> distribution(0.f,1.f);
        sites = init(w, h, size, rng, distribution);

        // image voronoisation
        colors = std::vector<Color>( size, Black() );
        std::vector<Color> old_colors( size, Black() );
        std::vector<vec2> gradients( size, vec2::zero() ), old_gradients( size, vec2::zero() );
        std::vector<float> areas( size, 0 );

        // gradient descent params
        const float delta0 = 0.02f * std::sqrt(float(w * w + h * h)); 
        const float sigma = 0.5f;

        // Nesterov optimistaion
        const float gamma = 0.9;

        for(int iter = 0; iter < max_iter; ++iter)
        {
            for(int i = 0; i < size; ++i)
            {
                float da = areas[i] > 0.005 ? 2.f / std::sqrt(areas[i]) : 1;
                sites[i] = sites[i] - gamma * da * normalize(old_gradients[i]);
            }

            // compute geometric Voronoi graph
            graph = build_graph( sites, w, h );

            // compute cell colors and areas
            evaluate_colors_areas( graph, image, colors, areas );

            // compute gradients
            evaluate_gradient( graph, colors, image, gradients );

            // move sites in the gradient direction
            const float exp = iter / (max_iter - iter);
            const float delta = delta0 * std::pow(sigma, exp); 

            for(int i = 0; i < size; ++i)
            {
                // Nesterov simple + optimisation w.r.t to cell area (small cells moves faster)
                float da = areas[i] > 0.005 ? 2.f / std::sqrt(areas[i]) : 1;
                sites[i] = sites[i] - delta * da * (gamma * normalize(old_gradients[i]) + normalize(gradients[i]));
                // sites[i] = sites[i] - delta * (gamma * normalize(old_gradients[i]) + normalize(gradients[i]));
                
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

            std::swap(gradients, old_gradients);
            // std::swap(colors, old_colors);
            std::fill(gradients.begin(), gradients.end(), vec2::zero());
            std::fill(colors.begin(), colors.end(), Black());
            std::fill(areas.begin(), areas.end(), 0);
        }

        // compute geometric Voronoi graph
        graph = build_graph( sites, w, h );

        // compute cell colors and areas
        evaluate_colors_areas( graph, image, colors, areas );
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
    std::vector<float> areas( sites.size(), 0 );
    evaluate_colors_areas( graph, image, colors, areas );

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

Image draw_cells_kd_manhattan( const std::vector<vec2>& sites, const std::vector<Color>& colors, const int w, const int h )
{
    typedef kdtree::Point<2> Point2;
    typedef kdtree::KDTree<Color, 2, 1> KDTree2Color;

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

Image draw_cells_kd_tchebychev( const std::vector<vec2>& sites, const std::vector<Color>& colors, const int w, const int h )
{
    typedef kdtree::Point<2> Point2;
    typedef kdtree::KDTree<Color, 2, std::numeric_limits<int>::max()> KDTree2Color;

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

Image draw_cells_kd_minkowski( const std::vector<vec2>& sites, const std::vector<Color>& colors, const int w, const int h )
{
    typedef kdtree::Point<2> Point2;
    typedef kdtree::KDTree<Color, 2, 4> KDTree2Color;

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

#endif
