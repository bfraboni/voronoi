#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <random>
#include <cstdlib>
#include <cassert>
#include <functional>

#include "voronoization.h"
#include "transport.h"

int main( int argc, char * argv[] )
{
    if(argc < 6 )
        printf("usage : voronoi <sites> <iteration> <nb_input> <input>");

    // parameters
    const int size = atoi(argv[1]);
    const int max_iter = atoi(argv[2]);
    const int dtype = argc > 5 ? atoi(argv[5]) : 2;

    // i/o info
    const int nb_images = atoi(argv[3]);
    std::vector<Image> images(nb_images); 
    for( int i = 0; i < nb_images; ++i )
        images[i] = read_image(argv[i+4]);

    std::string output_path = argc > 4+nb_images ? argv[4+nb_images] : "";

    // image voronoisation
    // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams, Nivoliers V, Lévy B, 2013
    std::vector<Voronoization> voronois;
    for( int i = 0; i < (int)images.size(); ++i )
    {
        voronois.push_back( Voronoization( images[i], size, max_iter ) );

        Image voronoi = draw_cells_kd( images[i], voronois.back().sites );
        Image graph = draw_graph( voronoi, voronois.back().sites );
        
        std::stringstream ss;
        ss << "voronoi-" << i << ".png";
        write_image(voronoi, ss.str().c_str());
        ss.str("");
        ss << "graph-" << i << ".png";
        write_image(graph, ss.str().c_str());
    }

    // smooth transition
    int frame = 0;
    for( int i = 0; i < (int)images.size()-1; ++i )
    {
        const int curr = i;
        const int next = i+1;

        typedef kdtree::Point<5> Point5;

        std::vector<Point5> points_curr( size );
        std::vector<Point5> points_next( size );

        const Voronoization&        v_curr = voronois[curr];
        const Voronoization&        v_next = voronois[next];

        cinekine::voronoi::Graph    graph_curr = build_graph( v_curr.sites, v_curr.w, v_curr.h );
        cinekine::voronoi::Graph    graph_next = build_graph( v_next.sites, v_next.w, v_next.h );
        
        std::vector<Color>          colors_curr = evaluate_colors( graph_curr, v_curr.image );
        std::vector<Color>          colors_next = evaluate_colors( graph_next, v_next.image );

        #pragma omp for
        for (int j = 0; j < size; ++j)
        {   
            Color color_curr = colors_curr[j];
            Color color_next = colors_next[j];

            // Color color_curr = rgb2hsv(colors_curr[j]);
            // color_curr.r /= 360.f;
            // Color color_next = rgb2hsv(colors_next[j]);
            // color_next.r /= 360.f;

            points_curr[j] = { v_curr.sites[j].x / v_curr.w, v_curr.sites[j].y / v_curr.h, color_curr.r, color_curr.g, color_curr.b };
            points_next[j] = { v_next.sites[j].x / v_next.w, v_next.sites[j].y / v_next.h, color_next.r, color_next.g, color_next.b };
        }

        int w = v_curr.w, h = v_curr.h, iter = 100;
        Transport<Point5> transport(points_curr, points_next, iter);

        auto draw = [output_path, size, w, h, &frame]( const std::vector<Point5>& points ) -> void 
        {
            std::vector<vec2> sites( size );
            std::vector<Color> colors( size );
            #pragma omp for
            for (int k = 0; k < size; ++k)
            {
                sites[k] = vec2(points[k][0] * w, points[k][1] * h); 
                colors[k] = Color( points[k][2], points[k][3], points[k][4] ); 
                // colors[k] = hsv2rgb(Color( points[k][2]*360.f, points[k][3], points[k][4] )); 
            }

            Image voronoi = draw_cells_kd( sites, colors, w, h );
            std::stringstream ss;
            ss << output_path << "smooth-" << std::setfill('0') << std::setw(3) << frame << ".png";
            write_image(voronoi, ss.str().c_str());
            frame++;
        };   

        // transport.transport<>( draw );
        transport.transport( );
        
        assert( (int)transport.tmap.size() == size );
        int frames = 100;
        float step = 1.f / (frames - 1.f);

        for( int j = 0; j < 5; ++j )
            draw( points_curr );

        std::vector<Point5> points( size );
        for( int j = 0; j < frames; ++j )
        {
            float t = j * step;
            int id= 0;
            for (const auto& pair : transport.tmap)
            {   
                points[id] = lerp( points_curr[pair.first], points_next[pair.second], t );
                ++id;
            }
            draw( points );
        }

        for( int j = 0; j < 5; ++j )
            draw( points_next );
    }

    return 0;
}