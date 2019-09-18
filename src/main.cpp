#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <random>
#include <cstdlib>
#include <cassert>
#include <functional>

#include "voronoization.h"

int main( int argc, char * argv[] )
{
    if(argc != 6 )
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

        typedef kdtree::Point<2> Point2;
        typedef kdtree::KDTree<Color, 2> KDTree2Color;

        std::vector<Point2> points_curr( size );
        const Voronoization& v_curr = voronois[curr];
        cinekine::voronoi::Graph graph_curr = build_graph( v_curr.sites, v_curr.w, v_curr.h );
        std::vector<Color> colors_curr = evaluate_colors( graph_curr, v_curr.image );

        #pragma omp for
        for (int j = 0; j < size; ++j)
            points_curr[j] = { v_curr.sites[j].x, v_curr.sites[j].y };

        KDTree2Color kdtree( points_curr, colors_curr );

        std::vector<Point2> points_next( size );
        const Voronoization& v_next = voronois[next];
        cinekine::voronoi::Graph graph_next = build_graph( v_next.sites, v_next.w, v_next.h );
        std::vector<Color> colors_next = evaluate_colors( graph_next, v_next.image );

        #pragma omp for
        for (int j = 0; j < size; ++j)
            points_next[j] = { v_next.sites[j].x, v_next.sites[j].y };

        std::vector<int> nearests( size );
        std::vector<Point2> npoints( size );
        std::vector<Color> ncolors( size );
        #pragma omp for
        for (int j = 0; j < size; ++j)
        {
            const Point2& query = points_next[j];
            int best = kdtree.nearest( query );
            nearests[j] = kdtree.nodes[best].data.id;
            npoints[j] = kdtree.nodes[best].data.p;
            ncolors[j] = kdtree.nodes[best].data.v;
        }

        int nb_frames = 20;
        float step = 1.f / (frame - 1.f);
        for(int j = 0; j < nb_frames; ++j)
        {
            float t = j * step;

            std::vector<vec2> sites( size );
            std::vector<Color> colors( size );
            #pragma omp for
            for (int k = 0; k < size; ++k)
            {
                Point2 p = lerp(npoints[k], points_next[k], t); 
                sites[k] = vec2(p[0], p[1]); 
                
                Color a = rgb2hsv(ncolors[k]);
                Color b = rgb2hsv(colors_next[k]);
                Color c = hsv2rgb(a + t * (b - a));

                if( t == 1 )
                    colors[k] = colors_next[k]; 
                else if( t == 0 )
                    colors[k] = ncolors[k]; 
                else
                    colors[k] = c; 
            }

            Image voronoi = draw_cells_kd( sites, colors, v_curr.w, v_curr.h );
            std::stringstream ss;
            ss << "smooth-" << std::setfill('0') << std::setw(3) << frame << ".png";
            write_image(voronoi, ss.str().c_str());

            frame++;
        }
    }

    return 0;
}