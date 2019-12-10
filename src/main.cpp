#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <random>
#include <cstdlib>
#include <cassert>
#include <functional>
#include <omp.h>

#include "colorspace.h"
#include "voronoization.h"
#include "transport.h"

int main( int argc, char * argv[] )
{
    if( argc < 5 )
    {
        printf("usage : voronoi <sites> <iteration> <nb_input> <input>\n");
        return 1;
    }

    omp_set_nested(1);

    // parameters
    const int size = atoi(argv[1]);
    const int max_iter = atoi(argv[2]);
    const int dtype = argc > 5 ? atoi(argv[5]) : 2;

    // load images + i/o info
    const int nb_images = atoi(argv[3]);
    std::vector<Image> images(nb_images); 
    #pragma omp parallel for
    for( int i = 0; i < nb_images; ++i )
        images[i] = read_image(argv[i+4]);

    const int w = images[0].width();
    const int h = images[0].height();

    std::string output_path = argc > 4+nb_images ? argv[4+nb_images] : "";

    // image voronoisation
    // cf:  Approximating Functions on a Mesh with Restricted Voronoï Diagrams, Nivoliers V, Lévy B, 2013
    std::vector<Voronoization> voronois(images.size());
    #pragma omp parallel for
    for( int i = 0; i < (int)images.size(); ++i )
    {
        printf("thread %d:\tprocessing %s\n", (int)omp_get_thread_num(), argv[i+4]);
        
        voronois[i] = Voronoization( images[i], size, max_iter );

        printf("thread %d:\tvoronoization done.\n", (int)omp_get_thread_num());

        Image voronoi = draw_cells_kd( voronois[i].sites, voronois[i].colors, images[i].width(), images[i].height() );
        Image graph = draw_graph( voronoi, voronois[i].sites );
        
        std::stringstream ss;
        ss << "voronoi-" << i << ".bmp";
        write_image(voronoi, ss.str().c_str());
        ss.str("");
        ss << "graph-" << i << ".bmp";
        write_image(graph, ss.str().c_str());
    }

    if( 1 )
    {
        // smooth transition interpolation
        const int fixed_frames = 15;
        const int moving_frames = 100;
        const int total_frames = 2 * fixed_frames + moving_frames;

        #pragma omp parallel for
        for( int i = 0; i < (int)images.size(); ++i )
        {
            const int curr = i;
            const int next = (i+1)%images.size();

            typedef kdtree::Point<5> Point5;

            std::vector<Point5> points_curr( size );
            std::vector<Point5> points_next( size );

            const Voronoization& v_curr = voronois[curr];
            const Voronoization& v_next = voronois[next];

            std::vector<Color>& colors_curr = voronois[curr].colors;
            std::vector<Color>& colors_next = voronois[next].colors;

            #pragma omp parallel for
            for (int j = 0; j < size; ++j)
            {   
                Color color_curr, color_next;
                
                // RGB space
                if( 0 )
                {
                    color_curr = colors_curr[j];
                    color_next = colors_next[j];
                }

                // HSV space
                if( 0 )
                {
                    color_curr = rgb2hsv(colors_curr[j]);
                    color_curr.r /= 360.f;
                    color_next = rgb2hsv(colors_next[j]);
                    color_next.r /= 360.f;
                }

                // YUV space
                if( 0 )
                {
                    color_curr = rgb2yuv(colors_curr[j]);
                    color_next = rgb2yuv(colors_next[j]);
                }

                // LAB space
                if( 1 )
                {
                    color_curr = rgb2lab(colors_curr[j]);
                    
                    color_curr.r /= 100.f;
                    
                    color_curr.g += 128.f; 
                    color_curr.g /= 256.f; 
                    
                    color_curr.b += 128.f; 
                    color_curr.b /= 256.f; 
                    
                    color_curr = clamp(color_curr);

                    color_next = rgb2lab(colors_next[j]);
                    
                    color_next.r /= 100.f;
                    
                    color_next.g += 128.f; 
                    color_next.g /= 256.f; 

                    color_next.b += 128.f; 
                    color_next.b /= 256.f; 
                    
                    color_next = clamp(color_next);
                }

                points_curr[j] = { v_curr.sites[j].x / w, v_curr.sites[j].y / h, color_curr.r, color_curr.g, color_curr.b };
                points_next[j] = { v_next.sites[j].x / w, v_next.sites[j].y / h, color_next.r, color_next.g, color_next.b };
            }

            auto draw = [output_path, size, w, h]( const std::vector<Point5>& points, int frame ) -> void 
            {
                std::vector<vec2> sites( size );
                std::vector<Color> colors( size );

                #pragma omp parallel for
                for (int k = 0; k < size; ++k)
                {
                    sites[k] = vec2(points[k][0] * w, points[k][1] * h); 

                    // RGB space
                    if( 0 )
                        colors[k] = Color( points[k][2], points[k][3], points[k][4] ); 
                    
                    // HSV space
                    if( 0 )
                        colors[k] = hsv2rgb(Color( points[k][2]*360.f, points[k][3], points[k][4] )); 

                    // YUV space
                    if( 0 )
                        colors[k] = yuv2rgb(Color( points[k][2], points[k][3], points[k][4] )); 

                    // LAB space
                    if( 1 )
                    {
                        colors[k] = lab2rgb(
                            Color( 
                                points[k][2]*100.f, 
                                points[k][3]*256.f - 128.f, 
                                points[k][4]*256.f - 128.f 
                            )); 
                    }
                }

                Image voronoi = draw_cells_kd( sites, colors, w, h );
                // Image voronoi = draw_cells_kd_manhattan( sites, colors, w, h );
                // Image voronoi = draw_cells_kd_minkowski( sites, colors, w, h );
                Image graph = draw_graph( voronoi, sites );
                std::stringstream ss;
                ss << output_path << "smooth-" << std::setfill('0') << std::setw(3) << frame << ".bmp";
                write_image(voronoi, ss.str().c_str());
                // write_image(graph, ss.str().c_str());
            };   
            
            const float step = 1.f / (moving_frames - 1.f);

            #pragma omp parallel for
            for( int j = 0; j < fixed_frames; ++j )
            {
                draw( points_curr, i * total_frames + j );
                printf("frame %d\n", i * total_frames + j);
            }

            // use result of transport
            Transport<Point5> transport(points_curr, points_next, max_iter);

            if( 0 )
            {
                // use OT to find point assignement and linearly interpolate
                // transport.transport( );
                // std::vector<Point5> points( size );
                // for( int j = 0; j < moving_frames; ++j )
                // {
                //     float t = j * step;
                //     int id= 0;
                //     for (const auto& pair : transport.tmap)
                //     {   
                //         points[id] = lerp( points_curr[pair.first], points_next[pair.second], t );
                //         ++id;
                //     }
                //     draw( points, i * total_frames + j + fixed_frames );
                // }
            }
            else
            {
                // use OT to find Wasserstein barycenters to direclty interpolate in transport space

                // #pragma omp parallel for
                std::vector<Point5> points = points_curr;
                for( int j = 0; j < moving_frames; ++j )
                {
                    float t = j * step;
                    std::vector<Point5> points_out = transport.wasserstein_barycenters( t, &points );
                    std::swap(points_out, points);
                    draw( points, i * total_frames + j + fixed_frames );
                    printf("frame %d\n", i * total_frames + j + fixed_frames);
                }
            }

            #pragma omp parallel for
            for( int j = 0; j < fixed_frames; ++j )
            {
                draw( points_next, i * total_frames + j + fixed_frames + moving_frames );
                printf("frame %d\n", i * total_frames + j + fixed_frames + moving_frames);
            }
        }
    }

    return 0;
}
