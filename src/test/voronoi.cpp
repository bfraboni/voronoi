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
        printf("usage : voronoi <sites> <iterations> <nb_input> <input> [path]");

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
    for( int i = 0; i < (int)images.size(); ++i )
    {
        Voronoization voronoization( images[i], size, max_iter );

        Image voronoi = draw_cells_kd( images[i], voronoization.sites );
        Image graph = draw_graph( voronoi, voronoization.sites );
        
        std::stringstream ss;
        ss << output_path << "voronoi-" << i << ".png";
        write_image(voronoi, ss.str().c_str());
        ss.str("");
        ss << output_path << "graph-" << i << ".png";
        write_image(graph, ss.str().c_str());
    }

    return 0;
}