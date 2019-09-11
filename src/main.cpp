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

typedef kdtree::Point<2> Point2;
typedef kdtree::KDTree<2> KDTree2;
typedef std::vector<Point2> Voronoi2;

int main( int argc, char * argv[] )
{
    if(argc < 4 || argc > 5)
        printf("usage : voronoi <sites> <input> <output> [metric]");

    int sites = atoi(argv[1]);
    Image image = read_image(argv[2]);
    int dtype = argc > 4 ? atoi(argv[4]) : 2;

    // image size
    int size = image.width() * image.height();

    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

    // generate sites position randomly
    std::vector<Point2> sites;
    for(int i = 0; i < sites; i++)
    {
        float y = distribution(rng) * image.width();
        float x = distribution(rng) * image.height();
        sites.push_back({x, y});
    } 

    // build kdtree for nearest neighbour search
    KDTree2 kdtree(sites);

    // output image
    Image out(image.width(), image.height());

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

    write_image(out, argv[3]);

    return 0;
}