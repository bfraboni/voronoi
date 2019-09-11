#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

#include "vec.h"
#include "color.h"
#include "image.h"
#include "image_io.h"


inline float square( float a ) { return a * a; }

inline float distance (const vec2& a, const vec2& b, const int p = 2 )
{
    // Manhattan
    if( p == 1 )
        return std::abs(a.x - b.x) + std::abs(a.y - b.y);
    // Euclidean
    else if( p == 2 )
        return std::sqrt(square(a.x - b.x) + square(a.y - b.y));
    // Tchebychev
    else if( p == -1 )
        return std::max(std::abs(a.x - b.x), std::abs(a.y - b.y));
    // Minkowski
    else if( p > 2 )
        return std::pow(std::pow(std::abs(a.x - b.x), p) + std::pow(std::abs(a.y - b.y), p), 1.f / p);
    // no distance
    else
        return 0;
}

struct Site
{
    vec2 position;
    Color color;
};

typedef std::vector<Site> Voronoi;

Image random_voronoisation(const Image& image, Voronoi& voronoi, int sites, int dtype = 2)
{
    // image size
    int size = image.width() * image.height();

    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

    // generate sites position
    for(int i = 0; i < sites; i++)
    {
        float y = distribution(rng) * image.width();
        float x = distribution(rng) * image.height();
        vec2 position(x, y);
        voronoi.push_back({position, Black()});
    } 

    // output image
    Image out(image.width(), image.height());

    // compute sites colors and output image
    #pragma omp for
    for( int i = 0; i < image.width(); i++ )
    for( int j = 0; j < image.height(); j++ )
    {   
        // find nearest site for this pixel
        int id = -1;
        float dmin = std::numeric_limits<float>::max();
        vec2 current(i,j);
        for( int k = 0; k < sites; k++ )
        {   
            Site& site = voronoi[k];
            float d = distance(site.position, current, dtype);
            if( d < dmin )
            {
                dmin = d;
                id = k;
            }
        }
        // printf("nearest : %d %d %d %f\n", i, j, id, dmin);
        
        // maj site mean
        assert(id >= 0);
        Site& nearest = voronoi[id];

        #pragma omp atomic
        nearest.color.r += image(i,j).r;
        #pragma omp atomic
        nearest.color.g += image(i,j).g;
        #pragma omp atomic
        nearest.color.b += image(i,j).b;
        #pragma omp atomic
        nearest.color.a += 1;

        // maj image
        out(i,j).r = id;
        out(i,j).g = dmin;
    }

    // compute new image
    #pragma omp for
    for( int i = 0; i < image.width(); i++ )
    for( int j = 0; j < image.height(); j++ )
    {
        Site& nearest = voronoi[out(i,j).r];
        if( nearest.color.a > 0 )
            out(i,j) = nearest.color / nearest.color.a;
    }

    return out;
}

int main( int argc, char * argv[] )
{
    int sites = atoi(argv[1]);
    Image image = read_image(argv[2]);
    int dtype = argc > 4 ? atoi(argv[4]) : 2;

    Voronoi voronoi;

    write_image(random_voronoisation(image, voronoi, sites, dtype), argv[3]);

    return 0;
}