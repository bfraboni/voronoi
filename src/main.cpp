#include <iostream>
#include <limits>
#include <random>

#include "vec.h"
#include "color.h"
#include "image.h"
#include "image_io.h"

inline float square( float a ) { return a * a; }

inline float distance (const vec2& a, const vec2& b, const float p = 2 )
{
    // Manhattan
    if( p == 1 )
        return std::abs(a.x - b.x) + std::abs(a.y - b.y);
    // Euclidean
    else if( p == 2 )
        return std::sqrt(square(a.x - b.x) + square(a.y - b.y));
    // Tchebychev
    else if( p == std::numeric_limits<float>::infinity() )
        return std::max(std::abs(a.x - b.x), std::abs(a.y - b.y));
    // Minkowski
    else if( p > 2 )
        return std::pow(std::pow(a.x - b.x, p) + std::pow(a.y - b.y, p), 1.f / p);
    
    else
        return 0;
}

struct Site
{
    vec2 position;
    Color color;
};

typedef std::vector<Site> Voronoi;

Image random_voronoisation(const Image& image, Voronoi& voronoi, int sites)
{
    // image size
    int size = image.width() * image.height();

    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<int> distribution(size);

    // generate sites position
    for(int i = 0; i < sites; i++)
    {
        int offset = distribution(rng);
        int y = offset / image.width();
        int x = offset % image.width();
        vec2 position(x, y);
        voronoi.push_back({position, Black()});
    } 

    // compute sites colors and output image
    #pragma omp for
    for( int i = 0; i < sites; i++ )
    {
        Site& site = voronoi[i];
        Color& color = site.color;
        const vec2& position = site.position;

        // for( )
    }

}

Image reconstruction(Voronoi& voronoi, float distance)
{

}

int main( int argc, char * argv[] )
{
    int sites = std::atoi(argv[1]);
    Image image = read_image(argv[2]);

    Voronoi voronoi;






    return 0;
}