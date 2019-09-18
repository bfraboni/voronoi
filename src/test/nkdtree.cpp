#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

#include "kdtree.h"

int main( int argc, char * argv[] )
{   
    typedef kdtree::Point<2> Point;
    typedef kdtree::KDTree<2, 2> KDTree;
    typedef kdtree::distance<2, 2> Distance;

    // random distribution of the points
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(-1024.f,1024.f);

    // generate sites position
    std::vector<Point> v;
    for(int i = 0; i < 10000; i++)
    {
        v.push_back(Point());
        for(int j = 0; j < 2; j++)
        {
            v.back()[j] = distribution(rng);
        }
    }

    // KDTree construction
    KDTree kdtree(v);

    // KDTree search
    // Point test = {12.4, 125.3, 954.5487, 519.76, -467.541};
    for( int px = -1024; px < 1024; px+=32)
    for( int py = -1024; py < 1024; py+=32)
    {
        Point test = {px, py};
        int id = kdtree.nearest( test );
        printf("\nnearest kdtree search:\t\t%d\t%f\n", id, Distance()(kdtree.nodes[id].point, test));

        // // Exhaustive search
        int best = -1;
        float dmin = std::numeric_limits<float>::max();
        for(int i = 0; i < (int)kdtree.nodes.size(); i++)
        {
            // printf("search %d\n", i);
            float d = Distance()(kdtree.nodes[i].point, test); 
            if( d < dmin )
            {
                best = i;
                dmin = d;
                // printf("\tbest %d %f\n", best, dmin);
            }
        }
        assert( id == best );
        printf("nearest exhaustive search:\t%d\t%f\n\n", best, dmin);
    }

    return 0;
}