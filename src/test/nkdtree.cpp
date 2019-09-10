#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

#include "kdtree.h"

int main( int argc, char * argv[] )
{   
    typedef kdtree::Point<10> Point;
    typedef kdtree::KDTree<10> KDTree;
    typename kdtree::distance<10,2> Euclidean;

    // random distribution of the points
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(-1024.f,1024.f);

    // generate sites position
    std::vector<Point> v;
    for(int i = 0; i < 200; i++)
    {
        v.push_back(Point());
        for(int j = 0; j < 10; j++)
        {
            v.back()[j] = distribution(rng);
        }
    }

    // KDTree construction
    KDTree kdtree(v);

    // KDTree search
    Point test = {12.4, 125.3, 954.5487, 519.76, -467.541, -986.23, 654.2, 247.99, -12.45, 1.2354};
    printf("nearest kdtree search: %d\n\n", kdtree.nearest( test, Euclidean ));

    // // Exhaustive search
    // int best = -1;
    // float dmin = std::numeric_limits<float>::max();
    // for(int i = 0; i < (int)kdtree.nodes.size(); i++)
    // {
    //     printf("search %d\n", i);
    //     float d = distance(kdtree.nodes[i].position, test); 
    //     if( d < dmin )
    //     {
    //         best = i;
    //         dmin = d;
    //         printf("\tbest %d %f\n", best, dmin);
    //     }
    // }
    // printf("nearest exhaustive search: %d\n\n", best);

    return 0;
}