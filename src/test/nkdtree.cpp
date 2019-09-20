#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

#include "kdtree.h"

int main( int argc, char * argv[] )
{   
    typedef kdtree::Point<2> Point;
    typedef kdtree::KDTree<float, 2, 2> KDTree;
    typedef kdtree::distance<2, 2> Distance;

    // random distribution of the points
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(-1024.f,1024.f);

    // generate sites position
    std::vector<Point> points;
    std::vector<float> values;
    for(int i = 0; i < 10000; i++)
    {
        points.push_back(Point());
        for(int j = 0; j < 2; j++)
        {
            points.back()[j] = distribution(rng);
        }
        values.push_back(distribution(rng));
    }

    // KDTree construction
    KDTree kdtree(points, values);

    // KDTree search
    // Point test = {12.4, 125.3, 954.5487, 519.76, -467.541};
    for( int px = -1024; px < 1024; px+=32)
    for( int py = -1024; py < 1024; py+=32)
    {
        Point test = {(float)px, (float)py};
        int id = kdtree.nearest( test );
        printf("\nnearest kdtree search:\t\t%d\t%f\n", id, Distance()(kdtree.nodes[id].data.p, test));

        // // Exhaustive search
        int best = -1;
        float dmin = std::numeric_limits<float>::max();
        for(int i = 0; i < (int)kdtree.nodes.size(); i++)
        {
            // printf("search %d\n", i);
            float d = Distance()(kdtree.nodes[i].data.p, test); 
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

    // KDTree KNN search
    for( int px = -1024; px < 1024; px+=32)
    for( int py = -1024; py < 1024; py+=32)
    {
        Point test = {(float)px, (float)py};
        std::vector<int> vec = kdtree.knearest( test, 10 );


        // // Exhaustive search
        std::vector< std::pair<float, int> > distances;
        for(int i = 0; i < (int)kdtree.nodes.size(); i++)
        {
            float d = Distance()(kdtree.nodes[i].data.p, test); 
            distances.push_back(std::make_pair(d, i));
        }
        std::sort(distances.begin(), distances.end(), [](const std::pair<float, int>& a, const std::pair<float, int>& b) -> bool {return a.first < b.first;});
        
        for( int i = 0; i < 10; ++i )
        {   
            printf("\n%d nearest\n", i);
            printf("kdtree search:\t\t%d\t%f\n", vec[i], Distance()(kdtree.nodes[vec[i]].data.p, test));
            printf("exhaustive search:\t%d\t%f\n\n", distances[i].second, distances[i].first);
            assert( vec[i] == distances[i].second );
        }

    }

    return 0;
}