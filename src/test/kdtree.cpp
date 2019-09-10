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

struct KDNode
{   
    // childs
    int left = -1, right = -1;

    // position
    vec2 position;

    // depth
    int depth;
};

struct KDTree
{
    std::vector<KDNode> nodes;
    std::vector<vec2>& sites;

    KDTree(std::vector<vec2>& sites): sites(sites)
    {
        build(0, sites.size(), 0);
        printf("nodes %d\n", (int) nodes.size());
        print(std::string(), (int) nodes.size()-1, false);
    }

    int build( const int begin, const int end, const int depth )
    {
        if( begin >= end ) return -1;

        int axis = depth % 2;
        int len = end - begin;
        int mid = begin + len / 2;

        auto compare = [axis](const vec2& a, const vec2& b) 
        {
            return axis ? a.x < b.x : a.y < b.y;
        };

        std::nth_element(sites.data() + begin, sites.data() + mid, sites.data() + end, compare);

        int left = build( begin, mid, depth + 1 );
        int right = build( mid + 1, end, depth + 1 );

        KDNode node;
        node.depth = depth;
        node.position = sites[mid];
        node.left = left;
        node.right = right;

        nodes.push_back(node);
        
        return (int) nodes.size() - 1;
    }

    int search( const vec2& position )
    {
        float dmin = std::numeric_limits<float>::max();
        int best = -1;
        search_recursive( position, (int) nodes.size() - 1, 0, dmin, best );
        return best;
    }

    void search_recursive(const vec2& position, const int id, const int depth, float& dmin, int& best) 
    {
        if( id < 0 ) return;

        printf("search %d\n", id);
        float d = distance(position, nodes[id].position);
        int axis = depth % 2;
        float dx = axis ? position.x - nodes[id].position.x : position.y - nodes[id].position.y;

        if(d < dmin) 
        {
            best = id;
            dmin = d;
            printf("\tbest %d %f\n", best, dmin);
        }

        int near = dx <= 0 ? nodes[id].left : nodes[id].right;
        int far = dx <= 0 ? nodes[id].right : nodes[id].left;

        search_recursive( position, near, depth+1, dmin, best );
        
        if( square(dx) >= dmin ) return;
        
        search_recursive(position, far, depth+1, dmin, best );
    }


    void print( const std::string& prefix, const int node, bool isleft )
    {
        if( node < 0 || node >= (int)nodes.size() ) return;
        std::cout << prefix;
        std::cout << (isleft ? "├──" : "└──" );
        printf("(%.3f,%.3f)\n", nodes[node].position.x, nodes[node].position.y);
        print( prefix + (isleft ? "│   " : "    "), nodes[node].left, true);
        print( prefix + (isleft ? "│   " : "    "), nodes[node].right, false);
    }
};

int main( int argc, char * argv[] )
{
    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

    // generate sites position
    std::vector<vec2> v;
    for(int i = 0; i < 200; i++)
    {
        v.push_back(vec2(distribution(rng), distribution(rng)));
        // printf("%f %f\n", v.back().x, v.back().y);
    }

    // KDTree construction
    KDTree kdtree(v);

    // KDTree search
    vec2 test(0.2,0.2);
    printf("nearest kdtree search: %d\n\n", kdtree.search( test ));

    // Exhaustive search
    int best = -1;
    float dmin = std::numeric_limits<float>::max();
    for(int i = 0; i < (int)kdtree.nodes.size(); i++)
    {
        printf("search %d\n", i);
        float d = distance(kdtree.nodes[i].position, test); 
        if( d < dmin )
        {
            best = i;
            dmin = d;
            printf("\tbest %d %f\n", best, dmin);
        }
    }
    printf("nearest exhaustive search: %d\n\n", best);



    return 0;
}