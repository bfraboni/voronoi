#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

#include "vec.h"
#include "color.h"
#include "image.h"
#include "image_io.h"

struct KDNode
{   
    // childs
    int left = -1, right = -1;

    // position
    vec2 position;
};

struct CompareX
{
    bool operator() (const KDNode& a, const KDNode& b ) { return a.position.x < b.position.x; }
};

struct CompareY
{
    bool operator() (const KDNode& a, const KDNode& b ) { return a.position.y < b.position.y; }
};

struct KDTree
{
    std::vector<KDNode> nodes;
    KDTree(){}
};

int main( int argc, char * argv[] )
{
    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(0.f,1.f);

    // generate sites position
    float y = distribution(rng);
    float x = distribution(rng);

    return 0;
}