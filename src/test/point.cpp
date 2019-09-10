#include <iostream>
#include <random>
#include "point.h"

int main(void)
{
    typedef kdtree::Point<1000> Point;
    typedef kdtree::distance<1000,1> Manhattan;
    typedef kdtree::distance<1000,2> Euclidean;
    typedef kdtree::distance<1000,std::numeric_limits<int>::max()> Tchebychev;
    typedef kdtree::distance<1000,6> Minkowski;

    // random distribution of the sites
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<float> distribution(-1024.f, 1024.f);

    Point a, b;

    for(int i = 0; i<1000; ++i)
    {
        a[i] = distribution(rng);
        b[i] = distribution(rng);
    }

    // std::cout << a << std::endl;
    // std::cout << b << std::endl;

    float man = Manhattan()(a,b); 
    float euc = Euclidean()(a,b); 
    float min = Minkowski()(a,b);
    float tch = Tchebychev()(a,b); 

    printf("%f %f %f %f\n", man, euc, min, tch); 

    return 0;
}