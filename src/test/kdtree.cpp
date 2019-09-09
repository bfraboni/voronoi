#include <iostream>
#include <limits>
#include <random>
#include <cstdlib>
#include <cassert>

#include "vec.h"
#include "color.h"
#include "image.h"
#include "image_io.h"

// template<int N>
// struct Point
// {
//     float coords[N];
// };

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
        // 1 On cherche à placer le point C dans l'arbre, on applique donc l'algorithme d'insertion.
        // 2 Une fois arrivé à une feuille, ce point est sauvegardé comme étant le « meilleur candidat pour l'instant ».
        // 3 On remonte l'arbre :
        //     1 Si le nœud courant est plus proche, il devient le meilleur candidat ;
        //     2 On détermine si une partie de la cellule située de l'autre côté de l'hyperplan délimité par le nœud peut contenir un point plus proche. Pour cela, on détermine l'intersection de l'hyperplan avec la boule de recherche3 :
        //         si la boule de recherche traverse l'hyperplan, on effectue une recherche en descendant la branche correspondante ;
        //         sinon, on élimine la branche et on continue à remonter l'arbre.
        // 4 On s'arrête à la racine.
    }

    void print( const std::string& prefix, const int node, bool isleft )
    {
        if( node < 0 || node >= (int)nodes.size() ) return;

        std::cout << prefix;
        std::cout << (isleft ? "├──" : "└──" );
        printf("(%.3f,%.3f)\n", nodes[node].position.x, nodes[node].position.y);
        // std::cout << "(" << nodes[node].position.x << "," << nodes[node].position.x << ")" << std::endl;
        
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
    for(int i = 0; i < 10; i++)
    {
        v.push_back(vec2(distribution(rng), distribution(rng)));
        // printf("%f %f\n", v.back().x, v.back().y);
    }

    // KDTree construction
    KDTree kdtree(v);


    return 0;
}