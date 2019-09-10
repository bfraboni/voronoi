#ifndef POINTN_H
#define POINTN_H

#include <cmath>
#include <limits>
#include <algorithm>
#include "point.h"

namespace kdtree
{
    template<int N>        
    struct KDNode
    {   
        Point<N> point;
        int left = -1, right = -1;
    };

    template<int N, int P>        
    struct KDTree
    {
        typedef KDNode<N> node_type;
        typedef Point<N> point_type;

        std::vector<node_type> nodes;

        KDTree(std::vector<point_type>& points)
        {
            build(points, 0, points.size(), 0);
            printf("nodes %d\n", (int) nodes.size());
            print(std::string(), (int) nodes.size()-1, false);
        }

        int build(std::vector<point_type>& points, const int begin, const int end, const int depth )
        {
            if( begin >= end ) return -1;

            int axis = depth % N;
            int len = end - begin;
            int mid = begin + len / 2;

            auto compare = [axis](const point_type& a, const point_type& b) 
            {
                return a[axis] < b[axis];
            };

            std::nth_element(points.data() + begin, points.data() + mid, points.data() + end, compare);

            int left = build( points, begin, mid, depth + 1 );
            int right = build( points, mid + 1, end, depth + 1 );

            node_type node;
            node.point = points[mid];
            node.left = left;
            node.right = right;

            nodes.push_back(node);
            
            return (int) nodes.size() - 1;
        }

        int nearest( const point_type& query )
        {
            float dmin = std::numeric_limits<float>::max();
            int best = 0;
            nearest_recursive( query, (int) nodes.size() - 1, 0, dmin, best );
            return best;
        }

        void nearest_recursive( const point_type& query, const int id, const int depth, float& dmin, int& best ) 
        {
            if( id < 0 ) return;

            // printf("search %d\n", id);
            float d = kdtree::distance<N,P>()(query, nodes[id].point);
            int axis = depth % N;
            float dx = query[axis] - nodes[id].point[axis];

            if(d < dmin) 
            {
                best = id;
                dmin = d;
            }

            int near = dx <= 0 ? nodes[id].left : nodes[id].right;
            int far = dx <= 0 ? nodes[id].right : nodes[id].left;

            nearest_recursive( query, near, depth+1, dmin, best );
            
            if( std::abs(dx) >= dmin ) return;
            
            nearest_recursive( query, far, depth+1, dmin, best );
        }

        void print( const std::string& prefix, const int node, bool isleft )
        {
            if( node < 0 || node >= (int)nodes.size() ) return;
            std::cout << prefix;
            std::cout << (isleft ? "├──" : "└──" );
            std::cout << node << " " << nodes[node].point;
            print( prefix + (isleft ? "│   " : "    "), nodes[node].left, true);
            print( prefix + (isleft ? "│   " : "    "), nodes[node].right, false);
        }
    };
}

#endif