#ifndef KDTREE_H
#define KDTREE_H

#include <cmath>
#include <limits>
#include <algorithm>
#include "point.h"

namespace kdtree
{
    template<typename T, int N>           
    struct KDData
    {   
        Point<N> p;
        T v;
        KDData() = default;
        KDData( const Point<N>& p, const T& v ) : p(p), v(v) {}
    };

    template<typename T, int N>        
    struct KDNode
    {   
        const KDData<T, N>& data;
        int left, right;

        KDNode( const KDData<T, N>& data, int left, int right ) : data(data), left(left), right(right) {}
    };

    template<typename T, int N, int P = 2>        
    struct KDTree
    {
        typedef T               value_type;
        typedef Point<N>        point_type;
        typedef KDData<T, N>    data_type;
        typedef KDNode<T, N>    node_type;

        std::vector<node_type> nodes;
        std::vector<data_type> data;

        KDTree( const std::vector<point_type>& points, const std::vector<value_type>& values ) : data(points.size())
        {
            #pragma omp for 
            for(int i = 0; i < (int)data.size(); ++i)
                data[i] = KDData<T, N>( points[i], values[i] );

            build(0, data.size(), 0);
        }

        int build( const int begin, const int end, const int depth )
        {
            if( begin >= end ) return -1;

            int axis = depth % N;
            int len = end - begin;
            int mid = begin + len / 2;

            auto compare = [axis](const data_type& a, const data_type& b) 
            {
                return a.p[axis] < b.p[axis];
            };

            std::nth_element(data.data() + begin, data.data() + mid, data.data() + end, compare);

            int left = build( begin, mid, depth + 1 );
            int right = build( mid + 1, end, depth + 1 );

            nodes.push_back( node_type( data[mid], left, right ) );
            
            return (int) nodes.size() - 1;
        }
        
        const data_type& nearest_data( const point_type& query )
        {
            float dmin = std::numeric_limits<float>::max();
            int best = 0;
            nearest( query, (int) nodes.size() - 1, 0, dmin, best );
            return nodes[best].data;
        }

        const value_type& nearest_value( const point_type& query )
        {
            float dmin = std::numeric_limits<float>::max();
            int best = 0;
            nearest( query, (int) nodes.size() - 1, 0, dmin, best );
            return nodes[best].data.v;
        }

        int nearest( const point_type& query )
        {
            float dmin = std::numeric_limits<float>::max();
            int best = 0;
            nearest( query, (int) nodes.size() - 1, 0, dmin, best );
            return best;
        }

        void nearest( const point_type& query, const int id, const int depth, float& dmin, int& best ) 
        {
            if( id < 0 ) return;

            float d = kdtree::distance<N,P>()(query, nodes[id].data.p);
            int axis = depth % N;
            float dx = query[axis] - nodes[id].data.p[axis];

            if(d < dmin) 
            {
                best = id;
                dmin = d;
            }

            int near = dx <= 0 ? nodes[id].left : nodes[id].right;
            int far = dx <= 0 ? nodes[id].right : nodes[id].left;

            nearest( query, near, depth+1, dmin, best );
            
            if( std::abs(dx) >= dmin ) return;
            
            nearest( query, far, depth+1, dmin, best );
        }

        void print( const std::string& prefix, const int node, bool isleft )
        {
            if( node < 0 || node >= (int)nodes.size() ) return;
            std::cout << prefix;
            std::cout << (isleft ? "├──" : "└──" );
            std::cout << node << " " << nodes[node].data.p;
            print( prefix + (isleft ? "│   " : "    "), nodes[node].left, true);
            print( prefix + (isleft ? "│   " : "    "), nodes[node].right, false);
        }
    };
}

#endif