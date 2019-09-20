#ifndef KDTREE_H
#define KDTREE_H

#include <cmath>
#include <limits>
#include <algorithm>
#include <queue>
#include "point.h"

namespace kdtree
{
    template<typename T, int N>           
    struct KDData
    {   
        Point<N> p;
        T v;
        int id;
        KDData() = default;
        KDData( const Point<N>& p, const T& v, int id ) : p(p), v(v), id(id) {}
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
                data[i] = KDData<T, N>( points[i], values[i], i );

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

        const point_type& nearest_point( const point_type& query )
        {
            float dmin = std::numeric_limits<float>::max();
            int best = 0;
            nearest( query, (int) nodes.size() - 1, 0, dmin, best );
            return nodes[best].data.p;
        }

        const int& nearest_index( const point_type& query )
        {
            float dmin = std::numeric_limits<float>::max();
            int best = 0;
            nearest( query, (int) nodes.size() - 1, 0, dmin, best );
            return nodes[best].data.id;
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

        template<class E, class Container = std::vector<T>, class Compare = std::less<typename Container::value_type>>
        class fixed_priority_queue : public std::priority_queue<E> 
        {
        public:
            explicit fixed_priority_queue(unsigned int size) : fixed_size(size) {}
            unsigned int max_size() const {return fixed_size;}
        private:
            const unsigned int fixed_size;
        };

        typedef std::pair<float, int> priority_type;
        struct priority_compare
        {
            bool operator()(const priority_type& a, const priority_type& b) const {return a.first < b.first;}
        };
        
        typedef fixed_priority_queue<priority_type, std::vector<priority_type>, priority_compare> queue_type;

        std::vector<int> knearest( const point_type& query, const int k ) 
        {
            queue_type queue(k);
            knearest( query, (int) nodes.size() - 1, 0, queue );
            std::vector<int> points;
            while ( !queue.empty() ) 
            { 
                points.push_back(queue.top().second);
                queue.pop();
            }
            std::reverse(points.begin(), points.end());
            return points;
        }

        void knearest( const point_type& query, const int id, const int depth, queue_type& queue )
        {
            if( id < 0 ) return;

            float d = kdtree::distance<N,P>()(query, nodes[id].data.p);
            int axis = depth % N;
            float dx = query[axis] - nodes[id].data.p[axis];

            // if point is nearer to the kth farthest, put point in queue
            if( queue.size() < queue.max_size() || d < queue.top().first ) 
            {
                queue.push( std::make_pair(d, id) );
                // keep k elements only
                if( queue.size() > queue.max_size() ) 
                    queue.pop(); 
            }

            int near = dx <= 0 ? nodes[id].left : nodes[id].right;
            int far = dx <= 0 ? nodes[id].right : nodes[id].left;

            knearest( query, near, depth+1, queue );
            
            if( std::abs(dx) >= queue.top().first && queue.size() >= queue.max_size() ) return;
            
            knearest( query, far, depth+1, queue );
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