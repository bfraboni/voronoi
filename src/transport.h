#ifndef TRANSPORT_H
#define TRANSPORT_H

#include <chrono>
#include <random>
#include "point.h"

template<typename point_type>
struct Transport
{
    std::vector<point_type> input1, input2;
    // number of slices, pointset size
    int m = 64, point_size, max_iter; 

    // random distribution of the sites
    std::mt19937 rng;
    std::normal_distribution<float> distribution;

    Transport( const std::vector<point_type>& input1, const std::vector<point_type>& input2, const int max_iter ) : 
        input1(input1), 
        input2(input2), 
        point_size(input1.size()), 
        max_iter(max_iter),
        rng(std::random_device()()),
        distribution(0.f, 1.f)
    {}

    bool converged()
    {
        return true;
    }

    point_type random_direction()
    {
        point_type p;
        for( int i = 0; i < p.size(); ++i )
            p[i] = distribution(rng);
        return p;
    }

    template<typename Function>
    void transport( Function f )
    {
        int iter = 0;

        std::vector<point_type> displacement( point_size, point_type::zero() ), old_displacement( point_size, point_type::zero() ); 
        std::vector< std::pair<float, int> > projections1( point_size ), projections2( point_size ); 
        float gamma = 0.9;

        while( iter < max_iter )
        // while( !converged() )
        {
            std::vector<point_type> copy = input1;
            #pragma omp for
            for( int j = 0; j < point_size; ++j ) 
                copy[j] = copy[j] + gamma * old_displacement[j];
            
            #pragma omp for
            for( int i = 0; i < m; ++i )
            {
                // random ND slice
                point_type d = normalize(random_direction());

                // project points on the slice
                for( int j = 0; j < point_size; ++j )
                {
                    projections1[j] = std::make_pair( dot( copy[j], d ), j );
                    projections2[j] = std::make_pair( dot( input2[j], d ), j );
                }

                // sort projections
                auto compare = [&](const std::pair<float, int>& a, const std::pair<float, int>& b) -> bool{return a.first < b.first;};
                std::sort( projections1.begin(), projections1.end(), compare );
                std::sort( projections2.begin(), projections2.end(), compare );

                // aggregate displacement from p1 to p2
                for( int j = 0; j < point_size; ++j ) 
                {   
                    float l = (projections2[j].first - projections1[j].first);
                    point_type v = d * l / float(m);

                    int id = projections1[j].second;

                    for( int k = 0; k < d.size(); ++k )
                    {
                        #pragma omp atomic
                        displacement[id][k] += v[k];
                    }
                }
            } 

            // compute distance
            float sum = 0;
            for( int j = 0; j < point_size; ++j ) 
                sum += length2(displacement[j]);
            sum = std::sqrt(sum / point_size);
            printf("iteration %d rmse %f\n", iter, sum);

            // move points of input1
            #pragma omp for
            for( int j = 0; j < point_size; ++j ) 
            {   
                input1[j] = input1[j] + float(iter) / float(max_iter - iter) * (gamma * old_displacement[j] + displacement[j]);

                // force color clamping
                for( int k = 2; k < input1[j].size(); ++k )
                    input1[j][k] = std::min(1.f, std::max(0.f, input1[j][k]));
            }

            // if( iter%100 == 0 )
            f( input1 );

            std::swap(displacement, old_displacement);
            std::fill(displacement.begin(), displacement.end(), point_type::zero());

            iter++;
        }
    }
};


#endif