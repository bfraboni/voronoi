#ifndef POINT_H
#define POINT_H

#include <cmath>
#include <limits>
#include <iostream>

namespace kdtree
{
    template<int N>
    struct Point
    {
        float data[N];
        int size() const {return N;}
        
        float& operator[](int i) {return data[i];}
        float operator[](int i) const {return data[i];}

        typedef float* iterator;
        typedef const float* const_iterator;
        iterator begin() {return data;}
        iterator end() {return data + N;}
        const_iterator begin() const {return data;}
        const_iterator end() const {return data + N;}

        friend std::ostream& operator<< (std::ostream& out, const Point& pt)
        {
            out << "(";
            for(auto it = pt.begin(); it != pt.end(); ++it)
                out << *it << " ";
            out << ")\n";
            return out;
        }
    };

    template <int N>
    bool operator==(const Point<N>& a, const Point<N>& b) 
    {
        return std::equal(a.begin(), a.end(), b.begin());
    }

    template <int N>
    bool operator!=(const Point<N>& a, const Point<N>& b) 
    {
        return !(a == b);
    }

    template <int N>
    Point<N> operator+(const Point<N>& a, const Point<N>& b) 
    {
        Point<N> pt;
        for(int i = 0; i < N; ++i)
                pt[i] = a[i] + b[i];
        return pt;
    }

    template <int N>
    Point<N> operator-(const Point<N>& a, const Point<N>& b) 
    {
        Point<N> pt;
        for(int i = 0; i < N; ++i)
                pt[i] = a[i] - b[i];
        return pt;
    }

    template <int N>
    Point<N> operator/(const Point<N>& a, const Point<N>& b) 
    {
        Point<N> pt;
        for(int i = 0; i < N; ++i)
                pt[i] = a[i] / b[i];
        return pt;
    }

    template <int N>
    Point<N> operator*(const Point<N>& a, const Point<N>& b) 
    {
        Point<N> pt;
        for(int i = 0; i < N; ++i)
                pt[i] = a[i] * b[i];
        return pt;
    }

    template <int N>
    float dot(const Point<N>& a, const Point<N>& b) 
    {
        float res = 0;
        for(int i = 0; i < N; ++i)
                res += a[i] * b[i];
        return res;
    }

    // P > 2 : Minkowski distance
    template<int N, int P>
    struct distance
    {
        float operator() ( const Point<N>& a, const Point<N>& b ) const
        {
            float sum = 0.0;
            for(int i = 0; i < N; ++i)
                sum += std::pow(std::abs(a[i] - b[i]), P);
            return std::pow(sum, 1.f / P);
        }
    };

    // P = 1 : Manhattan distance
    template<int N>
    struct distance <N, 1>
    {
        float operator() ( const Point<N>& a, const Point<N>& b ) const
        {
            float sum = 0.0;
            for(int i = 0; i < N; ++i)
                sum += std::abs(a[i] - b[i]);
            return sum;
        }
    };

    // P = 2 : Euclidean distance
    template<int N>
    struct distance <N, 2>
    {
        float operator() ( const Point<N>& a, const Point<N>& b ) const
        {
            float sum = 0.0;
            for(int i = 0; i < N; ++i)
                sum += std::abs(a[i] - b[i]) * std::abs(a[i] - b[i]);
            return std::sqrt(sum);
        }
    };

    // P = inf : Tchebychev distance
    template<int N>
    struct distance <N, std::numeric_limits<int>::max()>
    {
        float operator() ( const Point<N>& a, const Point<N>& b ) const
        {
            float vmax = std::numeric_limits<float>::min();
            for(int i = 0; i < N; ++i)
                if( std::abs(a[i] - b[i]) > vmax )
                    vmax = std::abs(a[i] - b[i]);
            return vmax;
        }
    };
}

#endif