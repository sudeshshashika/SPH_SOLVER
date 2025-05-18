/**used for the cell list acceleration structure**/
#pragma once
#include <iostream>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

class Cell{

public: 
    int x;
    CUDA_CALLABLE_MEMBER Cell( int x_ = -1 ): x{x_}{}
    CUDA_CALLABLE_MEMBER Cell( const Cell& src): x{src.x}{}
    CUDA_CALLABLE_MEMBER Cell &operator=( const Cell& src){
        if( this == &src)
            return *this;
        x = src.x;
        return *this;
    }
    CUDA_CALLABLE_MEMBER ~Cell(){}
};