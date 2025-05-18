/** used when allocating device bytes and number of threads per block and grid size**/
#pragma once 
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream>
#include <string>
#include <fstream>
#include "state.h"

class device{

public:
    int nCells;     
    int nBytes;    
    int cBytes;     
    int TPB;          
    int num_blocks;   
    int cNumBlocks;   
    
    void setThreadsPerBlock(int N){
        TPB = N;
    }
    void setNumBlocks(const state& src){
        num_blocks = (src.particleCount + TPB -1) / TPB;
    }
    void setcNumBlocks(){
        cNumBlocks = (nCells + TPB -1) / TPB;
    }
    void setBytes(const state& src){
        nCells = src.x_n * src.y_n * src.z_n;
        nBytes = src.particleCount * sizeof(Particle);
        cBytes = (nCells) * sizeof(Cell);
    }
    void printgpu(){
    std::cout << nCells           <<std::endl;
    std::cout << nBytes  <<std::endl;
    std::cout << cBytes    <<std::endl;
    std::cout << TPB         <<std::endl;
    std::cout << num_blocks  <<std::endl;
    std::cout << cNumBlocks    <<std::endl;
    }
};