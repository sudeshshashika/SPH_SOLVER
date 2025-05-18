#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <chrono>
#include <cuda_runtime.h>
#include "Particle.h"
#include "Cell.h"
#include "device.h"
#include "state.h"
#include "readWrite.h"
#include "kernels.h"

int main(){

    std::string path_par {"input/TANK3DDouble.par"};
    std::string path_in  {"input/TANK3DDouble.in"};
    std::string out_path {"output/"}; 
    double time{0};                                        
    int ITERATION{0};   
    state D; 
    device gpu;
    D.setState(path_par);
    D.setParticleCount(path_in);
    D.setDomainVolume(2.0); //1.0
    // D.setboundaryLimits(0,6,0,6,0,6); //2.5,3.5,0,6.0,2.5,6.0
    // D.setDampingFactors(0.15,0.15);   //0.000001
    gpu.setThreadsPerBlock(32);
    gpu.setBytes(D);
    gpu.setcNumBlocks();
    gpu.setNumBlocks(D);
    Particle* P = (Particle*)malloc(gpu.nBytes);
    Cell* C     = (Cell*)malloc(gpu.cBytes); 
    reading(path_in,P); 
    Particle * d_P = nullptr;
    Cell* d_C      = nullptr;  
    cudaMalloc(&d_C, gpu.cBytes);
    cudaMalloc(&d_P,gpu.nBytes); 
    cudaMemcpy(d_P , P, gpu.nBytes, cudaMemcpyHostToDevice);
    INIT_MASS <<< gpu.num_blocks, gpu.TPB >>> (d_P, D);
    cudaMemcpy(P, d_P, gpu.nBytes, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    auto START = std::chrono::system_clock::now();
    while( time <= D.time_end){
        WriteOut(P, out_path, ITERATION, D);
        WriteVTK(P, out_path, ITERATION, D);
        cudaMemcpy(d_P, P, gpu.nBytes, cudaMemcpyHostToDevice);
        cudaMemcpy(d_C, C, gpu.cBytes, cudaMemcpyHostToDevice);
        RESET            <<< gpu.num_blocks, gpu.TPB >>> (d_P, D);
        cellInit         <<< gpu.cNumBlocks, gpu.TPB >>> (d_C, gpu);
        MoveAllParticles <<< gpu.num_blocks, gpu.TPB >>> (d_C, d_P, D);
        CONTINUITY       <<< gpu.num_blocks, gpu.TPB >>> (d_C, d_P, D);
        MOMENTUM         <<< gpu.num_blocks, gpu.TPB >>> (d_C, d_P, D);
        UPDATE           <<< gpu.num_blocks, gpu.TPB >>> (d_P, D, time);     
        cudaMemcpy(P,d_P, gpu.nBytes, cudaMemcpyDeviceToHost);
        cudaMemcpy(C,d_C, gpu.cBytes, cudaMemcpyDeviceToHost);
        cudaDeviceSynchronize();
        time += D.timestep_length;
        ++ITERATION;
    }
    auto END   = std::chrono::system_clock::now();
    std::chrono::duration<double>elepsed_time = END-START;
    printf("RUN_TIME = %f sec. \n",elepsed_time.count());
    cudaFree(d_P);
    cudaFree(d_C);
    free(P);
    free(C);
    return 0;

}
