#pragma once 
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <limits>
#include "Particle.h"
#include "state.h"

//READING INPUT PARTICLE DATA
void reading(const std::string &path, Particle* Particles){
    std::ifstream infile3{path};
    int Number_of_particles;
    infile3 >> Number_of_particles;
    for(int i =0; i< Number_of_particles; ++i){
        infile3 >> Particles[i].m >> Particles[i].f_id >> Particles[i].d[0] >> Particles[i].d[1] >> Particles[i].d[2] 
                                                       >> Particles[i].v[0] >> Particles[i].v[1] >> Particles[i].v[2] ;
    }
    infile3.close();
}
//WRITE .OUT FILES
void WriteOut(Particle* P, std::string path, int ITERATION, const state D){
    
    std::string part_out_PATH  = path;                
    std::string file_extension = ".out";
    if( (ITERATION % D.part_out_freq)==0){
        std::ofstream outfile{part_out_PATH+D.part_out_name_base + std::to_string(ITERATION/D.part_out_freq)+ file_extension};
        outfile << D.particleCount <<"\n";
        for(int i = 0; i < D.particleCount; ++i){
            outfile<< std::fixed<< P[i].m << " " << P[i].d[0]  << " " << P[i].d[1]   << " " << P[i].d[2]
                                          << " " << P[i].v[0]  << " " << P[i].v[1]    << " " << P[i].v[2] <<"\n";
        }
        outfile.close();   
    }
}
//WRITE .VTK FILES
void WriteVTK(Particle* P, std::string path,int ITERATION, const state D){

    std::string part_out_PATH  = path;               
    std::string file_extension = ".vtk";
    if( (ITERATION % D.vtk_out_freq)==0){
        std::ofstream outfile{part_out_PATH+ D.vtk_out_name_base + std::to_string(ITERATION/D.vtk_out_freq)+ file_extension};
        outfile << "# vtk DataFile Version 4.0\n";
        outfile << "hesp visualization file\n";
        outfile << "ASCII\n";
        outfile << "DATASET UNSTRUCTURED_GRID\n";
        outfile << "POINTS " << D.particleCount << " double" <<"\n";
        for(int i = 0; i < D.particleCount; ++i){
            // if(P[i].f_id == -2){
            outfile<< std::fixed<< P[i].d[0] << " " << P[i].d[1]  << " " << P[i].d[2] << "\n";
            // }
        }
        outfile << "CELLS " << "0 " << "0\n";
        outfile << "CELL_TYPES " << "0\n";
        outfile << "POINT_DATA " << D.particleCount << "\n";
        outfile << "SCALARS rho double\n";
        outfile << "LOOKUP_TABLE default\n";
        for(int i = 0; i < D.particleCount; ++i){
            // if(P[i].f_id == -2){
            outfile<< std::fixed<< P[i].rho<< "\n";
            // }
        }
        outfile << "SCALARS pressure double\n";
        outfile << "LOOKUP_TABLE default\n";
        for(int i = 0; i < D.particleCount; ++i){
            outfile<< std::fixed<< D.k * (P[i].rho - D.rho_o)<< "\n";
        }
        outfile << "SCALARS flag int\n";
        outfile << "LOOKUP_TABLE default\n";
        for(int i = 0; i < D.particleCount; ++i){
            outfile<< std::fixed<< P[i].f_id<< "\n";
        }
        outfile << "VECTORS v double\n";
        for(int i = 0; i < D.particleCount; ++i){
            // if(P[i].f_id == -2){
            outfile<< std::fixed<< P[i].v[0] << " " << P[i].v[1]  << " " << P[i].v[2] << "\n";
            // }
        }
        outfile.close();
    }
}