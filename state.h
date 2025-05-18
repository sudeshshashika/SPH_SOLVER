/**State of the simulation data is degined here**/
#pragma once 

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include "Particle.h"
#include "Cell.h"

class state{

public:
       
    std::string  part_input_file;
    double       time_end;              
    double       timestep_length;      
    int          part_out_freq;         
    std::string  part_out_name_base;    
    int          vtk_out_freq;          
    std::string  vtk_out_name_base;                    
    double       g;                     // gravity                                       
    double       rho_o;                 // intial reference density
    double       mu;                    // dynamic viscosity
    double       k;                     // constatant in pressure calculation
    double       h;                     // smoothing length
    double       H;                     // initial hight when calculating hydrstatic pressure
    double       gamma;                 // constatn for state equation
    double       x_min;                 
    double       y_min;                 
    double       z_min;                 
    double       x_max;                 
    double       y_max;                 
    double       z_max;                 
    double       x_n;                   
    double       y_n;                   
    double       z_n;                   
    double       len_x;
    double       len_y;
    double       len_z;
    int          particleCount; 
    double       soundspeed;            // Soundsped when calculating state equation 
    double       domainVolume;
    double       x_neg;
    double       y_neg;
    double       z_neg;
    double       x_plus;
    double       y_plus;
    double       z_plus;
    double       vel_damp;    
    double       d_damp;          

void setState(std::string infile1Name){
    std::string temp; 
    std::fstream infile1(infile1Name);
    infile1 >> temp >> part_input_file;
    infile1 >> temp >> time_end;              
    infile1 >> temp >> timestep_length;      
    infile1 >> temp >> part_out_freq;         
    infile1 >> temp >> part_out_name_base;    
    infile1 >> temp >> vtk_out_freq;          
    infile1 >> temp >> vtk_out_name_base;          
    infile1 >> temp >> g;     
    infile1 >> temp >> rho_o;
    infile1 >> temp >> mu;
    infile1 >> temp >> k;
    infile1 >> temp >> h;
    infile1 >> temp >> H;
    infile1 >> temp >> gamma;
    infile1 >> temp >> x_min; 
    infile1 >> temp >> y_min; 
    infile1 >> temp >> z_min; 
    infile1 >> temp >> x_max; 
    infile1 >> temp >> y_max; 
    infile1 >> temp >> z_max; 
    infile1 >> temp >> x_n; 
    infile1 >> temp >> y_n; 
    infile1 >> temp >> z_n; 
    infile1.close();
    len_x      = abs(x_max - x_min) / x_n;
    len_y      = abs(y_max - y_min) / y_n;
    len_z      = abs(z_max - z_min) / z_n;
    soundspeed = sqrt(2*abs(g) * H) *10; 
}
void setParticleCount(std::string infile2Name){
    std::fstream infile2(infile2Name);
    infile2 >> particleCount;
    infile2.close();
    
}
void setDomainVolume(double voulume){
    domainVolume = voulume;
}
void setDampingFactors(double vel_damping, double pos_damping){
    vel_damp = vel_damping;
    d_damp   = pos_damping;
}
void setboundaryLimits(double x, double xx,double y,double yy,double z,double zz){
    x_neg  = x; 
    y_neg  = y;
    z_neg  = z;
    x_plus = xx;
    y_plus = yy;
    z_plus = zz;
}

void printState(){
    std::cout << part_input_file    << std::endl;        
    std::cout << time_end           << std::endl;  
    std::cout << timestep_length    << std::endl;  
    std::cout << part_out_freq      << std::endl;  
    std::cout << part_out_name_base << std::endl;  
    std::cout << vtk_out_freq       << std::endl;  
    std::cout << vtk_out_name_base  << std::endl;  
    std::cout << g                  << std::endl;  
    std::cout << rho_o              << std::endl;
    std::cout << mu                 << std::endl;
    std::cout << k                  << std::endl;
    std::cout << h                  << std::endl;
    std::cout << H                  << std::endl;
    std::cout << gamma              << std::endl;
    std::cout << x_min              << std::endl;  
    std::cout << y_min              << std::endl;  
    std::cout << z_min              << std::endl;  
    std::cout << x_max              << std::endl;  
    std::cout << y_max              << std::endl;  
    std::cout << z_max              << std::endl;  
    std::cout << x_n                << std::endl;  
    std::cout << y_n                << std::endl;  
    std::cout << z_n                << std::endl;  
    std::cout << len_x              << std::endl;
    std::cout << len_y              << std::endl;
    std::cout << len_z              << std::endl;
    std::cout << particleCount      << std::endl;
    std::cout << soundspeed         << std::endl;     
}
};

