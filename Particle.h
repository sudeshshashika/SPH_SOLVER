/**used to create the particle class**/
#pragma once 
#include <iostream>

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

class Particle{
    friend std::ostream &operator<<(std::ostream& os, const Particle &obj){
        os << obj.m << "     " << obj.rho << "     " << obj.f_id  << "     " <<  "     " << obj.d[0]    << " " << obj.d[1]    << " " << obj.d[2] 
                                                                                         << "     " << obj.v[0]    << " " << obj.v[1]    << " " << obj.v[2]
                                                                                         << "     " << obj.a_p[0]  << " " << obj.a_p[1]  << " " << obj.a_p[2]
                                                                                         << "     " << obj.a_v[0]  << " " << obj.a_v[1]  << " " << obj.a_v[2]
                                                                                         << "     " << obj.a_g[0]  << " " << obj.a_g[1]  << " " << obj.a_g[2]
                                                                                         << "     " << obj.a_v0[0] << " " << obj.a_v0[1] << " " << obj.a_v0[2];
        return os;
    }
public:
    double m;        //mass
    int    id;       //id
    int    f_id;     //fluid particle id f_id = -2 ---> fluid ; f_id = -3 ---> wall
    double rho;      //density
    int    p_gid;    //global id
    double rho_dt;   // density gradientS
    double d[3];     //position
    double v[3];     //velocity
    double a_p[3];   //pressure 
    double a_v[3];   //viscous
    double a_g[3];   //gravitational
    double a_v0[3];  //artificial viscousity
    double rf[3];    //repulsive force
    
    CUDA_CALLABLE_MEMBER Particle(double m=1.0, int id = 0,int f_id = -2, double rho = 0 , int p_gid = 0, double rho_dt = 0, double pressure= 0.0): m{m}, id{id},f_id{f_id},
                                                                                                                                                    rho{rho}, p_gid{p_gid}, rho_dt{rho_dt}  
    {
        for(int i = 0; i < 3; ++i){
            d[i]   = 0.0;
            v[i]   = 0.0;
            a_p[i] = 0.0;
            a_v[i] = 0.0;
            a_g[i] = 0.0;
            a_v0[i]= 0.0;
            rf[i]  = 0.0;
        }
    }
    CUDA_CALLABLE_MEMBER Particle(const Particle& src): m{src.m}, id{src.id},f_id{src.f_id},
                                                        rho{src.rho}, p_gid{src.p_gid},rho_dt{src.rho_dt}
    {
        for(int i =0; i < 3; ++i){
            d[i]   = 0.0;
            v[i]   = 0.0;
            a_p[i] = 0.0;
            a_v[i] = 0.0;
            a_g[i] = 0.0;
            a_v0[i]= 0.0;
            rf[i]  = 0.0;
        }
    }
    CUDA_CALLABLE_MEMBER Particle& operator=(const Particle& src){
        if(this == &src)
            return *this;
        m     = src.m;      
        id    = src.id;
        f_id  = src.f_id;    
        rho   = src.rho;     
        p_gid = src.p_gid;   
        rho_dt= src.rho_dt;
        for(int i =0; i < 3; ++i){
            d[i]   = 0.0;
            v[i]   = 0.0;
            a_p[i] = 0.0;
            a_v[i] = 0.0;
            a_g[i] = 0.0;
            a_v0[i]= 0.0;
            rf[i]  = 0.0;
        }
        return *this;           
    }
    CUDA_CALLABLE_MEMBER ~Particle(){}
};
