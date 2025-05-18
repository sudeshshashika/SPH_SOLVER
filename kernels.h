/*All the kernels define here**/
#pragma once 
#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream>
#include <string>
#include <fstream>
#include "Particle.h"
#include "Cell.h"
#include "state.h"
#include "device.h"

//DENSITY FROM PRESSURE
CUDA_CALLABLE_MEMBER double DensityfromState(Particle &A, double rho_o, double gama, double soundspeed, double gravity, double H){
    double B       = (rho_o * soundspeed * soundspeed) / gama;
    double val     = (rho_o * abs(gravity) * (H - A.d[1]) / B) + 1;
    double density = rho_o * ( pow(val, ( 1 / gama )) );
    return density;
}
//INITIALIZING MASS AND DENSITIES OF PARTICLES
__global__ void INIT_MASS(Particle* p , const state D){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        // double domainvolume = 1.0;
        if(p[idx].f_id == -2){
        p[idx].m  = D.rho_o * D.domainVolume / D.particleCount; 
        } 
        p[idx].rho  = D.rho_o;//DensityfromState(p[idx], D.rho_o, D.gamma, D.soundspeed, D.g, D.H);
    }
}
//RESET PARTICLE DATA
__global__ void RESET(Particle* p, const state D){//int N){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        p[idx].rho    = 0;
        p[idx].rho_dt = 0.0;
        for(int DIM = 0; DIM < 3; ++DIM){
            p[idx].a_p[DIM] = 0.0;
            p[idx].a_v[DIM] = 0.0;
        }
    }
}
//INITILIZING THE CELL LIST INDEX
__global__ void cellInit(Cell* c, const device gpu){
    long long idx = blockIdx.x * blockDim.x  + threadIdx.x;
    if(idx < gpu.nCells){
        c[idx].x = -1.0;
    }
}
//MOVE PARTICLES INDEXES WITH CELL INDEXES
__global__ void MoveAllParticles(Cell*c, Particle* P, const state D){
    long long idx = blockIdx.x * blockDim.x + threadIdx.x;   
    if(idx < D.particleCount){
        P[idx].id = idx;
        int x = (P[idx].d[0] - D.x_min) / D.len_x;
        int y = (P[idx].d[1] - D.y_min) / D.len_y;
        int z = (P[idx].d[2] - D.z_min) / D.len_z;
        int index = x * (D.y_n * D.z_n) + y * D.z_n + z;
        int val = atomicExch(&c[index].x , idx); 
        P[idx].id = val;
    }  
}
//DENSITY VARIATION OF ONE PARTICLE
CUDA_CALLABLE_MEMBER void continuity(Particle&  a, Particle& b, double h){
    double r2 = 0.0;
    double r_ij[3]{0.0,0.0,0.0};
    for(int dim = 0; dim < 3; ++dim){
        r_ij[dim] = a.d[dim] - b.d[dim];
        r2        = r2 + r_ij[dim] * r_ij[dim];
    }
    double h2 = h*h;
    double h9 = h*h*h*h*h*h*h*h*h;
    double wPoly6 = 315.0 / (64.0*3.14159*h9); 
    double r = sqrt(r2);
    if(r <  h){
        a.rho = a.rho + b.m * wPoly6 * (h2-(r*r))*(h2-(r*r))*(h2-(r*r)); 
    }
}
//DENSITY VARIATION OF ALL PARTICLES == CONTINUITY EQUATION
__global__ void CONTINUITY(Cell* c,Particle* p, const state D){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){  
        int X = (p[idx].d[0] - D.x_min) / D.len_x;
        int Y = (p[idx].d[1] - D.y_min) / D.len_y;
        int Z = (p[idx].d[2] - D.z_min) / D.len_z;
        for(int i = X-1; i <= X+1; ++i){
            for(int j = Y-1; j <= Y+1; ++j){
                for(int k = Z-1; k <= Z+1; ++k){
                    if(i>=0 && i<D.x_n && j>=0 && j<D.y_n && k>=0 && k<D.z_n){
                        int  ncell_ID            = k + (D.z_n * j ) + (D.z_n * D.y_n ) * i;
                        int HeadneighborParticle = c[ncell_ID].x;
                        while(HeadneighborParticle != -1){
                            if(idx != HeadneighborParticle){
                                continuity(p[idx], p[HeadneighborParticle],D.h); 
                            } 
                        HeadneighborParticle = p[HeadneighborParticle].id;
                        }
                    }
                }
            }
        }
    }
}
//PRESSURE FORCE AND VISCOUS FORCE CACLCULATION FOR ONE PARTICLE == MOMENTUM EQUATION
CUDA_CALLABLE_MEMBER void momentum(Particle& a, Particle& b, double h, double rho_o, double k, double mu){
    double pa  = k * (a.rho - rho_o);
    double pb  = k * (b.rho - rho_o);
    double r2  = 0.0;
    double r_ij[3]{0.0,0.0,0.0};
    for(int dim = 0; dim < 3; ++dim){
        r_ij[dim] = a.d[dim] - b.d[dim];
        r2        = r2 + r_ij[dim] * r_ij[dim];
    }
    double r  = sqrt(r2);
    double h6 = h*h*h*h*h*h;
    double wSpikyGradient = -45.0 / (3.14159 * h6);
    double wViscosityLaplacian = 45.0 / (3.14159 * h6);
    if(r <  h){
        for(int DIM = 0; DIM < 3; ++DIM){
            a.a_p[DIM] = a.a_p[DIM] - b.m * (0.5 * ( pa + pb ) / b.rho) * wSpikyGradient*(h-r)*(h-r)*(r_ij[DIM] / r);
            if(a.f_id == -2){//}
            a.a_v[DIM] = a.a_v[DIM] + (mu * b.m * (b.v[DIM] - a.v[DIM]) / b.rho) * wViscosityLaplacian * ( h - r);
            }
        }
    }
}
//TOTAL PRESSURE AND VISCOUS FORCE CALCULATION
__global__ void MOMENTUM(Cell*c,Particle* p, const state D){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        int X = (p[idx].d[0] - D.x_min) / D.len_x;
        int Y = (p[idx].d[1] - D.y_min) / D.len_y;
        int Z = (p[idx].d[2] - D.z_min) / D.len_z;
        for(int i = X-1; i <= X+1; ++i){
            for(int j = Y-1; j <= Y+1; ++j){
                for(int k = Z-1; k <= Z+1; ++k){
                    if(i>=0 && i<D.x_n && j>=0 && j<D.y_n && k>=0 && k<D.z_n){
                        int  ncell_ID = k + (D.z_n * j ) + (D.z_n * D.y_n ) * i;
                        int HeadneighborParticle = c[ncell_ID].x;
                        while(HeadneighborParticle != -1){
                            if(idx != HeadneighborParticle){
                                momentum(p[idx],p[HeadneighborParticle],D.h,D.rho_o,D.k,D.mu);
                            } 
                            HeadneighborParticle = p[HeadneighborParticle].id;
                        }
                    }
                }
            }
        }
    }
}
//REFLECTIVE BOUNDARY CONDITION WITH DAMPING FACTOR
CUDA_CALLABLE_MEMBER void boundarCondition(Particle& A, double x_plus,double x_neg, double y_plus,double y_neg, double z_plus,double z_neg, double vel_damp,double d_damp){
	if (A.d[0] > x_plus){ 
		A.v[0] = -A.v[0]*vel_damp; 
		A.d[0] = A.d[0] - d_damp * (A.d[0] - x_plus);
	}
	if(A.d[0] < x_neg){
		A.v[0] = -A.v[0]*vel_damp;
		A.d[0] =  A.d[0] - d_damp * (A.d[0] - x_neg);
	}

	if(A.d[1] > y_plus){
		A.v[1] = -A.v[1]*vel_damp;
		A.d[1] =A.d[1] - d_damp * (A.d[1] - y_plus);
	}

	if(A.d[1] < y_neg){
		A.v[1] = -A.v[1] *vel_damp;
		A.d[1] = A.d[1] - d_damp * (A.d[1] - y_neg);
	}

	if(A.d[2] > z_plus){
		A.v[2] = -A.v[2]*vel_damp;
		A.d[2] = A.d[2] - d_damp * (A.d[2] - z_plus);
	}
	if(A.d[2] < z_neg){
		A.v[2] = -A.v[2]*vel_damp;
		A.d[2] = A.d[2] - d_damp * (A.d[2] - z_neg);
	}
}
//UPDATE THE POSITION AND VELOCITIES
__global__ void UPDATE(Particle *p, const state D,double time){
    long long idx  = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx < D.particleCount){
        double g[3]{0.0, D.g, 0.0 };
        for(int dim = 0; dim< 3; ++dim){
            /**please comment out 192-197 if you use outline input when you want to visualize the outline wall**/
            if(p[idx].f_id == -2){
            p[idx].v[dim] = p[idx].v[dim] + D.timestep_length * (((p[idx].a_p[dim] + p[idx].a_v[dim]) / p[idx].rho) +g[dim]);
            p[idx].d[dim] = p[idx].d[dim] + D.timestep_length * p[idx].v[dim];
            //boundarCondition(p[idx],D.x_plus,D.x_neg,D.y_plus,D.y_neg,D.z_plus,D.z_neg,D.vel_damp,D.d_damp; //this is used for reflective and damped boundary condition uncomment if you wish to 
                                                                                                              // use this for reflective input files(Then you don't need static particles)
            }
        }
        /**Sloshing action is defined here.
         * intially tank is moved in Z direction with 0.5 m/s unitl time = 0.5 s,
         * then tank velocity is zero between 0.5s and 1.5 s,
         * then tank is moved in -Z direction with -0.5 m/s until time =2.5 s, 
         * then tank velocity becomes zero for the rest of the time**/
        if(time <= 0.5){
            if(p[idx].f_id == -3){
                for(int dim = 0; dim< 3; ++dim){
                    p[idx].v[2]   = 0.5;
                    p[idx].d[dim] = p[idx].d[dim] + D.timestep_length * p[idx].v[dim];
                }
            }
        }
        if(time > 0.5 && time <= 1.5){
            if(p[idx].f_id == -3){
                for(int dim = 0; dim< 3; ++dim){
                    p[idx].v[2]   = 0;
                    p[idx].d[dim] = p[idx].d[dim] + D.timestep_length * p[idx].v[dim];
                }
            }
        }
        if(time > 1.5 && time <= 2.5){
            if(p[idx].f_id == -3){
                for(int dim = 0; dim< 3; ++dim){
                    p[idx].v[2]   = -0.5;
                    p[idx].d[dim] = p[idx].d[dim] + D.timestep_length * p[idx].v[dim];
                }
            }
        }
        if(time > 2.5){
            if(p[idx].f_id == -3){
                for(int dim = 0; dim< 3; ++dim){
                    p[idx].v[2]   = 0;
                    p[idx].d[dim] = p[idx].d[dim] + D.timestep_length * p[idx].v[dim];
                }
            }
        }  

    }

}
