#ifndef __FUNCTIONS__
#define __FUNCTIONS__
//cublasXt.h
#if defined(GPU)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    // #include <thrust/universal_vector.h>
    #include "helper_cuda.h"
    #include "helper_string.h"
#endif
void alloc( GRID grid, RATE_STATE_LAW * RSL, FAULT * fault );
void alloc_RK43( GRID * grid, VARIABLE * k1, VARIABLE * k2, VARIABLE * k3, VARIABLE * k4, VARIABLE * k5, VARIABLE * k6, VARIABLE * y3, VARIABLE * y4, Runge_Kutta * f1, Runge_Kutta * f2, Runge_Kutta * f3, Runge_Kutta * f4, Runge_Kutta * f5, Runge_Kutta * f6 );
void cal_grid( GRID grid, COORD coord );
double minval( GRID * grid, ADAPTIVE_FACTOR AF );
double maxval( GRID * grid, FAULT fault, Runge_Kutta rk4, int i );
int power( int x, int y );
void fault_init( GRID grid, COORD coord, RATE_STATE_LAW RSL, FAULT fault, DISPLACEMENT U, PARAMS * params, STRUCTURE structure, MEDIUM medium );
void eval_stress_on_fault( GRID * grid, DISPLACEMENT U, PARAMS * params, FAULT fault, MEDIUM medium );
double partial1st( GRID * grid, DISPLACEMENT U, long long index );
void NRsearch( GRID * grid, FAULT fault, RATE_STATE_LAW RSL, VARIABLE y, Runge_Kutta rk4 );
void ode_init( GRID * grid, COORD * coord, DISPLACEMENT * U, DISPLACEMENT * U2, DISPLACEMENT * U3, DISPLACEMENT * U4, MEDIUM * medium, STRUCTURE * structure );
double interp_2d( double x[], double y[], double z[][4], int ni, int nj, double xi, double yi );
double interp2_2d( double x[], double y[], double z[][3], int ni, int nj, double xi, double yi );
void SetMedium( GRID grid, COORD coord, MEDIUM medium, STRUCTURE structure );
void ode_solve( int method, GRID * grid, DISPLACEMENT U, Real * u_host, Real * f_host, Real * u_dev, Real * f_dev, Real * mu_dev, Real * delta_res, cublasHandle_t handle );
void extend_crew( GRID * grid, DISPLACEMENT U );
void SOR( GRID * grid, DISPLACEMENT U );
void multigrid( GRID * grid, DISPLACEMENT U, Real * u_host, Real * f_host, Real * u_dev, Real * f_dev, Real * mu_dev, Real * delta_res, cublasHandle_t handle );
inline Real *allocMem( const int M, const int N, const int K );
inline Real *allocMem_dev(const int M, const int N, const int K);
void MGredblackgs_2D( Real *u, Real *f, const DIMINFO diminfo, const INNERS inners );
void smooth( Real *u, Real *f, Real *mu, const int level );
template <int phase>
__global__ void rbgs_kernel( Real * __restrict__ u, Real * __restrict__ f, Real * __restrict__ mu, const int level);
__global__ void residual_kernel( Real * __restrict__ r, Real * __restrict__ u, Real * __restrict__ f, Real * __restrict__ mu, const int level );
__global__ void restriction_kernel( Real * __restrict__ r2, Real * __restrict__ r, Real * __restrict__ e2, const int level );
__global__ void prolongate_kernel( Real * __restrict__ u, Real * __restrict__ e2, const int level );
__global__ void calcuDelta_dev( Real *u, Real *f, Real *mu, Real *res, const int level );
void vcycle( Real *u, Real *f, Real *mu, const int level );
void setBoundary( GRID * grid, DISPLACEMENT U, FAULT fault, RATE_STATE_LAW RSL, VARIABLE y, PARAMS * params, double ratio );
void cal_dpsi( GRID * grid, RATE_STATE_LAW RSL, FAULT fault, VARIABLE y, Runge_Kutta rk4 );
double computeError( GRID * grid, VARIABLE y3, VARIABLE y4 );
double computeStepSize( GRID * grid, PARAMS * params, const double totErr, double totTol, double errA[2] );
void RK43( int method, GRID * grid, PARAMS * params, DISPLACEMENT U, FAULT fault, MEDIUM medium, RATE_STATE_LAW RSL );
void setparams( PARAMS * params, GRID * grid );


#endif //__FUNCTIONS__
