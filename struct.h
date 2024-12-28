#ifndef __STRUCT__
#define __STRUCT__

#define Real double 
//#define Real float

#define PI Real(3.141592653589793)
#define LEND 1

#include <sys/time.h>

typedef struct PARAMS
{
	
        Real tmax;
        Real t_elapsed;
        int t_year;
        int t_month;
        int t_day;
        int t_hour;
        int t_min;
        int t_sec;
 
        Real dt;
        Real ini_dt;
        Real dt_ev; 
        Real dt_old;
        Real mindt;
        Real maxdt;
        Real totTol;
        

        int maxNumSteps;            
        int it_save;
        int it_save_jump;
        int io_unit;

}PARAMS;


typedef struct COORD
{	
	Real * x;
	Real * z;
}COORDINATE, COORD;



typedef struct GRID
{
        Real dh;

	int ni;
	int nk;

	int ni1;
	int ni2;
	int nk1;
	int nk2;

	int nx;
	int nz;

	int nx1;
	int nx2;
	int nz1;
	int nz2;

	int i0;
	int it;

}GRID;


typedef struct DISPLACEMENT
{
	Real * Uy; 
}DISPLACEMENT;


typedef struct DIMINFO {
    int M;
    int N;
    int K;
    int fin_coa;
    double DH;
} DIMINFO;


typedef struct INNERS {
    int istart;
    int iend;
    int jstart;
    int jend;
    int kstart;
    int kend;
} INNERS;

typedef struct FAULT
{ 
        Real * V_f;
        Real * Slip; 
	Real * Ts;
	Real * Tn;
        Real * Tau_n0;
        Real * SlipShift;
}FAULT;


typedef struct MEDIUM
{
	Real * lambda;
        Real * mu;
        Real * mu_half;
        Real * mu_equi;
}MEDIUM;

typedef struct STRUCTURE
{
	Real * Vs;
	Real * Vp;
	Real * rho;
}STRUCTURE;


typedef struct RATE_STATE_LAW
{
        Real * a;
        Real * b;
        Real * f0;
        Real * V0;
        Real * L;
        Real * eta;
        Real * fw;
        Real * Vw;
        Real * psi;
        Real * dpsi;
        Real * Vini;
}RATE_STATE_LAW;


typedef struct ADAPTIVE_FACTOR
{

        Real * A_xi;
        Real * B_xi;
        Real * chi;
        Real * xi;
        Real * xi_1;
        Real * xi_2;
        Real * dt_temp;

}ADAPTIVE_FACTOR;


typedef struct Runge_Kutta
{
        Real * V_f;
        Real * dpsi;
}Runge_Kutta;

typedef struct VARIABLE
{
        Real * Slip;
        Real * psi;
}VARIABLE;

inline double seconds()
{
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}



#endif //__STRUCT__
