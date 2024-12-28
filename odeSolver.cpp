#include "header.h"
const double relax_rbgs = 1.00;
const int n_levels = 6; // 0, 1, 2, 3--total 4; level 3 is the base solver
const int n_cycles = 32;

Real *res [n_levels - 1];
Real *res2[n_levels - 1];
Real *err2[n_levels - 1];
DIMINFO dimsinfo[n_levels];
INNERS innersinfo[n_levels];
__constant__ DIMINFO dimsinfo_dev[n_levels];
__constant__ INNERS innersinfo_dev[n_levels];

typedef void (*MGIterSolver_t)(Real *u, Real *f, const DIMINFO diminfo, const INNERS inners);
MGIterSolver_t MGIterSolver;


void ode_init( GRID * grid, COORD * coord, DISPLACEMENT * U, DISPLACEMENT * U2, DISPLACEMENT * U3, DISPLACEMENT * U4, MEDIUM * medium, STRUCTURE * structure )
{

	int i;
        int k;

        grid->nx = grid->ni + 2*LEND;
        grid->nz = grid->nk + 2*LEND;

        grid->nx1 = 1;
        grid->nx2 = grid->nx1 + grid->nx - 1;
        grid->nz1 = 1;
        grid->nz2 = grid->nz1 + grid->nz - 1;

        grid->ni1 = grid->nx1 + LEND;
        grid->ni2 = grid->ni1 + grid->ni - 1;
        grid->nk1 = grid->nz1 + LEND;
        grid->nk2 = grid->nk1 + grid->nk - 1;

        grid->i0 = grid->ni1;
        grid->it = 0;


        int nx = grid->nx;
        int nz = grid->nz;
        int ni = grid->ni;
        int nk = grid->nk;

        long long num1 = nx * nz;
        long long num2 = (2*nx-1) * (2*nz-1);

        Real * pCoord = NULL;
        Real * pUy = NULL; 
        Real * pU2y = NULL; 
        Real * pU3y = NULL; 
        Real * pU4y = NULL; 
        Real * pMedium = NULL;
        Real * pStructure = NULL;

        long long size_Coord = sizeof( Real ) * num1 * 2;
        long long size_Uy = sizeof( Real ) * num1;
        long long size_Medium = sizeof( Real ) * ( num1 * 2 + num2 * 2 );
        long long size_Structure = sizeof( Real ) * num1 * 3;

	pCoord = ( Real * )malloc( size_Coord );
        pUy = ( Real * )malloc( size_Uy );
        pU2y = ( Real * )malloc( size_Uy );
        pU3y = ( Real * )malloc( size_Uy );
        pU4y = ( Real * )malloc( size_Uy );
        pMedium = ( Real * )malloc( size_Medium );
        pStructure = ( Real * )malloc( size_Structure );

	memset(  pCoord, 0, size_Coord );
	memset(  pUy, 0, size_Uy );
	memset(  pU2y, 0, size_Uy );
	memset(  pU3y, 0, size_Uy );
	memset(  pU4y, 0, size_Uy );
	memset(  pMedium, 0, size_Medium );
	memset(  pStructure, 0, size_Structure );


        coord->x = pCoord;
        coord->z = pCoord + num1;

        U->Uy = pUy;
        U2->Uy = pU2y;
        U3->Uy = pU3y;
        U4->Uy = pU4y;

        medium->lambda = pMedium;
        medium->mu = pMedium + num1;
        medium->mu_half = pMedium + num1 * 2;
        medium->mu_equi = pMedium + num1 * 2 + num2;

        structure->Vs = pStructure;
        structure->Vp = pStructure + num1;
        structure->rho = pStructure + num1 * 2;
//~~~~~~~~~~~~~~~~~~~~~~input initial U~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//   
        Real * Uy0 = NULL;
        num2 = ni * nk;
        long long size_Uy0 = sizeof( Real ) * num2;
        Uy0 = ( Real * )malloc( size_Uy0 );
	memset(  Uy0, 0, size_Uy0 );
        FILE * fp;
        char fileName[256];
        sprintf( fileName, "/Nas/sunxu/code_2D_Cart_cpp/src/Basin88_800_order4.bin" );
        fp = fopen( fileName, "rb" );
        size_t sizeRead = fread( Uy0, sizeof( Real ), num2, fp );
        if( sizeRead != num2 )
        {
           printf( "Read error!\n" );
        }
        fclose( fp );
      
        int ni1 = grid->ni1;
        int ni2 = grid->ni2;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        long long index1 = 0;
        long long index2 = 0;

        for ( k = nk1 - 1; k < nk2; k ++ ) 
        {
           for ( i = ni1 - 1; i < ni2; i ++ )
           {
              index1 = i + k * nx;
              index2 = i - ni1 + 1 + ( k - nk1 + 1 ) * ni;
              U->Uy[index1] = Uy0[index2];
           }
        }  
        free(Uy0);

//~~~~~~~~~~~~~~~~~~~~~~input initial U~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//   

}

double interp_2d( double x[], double y[], double z[][4], int ni, int nj, double xi, double yi )
{
        int i, j, k;
        double zi;
        double Lx[ni], xt[ni], xb[ni];
        double Ly[nj], yt[nj], yb[nj];
        for ( k = 0; k < ni; k++ )
        {
           Lx[k] = 1.0;
        }
        for ( i = 0; i < ni; i++ )
        {
           for ( k = 0; k < ni; k++ )
           {
              xb[k] = x[k];
              xb[k] = x[i] - xb[k];
           }
           xb[i] = 1.0;
           for ( k = 0; k < ni; k++ )
           {
              xt[k] = x[k];
              xt[k] = xi - xt[k];
           }
           xt[i] = 1.0;
           for ( j = 0; j < ni; j++ )
           {
              Lx[i] = Lx[i]*xt[j]/xb[j];
           }
        }


        for ( k = 0; k < nj; k++ )
        {
           Ly[k] = 1.0;
        }
        for ( i = 0; i < nj; i++ )
        {
           for ( k = 0; k < nj; k++ )
           {
              yb[k] = y[k];
              yb[k] = y[i] - yb[k];
           }
           yb[i] = 1.0;
           for ( k = 0; k < nj; k++ )
           {
              yt[k] = y[k];
              yt[k] = yi - yt[k];
           }
           yt[i] = 1.0;
           for ( j = 0; j < nj; j++ )
           {
              Ly[i] = Ly[i]*yt[j]/yb[j];
           }
        }


        zi = 0.0;
        for ( j = 0; j < nj; j++ )
        {
           for ( i = 0; i < ni; i++ )
           {
              zi = zi + Lx[i]*Ly[j]*z[i][j];
           }
        }    

        return zi;
}

double interp2_2d( double x[], double y[], double z[][3], int ni, int nj, double xi, double yi )
{
        int i, j, k;
        double zi;
        double Lx[ni], xt[ni], xb[ni];
        double Ly[nj], yt[nj], yb[nj];
        for ( k = 0; k < ni; k++ )
        {
           Lx[k] = 1.0;
        }
        for ( i = 0; i < ni; i++ )
        {
           for ( k = 0; k < ni; k++ )
           {
              xb[k] = x[k];
              xb[k] = x[i] - xb[k];
           }
           xb[i] = 1.0;
           for ( k = 0; k < ni; k++ )
           {
              xt[k] = x[k];
              xt[k] = xi - xt[k];
           }
           xt[i] = 1.0;
           for ( j = 0; j < ni; j++ )
           {
              Lx[i] = Lx[i]*xt[j]/xb[j];
           }
        }


        for ( k = 0; k < nj; k++ )
        {
           Ly[k] = 1.0;
        }
        for ( i = 0; i < nj; i++ )
        {
           for ( k = 0; k < nj; k++ )
           {
              yb[k] = y[k];
              yb[k] = y[i] - yb[k];
           }
           yb[i] = 1.0;
           for ( k = 0; k < nj; k++ )
           {
              yt[k] = y[k];
              yt[k] = yi - yt[k];
           }
           yt[i] = 1.0;
           for ( j = 0; j < nj; j++ )
           {
              Ly[i] = Ly[i]*yt[j]/yb[j];
           }
        }


        zi = 0.0;
        for ( j = 0; j < nj; j++ )
        {
           for ( i = 0; i < ni; i++ )
           {
              zi = zi + Lx[i]*Ly[j]*z[i][j];
           }
        }    

        return zi;
}


 
void SetMedium( GRID grid, COORD coord, MEDIUM medium, STRUCTURE structure )
{
        int i;
        int k;
        int nx = grid.nx;
        int nz = grid.nz;
        int ni = grid.ni;
        int nk = grid.nk;
        int ni1 = grid.ni1;
        int ni2 = grid.ni2;
        int nk1 = grid.nk1; 
        int nk2 = grid.nk2; 
        Real temp;
        long long index = 0;
        long long index2 = 0;
        long long index3 = 0;
        long long index4 = 0;
        long long index5 = 0;
        long long index6 = 0;
        long long index7 = 0;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Basin Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

        Real mu_out, mu_in, rho_out, rho_in;
        Real W, D, c, r, rr, rw;
        
        mu_out = 36.0;
        mu_in = 8.0;
        rho_out = 2.8000;
        rho_in = 2.0000;
        W = 24.0;
        D = 8.0;
        c = 0.5 * W / D;
        rr = 0.25 * W * W;
        rw = 1.0 + 0.5 * W / D; 

        for ( k = nk1 - 1; k < nk2; k ++ ) 
        {
           for ( i = ni1 - 1; i < ni2; i ++ )
           {
              index = i + k * nx;
              r = coord.x[index] * coord.x[index] + c * c * coord.z[index] * coord.z[index];
              temp = (r - rr) / rw;
              medium.mu[index] = 0.5 * (mu_out-mu_in) * ( tanh(temp) + 1.0 ) + mu_in; 
              structure.rho[index] = 0.5 * (rho_out-rho_in) * ( tanh(temp) + 1.0 ) + rho_in;
              structure.Vs[index] = 3.4640; 
              structure.Vp[index] = 6.0000; 
           }
        }

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Basin Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FILE * fp;
        char fileName[256];
        sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/rho.bin" );
        fp = fopen( fileName, "wb" );
        fwrite( structure.rho, sizeof( Real ), nx*nz, fp );
        fclose( fp );

        sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/mu.bin" );
        fp = fopen( fileName, "wb" );
        fwrite( medium.mu, sizeof( Real ), nx*nz, fp );
        fclose( fp );
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       
        for ( k = nk1 - 1; k < nk2; k ++ ) 
        {
           for ( i = ni1 - 1; i < ni2; i ++ )
           {
              index = i + k * nx;
              index2 = i*2 + (k*2) * (nx*2-1);
              medium.mu_half[index2] = medium.mu[index];
              medium.mu_equi[index2] = medium.mu[index];
           }
        }
        for ( i = ni1 - 1; i < ni2; i ++ )
        {
           index = i*2 + 0 * (nx*2-1);
           index2 = i*2 + ((nk1-1) * 2) * (nx*2-1);
           medium.mu_half[index] = medium.mu_half[index2];
           medium.mu_equi[index] = medium.mu_half[index2];
        }
        for ( i = ni1 - 1; i < ni2; i ++ )
        {
           index = i*2 + ((nz-1) * 2) * (nx*2-1);
           index2 = i*2 + ((nk2-1) * 2) * (nx*2-1);
           medium.mu_half[index] = medium.mu_half[index2];
           medium.mu_equi[index] = medium.mu_half[index2];
        }

        for ( k = 0; k < nz; k ++ )
        {
           index = 0*2 + (k*2) * (nx*2-1);
           index2 = (ni1-1)*2 + (k*2) * (nx*2-1);
           medium.mu_half[index] = medium.mu_half[index2];
           medium.mu_equi[index] = medium.mu_half[index2];
        }

        for ( k = 0; k < nz; k ++ )
        {
           index = (nx-1)*2 + (k*2) * (nx*2-1);
           index2 = (ni2-1)*2 + (k*2) * (nx*2-1);
           medium.mu_half[index] = medium.mu_half[index2];
           medium.mu_equi[index] = medium.mu_half[index2];
        }


        mu_out = 36.0;
        mu_in = 8.0;
        for ( k = 0; k < nz; k ++ ) 
        {
           for ( i = ni1 - 1; i < ni2 - 1; i ++ )
           {
              index = i*2 + (k*2) * (nx*2-1);
              index2 = (i+1)*2 + (k*2) * (nx*2-1);
              index3 = (i*2+1) + (k*2) * (nx*2-1);
              temp = 0.50 * ( medium.mu_half[index] + medium.mu_half[index2] ); 
              medium.mu_half[index3] = 1.00 / temp; 
              medium.mu_half[index3] = temp; 
           }
        }
        for ( i = ni1 - 1; i < ni2; i ++ ) 
        {
           for ( k = 0; k < nz - 1; k ++ )
           {
              index = i*2 + (k*2) * (nx*2-1);
              index2 = i*2 + ((k+1) * 2) * (nx*2-1);
              index3 = i*2 + (k*2 + 1) * (nx*2-1);
              temp = 0.50 * ( medium.mu_half[index] + medium.mu_half[index2] ); 
              medium.mu_half[index3] = 1.00 / temp; 
              medium.mu_half[index3] = temp; 
           }
        }

        int NGRIDX, NGRIDY, nsamp;
        int mi, mj;
        double mu_interp, MI, MJ, mu0, mu_equal;
        double mu_matrix[3][4];
        double x_index[3] = {-1.0, 0.0, 1.0};
        double y_index[4] = {-1.5, -0.5, 0.5, 1.5};    
   
        double mu_matrix2[4][3];
        double x_index2[4] = {-1.5, -0.5, 0.5, 1.5};
        double y_index2[3] = {-1.0, 0.0, 1.0};       
   
        NGRIDX = 10;
        NGRIDY = 10;
        nsamp = (2*NGRIDX) * (2*NGRIDY);
 
        for ( k = nk1 - 1; k < nk2; k ++ ) 
        {
           for ( i = ni1; i < ni2; i ++ )
           {
              index = (i*2-1) + (k*2) * (nx*2-1);
              index2 = (i-2)*2 + ((k-1)*2) * (nx*2-1);
              index3 = (i-2)*2 + (k*2) * (nx*2-1);
              index4 = (i-2)*2 + ((k+1)*2) * (nx*2-1);
              mu0 = 0.0;
              mu_matrix[0][0] = medium.mu_half[index2];
              mu_matrix[0][1] = medium.mu_half[index2 + 2];
              mu_matrix[0][2] = medium.mu_half[index2 + 4];
              mu_matrix[0][3] = medium.mu_half[index2 + 6];
              mu_matrix[1][0] = medium.mu_half[index3];
	      mu_matrix[1][1] = medium.mu_half[index3 + 2];
              mu_matrix[1][2] = medium.mu_half[index3 + 4];
              mu_matrix[1][3] = medium.mu_half[index3 + 6];
              mu_matrix[2][0] = medium.mu_half[index4];
              mu_matrix[2][1] = medium.mu_half[index4 + 2];
              mu_matrix[2][2] = medium.mu_half[index4 + 4];
              mu_matrix[2][3] = medium.mu_half[index4 + 6];
              for ( mj = -1 * NGRIDY; mj <= NGRIDY; mj++ )
              {
                 if( mj == 0 )
                 {
                    continue;
                 } 
                 for ( mi = -1 * NGRIDX; mi <= NGRIDX; mi++ )
                 {
                    if( mi == 0 )
                    {
                       continue;
                    } 
                    MI = mi/( (NGRIDX+1.0)*2.0 );
                    MJ = mj/( (NGRIDY+1.0)*2.0 );
                    mu_interp = interp_2d( x_index, y_index, mu_matrix, 3, 4, MI, MJ );
                    mu0 = mu0 + 1.0/mu_interp;
                 }
              }
              medium.mu_equi[index] = nsamp / mu0;
           }
        }

        for ( k = nk1; k < nk2; k ++ ) 
        {
           for ( i = ni1; i < ni2 - 1; i ++ )
           {
              index = i*2 + (k*2-1) * (nx*2-1);
              index2 = (i-1)*2 + ((k-2)*2) * (nx*2-1);
              index3 = (i-1)*2 + ((k-1)*2) * (nx*2-1);
              index4 = (i-1)*2 + (k*2) * (nx*2-1);
              index5 = (i-1)*2 + ((k+1)*2) * (nx*2-1);
              mu0 = 0.0;
              mu_matrix2[0][0] = medium.mu_half[index2];
              mu_matrix2[0][1] = medium.mu_half[index2 + 2];
              mu_matrix2[0][2] = medium.mu_half[index2 + 4];
              mu_matrix2[1][0] = medium.mu_half[index3];
              mu_matrix2[1][1] = medium.mu_half[index3 + 2];
              mu_matrix2[1][2] = medium.mu_half[index3 + 4];
              mu_matrix2[2][0] = medium.mu_half[index4];
              mu_matrix2[2][1] = medium.mu_half[index4 + 2];
              mu_matrix2[2][2] = medium.mu_half[index4 + 4];
              mu_matrix2[3][0] = medium.mu_half[index5];
              mu_matrix2[3][1] = medium.mu_half[index5 + 2];
              mu_matrix2[3][2] = medium.mu_half[index5 + 4];
              for ( mj = -1 * NGRIDY; mj <= NGRIDY; mj++ )
              {
                 if( mj == 0 )
                 {
                    continue;
                 } 
                 for ( mi = -1 * NGRIDX; mi <= NGRIDX; mi++ )
                 {
                    if( mi == 0 )
                    {
                       continue;
                    } 
                    MI = mi/( (NGRIDX+1.0)*2.0 );
                    MJ = mj/( (NGRIDY+1.0)*2.0 );
                    mu_interp = interp2_2d( x_index2, y_index2, mu_matrix2, 4, 3, MI, MJ );
                    mu0 = mu0 + 1.0/mu_interp;
                 }
              }
              medium.mu_equi[index] = nsamp / mu0;
           }
        }


        for ( i = ni1; i < ni2 - 1; i ++ )
        {
           index = i*2 + ((nk1-1)*2-1) * (nx*2-1);
           index2 = (i-1)*2 + 0 * (nx*2-1);
           index3 = (i-1)*2 + 0 * (nx*2-1);
           index4 = (i-1)*2 + ((nk1-1)*2) * (nx*2-1);
           index5 = (i-1)*2 + (nk1*2) * (nx*2-1);
           mu0 = 0.0;
           mu_matrix2[0][0] = medium.mu_half[index2];
           mu_matrix2[0][1] = medium.mu_half[index2 + 2];
           mu_matrix2[0][2] = medium.mu_half[index2 + 4];
           mu_matrix2[1][0] = medium.mu_half[index3];
           mu_matrix2[1][1] = medium.mu_half[index3 + 2];
           mu_matrix2[1][2] = medium.mu_half[index3 + 4];
           mu_matrix2[2][0] = medium.mu_half[index4];
           mu_matrix2[2][1] = medium.mu_half[index4 + 2];
           mu_matrix2[2][2] = medium.mu_half[index4 + 4];
           mu_matrix2[3][0] = medium.mu_half[index5];
           mu_matrix2[3][1] = medium.mu_half[index5 + 2];
           mu_matrix2[3][2] = medium.mu_half[index5 + 4];
           for ( mj = -1 * NGRIDY; mj <= NGRIDY; mj++ )
           {
              if( mj == 0 )
              {
                 continue;
              } 
              for ( mi = -1 * NGRIDX; mi <= NGRIDX; mi++ )
              {
                 if( mi == 0 )
                 {
                    continue;
                 } 
                 MI = mi/( (NGRIDX+1.0)*2.0 );
                 MJ = mj/( (NGRIDY+1.0)*2.0 );
                 mu_interp = interp2_2d( x_index2, y_index2, mu_matrix2, 4, 3, MI, MJ );
                 mu0 = mu0 + 1.0/mu_interp;
              }
           }
           medium.mu_equi[index] = nsamp / mu0;
        }

        for ( i = ni1; i < ni2 - 1; i ++ )
        {
           index = i*2 + (nk2*2-1) * (nx*2-1);
           index2 = (i-1)*2 + ((nk2-2)*2) * (nx*2-1);
           index3 = (i-1)*2 + ((nk2-1)*2) * (nx*2-1);
           index4 = (i-1)*2 + (nk2*2) * (nx*2-1);
           index5 = (i-1)*2 + (nk2*2) * (nx*2-1);
           mu0 = 0.0;
           mu_matrix2[0][0] = medium.mu_half[index2];
           mu_matrix2[0][1] = medium.mu_half[index2 + 2];
           mu_matrix2[0][2] = medium.mu_half[index2 + 4];
           mu_matrix2[1][0] = medium.mu_half[index3];
           mu_matrix2[1][1] = medium.mu_half[index3 + 2];
           mu_matrix2[1][2] = medium.mu_half[index3 + 4];
           mu_matrix2[2][0] = medium.mu_half[index4];
           mu_matrix2[2][1] = medium.mu_half[index4 + 2];
           mu_matrix2[2][2] = medium.mu_half[index4 + 4];
           mu_matrix2[3][0] = medium.mu_half[index5];
           mu_matrix2[3][1] = medium.mu_half[index5 + 2];
           mu_matrix2[3][2] = medium.mu_half[index5 + 4];
           for ( mj = -1 * NGRIDY; mj <= NGRIDY; mj++ )
           {
              if( mj == 0 )
              {
                 continue;
              } 
              for ( mi = -1 * NGRIDX; mi <= NGRIDX; mi++ )
              {
                 if( mi == 0 )
                 {
                    continue;
                 } 
                 MI = mi/( (NGRIDX+1.0)*2.0 );
                 MJ = mj/( (NGRIDY+1.0)*2.0 );
                 mu_interp = interp2_2d( x_index2, y_index2, mu_matrix2, 4, 3, MI, MJ );
                 mu0 = mu0 + 1.0/mu_interp;
              }
           }
           medium.mu_equi[index] = nsamp / mu0;
        }



        //FILE * fp;
        //char fileName[256];
        sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/mu_equi.bin" );
        fp = fopen( fileName, "wb" );
        fwrite( medium.mu_equi, sizeof( Real ), (2*nx-1) * (2*nz-1), fp );
        fclose( fp );

 
}



void ode_solve( int method, GRID * grid, DISPLACEMENT U, Real * u_host, Real * f_host, Real * u_dev, Real * f_dev, Real * mu_dev, Real * delta_res, cublasHandle_t handle )
{
        multigrid( grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );  
        extend_crew( grid, U ); 
}


void extend_crew( GRID * grid, DISPLACEMENT U )
{       
        int i, k, n;        
        int nx = grid->nx;
        int nx1 = grid->nx1;
        int nx2 = grid->nx2;
        int nz1 = grid->nz1;
        int nz2 = grid->nz2;
        int ni1 = grid->ni1;
        int ni2 = grid->ni2;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        long long index1 = 0;
        long long index2 = 0;

        for ( k = nz1 - 1; k < nz2; k ++ )
        {
           for ( n = 1; n < LEND + 1; n ++ )
           {
              index1 = ni1 - 1 + k * nx;
              index2 = ni2 - 1 + k * nx;
              U.Uy[index1-n] = 2.0 * U.Uy[index1] - U.Uy[index1+n];
              U.Uy[index2+n] = 2.0 * U.Uy[index2] - U.Uy[index2-n];
           }
        }

        
        for ( n = 1; n < LEND + 1; n ++ )
        {
           for ( i = nx1 - 1; i < nx2; i ++ )
           {
              index1 = i + (nk1-1) * nx;
              index2 = i + (nk2-1) * nx;
              U.Uy[index1-n*nx] = U.Uy[index1+n*nx];
              U.Uy[index2+n*nx] = U.Uy[index2-n*nx];
           }
        }

}

void multigrid( GRID * grid, DISPLACEMENT U, Real * u_host, Real * f_host, Real * u_dev, Real * f_dev, Real * mu_dev, Real * delta_res, cublasHandle_t handle )
{

        int i, j;
        int M = 1;
        int N = grid->nk;
        int K = grid->ni; 
        int nx = grid->nx; 
        int ni1 = grid->ni1;
        int ni2 = grid->ni2;
        int nk1 = grid->nk1; 
        int nk2 = grid->nk2;
        size_t num = M * N * K;
        long long idx;
        long long index;
       

        for (j = nk1 - 1; j < nk2; j++) 
        {
           for (i = ni1 - 1; i < ni2; i++) 
           {
              index = i + j * nx;
              idx = i - ni1 + 1 + (j - nk1 + 1) * K;
              u_host[idx] = U.Uy[index];
           }
        }

        checkCudaErrors( cudaMemcpy(u_dev, u_host, sizeof(Real) * num, cudaMemcpyHostToDevice) );
        checkCudaErrors( cudaMemcpy(f_dev, f_host, sizeof(Real) * num, cudaMemcpyHostToDevice) );

        dim3 block(32, 16, 1);
        dim3 Grid;
        Grid.x = (innersinfo[0].kend - innersinfo[0].kstart + block.x - 1) / block.x;
        Grid.y = (innersinfo[0].jend - innersinfo[0].jstart + block.y - 1) / block.y;
        Grid.z = (innersinfo[0].iend - innersinfo[0].istart + block.z - 1) / block.z;

        int iter = 0;
        double eps = 1.0e-9;
        double delta = 1.0;
        std::vector<double> delta_vector;

        calcuDelta_dev <<<Grid, block>>> (u_dev, f_dev, mu_dev, delta_res, 0);
        checkCudaErrors( cublasDnrm2(handle, M*N*K, delta_res, 1, &delta) );
        delta_vector.push_back(delta);
        while (delta > eps) 
        {
           vcycle(u_dev, f_dev, mu_dev, 0);
           calcuDelta_dev <<<Grid, block>>> (u_dev, f_dev, mu_dev, delta_res, 0);
           checkCudaErrors( cublasDnrm2(handle, M*N*K, delta_res, 1, &delta) );
           delta_vector.push_back(delta);

           iter += (n_cycles) * (n_levels > 1 ? (2 * n_levels - 1) : 1);        
           fflush(stdout);
        }

        checkCudaErrors( cudaMemcpy(u_host, u_dev, sizeof(Real) * num, cudaMemcpyDeviceToHost) );
        for (j = nk1 - 1; j < nk2; j++) 
        {
           for (i = ni1 - 1; i < ni2; i++) 
           {
              index = i + j * nx;
              idx = i - ni1 + 1 + (j - nk1 + 1) * K;
              U.Uy[index] = u_host[idx];
           }
        }

        return;

}

inline Real *allocMem( const int M, const int N, const int K ) 
{
        size_t num = M * N * K;
        Real *ptr = ( Real * )malloc( sizeof(Real) * num );
        memset( ptr, 0, sizeof(Real) * num );
        return ptr;
}

inline Real *allocMem_dev(const int M, const int N, const int K) 
{
        size_t num = M * N * K;
        Real *ptr;
        checkCudaErrors( cudaMalloc((void **)&ptr, sizeof(Real) * num) );
        checkCudaErrors( cudaMemset(ptr, 0, sizeof(Real) * num) );
        return ptr;
}



void smooth( Real *u, Real *f, Real *mu, const int level )
{
        dim3 block(32, 16, 1);
        dim3 Grid;
        Grid.x = (innersinfo[level].kend - innersinfo[level].kstart + block.x - 1) / block.x;
        Grid.y = (innersinfo[level].jend - innersinfo[level].jstart + block.y - 1) / block.y;
        Grid.z = (innersinfo[level].iend - innersinfo[level].istart + block.z - 1) / block.z;
        int iter;
        for (iter = 0; iter < n_cycles; iter++) 
        {
           rbgs_kernel<0><<<Grid, block>>>(u, f, mu, level);
           rbgs_kernel<1><<<Grid, block>>>(u, f, mu, level);
        }

}

template <int phase>
__global__ void rbgs_kernel( Real * __restrict__ u, Real * __restrict__ f, Real * __restrict__ mu, const int level)
{
        int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int i = threadIdx.z + blockIdx.z * blockDim.z;

        DIMINFO diminfo = dimsinfo_dev[level];
        DIMINFO diminfo_0 = dimsinfo_dev[0];
        INNERS inners = innersinfo_dev[level];

        int Mdim = diminfo.M;
        int Ndim = diminfo.N;
        int Kdim = diminfo.K;
        int fin_coa = diminfo.fin_coa;
        int nx = diminfo_0.K + 2*LEND;
        double DH2 = diminfo.DH * diminfo.DH;
        long long index = 0;
        long long index1 = 0;
        long long index2 = 0;
        long long index3 = 0;
        long long index4 = 0;


        if ( (i + j + k) % 2 != phase ) return;

        if(level == 0)
        {  
           index =  (k+LEND)*2 + ((j+LEND)*2) * (nx*2-1); 
           index1 = ((k+LEND)*2 - 1) + ((j+LEND)*2) * (nx*2-1); 
           index2 = ((k+LEND)*2 + 1) + ((j+LEND)*2) * (nx*2-1); 
           index3 = (k+LEND)*2 + ((j+LEND)*2 - 1) * (nx*2-1); 
           index4 = (k+LEND)*2 + ((j+LEND)*2 + 1) * (nx*2-1);  
           if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) 
           {
              if (j == inners.jstart)
              {
                 u[INDEX(i, j, k)] = ( relax_rbgs / \
                                       ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) ) * \
                                    ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                      mu[index3] * u[INDEX(i, j+1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                                                                                  DH2 * f[INDEX(i, j, k)] ) + \
                                       (1.00 - relax_rbgs) * u[INDEX(i, j, k)];
              }
              else if (j == inners.jend - 1)
              {
                 u[INDEX(i, j, k)] = ( relax_rbgs / \
                                       ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) ) * \
                                    ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                      mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j-1, k)] - \
                                                                                  DH2 * f[INDEX(i, j, k)] ) + \
                                       (1.00 - relax_rbgs) * u[INDEX(i, j, k)];
              }
              else
              {
                 u[INDEX(i, j, k)] = ( relax_rbgs / \
                                       ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) ) * \
                                    ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                      mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                                                                                  DH2 * f[INDEX(i, j, k)] ) + \
                                       (1.00 - relax_rbgs) * u[INDEX(i, j, k)];
              }
           }
        }
        else
        {
           index =  (k*fin_coa+LEND)*2 + ((j*fin_coa+LEND)*2) * (nx*2-1); 
           index1 = ((k*2-1)*fin_coa/2+LEND)*2 + ((j*fin_coa+LEND)*2) * (nx*2-1); 
           index2 = ((k*2+1)*fin_coa/2+LEND)*2 + ((j*fin_coa+LEND)*2) * (nx*2-1); 
           index3 = (k*fin_coa+LEND)*2 + (((j*2-1)*fin_coa/2+LEND)*2) * (nx*2-1); 
           index4 = (k*fin_coa+LEND)*2 + (((j*2+1)*fin_coa/2+LEND)*2) * (nx*2-1); 
           if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) 
           {
              if (j == inners.jstart)
              {
                 u[INDEX(i, j, k)] = ( relax_rbgs / \
                                       ( mu[index1] + mu[index2] + mu[index4] + mu[index4] ) ) * \
                                    ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                      mu[index4] * u[INDEX(i, j+1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                                                                                  DH2 * f[INDEX(i, j, k)] ) + \
                                       (1.00 - relax_rbgs) * u[INDEX(i, j, k)];
              }
              else if (j == inners.jend - 1)
              {
                 u[INDEX(i, j, k)] = ( relax_rbgs / \
                                       ( mu[index1] + mu[index2] + mu[index3] + mu[index3] ) ) * \
                                    ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                      mu[index3] * u[INDEX(i, j-1, k)] + mu[index3] * u[INDEX(i, j-1, k)] - \
                                                                                  DH2 * f[INDEX(i, j, k)] ) + \
                                       (1.00 - relax_rbgs) * u[INDEX(i, j, k)];
              }
              else
              {
                 u[INDEX(i, j, k)] = ( relax_rbgs / \
                                       ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) ) * \
                                    ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                      mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                                                                                  DH2 * f[INDEX(i, j, k)] ) + \
                                       (1.00 - relax_rbgs) * u[INDEX(i, j, k)];
              }
           }
        }


}

__global__ void residual_kernel( Real * __restrict__ r, Real * __restrict__ u, Real * __restrict__ f, Real * __restrict__ mu, const int level )
{
        int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int i = threadIdx.z + blockIdx.z * blockDim.z;

        DIMINFO diminfo = dimsinfo_dev[level];
        DIMINFO diminfo_0 = dimsinfo_dev[0];
        INNERS inners = innersinfo_dev[level];

        int Mdim = diminfo.M;
        int Ndim = diminfo.N;
        int Kdim = diminfo.K;
        int fin_coa = diminfo.fin_coa;
        int nx = diminfo_0.K + 2*LEND;
        double DH2 = diminfo.DH * diminfo.DH;
        long long index = 0;
        long long index1 = 0;
        long long index2 = 0;
        long long index3 = 0;
        long long index4 = 0;

        if(level == 0)
        {
           index =  (k+LEND)*2 + ((j+LEND)*2) * (nx*2-1); 
           index1 = ((k+LEND)*2 - 1) + ((j+LEND)*2) * (nx*2-1); 
           index2 = ((k+LEND)*2 + 1) + ((j+LEND)*2) * (nx*2-1); 
           index3 = (k+LEND)*2 + ((j+LEND)*2 - 1) * (nx*2-1); 
           index4 = (k+LEND)*2 + ((j+LEND)*2 + 1) * (nx*2-1);  
           if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) 
           {
              if (j == inners.jstart)
              {
                 r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                   ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                     mu[index3] * u[INDEX(i, j+1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                   ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
              }
              else if (j == inners.jend - 1)
              {
                 r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                   ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                     mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j-1, k)] - \
                   ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
              }
              else
              {
                 r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                   ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                     mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                   ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
              }
           }
        }
        else
        {
           index =  (k*fin_coa+LEND)*2 + ((j*fin_coa+LEND)*2) * (nx*2-1); 
           index1 = ((k*2-1)*fin_coa/2+LEND)*2 + ((j*fin_coa+LEND)*2) * (nx*2-1); 
           index2 = ((k*2+1)*fin_coa/2+LEND)*2 + ((j*fin_coa+LEND)*2) * (nx*2-1); 
           index3 = (k*fin_coa+LEND)*2 + (((j*2-1)*fin_coa/2+LEND)*2) * (nx*2-1); 
           index4 = (k*fin_coa+LEND)*2 + (((j*2+1)*fin_coa/2+LEND)*2) * (nx*2-1); 
           if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) 
           {
              if (j == inners.jstart)
              {
                 r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                               ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                 mu[index4] * u[INDEX(i, j+1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                               ( mu[index1] + mu[index2] + mu[index4] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
              }
              else if (j == inners.jend - 1)
              {
                 r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                               ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                 mu[index3] * u[INDEX(i, j-1, k)] + mu[index3] * u[INDEX(i, j-1, k)] - \
                             ( mu[index1] + mu[index2] + mu[index3] + mu[index3] ) * u[INDEX(i, j, k)] ) / DH2;
              }
              else
              {
                 r[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                               ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                 mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                             ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
              }
           }
        }


}

__global__ void restriction_kernel( Real * __restrict__ r2, Real * __restrict__ r, Real * __restrict__ e2, const int level )
{
        int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int i = threadIdx.z + blockIdx.z * blockDim.z;
    
        DIMINFO diminfo_coar = dimsinfo_dev[level];
        DIMINFO diminfo_fine = dimsinfo_dev[level - 1];
    
        int Mdim_coar = diminfo_coar.M;
        int Ndim_coar = diminfo_coar.N;
        int Kdim_coar = diminfo_coar.K;
    
        int Mdim_fine = diminfo_fine.M;
        int Ndim_fine = diminfo_fine.N;
        int Kdim_fine = diminfo_fine.K;
    
    
        if (i >= 0 && i < Mdim_coar && j >= 0 && j < Ndim_coar && k >= 1 && k < Kdim_coar - 1) 
        {
           if (j == 0)
           {
              r2[INDEX_coar(i, j, k)] = 0.06250 * (4.0 * r[INDEX_fine(2*i, 2*j, 2*k)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j, 2*k-1)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j, 2*k+1)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j+1, 2*k)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j+1, 2*k)] + \
                                                         r[INDEX_fine(2*i, 2*j+1, 2*k-1)] + \
                                                         r[INDEX_fine(2*i, 2*j+1, 2*k+1)] + \
                                                         r[INDEX_fine(2*i, 2*j+1, 2*k-1)] + \
                                                         r[INDEX_fine(2*i, 2*j+1, 2*k+1)]);
              e2[INDEX_coar(i, j, k)] = 0.0;
           }
           else if (j == Ndim_coar - 1)
           {
              r2[INDEX_coar(i, j, k)] = 0.06250 * (4.0 * r[INDEX_fine(2*i, 2*j, 2*k)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j, 2*k-1)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j, 2*k+1)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j-1, 2*k)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j-1, 2*k)] + \
                                                         r[INDEX_fine(2*i, 2*j-1, 2*k-1)] + \
                                                         r[INDEX_fine(2*i, 2*j-1, 2*k+1)] + \
                                                         r[INDEX_fine(2*i, 2*j-1, 2*k-1)] + \
                                                         r[INDEX_fine(2*i, 2*j-1, 2*k+1)]);
              e2[INDEX_coar(i, j, k)] = 0.0;
           }
           else
           {
              r2[INDEX_coar(i, j, k)] = 0.06250 * (4.0 * r[INDEX_fine(2*i, 2*j, 2*k)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j, 2*k-1)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j, 2*k+1)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j-1, 2*k)] + \
                                                   2.0 * r[INDEX_fine(2*i, 2*j+1, 2*k)] + \
                                                         r[INDEX_fine(2*i, 2*j-1, 2*k-1)] + \
                                                         r[INDEX_fine(2*i, 2*j-1, 2*k+1)] + \
                                                         r[INDEX_fine(2*i, 2*j+1, 2*k-1)] + \
                                                         r[INDEX_fine(2*i, 2*j+1, 2*k+1)]);
              e2[INDEX_coar(i, j, k)] = 0.0;
           }
        }
    

}

__global__ void prolongate_kernel( Real * __restrict__ u, Real * __restrict__ e2, const int level )
{
        int k = threadIdx.x + blockIdx.x * blockDim.x;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int i = threadIdx.z + blockIdx.z * blockDim.z;
    
        DIMINFO diminfo_coar = dimsinfo_dev[level];
        DIMINFO diminfo_fine = dimsinfo_dev[level - 1];
    
        int Mdim_coar = diminfo_coar.M;
        int Ndim_coar = diminfo_coar.N;
        int Kdim_coar = diminfo_coar.K;
    
        int Mdim_fine = diminfo_fine.M;
        int Ndim_fine = diminfo_fine.N;
        int Kdim_fine = diminfo_fine.K;
    
        if (i >= 0 && i < Mdim_coar && j >= 0 && j < (Ndim_coar - 1) && k >= 0 && k < (Kdim_coar - 1)) 
        {
           u[INDEX_fine(2*i, 2*j, 2*k)]     += e2[INDEX_coar(i, j, k)];
           u[INDEX_fine(2*i, 2*j, 2*k+1)]   += 0.5 * (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j, k+1)]);
           u[INDEX_fine(2*i, 2*j+1, 2*k)]   += 0.5 * (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j+1, k)]);
           u[INDEX_fine(2*i, 2*j+1, 2*k+1)] += 0.25* (e2[INDEX_coar(i, j, k)] + e2[INDEX_coar(i, j, k+1)] + \
                                                      e2[INDEX_coar(i, j+1, k)] + e2[INDEX_coar(i, j+1, k+1)]);
        }

}

__global__ void calcuDelta_dev( Real *u, Real *f, Real *mu, Real *res, const int level )
{
        DIMINFO diminfo = dimsinfo_dev[level];
        DIMINFO diminfo_0 = dimsinfo_dev[0];
        INNERS inners = innersinfo_dev[level];
    
        int Mdim = diminfo.M;
        int Ndim = diminfo.N;
        int Kdim = diminfo.K;
        int fin_coa = diminfo.fin_coa;
        int nx = diminfo_0.K + 2*LEND;
        double DH2 = diminfo.DH * diminfo.DH;
        long long index = 0;
        long long index1 = 0;
        long long index2 = 0;
        long long index3 = 0;
        long long index4 = 0;
    
        int k = threadIdx.x + blockIdx.x * blockDim.x + 1;
        int j = threadIdx.y + blockIdx.y * blockDim.y;
        int i = threadIdx.z + blockIdx.z * blockDim.z;
    
    
        index =  (k+LEND)*2 + ((j+LEND)*2) * (nx*2-1); 
        index1 = ((k+LEND)*2 - 1) + ((j+LEND)*2) * (nx*2-1); 
        index2 = ((k+LEND)*2 + 1) + ((j+LEND)*2) * (nx*2-1); 
        index3 = (k+LEND)*2 + ((j+LEND)*2 - 1) * (nx*2-1); 
        index4 = (k+LEND)*2 + ((j+LEND)*2 + 1) * (nx*2-1);  
        if (i >= inners.istart && i < inners.iend && j >= inners.jstart && j < inners.jend && k >= inners.kstart && k < inners.kend) 
        {
           if (j == inners.jstart)
           {
              res[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                               ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                 mu[index3] * u[INDEX(i, j+1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                             ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
           }
           else if (j == inners.jend - 1)
           {
              res[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                               ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                 mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j-1, k)] - \
                             ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
           }
           else
           {
              res[INDEX(i, j, k)] = f[INDEX(i, j, k)] - \
                               ( mu[index1] * u[INDEX(i, j, k-1)] + mu[index2] * u[INDEX(i, j, k+1)] + \
                                 mu[index3] * u[INDEX(i, j-1, k)] + mu[index4] * u[INDEX(i, j+1, k)] - \
                             ( mu[index1] + mu[index2] + mu[index3] + mu[index4] ) * u[INDEX(i, j, k)] ) / DH2;
           }
        }


}

void vcycle( Real *u, Real *f, Real *mu, const int level )
{
        if (level >= n_levels - 1) 
        {
           smooth(u, f, mu, level);
           return;
        }
    
        dim3 block(32, 16, 1);
        dim3 Grid;
        dim3 Grid_1;
        dim3 Grid_2;
        Grid.x = (innersinfo[level].kend - innersinfo[level].kstart + block.x - 1) / block.x;
        Grid.y = (innersinfo[level].jend - innersinfo[level].jstart + block.y - 1) / block.y;
        Grid.z = (innersinfo[level].iend - innersinfo[level].istart + block.z - 1) / block.z;
    
        Real *r = res[level];
        Real *r2 = res2[level];
        Real *e2 = err2[level];
        // Pre-smoothing
        smooth(u, f, mu, level);
    
        // Compute residual error
        residual_kernel <<<Grid, block>>> (r, u, f, mu, level);
    
        // Restriction and fill zero
        Grid_1.x = (dimsinfo[level + 1].K - 2 + block.x - 1) / block.x;
        Grid_1.y = (dimsinfo[level + 1].N + block.y - 1) / block.y;
        Grid_1.z = (dimsinfo[level + 1].M + block.z - 1) / block.z;
        restriction_kernel <<<Grid_1, block>>> (r2, r, e2, level + 1);
    
        // Recursion vcycle
        vcycle(e2, r2, mu, level + 1);
    
        // Prolongation and correction
        Grid_2.x = (dimsinfo[level + 1].K - 1 + block.x - 1) / block.x;
        Grid_2.y = (dimsinfo[level + 1].N - 1 + block.y - 1) / block.y;
        Grid_2.z = (dimsinfo[level + 1].M + block.z - 1) / block.z;
        prolongate_kernel <<<Grid_2, block>>> (u, e2, level + 1);
    
        // Post-smoothing
        smooth(u, f, mu, level);
        // checkCudaErrors( cudaDeviceSynchronize() );
    
        return;

}


void setBoundary( GRID * grid, DISPLACEMENT U, FAULT fault, RATE_STATE_LAW RSL, VARIABLE y, PARAMS * params, double ratio )
{
        int k;
        int nx = grid->nx;
        int nz = grid->nz;
        int ni2 = grid->ni2;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        int i0 = grid->i0;
        long long index = 0;
        long long index2 = 0; 
        
        for ( k = nk1 - 1; k < nk2; k ++ )
        {
           index = k * nx + i0 - 1;
           U.Uy[index] = 0.5 * y.Slip[k];
        }

        for ( k = 0; k < nz; k ++ )
        {
           index = k * nx + ni2 - 1;
           U.Uy[index] = 0.50 * RSL.Vini[0] * ( params->t_elapsed + ratio * params->dt ) + fault.SlipShift[k];
        }

        for ( k = 0; k < nx; k ++ )
        {
           index = ( nk1 - 2 ) * nx + k;
           index2 = nk1 * nx + k;
           U.Uy[index] = U.Uy[index2];
        }

        for ( k = 0; k < nx; k ++ )   
        {
           index = nk2 * nx + k;
           index2 = ( nk2 - 2 ) * nx + k;
           U.Uy[index] = U.Uy[index2];
        }

}

void cal_dpsi( GRID * grid, RATE_STATE_LAW RSL, FAULT fault, VARIABLE y, Runge_Kutta rk4 )
{
        int k;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        double A;
        for ( k = nk1 - 1; k < nk2; k ++ )
        {
           A = exp( (RSL.f0[0] - y.psi[k]) / RSL.b[k] );
           rk4.dpsi[k] = 0.0;
           if( !isinf(A) && RSL.b[k]>1e-3 )
           {
              rk4.dpsi[k] = ( RSL.b[k] * RSL.V0[0] / RSL.L[0] ) * ( A - rk4.V_f[k] / RSL.V0[0] );
           }
        }

}




double computeError( GRID * grid, VARIABLE y3, VARIABLE y4 )
{
        int k;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        int nk = grid->nk;

        double err1 = 0.0;
        double err2 = 0.0;
        double totErr = 0.0;
        double temp1;
        double temp2;
  
        for ( k = nk1 - 1; k < nk2; k ++ )
        {
           temp1 = y3.Slip[k] - y4.Slip[k];
           temp2 = y3.psi[k] - y4.psi[k];
           err1 = err1 + temp1 * temp1;
           err2 = err2 + temp2 * temp2;
        }
        err1 = sqrt( err1 );
        err2 = sqrt( err2 );
  
        totErr = err1 / ( sqrt((double)nk) * 1.0 ) + err2 / ( sqrt((double)nk) * 1.0 );
        //1.0 is the scale;
       
        return totErr;
} 


double computeStepSize( GRID * grid, PARAMS * params, const double totErr, double totTol, double errA[2] )
{
        double ord = 4.0;                        //order, 阶数精度
        double kappa = 0.9;                        
        double stepRatio;
        double alpha;
        double beta;
        double gamma;
        alpha = 0.49/ord;
        beta = 0.34/ord;
        gamma = 0.1/ord;
        double deltaT;

        if( grid->it < 4 )
        {
           stepRatio = kappa * pow( totTol/totErr, 1.0/(1.0+ord) );
        }
        else
        {
           stepRatio = kappa * pow( totTol/totErr, alpha ) * pow( errA[0]/totTol, beta ) \
                             * pow( totTol/errA[1], gamma );
        }

        deltaT = stepRatio * params->dt;
        deltaT = MIN(params->dt * 5.0, deltaT);
        deltaT = MIN(params->maxdt, deltaT);
        deltaT = MAX(params->mindt, deltaT);

        return deltaT;         

}





void RK43( int method, GRID * grid, PARAMS * params, DISPLACEMENT U, FAULT fault, MEDIUM medium, RATE_STATE_LAW RSL )
{

        VARIABLE k1;
        VARIABLE k2;
        VARIABLE k3;
        VARIABLE k4;
        VARIABLE k5;
        VARIABLE k6;
        VARIABLE y3;
        VARIABLE y4;
        Runge_Kutta f1;
        Runge_Kutta f2;
        Runge_Kutta f3;
        Runge_Kutta f4;
        Runge_Kutta f5;
        Runge_Kutta f6;
        alloc_RK43( grid, &k1, &k2, &k3, &k4, &k5, &k6, &y3, &y4, &f1, &f2, &f3, &f4, &f5, &f6 );


        // coefficients
        double c2 = 1./2.;
        double c3 = 83./250.;
        double c4 = 31./50.;
        double c5 = 17./20.;
        double c6 = 1.;
       
        double b1 = 82889./524892.;
        double b3 = 15625./83664.;
        double b4 = 69875./102672.;
        double b5 = -2260./8211.;
        double b6 = 1./4.;
             
        double hb1 = 4586570599./29645900160.;
        double hb3 = 178811875./945068544.;
        double hb4 = 814220225./1159782912.;
        double hb5 = -3700637./11593932.;
        double hb6 = 61727./225920.;
        
        double a21 = 1./2.;                                                                          
    
        double a31 = 13861./62500.;                                                                  
        double a32 = 6889./62500.;
        
        double a41 = -116923316275./ 2393684061468.;
        double a42 = -2731218467317./15368042101831.;
        double a43 = 9408046702089./11113171139209.;
       
        double a51 = -451086348788./ 2902428689909.;
        double a52 = -2682348792572./7519795681897.;
        double a53 = 12662868775082./11960479115383.;
        double a54 = 3355817975965./11060851509271.;
       
        double a61 = 647845179188./3216320057751.;
        double a62 = 73281519250./8382639484533.;
        double a63 = 552539513391./3454668386233.;
        double a64 = 3354512671639./8306763924573.;
        double a65 = 4040./17871.;



        int i, j;
        int ni1 = grid->ni1;
        int ni2 = grid->ni2;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        int nx = grid->nx;
        int nz = grid->nz;
        int num1 = nx * nz;
        int attemptCount;
        long long num2 = nz;
        long long size_kf = sizeof( Real ) * num2 * 2;
        double totErr = 0.0;
        double totTol;
        totTol = params->totTol;
        double errA[2] = {0.0, 0.0}; 
        double newDeltaT;
        Real V_max, Ts_max;


        if( params->maxNumSteps == 0 )
        {
           printf( "ERROR: max number of steps == 0\n" );
           exit(0);
        }
 
        params->dt = params->ini_dt;
//~~~~~~~~~~~~~~~~~~~~~~~GPU Initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        int M = 1;
        int DIM = 2;
        int N = grid->nk;
        int K = grid->ni;
        int fin_coa; 
        size_t num = M * N * K;
        Real DH = grid->dh;

        Real * u_host = NULL;
        Real * f_host = NULL;
        u_host = allocMem(M, N, K);       
        f_host = allocMem(M, N, K);       
        Real *u_dev, *f_dev, *mu_dev;
        checkCudaErrors( cudaMalloc(&u_dev, sizeof(Real) * num) );
        checkCudaErrors( cudaMalloc(&f_dev, sizeof(Real) * num) );
        size_t num_mu = (2*nx-1) * (2*nz-1);
        checkCudaErrors( cudaMalloc(&mu_dev, sizeof(Real) * num_mu) );
        checkCudaErrors( cudaMemcpy(mu_dev, medium.mu_equi, sizeof(Real) * num_mu, cudaMemcpyHostToDevice) );
         
        fin_coa = power(2, 0);
        dimsinfo[0] = {M, N, K, fin_coa, DH};
        int m = M;
        int n = N;
        int k = K;
        double dh = DH;
        // allocate memory of residual ans error for different level
        for (i = 0; i < n_levels - 1; i++) 
        {
           res[i]  = allocMem_dev(m, n, k);
           n  = (n-1)/2 + 1;
           k  = (k-1)/2 + 1;
           printf("n = %d\n",n);
           printf("k = %d\n",k);
           dh *= 2;
           res2[i] = allocMem_dev(m, n, k);
           err2[i] = allocMem_dev(m, n, k);
           fin_coa = power(2, i+1);
           printf("fin_coa = %d\n",fin_coa);
           printf("dh = %f\n",dh);
           dimsinfo[i + 1] = {m, n, k, fin_coa, dh};
        }

        switch (DIM) 
        {
           case (2):
              for (i = 0; i < n_levels; i++) 
              {
                 innersinfo[i].istart = 0;     innersinfo[i].iend = 1;
                 innersinfo[i].jstart = 0;     innersinfo[i].jend = dimsinfo[i].N;
                 innersinfo[i].kstart = 1;     innersinfo[i].kend = dimsinfo[i].K - 1;
              }
              break;
        }

        checkCudaErrors( cudaMemcpyToSymbol(dimsinfo_dev, dimsinfo, sizeof(DIMINFO)*n_levels) );
        checkCudaErrors( cudaMemcpyToSymbol(innersinfo_dev, innersinfo, sizeof(INNERS)*n_levels) );
        cublasHandle_t handle;
        checkCudaErrors( cublasCreate(&handle) );
        
        Real *delta_res;
        delta_res = allocMem_dev(M, N, K);
//~~~~~~~~~~~~~~~~~~~~~~~GPU Initialization~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

        for ( k = nk1 - 1; k < nk2; k ++ )
        {
           k1.Slip[k] = fault.Slip[k];
           k1.psi[k]  = RSL.psi[k];
        }
        setBoundary( grid, U, fault, RSL, k1, params, 0.0 );
        ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
        eval_stress_on_fault( grid, U, params, fault, medium );
        NRsearch( grid, fault, RSL, k1, f1 ); 
        cal_dpsi( grid, RSL, fault, k1, f1 );
        V_max = maxval( grid, fault, f1, 0 );
        Ts_max = maxval( grid, fault, f1, 1);
        printf( "Basin88_800 it=%d, dt=%.15lf, V_max=%.15e, Ts_max=%.15lf \n", grid->it, params->dt, V_max, Ts_max );
//~~~~~~~~~~~~~~~~~~~~~~~~data output~~~~~~~~~~~~~~~~~~~~~~~
        FILE * fp;
        char fileName[256];
        sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/output_fault_0.bin" );
        fp = fopen( fileName, "wb" );
        fwrite( &params->t_elapsed, sizeof( Real ), 1, fp );
        fwrite( f1.V_f, sizeof( Real ), nz, fp );
        fwrite( fault.Slip, sizeof( Real ), nz, fp );
        fwrite( fault.Ts, sizeof( Real ), nz, fp );
        fwrite( RSL.psi, sizeof( Real ), nz, fp );
        fwrite( U.Uy, sizeof( Real ), nx*nz, fp );
        fwrite( fault.Tn, sizeof( Real ), nz, fp );
        fwrite( RSL.a, sizeof( Real ), nz, fp );
        fwrite( RSL.b, sizeof( Real ), nz, fp );
        fwrite( &params->dt, sizeof( Real ), 1, fp );
        fclose( fp );
        params->it_save = params->it_save + 1;
//~~~~~~~~~~~~~~~~~~~~~~~~data output~~~~~~~~~~~~~~~~~~~~~~~



        while( grid->it < params->maxNumSteps && params->t_elapsed < params->tmax )
        {
           grid->it = grid->it + 1;
           attemptCount = 0;
           while( attemptCount < 100 )
           {
              attemptCount = attemptCount + 1;
              if( attemptCount >= 100)
              {
                 printf( "RK43 WARNING: maximum number of attempts reached\n" );
              }
              if ( params->t_elapsed + params->dt > params->tmax )  
              {
                 params->dt = params->tmax - params->t_elapsed;
              }
              memset( k2.Slip, 0, size_kf );      
              memset( k3.Slip, 0, size_kf );      
              memset( k4.Slip, 0, size_kf );      
              memset( k5.Slip, 0, size_kf );      
              memset( k6.Slip, 0, size_kf );      
              memset( y3.Slip, 0, size_kf );      
              memset( y4.Slip, 0, size_kf );      
              memset( f2.V_f, 0, size_kf );      
              memset( f3.V_f, 0, size_kf );      
              memset( f4.V_f, 0, size_kf );      
              memset( f5.V_f, 0, size_kf );      
              memset( f6.V_f, 0, size_kf );      


              for ( k = nk1 - 1; k < nk2; k ++ )
              {
                 k2.Slip[k] = fault.Slip[k] + a21 * params->dt * f1.V_f[k];
                 k2.psi[k] = RSL.psi[k] + a21 * params->dt * f1.dpsi[k];
              }
              setBoundary( grid, U, fault, RSL, k2, params, c2 );
              ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
              eval_stress_on_fault( grid, U, params, fault, medium );
              NRsearch( grid, fault, RSL, k2, f2 ); 
              cal_dpsi( grid, RSL, fault, k2, f2 );

              for ( k = nk1 - 1; k < nk2; k ++ )
              {
                 k3.Slip[k] = fault.Slip[k] + a31 * params->dt * f1.V_f[k];
                 k3.Slip[k] = k3.Slip[k] + a32 * params->dt * f2.V_f[k];
                 k3.psi[k] = RSL.psi[k] + a31 * params->dt * f1.dpsi[k];
                 k3.psi[k] = k3.psi[k] + a32 * params->dt * f2.dpsi[k];
              }
              setBoundary( grid, U, fault, RSL, k3, params, c3 );
              ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
              eval_stress_on_fault( grid, U, params, fault, medium );
              NRsearch( grid, fault, RSL, k3, f3 ); 
              cal_dpsi( grid, RSL, fault, k3, f3 );


              for ( k = nk1 - 1; k < nk2; k ++ )
              {
                 k4.Slip[k] = fault.Slip[k] + a41 * params->dt * f1.V_f[k];
                 k4.Slip[k] = k4.Slip[k] + a42 * params->dt * f2.V_f[k];
                 k4.Slip[k] = k4.Slip[k] + a43 * params->dt * f3.V_f[k];
                 k4.psi[k] = RSL.psi[k] + a41 * params->dt * f1.dpsi[k];
                 k4.psi[k] = k4.psi[k] + a42 * params->dt * f2.dpsi[k];
                 k4.psi[k] = k4.psi[k] + a43 * params->dt * f3.dpsi[k];
              }
              setBoundary( grid, U, fault, RSL, k4, params, c4 );
              ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
              eval_stress_on_fault( grid, U, params, fault, medium );
              NRsearch( grid, fault, RSL, k4, f4 ); 
              cal_dpsi( grid, RSL, fault, k4, f4 );


              for ( k = nk1 - 1; k < nk2; k ++ )
              {
                 k5.Slip[k] = fault.Slip[k] + a51 * params->dt * f1.V_f[k];
                 k5.Slip[k] = k5.Slip[k] + a52 * params->dt * f2.V_f[k];
                 k5.Slip[k] = k5.Slip[k] + a53 * params->dt * f3.V_f[k];
                 k5.Slip[k] = k5.Slip[k] + a54 * params->dt * f4.V_f[k];
                 k5.psi[k] = RSL.psi[k] + a51 * params->dt * f1.dpsi[k];
                 k5.psi[k] = k5.psi[k] + a52 * params->dt * f2.dpsi[k];
                 k5.psi[k] = k5.psi[k] + a53 * params->dt * f3.dpsi[k];
                 k5.psi[k] = k5.psi[k] + a54 * params->dt * f4.dpsi[k];
              }
              setBoundary( grid, U, fault, RSL, k5, params, c5 );
              ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
              eval_stress_on_fault( grid, U, params, fault, medium );
              NRsearch( grid, fault, RSL, k5, f5 ); 
              cal_dpsi( grid, RSL, fault, k5, f5 );

              for ( k = nk1 - 1; k < nk2; k ++ )
              {
                 k6.Slip[k] = fault.Slip[k] + a61 * params->dt * f1.V_f[k];
                 k6.Slip[k] = k6.Slip[k] + a62 * params->dt * f2.V_f[k];
                 k6.Slip[k] = k6.Slip[k] + a63 * params->dt * f3.V_f[k];
                 k6.Slip[k] = k6.Slip[k] + a64 * params->dt * f4.V_f[k];
                 k6.Slip[k] = k6.Slip[k] + a65 * params->dt * f5.V_f[k];
                 k6.psi[k] = RSL.psi[k] + a61 * params->dt * f1.dpsi[k];
                 k6.psi[k] = k6.psi[k] + a62 * params->dt * f2.dpsi[k];
                 k6.psi[k] = k6.psi[k] + a63 * params->dt * f3.dpsi[k];
                 k6.psi[k] = k6.psi[k] + a64 * params->dt * f4.dpsi[k];
                 k6.psi[k] = k6.psi[k] + a65 * params->dt * f5.dpsi[k];
              }
              setBoundary( grid, U, fault, RSL, k6, params, c6 );
              ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
              eval_stress_on_fault( grid, U, params, fault, medium );
              NRsearch( grid, fault, RSL, k6, f6 ); 
              cal_dpsi( grid, RSL, fault, k6, f6 );


              for ( k = nk1 - 1; k < nk2; k ++ )
              {
                 y3.Slip[k] = fault.Slip[k] + hb1 * params->dt * f1.V_f[k];
                 y3.Slip[k] = y3.Slip[k] + hb3 * params->dt * f3.V_f[k];
                 y3.Slip[k] = y3.Slip[k] + hb4 * params->dt * f4.V_f[k];
                 y3.Slip[k] = y3.Slip[k] + hb5 * params->dt * f5.V_f[k];
                 y3.Slip[k] = y3.Slip[k] + hb6 * params->dt * f6.V_f[k];
                 y3.psi[k] = RSL.psi[k] + hb1 * params->dt * f1.dpsi[k];
                 y3.psi[k] = y3.psi[k] + hb3 * params->dt * f3.dpsi[k];
                 y3.psi[k] = y3.psi[k] + hb4 * params->dt * f4.dpsi[k];
                 y3.psi[k] = y3.psi[k] + hb5 * params->dt * f5.dpsi[k];
                 y3.psi[k] = y3.psi[k] + hb6 * params->dt * f6.dpsi[k];


                 y4.Slip[k] = fault.Slip[k] + b1 * params->dt * f1.V_f[k];
                 y4.Slip[k] = y4.Slip[k] + b3 * params->dt * f3.V_f[k];
                 y4.Slip[k] = y4.Slip[k] + b4 * params->dt * f4.V_f[k];
                 y4.Slip[k] = y4.Slip[k] + b5 * params->dt * f5.V_f[k];
                 y4.Slip[k] = y4.Slip[k] + b6 * params->dt * f6.V_f[k];
                 y4.psi[k] = RSL.psi[k] + b1 * params->dt * f1.dpsi[k];
                 y4.psi[k] = y4.psi[k] + b3 * params->dt * f3.dpsi[k];
                 y4.psi[k] = y4.psi[k] + b4 * params->dt * f4.dpsi[k];
                 y4.psi[k] = y4.psi[k] + b5 * params->dt * f5.dpsi[k];
                 y4.psi[k] = y4.psi[k] + b6 * params->dt * f6.dpsi[k];
              }

              totErr = computeError( grid, y3, y4 );
              if( totErr < totTol )
              {
                 break;
              }
              
              params->dt = computeStepSize( grid, params, totErr, totTol, errA );
              if( params->mindt == params->dt ) { break; } 

           }
           params->t_elapsed = params->t_elapsed + params->dt;
           for ( k = nk1 - 1; k < nk2; k ++ )
           {
              fault.Slip[k] = y4.Slip[k];
              k1.Slip[k] = y4.Slip[k];
              f1.V_f[k] = 0.0;
              RSL.psi[k] = y4.psi[k];
              k1.psi[k] = y4.psi[k];
              f1.dpsi[k] = 0.0;
           }
           setBoundary( grid, U, fault, RSL, k1, params, 0.0 );
           ode_solve( method, grid, U, u_host, f_host, u_dev, f_dev, mu_dev, delta_res, handle );
           eval_stress_on_fault( grid, U, params, fault, medium );
           NRsearch( grid, fault, RSL, k1, f1 ); 
           cal_dpsi( grid, RSL, fault, k1, f1 );

           if( totErr != 0.0 )
           {
              newDeltaT = computeStepSize( grid, params, totErr, totTol, errA );
           }
           else
           {
              newDeltaT = params->dt;
           }

           errA[1] = errA[0];
           errA[0] = totErr;
           params->dt = newDeltaT;

           V_max = maxval( grid, fault, f1, 0 );
           Ts_max = maxval( grid, fault, f1, 1);
           printf( "Basin88_800 it=%d, dt=%.15lf, V_max=%.15e, Ts_max=%.15lf \n", grid->it, params->dt, V_max, Ts_max );
//~~~~~~~~~~~~~~~~~~~~~~~~data output~~~~~~~~~~~~~~~~~~~~~~~
           if( grid->it % params->it_save_jump == 0 )
           {
              FILE * fp;
              char fileName[256];
              sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/output_fault_%d.bin", params->it_save );
              fp = fopen( fileName, "wb" );
              fwrite( &params->t_elapsed, sizeof( Real ), 1, fp );
              fwrite( f1.V_f, sizeof( Real ), nz, fp );
              fwrite( fault.Slip, sizeof( Real ), nz, fp );
              fwrite( fault.Ts, sizeof( Real ), nz, fp );
              fwrite( RSL.psi, sizeof( Real ), nz, fp );
              fwrite( &params->dt, sizeof( Real ), 1, fp );
              fclose( fp );
              params->it_save = params->it_save + 1;
           }
//~~~~~~~~~~~~~~~~~~~~~~~~data output~~~~~~~~~~~~~~~~~~~~~~~

        } 

//~~~~~~~~~~~~~~~~~~~~~~~GPU Destroy~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
        checkCudaErrors( cudaFree(delta_res) );
        for (i = 0; i < n_levels - 1; i++) 
        {
           checkCudaErrors( cudaFree(res[i] ) );
           checkCudaErrors( cudaFree(res2[i]) );
           checkCudaErrors( cudaFree(err2[i]) );
        }
        
        checkCudaErrors( cudaFree(u_dev) );
        checkCudaErrors( cudaFree(f_dev) );
        checkCudaErrors( cudaFree(mu_dev) );
        checkCudaErrors( cublasDestroy(handle) );
        free(u_host);
        free(f_host);
//~~~~~~~~~~~~~~~~~~~~~~~GPU Destroy~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

} 

 
