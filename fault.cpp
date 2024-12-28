#include "header.h"

double minval( GRID * grid, ADAPTIVE_FACTOR AF )
{
        int k;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        Real smallest;
        smallest = AF.dt_temp[nk1-1];

        for ( k = nk1; k < nk2; k ++ )
        {
           if( smallest > AF.dt_temp[k] )
           {
              smallest = AF.dt_temp[k];
           }
        }
        return smallest;   
}


double maxval( GRID * grid, FAULT fault, Runge_Kutta rk4, int i )
{
        int k;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        Real max;
 
        if( i == 0 )
        {
           max = rk4.V_f[nk1-1];
           for ( k = nk1; k < nk2; k ++ )
           {
              if( max < rk4.V_f[k] )
              {
                 max = rk4.V_f[k];
              }
           }
        }
        else
        {
           max = fault.Ts[nk1-1];
           for ( k = nk1; k < nk2; k ++ )
           {
              if( max < fault.Ts[k] )
              {
                 max = fault.Ts[k];
              }
           }
        }
        return max;   
}

int power( int x, int y )
{
        int product = 1;
        if (y == 0)
        {
           product = 1;
        }
        else
        {
           for(int i = 1; i <= y; i++)
           {
              product *= x;
           }
        }
        return product;
}



void fault_init( GRID grid, COORD coord, RATE_STATE_LAW RSL, FAULT fault, DISPLACEMENT U, PARAMS * params, STRUCTURE structure, MEDIUM medium )
{
        int k;
        int nx = grid.nx;
        int nz = grid.nz;
        int ni2 = grid.ni2;
        int nk1 = grid.nk1;
        int nk2 = grid.nk2;
        int i0 = grid.i0;
        long long index = 0;
        Real b0, b_max, depth;

        RSL.f0[0] = 0.6;
        RSL.V0[0] = 1e-6;
        RSL.L[0] = 0.008;
        RSL.fw[0] = 0.2;
        RSL.Vini[0] = 1e-9;

        b0 = 0.020;
        b_max = 0.00;

        for ( k = nk1 - 1; k < nk2; k ++ )
        {
           index = k * nx + 1;
           depth = -coord.z[index];
           if( depth <= 12.0000 )
           {
              RSL.b[k] = b0;
           }
           else if( depth < 17.0000 )
           {
              RSL.b[k] = b0 + ( b_max - b0 ) * ( depth - 12.0000 ) / 5.0000;
           }
           else
           {
              RSL.b[k] = b_max;
           }
 
           RSL.a[k] = 0.015;
           fault.V_f[k] = RSL.Vini[0];

           fault.Tau_n0[k] = -50.0;

//compute initial psi       
           RSL.psi[k] = log( RSL.V0[0]/RSL.Vini[0] );
           RSL.psi[k] = RSL.f0[0] + RSL.b[k]*RSL.psi[k]; 
           RSL.eta[k] = structure.rho[index] * medium.mu[index]; 
           RSL.eta[k] = 0.5 * sqrt( RSL.eta[k] );
        }

        for ( k = nk1 - 1; k < nk2; k ++ )
        {
           index = k * nx + i0 - 1;
           fault.Slip[k] = 2.0 * U.Uy[index];
        }

        for ( k = 0; k < nz; k ++ )
        {
           index = k * nx + ni2 - 1;
           fault.SlipShift[k] = U.Uy[index];
        }

}



void eval_stress_on_fault( GRID * grid, DISPLACEMENT U, PARAMS * params, FAULT fault, MEDIUM medium )
{
        int k;
        Real du;
        int nx = grid->nx; 
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        int i0 = grid->i0;
        int i1 = grid->ni2;
        long long index = 0;
        for ( k = nk1 -1 ; k < nk2; k ++ )
        {
           index = i0 - 1 + k * nx;
           du = partial1st( grid, U, index );
           fault.Ts[k] = medium.mu[index] * du * 1.0e3;
           fault.Tn[k] = fault.Tau_n0[k];
        }
}

double partial1st( GRID * grid, DISPLACEMENT U, long long index )
{
        Real du;
        Real dh = grid->dh;
        du = (-24.0/(17.0*dh)) * U.Uy[index] + (59.0/(34.0*dh)) * U.Uy[index+1] + (-4.0/(17.0*dh)) * U.Uy[index+2] + (-3.0/(34.0*dh)) * U.Uy[index+3];
        return du;
}




void NRsearch( GRID * grid, FAULT fault, RATE_STATE_LAW RSL, VARIABLE y, Runge_Kutta rk4 )
{
        int k, iter, iter_max;
        Real jacvec, Tau_s;
        Real fa, fb, fy, df, d, var_new, var_old, err;
        int flag_nan;
        Real friction_update, RSL_S;
        Real temp;
        int nk1 = grid->nk1;
        int nk2 = grid->nk2;
        
        for ( k = nk1 -1 ; k < nk2; k ++ )
        {
           flag_nan = 0;
           LOOP: iter = 0;
           var_old = 0.5 * rk4.V_f[k] / RSL.V0[0] * exp( y.psi[k]/RSL.a[k] );
           var_old = asinh( var_old );
         
           fa = RSL.a[k] * abs( fault.Tn[k] );
           fb = fault.Ts[k];

           iter_max = 100000;
           while( iter < iter_max )
           {
              fy = 2.0*RSL.eta[k] * RSL.V0[0]*exp( -y.psi[k]/RSL.a[k] ) * sinh( var_old ) + fa*var_old - fb;
              df = 2.0*RSL.eta[k] * RSL.V0[0]*exp( -y.psi[k]/RSL.a[k] ) * cosh( var_old ) + fa;
              d = -fy / ( df + 1e-100 );
              var_new = var_old + d;
              err = abs(d) / ( abs(var_old) + 1e-100 );
             
              iter = iter + 1;
              if( flag_nan )
              {
                 printf( "check search:%d, %f, %f \n", iter, var_new, var_old );
              }
              var_old = var_new;
              if ( err < 1e-10 ) break;
           }

           friction_update = RSL.a[k] * var_new;
           Tau_s = -fault.Tn[k] * friction_update;
           
           if( isnan( friction_update ) )
           {
              printf( "%d, %f, %f \n", k, friction_update, Tau_s );
              if( flag_nan )
              {
                 printf( "k=%d \n", k );
                 printf( "NaN, error for friction\n" );
                 exit(0);
              }
              else
              {
                 flag_nan = 1;
                 goto LOOP;
              }
           }

           rk4.V_f[k] = 2.0*RSL.V0[0]*exp( -y.psi[k]/RSL.a[k] );
           rk4.V_f[k] = rk4.V_f[k] * sinh(var_new);
	                 
        }

}












