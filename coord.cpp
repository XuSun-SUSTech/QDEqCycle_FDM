#include "header.h"

void cal_grid( GRID grid, COORD coord )
{
        int i;
        int k;
        int nx = grid.nx;
        int nz = grid.nz; 
        int nk2 = grid.nk2;
        Real dh = grid.dh;  
        long long index = 0;

        for ( k = 0; k < nz; k ++ ) 
        {
           for ( i = 0; i < nx; i ++ )
           {
              index = i + k * nx;
              coord.x[index] = (i-1) * dh * 1.0e-3;
           }
        } 

        for ( k = 0; k < nz; k ++ ) 
        {
           for ( i = 0; i < nx; i ++ )
           {
              index = i + k * nx;
              coord.z[index] = (k+1-nk2) * dh * 1.0e-3;
           }
        }
//*
        FILE * fp;
        char fileName[256];
        sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/coord_x.bin" );
        fp = fopen( fileName, "wb" ); 
        fwrite( coord.x, sizeof( Real ), nx * nz, fp );
        fclose( fp ); 

        sprintf( fileName, "/Nas/sunxu/data_test/GJI/MG/Basin88_800/coord_z.bin" );
        fp = fopen( fileName, "wb" ); 
        fwrite( coord.z, sizeof( Real ), nx * nz, fp );
        fclose( fp ); 
//*/
} 
