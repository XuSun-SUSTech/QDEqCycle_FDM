#include "header.h"

void alloc( GRID grid, RATE_STATE_LAW * RSL, FAULT * fault ) 
{

        int nz = grid.nz;

        long long num = nz;

        Real * pRSL = NULL;
        Real * pfault = NULL; 

        long long size_RSL = sizeof( Real ) * ( num * 6 + 5 );
        long long size_fault = sizeof( Real ) * num * 6;

	pRSL   = ( Real * )malloc( size_RSL );
        pfault = ( Real * )malloc( size_fault );

	memset(  pRSL, 0, size_RSL );
	memset(  pfault, 0, size_fault );

        RSL->a = pRSL;
        RSL->b = pRSL + num;
        RSL->Vw = pRSL + num * 2;
        RSL->psi = pRSL + num * 3;
        RSL->dpsi = pRSL + num * 4;
        RSL->eta = pRSL + num * 5; 
        RSL->f0 = pRSL + num * 6;
        RSL->V0 = pRSL + num * 6 + 1;
        RSL->L = pRSL + num * 6 + 2;
        RSL->fw = pRSL + num * 6 + 3;
        RSL->Vini = pRSL + num * 6 + 4;

        fault->V_f = pfault;
        fault->Slip = pfault + num;
        fault->Ts = pfault + num * 2;
        fault->Tn = pfault + num * 3;
        fault->Tau_n0 = pfault + num * 4;
        fault->SlipShift = pfault + num * 5;


}

void alloc_RK43( GRID * grid, VARIABLE * k1, VARIABLE * k2, VARIABLE * k3, VARIABLE * k4, VARIABLE * k5, VARIABLE * k6, VARIABLE * y3, VARIABLE * y4, Runge_Kutta * f1, Runge_Kutta * f2, Runge_Kutta * f3, Runge_Kutta * f4, Runge_Kutta * f5, Runge_Kutta * f6 ) 
{

        int nz = grid->nz;

        long long num = nz;


        Real * pk1 = NULL; 
        Real * pk2 = NULL; 
        Real * pk3 = NULL; 
        Real * pk4 = NULL; 
        Real * pk5 = NULL; 
        Real * pk6 = NULL; 
        Real * py3 = NULL; 
        Real * py4 = NULL; 
        Real * pf1 = NULL; 
        Real * pf2 = NULL; 
        Real * pf3 = NULL; 
        Real * pf4 = NULL; 
        Real * pf5 = NULL; 
        Real * pf6 = NULL; 

        long long size_k = sizeof( Real ) * num * 2;
        long long size_f = sizeof( Real ) * num * 2;


        pk1 = ( Real * )malloc( size_k );
        pk2 = ( Real * )malloc( size_k );
        pk3 = ( Real * )malloc( size_k );
        pk4 = ( Real * )malloc( size_k );
        pk5 = ( Real * )malloc( size_k );
        pk6 = ( Real * )malloc( size_k );
        py3 = ( Real * )malloc( size_k );
        py4 = ( Real * )malloc( size_k );
        pf1 = ( Real * )malloc( size_f );
        pf2 = ( Real * )malloc( size_f );
        pf3 = ( Real * )malloc( size_f );
        pf4 = ( Real * )malloc( size_f );
        pf5 = ( Real * )malloc( size_f );
        pf6 = ( Real * )malloc( size_f );

	memset(  pk1, 0, size_k );
	memset(  pk2, 0, size_k );
	memset(  pk3, 0, size_k );
	memset(  pk4, 0, size_k );
	memset(  pk5, 0, size_k );
	memset(  pk6, 0, size_k );
	memset(  py3, 0, size_k );
	memset(  py4, 0, size_k );
	memset(  pf1, 0, size_f );
	memset(  pf2, 0, size_f );
	memset(  pf3, 0, size_f );
	memset(  pf4, 0, size_f );
	memset(  pf5, 0, size_f );
	memset(  pf6, 0, size_f );


        k1->Slip = pk1;
        k1->psi = pk1 + num;

        k2->Slip = pk2;
        k2->psi = pk2 + num;

        k3->Slip = pk3;
        k3->psi = pk3 + num;

        k4->Slip = pk4;
        k4->psi = pk4 + num;

        k5->Slip = pk5;
        k5->psi = pk5 + num;

        k6->Slip = pk6;
        k6->psi = pk6 + num;

        y3->Slip = py3;
        y3->psi = py3 + num;

        y4->Slip = py4;
        y4->psi = py4 + num;

        f1->V_f = pf1;
        f1->dpsi = pf1 + num;

        f2->V_f = pf2;
        f2->dpsi = pf2 + num;

        f3->V_f = pf3;
        f3->dpsi = pf3 + num;

        f4->V_f = pf4;
        f4->dpsi = pf4 + num;

        f5->V_f = pf5;
        f5->dpsi = pf5 + num;

        f6->V_f = pf6;
        f6->dpsi = pf6 + num;
}

