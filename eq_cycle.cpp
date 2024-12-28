#include "header.h"


int main()
{

	PARAMS params;
	GRID grid = { 0 };
	STRUCTURE structure;	
        MEDIUM medium;
	COORD coord;
        DISPLACEMENT U;
        DISPLACEMENT U2;
        DISPLACEMENT U3;
        DISPLACEMENT U4;
        FAULT fault;
        RATE_STATE_LAW RSL;

        int ODE_method = 1;

	setparams( &params, &grid );	
        	
        ode_init( &grid, &coord, &U, &U2, &U3, &U4, &medium, &structure );
        extend_crew( &grid, U );
        cal_grid( grid, coord );
        SetMedium( grid, coord, medium, structure );
        alloc( grid, &RSL, &fault );
        fault_init( grid, coord, RSL, fault, U, &params, structure, medium );

        
        RK43( ODE_method, &grid, &params, U, fault, medium, RSL );

	return 0;
}
