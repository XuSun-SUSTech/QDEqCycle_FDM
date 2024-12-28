#include "header.h"

void setparams( PARAMS * params, GRID * grid )
{  
	params->t_year = 3000;
	params->t_month = 0;
        params->t_day = 0;
        params->t_hour = 0;
        params->t_min = 1;
        params->t_sec = 0;
        params->tmax = params->t_sec + 60.0*(params->t_min + 60.0*(params->t_hour + 24.0*( 
                       params->t_day + 30.0*(params->t_month + 12.0*params->t_year))));
        params->t_elapsed = 0.0;
        params->ini_dt = 1e-3;
        params->mindt = 1e-5;
        params->maxdt = 1e8;
        params->totTol = 1e-7;
        params->maxNumSteps = 1e8;


        grid->ni = 801;
        grid->nk = 801;
        grid->dh = 30.0;
        grid->it = 0;       

        params->it_save = 0;
        params->it_save_jump = 2;
        params->io_unit = 20;
        
}
