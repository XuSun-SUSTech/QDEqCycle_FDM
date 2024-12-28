#ifndef HEADER_H
#define HEADER_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
# include <string.h>
#include <math.h>
#include "struct.h"
#include "functions.h"
#include <sys/time.h>

#include <vector>
#include <functional>


#if defined(GPU)
    #include <cuda_runtime.h>
    #include <cublas_v2.h>
    // #include <thrust/universal_vector.h>
    #include "helper_cuda.h"
    #include "helper_string.h"
#endif


#define INDEX(i, j, k) ( (k) + (j) * (Kdim) + (i) * (Kdim) * (Ndim) )
#define INDEX_fine(i, j, k) ( (k) + (j) * (Kdim_fine) + (i) * (Kdim_fine) * (Ndim_fine) )
#define INDEX_coar(i, j, k) ( (k) + (j) * (Kdim_coar) + (i) * (Kdim_coar) * (Ndim_coar) )
#define MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define MIN(a,b) (((a) <= (b)) ? (a) : (b))



#if defined(GPU)
extern __constant__ int dim_dev[3];
extern __constant__ int inn_dev[6];
extern cublasHandle_t handle;
#endif

#define KHZ 4


#endif // !HEADER_H
