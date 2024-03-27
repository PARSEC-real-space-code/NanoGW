#include <hip/hip_runtime.h>

__global__ void hadamard_prod_mat_sym_(double *d_PsiV, double *d_PsiC, 
    double *d_Cmtrx, int *d_lcrep, int ldn, int iv, int nv, int nc,
    int maxnc_sym, int lcrep_bound1, int lcrep_bound2)
{
    size_t index = blockIdx.x * blockDim.x + threadIdx.x; 
    size_t stride = blockDim.x * gridDim.x;
    size_t nmat, jr, ic;
    nmat = ldn * (lcrep_bound2 - lcrep_bound1 + 1);
    for (size_t i = index; i < nmat; i += stride)
    {
       jr = i % ldn; 
       ic = d_lcrep[i / ldn + lcrep_bound1 - 1] - 1;
       d_Cmtrx[i] = d_PsiV[iv*ldn+jr] * d_PsiC[ic*ldn+jr];
    }
}

extern "C"
{
    void hadamard_prod_mat_sym( double **d_PsiV, double **d_PsiC, 
                 double **d_Cmtrx, int **d_lcrep, int ldn, int iv, 
                 int nv, int nc, int maxnc_sym, 
                 int lcrep_bound1, int lcrep_bound2 )
    {
        hipLaunchKernelGGL((hadamard_prod_mat_sym_), 
         dim3(256), dim3(256), 0, 0,
        *d_PsiV, *d_PsiC, *d_Cmtrx, *d_lcrep, ldn, iv, 
         nv, nc, maxnc_sym, lcrep_bound1, lcrep_bound2);
    }
}

