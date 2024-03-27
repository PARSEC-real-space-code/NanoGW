#include <hip/hip_runtime.h>

__global__ void hadamard_prod_mat_(double *d_PsiV, double *d_PsiC, 
    double *d_Cmtrx, int ldn, int iv, int nv, int nc)
{
    size_t index = blockIdx.x * blockDim.x + threadIdx.x; 
    size_t stride = blockDim.x * gridDim.x;
    size_t nmat, jv, jc;
    nmat = ldn * nc;
    for (size_t i = index; i < nmat; i += stride)
    {
       jv = i % ldn; 
       jc = i;
       d_Cmtrx[i] = d_PsiV[iv*ldn+jv] * d_PsiC[jc];  
    }
}

extern "C"
{
    void hadamard_prod_mat(double **d_PsiV, double **d_PsiC, 
                 double **d_Cmtrx, int ldn, int iv, int nv, int nc)
    {
        hipLaunchKernelGGL((hadamard_prod_mat_), dim3(256), dim3(256), 0, 0,
                 *d_PsiV, *d_PsiC, *d_Cmtrx, ldn, iv, nv, nc);
    }
}
