#ifndef _CUDACOMP_H
#define _CUDACOMP_H


#ifdef HAVE_CUDA

#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <device_types.h>
#include <pthread.h>

 #endif

#ifdef HAVE_CUDA

// data passed to each thread
typedef struct
{
    int thread_no;
    long numl0;
    int cindex; // computation index
    int *status; // where to white status
} THDATA;



/** \brief This structure holds the GPU computation setup for matrix multiplication
 *
 *  By declaring an array of these structures,
 * several parallel computations can be executed
 *
 */

typedef struct
{
    int init; /// 1 if initialized
    int *refWFSinit; /// reference init
    int alloc; /// 1 if memory has been allocated
    long CM_ID;
    long CM_cnt;
    long timerID;
    
    int M;
    int N;


    // synchronization
    int sem; // if sem = 1, wait for semaphore to perform computation
    int gpuinit;
    // one semaphore per thread
    sem_t **semptr1; // command to start matrix multiplication (input)
    sem_t **semptr2; // memory transfer to device completed (output)
    sem_t **semptr3; // computation done (output)
    sem_t **semptr4; // command to start transfer to host (input)
    sem_t **semptr5; // output transfer to host completed (output)

    // computer memory (host)
    float *cMat;
    float **cMat_part;
    float *wfsVec;
    float **wfsVec_part;
    float *wfsRef;
    float **wfsRef_part;
    float *dmVec;
    float *dmVecTMP;
    float **dmVec_part;
    float **dmRef_part;

    // GPU memory (device)
    float **d_cMat;
    float **d_wfsVec;
    float **d_dmVec;
    float **d_wfsRef;
    float **d_dmRef;

    // threads
    THDATA *thdata;
    int *iret;
    pthread_t *threadarray;
    int NBstreams;
    cudaStream_t *stream;
    cublasHandle_t *handle;

    // splitting limits
    int *Nsize;
    int *Noffset;

    int *GPUdevice;

    int orientation;
    long IDout;


} GPUMATMULTCONF;
#endif



int init_cudacomp();

int CUDACOMP_init();

#ifdef HAVE_CUDA
void matrixMulCPU(float *cMat, float *wfsVec, float *dmVec, int M, int N);
int GPUloadCmat(int index);
int GPU_loop_MultMat_setup(int index, char *IDcontrM_name, char *IDwfsim_name, char *IDoutdmmodes_name, long NBGPUs, int *GPUdevice, int orientation, int USEsem, int initWFSref, long loopnb);
int GPU_loop_MultMat_execute(int index, int *status, int *GPUstatus, float alpha, float beta, int timing);
int GPU_loop_MultMat_free(int index);

#ifdef HAVE_MAGMA
int CUDACOMP_magma_compute_SVDpseudoInverse(char *ID_Rmatrix_name, char *ID_Cmatrix_name, double SVDeps, long MaxNBmodes, char *ID_VTmatrix_name);
#endif

void *compute_function( void *ptr );
int CUDACOMP_Coeff2Map_Loop(char *IDmodes_name, char *IDcoeff_name, int GPUindex, char *IDoutmap_name, int offsetmode, char *IDoffset_name);

int CUDACOMP_extractModesLoop(char *in_stream, char *intot_stream, char *IDmodes_name, char *IDrefin_name, char *IDrefout_name, char *IDmodes_val_name, int GPUindex, int PROCESS, int TRACEMODE, int MODENORM, int insem, int axmode);

int GPUcomp_test(long NBact, long NBmodes, long WFSsize, long GPUcnt);

#endif


#endif
