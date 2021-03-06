/*
 * mxm_gpu.cu
 *  @author azamat, mmin
 *  @since  July 13, 2012
 */

#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

#define KERNEL  1
#define TILE   8 //autotune-able
#define VERBOSE 1
#define CUBLAS 0

#if VERBOSE
int dbg=1;
#else
int dbg=0;
#endif

#if CUBLAS
#include <cublas_v2.h>
cublasHandle_t cublas_h;
cublasStatus_t stat = cublasCreate(&cublas_h);
#endif

static int once=0;
cudaEvent_t tstart, tstop, start, stop;
float kern=0.0f, xfer=0.0f;

int mpirank=0, devid=0;

#define onceMallocMemcpy(x,dbg) do{                                \
  if ((x)->sync&0x1) {                                             \
    cudaMalloc(&(x)->dev,(x)->sz);                                 \
    if(dbg){                                                       \
      printf("r%d.d%d cudaMalloc'ed:     %s, %d B\n",mpirank,devid,(x)->vname,(x)->sz);  \
    }                                                              \
    (x)->sync^=0x1;                                                \
  }                                                                \
  if ((x)->sync&0x2) {                                             \
    cudaEventRecord(tstart,0);                                     \
    cudaMemcpy((x)->dev,(x)->host,(x)->sz,cudaMemcpyHostToDevice); \
    cudaEventRecord(tstop,0); cudaEventSynchronize(tstop);         \
    if(dbg){                                                       \
      cudaEventElapsedTime(&xfer,tstart,tstop);                    \
      printf("r%d.d%d cudaMemcpy'ed H2D: %s, %d B, %f ms, %.2f MB/s\n",mpirank,devid,(x)->vname,(x)->sz,xfer,(1e3f*(x)->sz)/(xfer*(1<<20)));  \
    }                                                              \
    (x)->sync^=0x2;                                                \
  }                                                                \
}while(0)
#define onceMemcpyFree(x,dbg) do{                                  \
  if ((x)->sync&0x4) {                                             \
    cudaEventRecord(tstart,0);                                     \
    cudaMemcpy((x)->host,(x)->dev,(x)->sz,cudaMemcpyDeviceToHost); \
    cudaEventRecord(tstop,0); cudaEventSynchronize(tstop);         \
    if(dbg){                                                       \
      cudaEventElapsedTime(&xfer,tstart,tstop);                    \
      printf("r%d.d%d cudaMemcpy'ed D2H: %s, %d B, %f ms, %.2f MB/s\n",mpirank,devid,(x)->vname,(x)->sz,xfer,(1e3f*(x)->sz)/(xfer*(1<<20)));  \
    }                                                              \
    (x)->sync^=0x4;                                                \
  }                                                                \
  if ((x)->sync&0x8) {                                             \
    cudaFree((x)->dev);                                            \
    if(dbg){                                                       \
      printf("r%d.d%d cudaFree'ed:       %s\n",mpirank,devid,(x)->vname);                \
    }                                                              \
    (x)->sync^=0x8;                                                \
  }                                                                \
}while(0)


//=============================================================================
extern "C" {
  struct memptr {
    int sync; //sync flags: 0x1->allocate, 0x2->copy H2D, 0x4->copy D2H, 0x8->deallocate
    int sz;
    double* host;
    double* dev;
    char* vname;
  };
  typedef struct memptr memptr_t;
  void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p);
  void local_grad3_gpu_(
    memptr_t *u1r, memptr_t *u1s, memptr_t *u1t,
    memptr_t *u2r, memptr_t *u2s, memptr_t *u2t,
    memptr_t *u3r, memptr_t *u3s, memptr_t *u3t,
    memptr_t *u1 , memptr_t *u2 , memptr_t *u3 ,
    memptr_t *mp_d, memptr_t *mp_dt,
    int *n, int *nelts, int *lpts1, int *rank);
  void curl_gpu_(
    memptr_t *u1r, memptr_t *u1s, memptr_t *u1t,
    memptr_t *u2r, memptr_t *u2s, memptr_t *u2t,
    memptr_t *u3r, memptr_t *u3s, memptr_t *u3t,
    memptr_t *rxmn,memptr_t *sxmn,memptr_t *txmn,
    memptr_t *rymn,memptr_t *symn,memptr_t *tymn,
    memptr_t *rzmn,memptr_t *szmn,memptr_t *tzmn,
    memptr_t *w1,  memptr_t *w2,  memptr_t *w3,
    memptr_t *w3mn,
    int *nxyz, int *nelts, int *lpts1);
}


//=============================================================================
// basic curl kernel impl
// this source: 44 registers/thread, 200 bytes cmem[0]
// can improve bandwidth: ~83% of peak (140 GB/s) due gmem cache misses; global_load_miss/inst_issued=26%
__global__ void curl_k(
    const double* __restrict__ rxmn,const double* __restrict__ rymn,const double* __restrict__ rzmn,
    const double* __restrict__ sxmn,const double* __restrict__ symn,const double* __restrict__ szmn,
    const double* __restrict__ txmn,const double* __restrict__ tymn,const double* __restrict__ tzmn,
    const double* __restrict__ u1r, const double* __restrict__ u1s, const double* __restrict__ u1t,
    const double* __restrict__ u2r, const double* __restrict__ u2s, const double* __restrict__ u2t,
    const double* __restrict__ u3r, const double* __restrict__ u3s, const double* __restrict__ u3t,
    const double* __restrict__ w3mn,const int lpts1, 
    double* const __restrict__ w1
  ){
  const int k=blockIdx.x*blockDim.x+threadIdx.x;
  const double w3mk=w3mn[threadIdx.x];
  double* const __restrict__ w2 = &w1[lpts1];
  double* __restrict__ const w3 = &w2[lpts1];

  w1[k]= w3mk*u3r[k]*rymn[k]
       + w3mk*u3s[k]*symn[k]
       + w3mk*u3t[k]*tymn[k]
       - w3mk*u2r[k]*rzmn[k]
       - w3mk*u2s[k]*szmn[k]
       - w3mk*u2t[k]*tzmn[k];

  w2[k]= w3mk*u1r[k]*rzmn[k]
       + w3mk*u1s[k]*szmn[k]
       + w3mk*u1t[k]*tzmn[k]
       - w3mk*u3r[k]*rxmn[k]
       - w3mk*u3s[k]*sxmn[k]
       - w3mk*u3t[k]*txmn[k];

  w3[k]= w3mk*u2r[k]*rxmn[k]
       + w3mk*u2s[k]*sxmn[k]
       + w3mk*u2t[k]*txmn[k]
       - w3mk*u1r[k]*rymn[k]
       - w3mk*u1s[k]*symn[k]
       - w3mk*u1t[k]*tymn[k];
}


//=============================================================================
// basic multi-mxm impl
__global__ void mxm_vanilla(const double* __restrict__ a, const int m,
                            const double* __restrict__ b, const int n,
                            double* __restrict__ c, const int p,
                            const int nelts, const int ldims){
  const int row=blockIdx.y*blockDim.y+threadIdx.y;
  const int col=blockIdx.x*blockDim.x+threadIdx.x;
  if(row<m && col<p){ //eliminate out-of-bounds threads
    double s;
    int lda=( ldims&0x1)    *m*n    //if a's bit (0x1) is set, its leading dim is of size m*n 
      , ldb=((ldims&0x2)>>1)*n*p
      , ldc=((ldims&0x4)>>2)*m*p
      , ldi=((ldims&0x8)>>3)*m*n*p; //for inner dimensions
    if(ldims<8){ //no inner iterations
      for(int e=0; e<nelts; e++){
        s=0.0;
        for(int k=0; k<n; k++){
          s+=a[e*lda+k*m+row]*b[e*ldb+col*n+k];
        }
        c[e*ldc+col*m+row]=s;
      }
    }else{
      for(int e=0; e<nelts; e++){
        for(int i=0; i<m; i++){
          s=0.0;
          for(int k=0; k<n; k++){
            s+=a[e*ldi+i*lda+k*m+row]*b[col*n+k];
          }
          c[e*ldi+i*ldc+col*m+row]=s;
        }
      }
    }
  }
}


//=============================================================================
// mxm: R = D * U
__global__ void mxmr_any(
    const double* __restrict__ a, const int m,
    const double* __restrict__ b, const int n,
          double* __restrict__ c, const int p){
  const int col=threadIdx.z*blockDim.y+threadIdx.y;
  register double s=0.0;
  #pragma unroll
  for(int k=0; k<n; k++){
    s+=a[threadIdx.x+m*k]*b[blockIdx.x*n*p+col*n+k];
  }
  c[blockIdx.x*m*p+col*m+threadIdx.x]=s;
}

// bandwidth: ~44% of peak (75 GB/s) 
__global__ void mxmr8(
    const double* __restrict__ a,
    const double* __restrict__ b,
          double* __restrict__ c){
  __shared__ double as[64], bs[512];
  const int col=8*threadIdx.z+threadIdx.y;
  as[8*threadIdx.y+threadIdx.x]=a[8*threadIdx.y+threadIdx.x];
  bs[8*col+threadIdx.x]=b[512*blockIdx.x+8*col+threadIdx.x];
  __syncthreads();
  register double s=0.0;
  #pragma unroll 8
  for(int k=0; k<8; k++){
    s+=as[8*k+threadIdx.x]*bs[8*col+k];
  }
  c[512*blockIdx.x+8*col+threadIdx.x]=s;
}


//=============================================================================
// mxm: S = U * D'
__global__ void mxms_any(
    const double* __restrict__ a, const int m,
    const double* __restrict__ b, const int n,
          double* __restrict__ c, const int p){
  const int col=threadIdx.z*blockDim.y+threadIdx.y;
  register double s=0.0;
  #pragma unroll
  for(int k=0; k<n; k++){
    s+=a[blockIdx.x*m*n*p+threadIdx.z*m*n+m*k+threadIdx.x]*b[threadIdx.y*n+k];
  }
  c[blockIdx.x*m*n*p+col*m+threadIdx.x]=s;
}

// bandwidth: ~44% of peak (74 GB/s) 
__global__ void mxms8(
    const double* __restrict__ a,
    const double* __restrict__ b,
          double* __restrict__ c){
  __shared__ double as[512], bs[64];
  const int col=8*threadIdx.z+threadIdx.y;
  as[8*col+threadIdx.x]=a[512*blockIdx.x+8*col+threadIdx.x];
  bs[8*threadIdx.y+threadIdx.x]=b[8*threadIdx.y+threadIdx.x];
  __syncthreads();
  register double s=0.0;
  #pragma unroll 8
  for(int k=0; k<8; k++){
    s+=as[64*threadIdx.z+8*k+threadIdx.x]*bs[8*threadIdx.y+k];
  }
  c[512*blockIdx.x+8*col+threadIdx.x]=s;
}


//=============================================================================
// mxm: T = U * D'
__global__ void mxmt_any(
    const double* __restrict__ a, const int m,
    const double* __restrict__ b, const int n,
          double* __restrict__ c, const int p){
  const int row=threadIdx.z*blockDim.x+threadIdx.x;
  register double s=0.0;
  #pragma unroll
  for(int k=0; k<n; k++){
    s+=a[blockIdx.x*m*n+row+k*m]*b[threadIdx.y*n+k];
  }
  c[blockIdx.x*m*p+threadIdx.y*m+row]=s;
}

__global__ void mxmt8(
    const double* __restrict__ a,
    const double* __restrict__ b,
          double* __restrict__ c){
  __shared__ double as[512], bs[64];
  const int row=8*threadIdx.z+threadIdx.x;
  as[64*threadIdx.y+row]=a[512*blockIdx.x+64*threadIdx.y+row];
  bs[ 8*threadIdx.y+threadIdx.x]=b[8*threadIdx.y+threadIdx.x];
  __syncthreads();
  register double s=0.0;
  #pragma unroll 8
  for(int k=0; k<8; k++){
    s+=as[8*threadIdx.z+threadIdx.x+k*64]*bs[8*threadIdx.y+k];
  }
  c[512*blockIdx.x+64*threadIdx.y+row]=s;
}


//=============================================================================
// mxm with 1D arrays
__global__ void mxm_1d(double* a, const int m, double* b, const int n, double* c, const int p){
  const int i=blockIdx.x*blockDim.x+threadIdx.x;
  if (i<m){
    for(int k=0; k<p; k++){
      double s=0.0;
      for(int j=0; j<n; j++){
        s+=a[j*m+i]*b[k*n+j];
      }
      c[k*m+i]=s;
    }
  }
}


// mxm with 2D arrays
__global__ void mxm_shared(double* a, const int m, double* b, const int n, double* c, const int p){
  __shared__ double as[TILE][TILE];
  __shared__ double bs[TILE][TILE];
  int bx=blockIdx.x, by=blockIdx.y, tx=threadIdx.x, ty=threadIdx.y;
  const int row=by*TILE+ty;
  const int col=bx*TILE+tx;
  double s=0.0;
  for(int t=0;t<m/TILE;t++){
    as[ty][tx]=a[col*m+t*TILE+tx];
    bs[ty][tx]=b[col*n+t*TILE+ty];
    __syncthreads();
    for(int k=0; k<TILE; k++){
      s+=as[ty][k]*bs[k][tx];
    }
    __syncthreads();
    c[col*m+row]=s;
  }
}


// globally-visible basic mxm implementation for small matrices
void mxm_std_gpu_(double* a, int* m, double* b, int* n, double* c, int* p){
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=*m*(*n)*sizeof(double)
    , sizeofB=*n*(*p)*sizeof(double)
    , sizeofC=*m*(*p)*sizeof(double);
  /*malloc and memcopy data H2D*/
  cudaMalloc(&dev_a,sizeofA);
  cudaMalloc(&dev_b,sizeofB);
  cudaMalloc(&dev_c,sizeofC);
  cudaMemcpy(dev_a,a,sizeofA,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b,b,sizeofB,cudaMemcpyHostToDevice);
  /*thread dimensions*/
  dim3 dimBlock, dimGrid;
#if KERNEL==1
  dimBlock.x=TILE; dimGrid.x=(*p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(*m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p,1,0);
#elif KERNEL==2
  dimBlock.x=TILE; dimGrid.x=(*m+dimBlock.x-1)/dimBlock.x;
  mxm_1d<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p);
#else
  dimBlock.x=TILE; dimGrid.x=(*p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(*m+dimBlock.y-1)/dimBlock.y;
  mxm_shared<<<dimGrid,dimBlock>>>(dev_a,*m,dev_b,*n,dev_c,*p);
#endif
  /*memcopy D2H*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
}


// sets up the aggregated mxm kernel launch
void mxm_gpu2(double* a, int as, int m
             ,double* b, int bs, int n
             ,double* c, int cs, int p
             ,int nelts, int mask, int dev){
  cudaSetDevice(dev);
  /*device variables*/
  double *dev_a, *dev_b, *dev_c;
  int sizeofA=as*sizeof(double)
    , sizeofB=bs*sizeof(double)
    , sizeofC=cs*sizeof(double);
  /*malloc and memcopy H2D*/
  cudaMalloc(&dev_a,sizeofA);
  cudaMalloc(&dev_b,sizeofB);
  cudaMalloc(&dev_c,sizeofC);
  cudaMemcpy(dev_a,a,sizeofA,cudaMemcpyHostToDevice);
  cudaMemcpy(dev_b,b,sizeofB,cudaMemcpyHostToDevice);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(dev_a,m, dev_b,n, dev_c,p, nelts,mask);
  /*memcopy D2H*/
  cudaMemcpy(c,dev_c,sizeofC,cudaMemcpyDeviceToHost);
  cudaFree(dev_a);
  cudaFree(dev_b);
  cudaFree(dev_c);
}

//=============================================================================
// sets up the aggregated mxm kernel launch
void mxm_gpu_agg(memptr_t *a, int m
                ,memptr_t *b, int n
                ,memptr_t *c, int p
                ,int nelts, int mask, int dev){
  cudaSetDevice(dev);
  /*malloc and memcopy H2D*/
  onceMallocMemcpy(a,dbg);
  onceMallocMemcpy(b,dbg);
  onceMallocMemcpy(c,dbg);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=TILE; dimGrid.x=(p+dimBlock.x-1)/dimBlock.x;
  dimBlock.y=TILE; dimGrid.y=(m+dimBlock.y-1)/dimBlock.y;
  mxm_vanilla<<<dimGrid,dimBlock>>>(a->dev,m, b->dev,n, c->dev,p, nelts,mask);
  /*memcopy D2H and dealloc*/
  onceMemcpyFree(a,dbg);
  onceMemcpyFree(b,dbg);
  onceMemcpyFree(c,dbg);
}


//=============================================================================
/**
 * Performs aggregated mxm for all elements at once.
 *
 * foreach e in 0..nelts
 *   u@r_{NxN^2} = d_{NxN} * u@_{NxN^2}^{e} // here @ is either 1, 2 or 3
 *   foreach k in 0..N
 *     u@s_{NxN}^{k} = u@_{NxN}^{k,e} * dt_{NxN}
 *   u@t_{N^2xN} = u@_{N^2xN}^{e} * dt_{NxN}
 */
void local_grad3_gpu_(memptr_t *u1r, memptr_t *u1s, memptr_t *u1t,  
                      memptr_t *u2r, memptr_t *u2s, memptr_t *u2t,  
                      memptr_t *u3r, memptr_t *u3s, memptr_t *u3t,  
                      memptr_t *u1 , memptr_t *u2 , memptr_t *u3 ,  
                      memptr_t *d,   memptr_t *dt,
                      int *n, int *nelts, int *lpts1, int *rank){
  int n1=*n, n2=n1*n1, n3=n1*n2, ne=*nelts;
  float gbytes = 1e3f*((2*ne*n3+n2)*3*8.0f)/(1<<30);
  float gflops = 1e3f*2*n3*n1*ne*3/(1<<30);

  // select the device
  int devs = 0;
  cudaGetDeviceCount(&devs);
  if (devs==1) {
    devid = 0;
  } else {
    devid = *rank%2;
  }
  cudaSetDevice(devid);
  mpirank=*rank;

  if (!once) {
    d->vname   = "d";
    dt->vname  = "dt";
    u1r->vname = "u1r";
    u1s->vname = "u1s";
    u1t->vname = "u1t";
    u2r->vname = "u2r";
    u2s->vname = "u2s";
    u2t->vname = "u2t";
    u3r->vname = "u3r";
    u3s->vname = "u3s";
    u3t->vname = "u3t";
    u1->vname  = "u1";
    u2->vname  = "u2";
    u3->vname  = "u3";
    cudaEventCreate(&tstart); cudaEventCreate(&tstop);
    cudaEventCreate(&start);  cudaEventCreate(&stop);
  }

  onceMallocMemcpy(d,  dbg);
  onceMallocMemcpy(u1r,dbg);
  onceMallocMemcpy(u2r,dbg);
  onceMallocMemcpy(u3r,dbg);
  // u1,u2,u3 are contiguous, do a single transfer
  u1->sz=*lpts1*3*sizeof(double);
  onceMallocMemcpy(u1, dbg);
  u2->dev=u1->dev+(*lpts1);
  u3->dev=u2->dev+(*lpts1);

  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;

  cudaEventRecord(start,0);

#if CUBLAS
  const double alpha = 1.0;
  const double beta  = 0.0;
  int inci, incj;
  for(int i=0; i<ne; i++){
    inci = i*n3;
    for(int j=0; j<*n; j++){
      incj = j*n2;
      cublasDgemm(cublas_h, CUBLAS_OP_N, CUBLAS_OP_N, *n,*n,*n, &alpha,
        d->dev,*n, u1->dev+inci+incj,*n, &beta, u1r->dev+inci+incj,*n);

      cublasDgemm(cublas_h, CUBLAS_OP_N, CUBLAS_OP_N, *n,*n,*n, &alpha,
        d->dev,*n, u2->dev+inci+incj,*n, &beta, u2r->dev+inci+incj,*n);

      cublasDgemm(cublas_h, CUBLAS_OP_N, CUBLAS_OP_N, *n,*n,*n, &alpha,
        d->dev,*n, u3->dev+inci+incj,*n, &beta, u3r->dev+inci+incj,*n);
    }
  } // this gets 0.19 GB/s
#else
  /* D_{NxN} * U_{NxN^2} = R_{NxN^2} foreach e */
  dimBlock.x=*n; dimBlock.y=*n, dimBlock.z=*n;
  dimGrid.x=ne;  dimGrid.y=1;   dimGrid.z=1;
  if (*n==8){
    mxmr8<<<dimGrid,dimBlock>>>(d->dev, u1->dev, u1r->dev);
    mxmr8<<<dimGrid,dimBlock>>>(d->dev, u2->dev, u2r->dev);
    mxmr8<<<dimGrid,dimBlock>>>(d->dev, u3->dev, u3r->dev);
  }else{
    // todo: dispatch to other specialized mxmr kernels
    mxmr_any<<<dimGrid,dimBlock>>>(d->dev,*n, u1->dev,*n, u1r->dev,n2);
    mxmr_any<<<dimGrid,dimBlock>>>(d->dev,*n, u2->dev,*n, u2r->dev,n2);
    mxmr_any<<<dimGrid,dimBlock>>>(d->dev,*n, u3->dev,*n, u3r->dev,n2);
  }
#endif
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("r%d.d%d r kernel time:     %f ms, eff.bw: %f GB/s, perf: %f GFlop/s\n",mpirank,devid,kern,gbytes/kern,gflops/kern);
  }

  onceMallocMemcpy(dt, dbg);
  onceMallocMemcpy(u1s,dbg);
  onceMallocMemcpy(u2s,dbg);
  onceMallocMemcpy(u3s,dbg);
  /* U_{NxN} * D'_{NxN} = S_{NxN} foreach e,k */
  dimBlock.x=*n; dimBlock.y=*n; dimBlock.z=*n;
  dimGrid.x=ne;  dimGrid.y=1;   dimGrid.z=1;
  cudaEventRecord(start,0);
  if (*n==8){
    mxms8<<<dimGrid,dimBlock>>>(u1->dev, dt->dev, u1s->dev);
    mxms8<<<dimGrid,dimBlock>>>(u2->dev, dt->dev, u2s->dev);
    mxms8<<<dimGrid,dimBlock>>>(u3->dev, dt->dev, u3s->dev);
  }else{
    // todo: dispatch to other specialized mxms kernels
    mxms_any<<<dimGrid,dimBlock>>>(u1->dev,*n, dt->dev,*n, u1s->dev,*n);
    mxms_any<<<dimGrid,dimBlock>>>(u2->dev,*n, dt->dev,*n, u2s->dev,*n);
    mxms_any<<<dimGrid,dimBlock>>>(u3->dev,*n, dt->dev,*n, u3s->dev,*n);
  }
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("r%d.d%d s kernel time:     %f ms, eff.bw: %f GB/s, perf: %f GFlop/s\n",mpirank,devid,kern,gbytes/kern,gflops/kern);
  }

  onceMallocMemcpy(u1t,dbg);
  onceMallocMemcpy(u2t,dbg);
  onceMallocMemcpy(u3t,dbg);
  /* U_{N^2xN} * D'_{NxN} = T_{N^2xN} foreach e */
  dimBlock.x=*n; dimBlock.y=*n; dimBlock.z=*n;
  dimGrid.x=ne;  dimGrid.y=1;   dimGrid.z=1;
  cudaEventRecord(start,0);
  if (*n==8){
    mxmt8<<<dimGrid,dimBlock>>>(u1->dev, dt->dev, u1t->dev);
    mxmt8<<<dimGrid,dimBlock>>>(u2->dev, dt->dev, u2t->dev);
    mxmt8<<<dimGrid,dimBlock>>>(u3->dev, dt->dev, u3t->dev);
  }else{
    // todo: dispatch to other specialized mxmt kernels
    mxmt_any<<<dimGrid,dimBlock>>>(u1->dev,n2, dt->dev,*n, u1t->dev,*n);
    mxmt_any<<<dimGrid,dimBlock>>>(u2->dev,n2, dt->dev,*n, u2t->dev,*n);
    mxmt_any<<<dimGrid,dimBlock>>>(u3->dev,n2, dt->dev,*n, u3t->dev,*n);
  }
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("r%d.d%d t kernel time:     %f ms, eff.bw: %f GB/s, perf: %f GFlop/s\n",mpirank,devid,kern,gbytes/kern,gflops/kern);
  }

  // nothing to copy D2H or to free
  //cudaDeviceSynchronize();
}

//=============================================================================
// Sets up the curl kernel
void curl_gpu_(memptr_t *u1r,  memptr_t *u1s,  memptr_t *u1t,
               memptr_t *u2r,  memptr_t *u2s,  memptr_t *u2t,
               memptr_t *u3r,  memptr_t *u3s,  memptr_t *u3t,
               memptr_t *rxmn, memptr_t *sxmn, memptr_t *txmn,
               memptr_t *rymn, memptr_t *symn, memptr_t *tymn,
               memptr_t *rzmn, memptr_t *szmn, memptr_t *tzmn,
               memptr_t *w1,   memptr_t *w2,   memptr_t *w3,
               memptr_t *w3mn, int *nxyz, int *nelts, int *lpts1){
  int n3=*nxyz, npts=*nelts*n3;
  float gbytes = 1e3f*((n3+21*npts)*8.0f)/(1<<30);
  float gflops = 1e3f*(51*npts)/(1<<30);
  if (!once){
    rxmn->vname="rxmn"; sxmn->vname="sxmn"; txmn->vname="txmn";
    rymn->vname="rymn"; symn->vname="symn"; tymn->vname="tymn";
    rzmn->vname="rzmn"; szmn->vname="szmn"; tzmn->vname="tzmn";
    w3mn->vname="w3mn";
    w1->vname="w1"; w2->vname="w2"; w3->vname="w3";
    once=1;
  }
  /*malloc and memcopy H2D*/
  onceMallocMemcpy(rxmn,dbg);
  onceMallocMemcpy(rymn,dbg);
  onceMallocMemcpy(rzmn,dbg);
  onceMallocMemcpy(sxmn,dbg);
  onceMallocMemcpy(symn,dbg);
  onceMallocMemcpy(szmn,dbg);
  onceMallocMemcpy(txmn,dbg);
  onceMallocMemcpy(tymn,dbg);
  onceMallocMemcpy(tzmn,dbg);
  onceMallocMemcpy(w3mn,dbg);
  onceMallocMemcpy(u1r, dbg);
  onceMallocMemcpy(u1s, dbg);
  onceMallocMemcpy(u1t, dbg);
  onceMallocMemcpy(u2r, dbg);
  onceMallocMemcpy(u2s, dbg);
  onceMallocMemcpy(u2t, dbg);
  onceMallocMemcpy(u3r, dbg);
  onceMallocMemcpy(u3s, dbg);
  onceMallocMemcpy(u3t, dbg);
  // w1,w2,w3 are contiguous, do a single transfer
  w1->sz=*lpts1*3*sizeof(double);
  onceMallocMemcpy(w1,  dbg);
  /*thread grid dimensions*/
  dim3 dimBlock, dimGrid;
  dimBlock.x=*nxyz; dimGrid.x=*nelts;
  cudaEventRecord(start,0);
  curl_k<<<dimGrid,dimBlock>>>(
    rxmn->dev,rymn->dev,rzmn->dev,
    sxmn->dev,symn->dev,szmn->dev,
    txmn->dev,tymn->dev,tzmn->dev,
    u1r->dev, u1s->dev, u1t->dev,
    u2r->dev, u2s->dev, u2t->dev,
    u3r->dev, u3s->dev, u3t->dev,
    w3mn->dev, *lpts1,
    w1->dev
  );
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&kern,start,stop);
  if(dbg){
    printf("r%d.d%d curl kernel time:  %f ms, eff.bw: %f GB/s, perf: %f GFlop/s\n",mpirank,devid,kern,gbytes/kern,gflops/kern);
  }
  onceMemcpyFree(rxmn,dbg);
  onceMemcpyFree(rymn,dbg);
  onceMemcpyFree(rzmn,dbg);
  onceMemcpyFree(sxmn,dbg);
  onceMemcpyFree(symn,dbg);
  onceMemcpyFree(szmn,dbg);
  onceMemcpyFree(txmn,dbg);
  onceMemcpyFree(tymn,dbg);
  onceMemcpyFree(tzmn,dbg);
  onceMemcpyFree(w3mn,dbg);
  onceMemcpyFree(u1r, dbg);
  onceMemcpyFree(u1s, dbg);
  onceMemcpyFree(u1t, dbg);
  onceMemcpyFree(u2r, dbg);
  onceMemcpyFree(u2s, dbg);
  onceMemcpyFree(u2t, dbg);
  onceMemcpyFree(u3r, dbg);
  onceMemcpyFree(u3s, dbg);
  onceMemcpyFree(u3t, dbg);
  onceMemcpyFree(w1,  dbg);
}

