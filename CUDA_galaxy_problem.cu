// on dione, first load the cuda module
//    module load cuda
//
// compile your program with
//    nvcc -O3 -arch=sm_70 --ptxas-options=-v -o galaxy galaxy_cuda.cu -lm
//
// run your program with 
//    srun -p gpu -c 1 --mem=10G ./galaxy_cuda RealGalaxies_100k_arcmin.dat SyntheticGalaxies_100k_arcmin.dat omega.out



#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

int    NoofReal;
int    NoofRand;
float *real_rasc, *real_decl;
float *rand_rasc, *rand_decl;

int *histogramDR, *histogramDD, *histogramRR;

long int CPUMemory = 0L;
long int GPUMemory = 0L;

const int totaldegrees = 360;
const int binsperdegree = 4;

const int N = 32; //GPU kernels are launched with N*N threads in a block
                    //Values of N lower than 19 do not work, as a minimum of 360 threads per block are required for copying the local histogram into global memory

// put here your GPU kernel(s) to calculate the histograms

__shared__ int localHist[totaldegrees*binsperdegree];
__global__ void fillHistogram(float *rasc_1, float *dec_1, float *rasc_2, float *dec_2, int *histogram) {

    // 0-initialize temporary histogram
    int tid = threadIdx.x + threadIdx.y * N;
    if (tid <= totaldegrees*binsperdegree) {
            localHist[tid] = 0;
    }
    __syncthreads();

    //do the calculations
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i < 100000 && j < 100000) {
        int bin;
        float innerexpr, result;
        innerexpr = sinf(dec_1[i]) * sinf(dec_2[j])+ cosf(dec_1[i]) * cosf(dec_2[j]) * cosf(rasc_1[i] - rasc_2[j]);
        if (innerexpr > 1.0f) {
            innerexpr = 1.0f;
        }
        result = acosf(innerexpr) *180/M_PI;   //final result, in degrees
        bin = (int)(result/0.25);   //the histogram bin the result maps to
        if (bin <= totaldegrees*binsperdegree) {
            atomicAdd_block(&localHist[bin], 1);
        }
    }
    __syncthreads();
    
    //flush results into global memory
    if (tid <= totaldegrees*binsperdegree) {
            atomicAdd(&(histogram[tid]), localHist[tid]);
    }
}

int main(int argc, char *argv[])
{
   int    readdata(char *argv1, char *argv2);
   long int histogramDRsum, histogramDDsum, histogramRRsum;
   double walltime;
   struct timeval _ttime;
   struct timezone _tzone;
   int getDevice(void);

   FILE *outfil;

   if ( argc != 4 ) {printf("Usage: a.out real_data random_data output_data\n");return(-1);}

   gettimeofday(&_ttime, &_tzone);
   walltime = (double)_ttime.tv_sec + (double)_ttime.tv_usec/1000000.;

   if ( readdata(argv[1], argv[2]) != 0 ) return(-1);

// For your entertainment: some performance parameters of the GPU you are running your programs on!
   //if ( getDevice() != 0 ) return(-1);
    histogramDR = (int *)calloc(totaldegrees*binsperdegree,sizeof(int));
    histogramDD = (int *)calloc(totaldegrees*binsperdegree,sizeof(int));
    histogramRR = (int *)calloc(totaldegrees*binsperdegree,sizeof(int));
    CPUMemory += 3*1L*(totaldegrees*binsperdegree)*sizeof(int);
   // input data is available in the arrays float real_rasc[], real_decl[], rand_rasc[], rand_decl[];
   // allocate memory on the GPU for input data and histograms
   // and initialize the data on GPU by copying the real and rand data to the GPU    
        
    //allocating and copying data arrays
    float *d_real_rasc; cudaMalloc(&d_real_rasc, (size_t)(sizeof(float)*100000));
    float *d_real_decl; cudaMalloc(&d_real_decl, (size_t)(sizeof(float)*100000));
    float *d_rand_rasc; cudaMalloc(&d_rand_rasc, (size_t)(sizeof(float)*100000));
    float *d_rand_decl; cudaMalloc(&d_rand_decl, (size_t)(sizeof(float)*100000));
    GPUMemory += 4*sizeof(float)*NoofReal;
        
    cudaMemcpy(d_real_rasc, real_rasc, (size_t)sizeof(float)*100000, cudaMemcpyHostToDevice);
    cudaMemcpy(d_real_decl, real_decl, (size_t)sizeof(float)*100000, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rand_rasc, rand_rasc, (size_t)sizeof(float)*100000, cudaMemcpyHostToDevice);
    cudaMemcpy(d_rand_decl, rand_decl, (size_t)sizeof(float)*100000, cudaMemcpyHostToDevice);
    
    //allocate histogram on device
    int *d_hist; cudaMalloc(&d_hist, (size_t)(totaldegrees*binsperdegree * sizeof(int)));
    GPUMemory += (totaldegrees*binsperdegree)*sizeof(int);

    // call the GPU kernel(s) that fill the three histograms
    dim3 threadsInBlock = dim3(N, N);
    dim3 blocksInGrid = dim3((NoofReal+N-1)/N, (NoofReal+N-1)/N);
    
    //histogram DD
    cudaMemset(d_hist, 0, sizeof(int)*totaldegrees*binsperdegree);
    fillHistogram<<<blocksInGrid, threadsInBlock, sizeof(int)*totaldegrees*binsperdegree>>>(d_real_rasc, d_real_decl, d_real_rasc, d_real_decl, d_hist);
	cudaDeviceSynchronize();
	cudaMemcpy(histogramDD, d_hist, sizeof(int)*totaldegrees*binsperdegree, cudaMemcpyDeviceToHost);
        
    //histogram DR
    cudaMemset(d_hist, 0, sizeof(int)*totaldegrees*binsperdegree);
    fillHistogram<<<blocksInGrid, threadsInBlock, sizeof(int)*totaldegrees*binsperdegree>>>(d_real_rasc, d_real_decl, d_rand_rasc, d_rand_decl, d_hist);
	cudaDeviceSynchronize();
	cudaMemcpy(histogramDR, d_hist, sizeof(int)*totaldegrees*binsperdegree, cudaMemcpyDeviceToHost);
        
    //histogram RR
    cudaMemset(d_hist, 0, sizeof(int)*totaldegrees*binsperdegree);
    fillHistogram<<<blocksInGrid, threadsInBlock, sizeof(int)*totaldegrees*binsperdegree>>>(d_rand_rasc, d_rand_decl, d_rand_rasc, d_rand_decl, d_hist);
	cudaDeviceSynchronize();
	cudaMemcpy(histogramRR, d_hist, sizeof(int)*totaldegrees*binsperdegree, cudaMemcpyDeviceToHost);
    
    
    //all done on GPU, free up memory
    cudaFree(d_real_rasc);
    cudaFree(d_real_decl);
    cudaFree(d_rand_rasc);
    cudaFree(d_rand_decl);
    
    cudaFree(d_hist);

// checking to see if your histograms have the right number of entries
   histogramDRsum = 0L;
   for ( int i = 0; i < binsperdegree*totaldegrees;++i ) histogramDRsum += histogramDR[i];
   printf("   DR histogram sum = %ld\n",histogramDRsum);
   //if ( histogramDRsum != 10000000000L ) {printf("   Incorrect histogram sum, exiting..\n");return(0);}

   histogramDDsum = 0L;
   for ( int i = 0; i < binsperdegree*totaldegrees;++i )
        histogramDDsum += histogramDD[i];
   printf("   DD histogram sum = %ld\n",histogramDDsum);
   //if ( histogramDDsum != 10000000000L ) {printf("   Incorrect histogram sum, exiting..\n");return(0);}

   histogramRRsum = 0L;
   for ( int i = 0; i < binsperdegree*totaldegrees;++i )
        histogramRRsum += histogramRR[i];
   printf("   RR histogram sum = %ld\n",histogramRRsum);
   if ( histogramRRsum != 10000000000L ) {printf("   Incorrect histogram sum, exiting..\n");return(0);}


   printf("   Omega values:");

   outfil = fopen(argv[3],"w");
   if ( outfil == NULL ) {printf("Cannot open output file %s\n",argv[3]);return(-1);}
   fprintf(outfil,"bin start\tomega\t        hist_DD\t        hist_DR\t        hist_RR\n");
   for ( int i = 0; i < binsperdegree*totaldegrees; ++i )
       {
       if ( histogramRR[i] > 0 )
          {
          double omega =  (histogramDD[i]-2*histogramDR[i]+histogramRR[i])/((double)(histogramRR[i]));

          fprintf(outfil,"%6.3f\t%15lf\t%15ld\t%15ld\t%15ld\n",((float)i)/binsperdegree, omega,
             histogramDD[i], histogramDR[i], histogramRR[i]);
          if ( i < 5 ) printf("   %6.4lf",omega);
          }
       else
          if ( i < 5 ) printf("         ");
       }

   printf("\n");

   fclose(outfil);

   printf("   Results written to file %s\n",argv[3]);
   printf("   CPU memory allocated  = %.2lf MB\n",CPUMemory/1000000.0);
   printf("   GPU memory allocated  = %.2lf MB\n",GPUMemory/1000000.0);

   gettimeofday(&_ttime, &_tzone);
   walltime = (double)(_ttime.tv_sec) + (double)(_ttime.tv_usec/1000000.0) - walltime;

   printf("   Total wall clock time = %.2lf s\n", walltime);

   return(0);
}

int readdata(char *argv1, char *argv2)
{
  int    i,linecount;
  char   inbuf[80];
  double ra, dec, dpi;
  FILE  *infil;
                                         
  printf("   Assuming data is in arc minutes!\n");
                          // phi   = ra/60.0 * dpi/180.0;
                          // theta = (90.0-dec/60.0)*dpi/180.0;
                          // otherwise use 
                          // phi   = ra * dpi/180.0;
                          // theta = (90.0-dec)*dpi/180.0;

  dpi = acos(-1.0);
  infil = fopen(argv1,"r");
  if ( infil == NULL ) {printf("Cannot open input file %s\n",argv1);return(-1);}

  linecount =0;
  while ( fgets(inbuf,80,infil) != NULL ) ++linecount;
  rewind(infil);

  printf("   %s contains %d galaxies\n",argv1, linecount-1);

  NoofReal = linecount-1;

  if ( NoofReal != 100000 ) {printf("Incorrect number of galaxies\n");return(1);}

  real_rasc = (float *)calloc(NoofReal,sizeof(float));
  real_decl = (float *)calloc(NoofReal,sizeof(float));
  CPUMemory += 2L*NoofReal*sizeof(float);

  fgets(inbuf,80,infil);
  sscanf(inbuf,"%d",&linecount);
  if ( linecount != 100000 ) {printf("Incorrect number of galaxies\n");return(1);}

  i = 0;
  while ( fgets(inbuf,80,infil) != NULL )
      {
      if ( sscanf(inbuf,"%lf %lf",&ra,&dec) != 2 ) 
         {
         printf("   Cannot read line %d in %s\n",i+1,argv1);
         fclose(infil);
         return(-1);
         }
      real_rasc[i] = (float)( ra/60.0*dpi/180.0);
      real_decl[i] = (float)(dec/60.0*dpi/180.0);
      ++i;
      }

  fclose(infil);

  if ( i != NoofReal ) 
      {
      printf("   Cannot read %s correctly\n",argv1);
      return(-1);
      }

  infil = fopen(argv2,"r");
  if ( infil == NULL ) {printf("Cannot open input file %s\n",argv2);return(-1);}

  linecount =0;
  while ( fgets(inbuf,80,infil) != NULL ) ++linecount;
  rewind(infil);

  printf("   %s contains %d galaxies\n",argv2, linecount-1);

  NoofRand = linecount-1;
  if ( NoofRand != 100000 ) {printf("Incorrect number of random galaxies\n");return(1);}

  rand_rasc = (float *)calloc(NoofRand,sizeof(float));
  rand_decl = (float *)calloc(NoofRand,sizeof(float));
  CPUMemory += 2L*NoofRand*sizeof(float);

  fgets(inbuf,80,infil);
  sscanf(inbuf,"%d",&linecount);
  if ( linecount != 100000 ) {printf("Incorrect number of random galaxies\n");return(1);}

  i =0;
  while ( fgets(inbuf,80,infil) != NULL )
      {
      if ( sscanf(inbuf,"%lf %lf",&ra,&dec) != 2 ) 
         {
         printf("   Cannot read line %d in %s\n",i+1,argv2);
         fclose(infil);
         return(-1);
         }
      rand_rasc[i] = (float)( ra/60.0*dpi/180.0);
      rand_decl[i] = (float)(dec/60.0*dpi/180.0);
      ++i;
      }

  fclose(infil);

  if ( i != NoofReal ) 
      {
      printf("   Cannot read %s correctly\n",argv2);
      return(-1);
      }

  return(0);
}




int getDevice(void)
{

  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  printf("   Found %d CUDA devices\n",deviceCount);
  if ( deviceCount < 0 || deviceCount > 128 ) return(-1);
  int device;
  for (device = 0; device < deviceCount; ++device) {
       cudaDeviceProp deviceProp;
       cudaGetDeviceProperties(&deviceProp, device);
       printf("      Device %s                  device %d\n", deviceProp.name,device);
       printf("         compute capability           =         %d.%d\n", deviceProp.major, deviceProp.minor);
       printf("         totalGlobalMemory            =        %.2lf GB\n", deviceProp.totalGlobalMem/1000000000.0);
       printf("         l2CacheSize                  =    %8d B\n", deviceProp.l2CacheSize);
       printf("         regsPerBlock                 =    %8d\n", deviceProp.regsPerBlock);
       printf("         multiProcessorCount          =    %8d\n", deviceProp.multiProcessorCount);
       printf("         maxThreadsPerMultiprocessor  =    %8d\n", deviceProp.maxThreadsPerMultiProcessor);
       printf("         sharedMemPerBlock            =    %8d B\n", (int)deviceProp.sharedMemPerBlock);
       printf("         warpSize                     =    %8d\n", deviceProp.warpSize);
       printf("         clockRate                    =    %8.2lf MHz\n", deviceProp.clockRate/1000.0);
       printf("         maxThreadsPerBlock           =    %8d\n", deviceProp.maxThreadsPerBlock);
       printf("         asyncEngineCount             =    %8d\n", deviceProp.asyncEngineCount);
       printf("         f to lf performance ratio    =    %8d\n", deviceProp.singleToDoublePrecisionPerfRatio);
       printf("         maxGridSize                  =    %d x %d x %d\n",
                          deviceProp.maxGridSize[0], deviceProp.maxGridSize[1], deviceProp.maxGridSize[2]);
       printf("         maxThreadsDim                =    %d x %d x %d\n",
                          deviceProp.maxThreadsDim[0], deviceProp.maxThreadsDim[1], deviceProp.maxThreadsDim[2]);
       printf("         concurrentKernels            =    ");
       if(deviceProp.concurrentKernels==1) printf("     yes\n"); else printf("    no\n");
       printf("         deviceOverlap                =    %8d\n", deviceProp.deviceOverlap);
       if(deviceProp.deviceOverlap == 1)
       printf("            Concurrently copy memory/execute kernel\n");
       }

    cudaSetDevice(0);
    cudaGetDevice(&device);
    if ( device != 0 ) printf("   Unable to set device 0, using %d instead",device);
    else printf("   Using CUDA device %d\n\n", device);

return(0);
}

