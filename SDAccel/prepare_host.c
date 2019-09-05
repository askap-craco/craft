/*******************************************************************************
 ** HOST Code
 *******************************************************************************/

#include "prepare_host.h"
#include "util_sdaccel.h"
#include "prepare.h"

int main(int argc, char* argv[])
{
  /* Prepare the data arraies */
  cl_uint nsamp_raw;
  cl_uint npol_in;
  cl_uint npol_out;
  cl_uint ndata_in;
  cl_uint ndata_out;
  cl_uint nsamp_average;
  cl_uint npol_average;
  cl_uint npol_sky;
  cl_uint npol_cal;
  cl_uint ndata_average;
  cl_uint ndata_sky;
  cl_uint ndata_cal;
  
  nsamp_raw = NTIME_PER_BUFBLOCK * NSAMP_PER_TIME;
  npol_in = nsamp_raw * 2;
  npol_out = nsamp_raw;
  ndata_in = npol_in * 2;
  ndata_out = npol_out * 2;

  nsamp_average = NSAMP_PER_TIME;
  npol_average = nsamp_average * 2;
  ndata_average = npol_average * 2;

  npol_sky = nsamp_average;
  ndata_sky = npol_sky * 2;
  
  npol_cal = nsamp_average * 2;
  ndata_cal = npol_cal * 2;
  
  cl_float *source_in = NULL;
  cl_float *source_out = NULL; 
  cl_float *source_cal = NULL;
  cl_float *source_sky = NULL;
  cl_float *source_average = NULL;

  source_in = (cl_float *)aligned_alloc(MEM_ALIGNMENT, ndata_in * sizeof(cl_float));
  source_out = (cl_float *)aligned_alloc(MEM_ALIGNMENT, ndata_out * sizeof(cl_float));
  source_cal = (cl_float *)aligned_alloc(MEM_ALIGNMENT, ndata_cal * sizeof(cl_float));
  source_sky = (cl_float *)aligned_alloc(MEM_ALIGNMENT, ndata_sky * sizeof(cl_float));
  source_average = (cl_float *)aligned_alloc(MEM_ALIGNMENT, ndata_average * sizeof(cl_float));
  printf("INFO: %f MB memory used in total\n",
	  (ndata_in + ndata_out + ndata_cal + ndata_sky + ndata_average) *
	  sizeof(cl_float)/(1024.*1024.));
  printf("INFO: %f MB memory used for raw input\n",
	  ndata_in * sizeof(cl_float)/(1024.*1024.));  
  printf("INFO: %f MB memory used for raw output\n",
	  ndata_out * sizeof(cl_float)/(1024.*1024.));  
  printf("INFO: %f MB memory used for average output\n",
	  ndata_average * sizeof(cl_float)/(1024.*1024.));  
  printf("INFO: %f MB memory used for calibration input\n",
	  ndata_cal * sizeof(cl_float)/(1024.*1024.));
  printf("INFO: %f MB memory used for sky model\n",
	  ndata_sky * sizeof(cl_float)/(1024.*1024.));

  ///* Setup device */
  //cl_int err;                        
  //cl_device_id device_id;
  //if (argc != 3) {
  //  printf("Usage: %s xclbin\n", argv[0]);
  //  return EXIT_FAILURE;    
  //}	
  //const char* target_device_name = argv[2];
  //device_id = get_device_id(target_device_name);  
  
  /* Do the calculation */
  srand(time(NULL));
  cl_uint i;
  cl_uint j;
  cl_uint k;
  cl_float elapsed_time;
  struct timespec start_average;
  struct timespec stop_average;
  struct timespec start_loop;
  struct timespec stop_loop;
  for(i=0; i<NAVEAGE; i++){
    memset(source_average, 0x00, ndata_average * sizeof(cl_float)); // Get memory reset for the average 
    
    /* Prepare the sky model and calibration data */
    for(j=0; j<nsamp_average; j++){
      source_sky[2*j] = rand() % DATA_SIZE;
      source_sky[2*j + 1] = rand() % DATA_SIZE;

      source_cal[4*j] = rand() % DATA_SIZE;
      source_cal[4*j + 1] = rand() % DATA_SIZE;
      source_cal[4*j + 2] = rand() % DATA_SIZE;
      source_cal[4*j + 3] = rand() % DATA_SIZE;
    }

    clock_gettime(CLOCK_REALTIME, &start_average);
    for(j=0; j<NBUFBLOCK_AVERAGE; j++){
      clock_gettime(CLOCK_REALTIME, &start_loop);
      /* Prepare the input raw data */
      printf("INFO: Running the %d average and %d loop, ", i, j);
      for(k=0; k<nsamp_raw; k++){
	source_in[4*j] = rand() % DATA_SIZE;
	source_in[4*j + 1] = rand() % DATA_SIZE;
	source_in[4*j + 2] = rand() % DATA_SIZE;
	source_in[4*j + 3] = rand() % DATA_SIZE;
      }
      
      /* Calculate on host */
      prepare(source_in, source_cal, source_sky, source_out, source_average);
      
      clock_gettime(CLOCK_REALTIME, &stop_loop);
      elapsed_time = (stop_loop.tv_sec - start_loop.tv_sec) + (stop_loop.tv_nsec - start_loop.tv_nsec)/1.0E9L;
      printf("INFO: elapsed_time of one loop is %f seconds\n", elapsed_time);
    }
    clock_gettime(CLOCK_REALTIME, &stop_average);
    elapsed_time = (stop_average.tv_sec - start_average.tv_sec) + (stop_average.tv_nsec - start_average.tv_nsec)/1.0E9L;
    printf("INFO: elapsed_time of one average is %f seconds\n", elapsed_time);
  }
  
  /* Free memory */
  free(source_in);
  free(source_out);
  free(source_cal);
  free(source_sky);
  free(source_average);
  
  return EXIT_SUCCESS;
}
