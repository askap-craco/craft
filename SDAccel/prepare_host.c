/*******************************************************************************
 ** HOST Code
 *******************************************************************************/

#include "prepare_host.h"
#include "util_sdaccel.h"
#include "prepare.h"

int main(int argc, char* argv[])
{
  /* Prepare the data arraies */
  cl_uint nsamp_raw, npol_in, npol_out, ndata_in, ndata_out;
  cl_uint nsamp_average, npol_average, npol_sky, npol_cal, ndata_average, ndata_sky, ndata_cal;
  
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
  
  cl_float *in = NULL;
  cl_float *out = NULL; 
  cl_float *cal = NULL;
  cl_float *sky = NULL;
  cl_float *average = NULL;

  in = (cl_float *)malloc(ndata_in * sizeof(cl_float));
  out = (cl_float *)malloc(ndata_out * sizeof(cl_float));
  cal = (cl_float *)malloc(ndata_cal * sizeof(cl_float));
  sky = (cl_float *)malloc(ndata_sky * sizeof(cl_float));
  average = (cl_float *)malloc(ndata_average * sizeof(cl_float));
  fprintf(stdout, "%f MB memory used in total\n",
	  (ndata_in + ndata_out + ndata_cal + ndata_sky + ndata_average) *
	  sizeof(cl_float)/(1024.*1024.));
  fprintf(stdout, "%f MB memory used for raw input\n",
	  ndata_in * sizeof(cl_float)/(1024.*1024.));  
  fprintf(stdout, "%f MB memory used for raw output\n",
	  ndata_out * sizeof(cl_float)/(1024.*1024.));  
  fprintf(stdout, "%f MB memory used for average output\n",
	  ndata_average * sizeof(cl_float)/(1024.*1024.));  
  fprintf(stdout, "%f MB memory used for calibration input\n",
	  ndata_cal * sizeof(cl_float)/(1024.*1024.));
  fprintf(stdout, "%f MB memory used for sky model\n",
	  ndata_sky * sizeof(cl_float)/(1024.*1024.));
  
  fflush(stdout);
  
  /* Do the calculation */
  srand(time(NULL));
  cl_uint i, j, k;
  cl_float elapsed_time;
  struct timespec start_average, stop_average;
  struct timespec start_loop, stop_loop;
  for(i = 0; i < NAVEAGE; i++){
    memset(average, 0x00, ndata_average * sizeof(cl_float)); // Get memory reset for the average 
    
    /* Prepare the sky model and calibration data */
    for(j = 0; j < nsamp_average; j++){
      sky[2*j] = rand() % DATA_SIZE;
      sky[2*j + 1] = rand() % DATA_SIZE;

      cal[4*j] = rand() % DATA_SIZE;
      cal[4*j + 1] = rand() % DATA_SIZE;
      cal[4*j + 2] = rand() % DATA_SIZE;
      cal[4*j + 3] = rand() % DATA_SIZE;
    }

    clock_gettime(CLOCK_REALTIME, &start_average);
    for(j = 0; j < NBUFBLOCK_AVERAGE; j++){
      clock_gettime(CLOCK_REALTIME, &start_loop);
      /* Prepare the input raw data */
      fprintf(stdout, "Running the %d average and %d loop, ", i, j);
      for(k = 0; k < nsamp_raw; k++){
	in[4*j] = rand() % DATA_SIZE;
	in[4*j + 1] = rand() % DATA_SIZE;
	in[4*j + 2] = rand() % DATA_SIZE;
	in[4*j + 3] = rand() % DATA_SIZE;
      }
      
      /* Calculate on host */
      prepare(in, cal, sky, out, average);
      
      clock_gettime(CLOCK_REALTIME, &stop_loop);
      elapsed_time = (stop_loop.tv_sec - start_loop.tv_sec) + (stop_loop.tv_nsec - start_loop.tv_nsec)/1.0E9L;
      fprintf(stdout, "elapsed_time of one loop is %f seconds\n", elapsed_time);
      fflush(stdout);
    }
    clock_gettime(CLOCK_REALTIME, &stop_average);
    elapsed_time = (stop_average.tv_sec - start_average.tv_sec) + (stop_average.tv_nsec - start_average.tv_nsec)/1.0E9L;
    fprintf(stdout, "elapsed_time of one average is %f seconds\n", elapsed_time);
    fflush(stdout);
  }
  
  /* Free memory */
  free(in);
  free(out);
  free(cal);
  free(sky);
  free(average);
  
  return EXIT_SUCCESS;
}
