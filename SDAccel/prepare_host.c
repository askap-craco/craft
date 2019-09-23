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
  for(i=0; i<NAVERAGE; i++){
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

  /* Do the calculation on card */
  if (argc != 3) {
    printf("Usage: %s xclbin\n", argv[0]);
    return EXIT_FAILURE;
  }	
  const char* target_device_name = argv[2];
  
  cl_device_id device_id;            
  cl_context context;                
  cl_command_queue commands;         
  cl_program program;                
  cl_kernel kernel_prepare;       
 
  cl_mem buffer_in;                 
  cl_mem buffer_cal;
  cl_mem buffer_sky;
  cl_mem buffer_out;
  cl_mem buffer_average;
  cl_mem pt[5];
  cl_int err;
  
  device_id = get_device_id(target_device_name);
  context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
  if (!context) {
    printf("Error: Failed to create a compute context!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
  if (!commands) {
    printf("Error: Failed to create a command commands!\n");
    printf("Error: code %i\n",err);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  
  cl_int status;
  unsigned char *kernelbinary;
  char *xclbin = argv[1];
  printf("INFO: loading xclbin %s\n", xclbin);
  cl_uint n_i0 = load_file_to_memory(xclbin, (char **) &kernelbinary);
  if (n_i0 < 0) {
    printf("failed to load kernel from xclbin: %s\n", xclbin);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }

  size_t n0 = n_i0;
  program = clCreateProgramWithBinary(context, 1, &device_id, &n0,
				      (const unsigned char **) &kernelbinary, &status, &err);
  free(kernelbinary);
  
  if ((!program) || (err!=CL_SUCCESS)) {
    printf("Error: Failed to create compute program from binary %d!\n", err);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  
  err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  if (err != CL_SUCCESS) {
    size_t len;
    char buffer[2048];

    printf("Error: Failed to build program executable!\n");
    clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
    printf("%s\n", buffer);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  
  kernel_prepare = clCreateKernel(program, "prepare", &err);
  if (!kernel_prepare || err != CL_SUCCESS) {
    printf("Error: Failed to create compute kernel_vector_add!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  buffer_in = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_float) * ndata_in, source_in, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }
  buffer_sky = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_float) * ndata_sky, source_sky, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }
  buffer_cal = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_float) * ndata_cal, source_cal, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }

  buffer_out = clCreateBuffer(context,  CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_float) * ndata_out, source_out, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }  
  buffer_average = clCreateBuffer(context,  CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_float) * ndata_average, source_average, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }
  
  if (!(buffer_in&&buffer_out&&buffer_sky&&buffer_cal&&buffer_average)) {
    printf("Error: Failed to allocate device memory!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }

  pt[0] = buffer_in;
  pt[1] = buffer_cal;
  pt[2] = buffer_sky;
  pt[3] = buffer_out;
  pt[4] = buffer_average;

  struct timespec start_kernel;
  struct timespec stop_kernel;

  clock_gettime(CLOCK_REALTIME, &start_kernel);
  err = clEnqueueMigrateMemObjects(commands,(cl_uint)3,pt, 0 ,0,NULL, NULL);
  
  err = 0;
  err |= clSetKernelArg(kernel_prepare, 0, sizeof(cl_mem), &buffer_in); 
  err |= clSetKernelArg(kernel_prepare, 1, sizeof(cl_mem), &buffer_cal); 
  err |= clSetKernelArg(kernel_prepare, 2, sizeof(cl_mem), &buffer_sky);
  err |= clSetKernelArg(kernel_prepare, 3, sizeof(cl_mem), &buffer_out); 
  err |= clSetKernelArg(kernel_prepare, 4, sizeof(cl_mem), &buffer_average);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel_vector_add arguments! %d\n", err);
    printf("Test failed\n");
 
  }
  
  err = clEnqueueTask(commands, kernel_prepare, 0, NULL, NULL);
  if (err) {
    printf("Error: Failed to execute kernel! %d\n", err);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  
  err = 0;
  err |= clEnqueueMigrateMemObjects(commands,(cl_uint)2,&pt[3], CL_MIGRATE_MEM_OBJECT_HOST,0,NULL, NULL);
  
  if (err != CL_SUCCESS) {
    printf("Error: Failed to write to source array: %d!\n", err);
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  
  err = clFinish(commands);
  clock_gettime(CLOCK_REALTIME, &stop_kernel);
  elapsed_time = (stop_kernel.tv_sec - start_kernel.tv_sec) + (stop_kernel.tv_nsec - start_kernel.tv_nsec)/1.0E9L;
  fprintf(stdout, "Elapsed time of kernel is %f seconds \n", elapsed_time);
  
  /* Free memory */
  clReleaseMemObject(buffer_in);
  clReleaseMemObject(buffer_sky);
  clReleaseMemObject(buffer_cal);
  clReleaseMemObject(buffer_out);
  clReleaseMemObject(buffer_average);
  
  free(source_in);
  free(source_out);
  free(source_cal);
  free(source_sky);
  free(source_average);
  
  clReleaseProgram(program);
  clReleaseKernel(kernel_prepare);
  clReleaseCommandQueue(commands);
  clReleaseContext(context);
  
  return EXIT_SUCCESS;
}
