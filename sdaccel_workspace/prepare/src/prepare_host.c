/*******************************************************************************
 ** HOST Code
 *******************************************************************************/

#include "prepare_host.h"
#include "util_sdaccel.h"
#include "prepare.h"

int main(int argc, char* argv[])
{
  /* Prepare the data arraies */
  cl_uint nsamp1;
  cl_uint npol_in;
  cl_uint npol_out;
  cl_uint nsamp2;
  cl_uint npol_average;
  cl_uint npol_sky;
  cl_uint npol_cal;

  nsamp1       = NSAMP_PER_TIME;
  npol_average = nsamp1*2;
  npol_cal     = nsamp1*2;
  npol_sky     = nsamp1;
  
  nsamp2   = NTIME_PER_BUFBLOCK*nsamp1;
  npol_in  = nsamp2*2;
  npol_out = nsamp2;

  cl_complex *source_in = NULL;
  cl_complex *source_sw_out = NULL;
  cl_complex *source_hw_out = NULL; 
  cl_complex *source_cal = NULL;
  cl_complex *source_sky = NULL;
  cl_complex *source_sw_average = NULL;
  cl_complex *source_hw_average = NULL;

  source_in     = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_in*sizeof(cl_complex));
  source_sw_out = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_out*sizeof(cl_complex));
  source_hw_out = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_out*sizeof(cl_complex));
  source_cal    = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_cal*sizeof(cl_complex));
  source_sky    = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_sky*sizeof(cl_complex));
  source_sw_average = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_average*sizeof(cl_complex));
  source_hw_average = (cl_complex *)aligned_alloc(MEM_ALIGNMENT, npol_average*sizeof(cl_complex));
  printf("INFO: %f MB memory used on host in total\n",
	  (npol_in + 2*npol_out + npol_cal + npol_sky + 2*npol_average) *
	  sizeof(cl_complex)/(1024.*1024.));
  printf("INFO: %f MB memory used on device in total\n",
	  (npol_in + npol_out + npol_cal + npol_sky + npol_average) *
	  sizeof(cl_complex)/(1024.*1024.));
  printf("INFO: %f MB memory used on device for raw input\n",
	  npol_in*sizeof(cl_complex)/(1024.*1024.));  
  printf("INFO: %f MB memory used on device for raw output\n",
	  npol_out*sizeof(cl_complex)/(1024.*1024.));  
  printf("INFO: %f MB memory used on device for average output\n",
	  npol_average*sizeof(cl_complex)/(1024.*1024.));  
  printf("INFO: %f MB memory used on device for calibration input\n",
	  npol_cal*sizeof(cl_complex)/(1024.*1024.));
  printf("INFO: %f MB memory used on device for sky model\n",
	  npol_sky*sizeof(cl_complex)/(1024.*1024.));

  /* Do the calculation */
  srand(time(NULL));
  cl_uint i;
  cl_float elapsed_time;
  struct timespec start_host;
  struct timespec stop_host;
    
  /* Prepare the sky model and calibration data */
  for(i=0; i<npol_sky; i++) {
    source_sky[i].real((float)(rand() % DATA_RANGE_RAND));
    source_sky[i].imag((float)(rand() % DATA_RANGE_RAND));
  }
  for(i=0; i<npol_cal; i++){
    source_cal[i].real((float)(rand() % DATA_RANGE_RAND));
    source_cal[i].imag((float)(rand() % DATA_RANGE_RAND));
  }
  
  /* Prepare the input raw data */
  for(i=0; i<npol_in; i++){
    source_in[i].real((float)(rand() % DATA_RANGE_RAND));
    source_in[i].imag((float)(rand() % DATA_RANGE_RAND));
  }
  
  /* Calculate on host */
  memset(source_sw_average, 0x00, npol_average*sizeof(cl_complex)); // Get memory reset for the average
  memset(source_hw_average, 0x00, npol_average*sizeof(cl_complex)); // Get memory reset for the average
  memset(source_sw_out, 0x00, npol_out*sizeof(cl_complex)); // Get memory reset for the out
  memset(source_hw_out, 0x00, npol_out*sizeof(cl_complex)); // Get memory reset for the out
  clock_gettime(CLOCK_REALTIME, &start_host);
  prepare(source_in, source_cal, source_sky, source_sw_out, source_sw_average);
  clock_gettime(CLOCK_REALTIME, &stop_host);
  elapsed_time = (stop_host.tv_sec - start_host.tv_sec) + (stop_host.tv_nsec - start_host.tv_nsec)/1.0E9L;
  printf("INFO: elapsed_time of one loop is %E seconds\n", elapsed_time);
  
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
  
  kernel_prepare = clCreateKernel(program, "prepare_knl", &err);
  if (!kernel_prepare || err != CL_SUCCESS) {
    printf("Error: Failed to create compute kernel_vector_add!\n");
    printf("Test failed\n");
    return EXIT_FAILURE;
  }
  buffer_in = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_complex)*npol_in, source_in, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }
  buffer_sky = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_complex)*npol_sky, source_sky, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }
  buffer_cal = clCreateBuffer(context,  CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_complex)*npol_cal, source_cal, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }

  buffer_out = clCreateBuffer(context,  CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_complex)*npol_out, source_hw_out, &err);
  if (err != CL_SUCCESS) {
    std::cout << "Return code for clCreateBuffer - in1" << err << std::endl;
  }  
  buffer_average = clCreateBuffer(context,  CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,  sizeof(cl_complex)*npol_average, source_hw_average, &err);
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

  struct timespec start_device;
  struct timespec stop_device;

  clock_gettime(CLOCK_REALTIME, &start_device);
  err = clEnqueueMigrateMemObjects(commands,(cl_uint)3,pt, 0 ,0,NULL, NULL);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel_prepare arguments! %d\n", err);
    printf("Test failed\n");
  }
  
  err = 0;
  err |= clSetKernelArg(kernel_prepare, 0, sizeof(cl_mem), &buffer_in); 
  err |= clSetKernelArg(kernel_prepare, 1, sizeof(cl_mem), &buffer_cal); 
  err |= clSetKernelArg(kernel_prepare, 2, sizeof(cl_mem), &buffer_sky);
  err |= clSetKernelArg(kernel_prepare, 3, sizeof(cl_mem), &buffer_out); 
  err |= clSetKernelArg(kernel_prepare, 4, sizeof(cl_mem), &buffer_average);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to set kernel_prepare arguments! %d\n", err);
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
  clock_gettime(CLOCK_REALTIME, &stop_device);
  elapsed_time = (stop_device.tv_sec - start_device.tv_sec) + (stop_device.tv_nsec - start_device.tv_nsec)/1.0E9L;
  fprintf(stdout, "Elapsed time of kernel is %E seconds \n", elapsed_time);

  /*
  for(i=0;i<npol_average;i++){
    if(source_sw_average[i]!=source_hw_average[i])
      std::cout << "Mismatch on average: " <<i << '\t' << source_sw_average[i]<< '\t'<< source_hw_average[i] << '\n';
  }
  
  for(i=0;i<npol_out;i++){
    if(source_sw_out[i]!=source_hw_out[i]){
      std::cout << "Mismatch on out: " <<i << '\t'<< source_sw_out[i]<< '\t'<< source_hw_out[i] << '\n';
      //break;
    }
  }
  */
  
  /* Free memory */
  clReleaseMemObject(buffer_in);
  clReleaseMemObject(buffer_sky);
  clReleaseMemObject(buffer_cal);
  clReleaseMemObject(buffer_out);
  clReleaseMemObject(buffer_average);
  
  free(source_in);
  free(source_sw_out);
  free(source_hw_out);
  free(source_cal);
  free(source_sky);
  free(source_sw_average);
  free(source_hw_average);
  
  clReleaseProgram(program);
  clReleaseKernel(kernel_prepare);
  clReleaseCommandQueue(commands);
  clReleaseContext(context);
  
  return EXIT_SUCCESS;
}
