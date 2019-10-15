/*******************************************************************************
** SDACCEL UTIL Code FILE
*******************************************************************************/

#include "util_sdaccel.h"

cl_uint load_file_to_memory(const char *filename, char **result)
{
	cl_uint size = 0;
	FILE *f = fopen(filename, "rb");
	if (f == NULL) {
        *result = NULL;
        return -1; // -1 means file opening fail
    }
    fseek(f, 0, SEEK_END);
    size = ftell(f);
    fseek(f, 0, SEEK_SET);
    *result = (char *)malloc(size+1);
    if (size != fread(*result, sizeof(char), size, f)) {
        free(*result);
        return -2; // -2 means file reading fail
    }
    fclose(f);
    (*result)[size] = 0;
    return size;
}

cl_device_id get_device_id(const char* target_device_name)
{
  cl_platform_id platforms[16];       
  cl_platform_id platform_id;
  cl_uint platform_count;
  cl_uint platform_found = 0;
  
  cl_device_id devices[16];
  cl_device_id device_id;
  cl_uint device_found = 0;

  char cl_platform_vendor[1001];
  
  cl_uint num_devices;
  char cl_device_name[1001];
  cl_int err;
  cl_uint iplat;
  cl_uint i;
  
  // ------------------------------------------------------------------------------------
  // Step 1: Get All PLATFORMS, then search for Target_Platform_Vendor (CL_PLATFORM_VENDOR)
  // ------------------------------------------------------------------------------------
	
  // Get the number of platforms
  // ..................................................
	  
  err = clGetPlatformIDs(16, platforms, &platform_count);
  if (err != CL_SUCCESS) {
    printf("Error: Failed to find an OpenCL platform!\n");
    printf("Test failed\n");
    exit(EXIT_FAILURE);
  }

  printf("INFO: Found %d platforms\n", platform_count);
 
  // ....................................................................................
  // step 1:  Search for Platform (ex: Xilinx) using: CL_PLATFORM_VENDOR = Target_Platform_Vendor
  // Check if the current platform matches Target_Platform_Vendor
  // ....................................................................................
 	
  for (iplat=0; iplat<platform_count; iplat++) {
    err = clGetPlatformInfo(platforms[iplat], CL_PLATFORM_VENDOR, 1000, (void *)cl_platform_vendor,NULL);
    if (err != CL_SUCCESS) {
      printf("Error: clGetPlatformInfo(CL_PLATFORM_VENDOR) failed!\n");
      printf("Test failed\n");
      exit(EXIT_FAILURE);
    }
    if (strcmp(cl_platform_vendor, "Xilinx") == 0) {
      printf("INFO: Selected platform %d from %s\n", iplat, cl_platform_vendor);
      platform_id = platforms[iplat];
      platform_found = 1;
    }
  }
  if (!platform_found) {
    printf("ERROR: Platform Xilinx not found. Exit.\n");
    exit(EXIT_FAILURE);
  }
  
  // ------------------------------------------------------------------------------------  
  // Step 1:  Get All Devices for selected platform Target_Platform_ID
  //            then search for Xilinx platform (CL_DEVICE_TYPE_ACCELERATOR = Target_Device_Name)
  // ------------------------------------------------------------------------------------

	
  err = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_ACCELERATOR, 16, devices, &num_devices);
  printf("INFO: Found %d devices\n", num_devices);
  if (err != CL_SUCCESS) {
    printf("ERROR: Failed to create a device group!\n");
    printf("ERROR: Test failed\n");
    exit(EXIT_FAILURE);
  }
  // ------------------------------------------------------------------------------------  
  // Step 1:  Search for CL_DEVICE_NAME = Target_Device_Name
  // ............................................................................
	
  for (i=0; i<num_devices; i++) {
    err = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, 1024, cl_device_name, 0);
    if (err != CL_SUCCESS) {
      printf("Error: Failed to get device name for device %d!\n", i);
      printf("Test failed\n");
      exit(EXIT_FAILURE);
    }
    printf("CL_DEVICE_NAME %s\n", cl_device_name);

    // ............................................................................  
    // Step 1: Check if the current device matches Target_Device_Name
    // ............................................................................ 

    if(strcmp(cl_device_name, target_device_name) == 0) {
      device_id = devices[i];
      device_found = 1;
      printf("Selected %s as the target device\n", cl_device_name);
    }
  }
  if(!device_found) {
    printf("ERROR: Failed to get device %s. EXit.\n", cl_device_name);
    exit(EXIT_FAILURE);
  }

  return device_id;
}
