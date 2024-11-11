// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <iostream>
#include <cuda_runtime.h>

int main (int /*argc*/, char** /*argv*/)
{
  //
  // get number of devices in the system, and perform proper error checking:
  // all CUDA API functions return a cudaError_t "object"
  //
  int numDevices(-1);
  cudaError_t error = cudaGetDeviceCount(&numDevices);
  if (error != cudaSuccess)
  {
    std::cerr << "CUDA ERROR (cudaGetDeviceCount): " << cudaGetErrorString(error) << ". Exiting..." << "\n";
    return 1;
  }
  std::cout << "Number of devices in this machine: " << numDevices << "\n" << "\n";

  //
  // manually treat the case that this is run on a system without a CUDA-capable
  // GPU that has the toolkit installed (like our servers, which share the software
  // stack even though not all of them have a GPU)
  //
  if (numDevices == 0)
  {
    std::cerr << "No CUDA-capable device found." << "\n";
    return 2;
  }

  //
  // loop over all devices and query (some of) their properties
  //
  for (int idevice(0); idevice<numDevices; ++idevice)
  {
    // set CUDA device
    error = cudaSetDevice(idevice);
    if (error != cudaSuccess)
    {
      std::cerr << "CUDA ERROR (cudaSetDevice): " << cudaGetErrorString(error) << ". Exiting..." << "\n";
      return 3;
    }
    // get device properties
    cudaDeviceProp prop;
    error = cudaGetDeviceProperties (&prop, idevice);
    if (error != cudaSuccess)
    {
      std::cerr << "CUDA ERROR (cudaGetDeviceProperties): " << cudaGetErrorString(error) << ". Exiting..." << "\n";
      return 4;
    }
    // print out device name and compute capabilities
    std::cout << "Device " << idevice << ": " << prop.name;
    std::cout << " (cc " << prop.major << "." << prop.minor << ")" << "\n";
    std::cout << "Managed Memory Support: " << ((prop.managedMemory == 1) ? "yes" : "no") << "\n";
    std::cout << "Important: configure all future builds with \" --cuda_arch=sm_";
    std::cout << prop.major << prop.minor << "\" on this machine." << "\n";

    if (idevice + 1 < numDevices)
      std::cout<<"\n";
  }
  return 0;
}
