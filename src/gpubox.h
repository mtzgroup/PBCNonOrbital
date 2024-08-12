#ifndef GPUBOX_H_
#define GPUBOX_H_

#include <cuda_runtime.h>

#include <string.h>

class GPUBox {
public:
    static int NumGPUstreams() { return 1; }
    static GPUBox* Checkout() { return new GPUBox(); }
    static void Checkin(GPUBox* gpu) { delete gpu; }

    void* cpuAlloc(size_t size) const { return malloc(size); }
    void* gpuAlloc(size_t size) const { char* p = NULL; cudaMalloc(&p, size); return p; }

    void cpuFree(void* p) const { free(p); }
    void gpuFree(void* p) const { cudaFree(p); }
};

#endif