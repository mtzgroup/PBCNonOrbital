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

    size_t MaxAlloc() const {
        int i_device; cudaGetDevice(&i_device);
        struct cudaDeviceProp gpu_properties; cudaGetDeviceProperties(&gpu_properties, i_device);
        return gpu_properties.totalGlobalMem;
    }
    size_t Align() const { return 0; }
};

class VecPool {
private:
    int n;
    double* array;
public:
    VecPool(const int n_thread, const int array_size) { array = new double[array_size]; n = array_size; }
    ~VecPool() { delete[] array; }

    void zero() { memset(array, 0, n * sizeof(double)); }
    void reduce(double* target) const { for (int i = 0; i < n; i++) target[i] += array[i]; }

    double* checkout() const { return array; }
    void checkin(double* target) const {}
};

#endif