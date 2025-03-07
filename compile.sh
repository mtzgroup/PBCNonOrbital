
nvcc -Xcompiler -fPIC src/periodic_becke/periodic_becke_kernel.cu -c -o periodic_becke_kernel.o &&

g++ -fPIC -shared src/periodic_cutoff.cpp src/periodic_lattice.cpp src/periodic_nuclear_repulsion.cpp src/helper.cpp src/periodic_grid.cpp src/periodic_kernel_data.cpp src/periodic_becke/periodic_becke.cpp src/periodic_becke/periodic_becke_weight.cpp src/periodic_becke/periodic_becke_gradient.cpp periodic_becke_kernel.o -o periodic_nuclear_repulsion.so &&

# g++ example/benzene_nuclear_repulsion.cpp periodic_nuclear_repulsion.so -lblas -lcudart -o test.exe
g++ example/lif_becke_weight.cpp periodic_nuclear_repulsion.so -lblas -lcudart -o test.exe
