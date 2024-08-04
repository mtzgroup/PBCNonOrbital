
g++ -fPIC -shared src/periodic_cutoff.cpp src/periodic_lattice.cpp src/periodic_nuclear_repulsion.cpp src/helper.cpp -o periodic_nuclear_repulsion.so

g++ example/benzene.cpp periodic_nuclear_repulsion.so -lblas -o test.exe
