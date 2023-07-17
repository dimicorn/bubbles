## Boost library installation guide:
```
$ wget -O boost_1_78_0.tar.gz https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz   
$ tar xzvf boost_1_78_0.tar.gz
```
Moving `/boost` folder from `/boost_1_78_0` to `/usr/local/include`
```
$ mv /path/to/boost_1_78_0/boost /usr/local/include
```
If you do not have access to `/usr/local/include`:
```
$ cd path/to/boost_1_78_0  
$ ./bootstrap.sh  
$ ./b2 install
```
Check that `/data` folder exists in the cloned repository folder!!!  
That is the directory, where the output of the program will be stored.

## MPI library installation guide:
```
$ wget https://download.open-mpi.org/release/open-mpi/v2.0/openmpi-2.0.4.tar.gz
$ tar -xf openmpi-2.0.4.tar.gz
$ cd openmpi-2.0.4/
$ ./configure --prefix=/usr/local
$ make all
$ make install
```
## To run the program: 
### Single thread:
```
$ g++ main.cpp class.cpp -o out
$ ./out
```
### Multiple threads (using MPI):
```
$ mpic++ -I path/to/boost_1_78_0 main2.cpp class2.cpp -o out  
$ mpirun --use-hwthread-cpus -np number_of_threads ./out
```