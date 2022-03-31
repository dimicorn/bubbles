### Program uses Boost library, installation guide:
$ wget -O boost_1_78_0.tar.gz https://boostorg.jfrog.io/artifactory/main/release/1.78.0/source/boost_1_78_0.tar.gz   
$ tar xzvf boost_1_78_0.tar.gz  
#### Moving boost folder from boost_1_78_0 to /usr/local/include
$ mv /path/to/boost_1_78_0/boost /usr/local/include  
#### If you do not have access to /usr/local/include:
$ cd path/to/boost_1_78_0  
$ ./bootstrap.sh  
$ ./b2 install  
#### Check that /data folder exists in the cloned repository folder!!!  
#### That is the directory, where the output of the program will be stored.  
### To run the program (single thread):
$ g++ main.cpp class.cpp -o out  
$ ./out  
### Running with multiple threads (using MPI):
$ mpic++ -I path/to/boost_1_78_0 main2.cpp class2.cpp -o out  
$ mpirun --use-hwthread-cpus -np number_of_threads ./out
