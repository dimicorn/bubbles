# Prerequisites
1. Boost
2. Open-MPI
3. Change in compile.sh script paths to these libraries

# Using MPI
```bash
./compile.sh main.cpp bubbles.cpp -p
```
## Runge-Kutta method (default)
```
mpirun -n number_of_threads ./out 
```
## Rosenbrock method
```
mpirun -n number_of_threads ./out -r
```

# Single thread
```bash
./compile.sh main.cpp bubbles.cpp
```
## Runge-Kutta method (default)
```
./out
```
## Rosenbrock method
```
./out -r
```