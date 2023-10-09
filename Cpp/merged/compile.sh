#!/bin/bash

flag='-p'

if [[ -z "$1" && -z "$2" ]]; then
	echo "No files for compilation"
	exit 1
elif [[ -z "$1" || -z "$2" ]]; then
	echo "No file for compilation"
	exit 1
fi

if [ -z "$3" ]; then
	g++ $1 $2 -o out -std=c++14 -O3 -Wno-deprecated-declarations -I/opt/homebrew/Cellar/open-mpi/4.1.5/include -I/opt/homebrew/Cellar/boost/1.82.0_1/include -L/opt/homebrew/Cellar/open-mpi/4.1.5/lib -lmpi
elif [ "$3" = $flag ]; then
	mpic++ $1 $2 -o out -std=c++14 -O3 -Wno-deprecated-declarations -I/opt/homebrew/Cellar/open-mpi/4.1.5/include -I/opt/homebrew/Cellar/boost/1.82.0_1/include -L/opt/homebrew/Cellar/open-mpi/4.1.5/lib
fi
