#!/bin/bash
g++ $1 $2 -I/opt/homebrew/Cellar/boost/1.82.0_1/include -std=c++14 -o out -w
