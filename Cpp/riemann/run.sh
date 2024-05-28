#!/bin/bash
g++ $1 $2 -std=c++17 -o out -w && ./out $3
