#!/bin/bash


domlen=1.0
diaph=0.5
cells=300
gamma=2.0
timeout=0.005 

dls=(0.001 0.01 1.0 10 100)
uls=(12.0 3.0 2.0 0.6 0.3)
pls=(0.18 2 0.02)

dr=1.0
ur=0.0
pr=0.02

mpa=1

for dl in ${dls[@]}; do
    for ul in ${uls[@]}; do
        for pl in ${pls[@]}; do
            echo $domlen $diaph $cells $gamma $timeout $dl $ul $pl $dr $ur $pr $mpa > shang.ini
            ./run.sh main.cpp riemann_solver.cpp shang.ini
            python3 plot.py exact.out $gamma $dl $ul $pl
        done
    done
done
