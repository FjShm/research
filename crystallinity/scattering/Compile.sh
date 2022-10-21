#!/bin/bash
# please rewrite include and lib path to yours

g++ -fopenmp -I/home/fujii/local/include/ -L/home/fujii/local/lib64 scattering.cpp -lyaml-cpp -o Scattering.o
