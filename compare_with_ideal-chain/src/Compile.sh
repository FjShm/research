#!/bin/bash
# plese rewrite include and lib path to yours

g++ -I/home/fujii/local/include/ -L/home/fujii/local/lib64 compare_with_ideal_chain.cpp -lyaml-cpp -o Compare_with_ideal_chain.o
