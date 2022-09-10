#!/bin/bash
# plese rewrite include and lib path to yours

g++ -I/home/fujii/local/include/ -L/home/fujii/local/lib64 crystallinity-C.cpp -lyaml-cpp -o Crystallinity-C.o
