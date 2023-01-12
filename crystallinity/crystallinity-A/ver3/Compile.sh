#!/bin/bash
# please rewrite include and lib path to yours

g++ -std=c++17 -I/home/fujii/local/include/ -L/home/fujii/local/lib64 crystallinity-A.cpp -lyaml-cpp -o Crystallinity-A.o
