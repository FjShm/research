#!/bin/bash
# plese rewrite include and lib path to yours

g++ -I/home/fujii/local/include/ -L/home/fujii/local/lib64 rotate_dump.cpp -lyaml-cpp -o RotateDump.o
