#!/bin/bash
# please rewrite include and lib path to yours

g++ -I~/local/include/ -L~/local/lib64 r2n.cpp -lyaml-cpp -o R2n.o
