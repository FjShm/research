#!/bin/bash
# plese rewrite include and lib path to yours

g++ -std=c++17 -I/home/fujii/local/include/ -L/home/fujii/local/lib64 ete_auto_corr.cpp -lyaml-cpp -o Ete_auto_corr.o
