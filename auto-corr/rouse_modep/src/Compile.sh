#!/bin/bash
# please rewrite include and lib path to yours

g++ -std=c++17 -fconcepts -I/home/fujii/local/include/ -L/home/fujii/local/lib64 rouse_modep.cpp -lyaml-cpp -o Rouse_modep.o
