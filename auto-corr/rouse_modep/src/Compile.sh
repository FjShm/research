#!/bin/bash
# please rewrite include and lib path to yours

g++ -std=c++17 -fconcepts -I/home/shoma/local/include/ -L/home/shoma/local/lib rouse_modep.cpp -lyaml-cpp -o Rouse_modep.o
