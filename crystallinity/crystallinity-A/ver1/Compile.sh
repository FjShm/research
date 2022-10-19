#!/bin/bash
# please rewrite include and lib path to yours

#g++ -I/home/fujii/local/include/ -L/home/fujii/local/lib64 crystallinity-A.cpp -lyaml-cpp -o Crystallinity-A.o
g++ -I/home/fujii/local/include/ -L/home/fujii/local/lib64 show_crystallinity.cpp -lyaml-cpp -o Show_crystallinity.o
