#!/bin/bash
# please rewrite include and lib path to yours
#
# About compiler
# $ g++ --version
#   > g++ (GCC) 9.3.1 20200408 (Red Hat 9.3.1-2)
#   > Copyright (C) 2019 Free Software Foundation, Inc.
#   > This is free software; see the source for copying conditions.  There is NO
#   > warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


g++ -std=c++17 -I~/local/include/ -L~/local/lib64 cluster.cpp -lyaml-cpp -o Cluster.o

