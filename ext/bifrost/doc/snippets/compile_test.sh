clang++ -Weverything -fexceptions -O3 -std=c++11 -march=native -Wno-c++98-compat -c test.cpp -o test.o
clang++ -o test test.o -s -O3 -lbifrost -pthread -lz
