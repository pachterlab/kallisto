#!/bin/sh

# this script cleans up an ec_vecs.kal file so that it can be read into
# kallisto c++

sed -e "s/,//g" -e "s/\|//g" $1
