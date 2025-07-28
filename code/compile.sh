#!/bin/bash

filename=$1

g++ -I/home/jangh/eigen-3.4.0 -o ${filename} ${filename}.cc -O2 -std=c++17

echo "Compling is done!"
