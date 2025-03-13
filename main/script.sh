#!/bin/sh
make clean
make
echo alpha-n512-1.txt beta-n512-1.txt norm-n512-1.txt 1 513 12 0.0 10.0 | ./exec 

