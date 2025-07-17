#!/bin/sh
make clean
make
echo alpha-n16-1.txt beta-n16-1.txt norm-n512-1.txt 1 17 1 | ./exec 

