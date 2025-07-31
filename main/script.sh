#!/bin/sh
make clean
make
echo alpha-n11-1.txt beta-n11-1.txt norm-n11-1.txt 1 11 1 0.0 10.0 | ./exec
echo alpha-n11-2.txt beta-n11-2.txt norm-n11-2.txt 1 11 2 0.0 10.0 | ./exec
echo alpha-n11-3.txt beta-n11-3.txt norm-n11-3.txt 1 11 3 0.0 10.0 | ./exec
echo alpha-n11-4.txt beta-n11-4.txt norm-n11-4.txt 1 11 4 0.0 10.0 | ./exec
echo alpha-n11-5.txt beta-n11-5.txt norm-n11-5.txt 1 11 5 0.0 10.0 | ./exec
echo alpha-n11-6.txt beta-n11-6.txt norm-n11-6.txt 1 11 6 0.0 10.0 | ./exec
echo alpha-n11-7.txt beta-n11-7.txt norm-n11-7.txt 1 11 7 0.0 10.0 | ./exec
echo alpha-n11-8.txt beta-n11-8.txt norm-n11-8.txt 1 11 8 0.0 10.0 | ./exec
echo alpha-n11-9.txt beta-n11-9.txt norm-n11-9.txt 1 11 9 0.0 10.0 | ./exec
echo alpha-n11-10.txt beta-n11-10.txt norm-n11-10.txt 1 11 10 0.0 10.0 | ./exec
echo alpha-n11-11.txt beta-n11-11.txt norm-n11-11.txt 1 11 11 0.0 10.0 | ./exec
