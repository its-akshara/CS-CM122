#! /bin/bash

if [ "$1" = "1" ]
then
    python3 basic_aligner.py -g practice_W_1/ref_practice_W_1_chr_1.txt -r practice_W_1/reads_practice_W_1_chr_1.txt -o test_output.txt -t practice_W_1_chr_1
else
    python3 basic_aligner.py -g hw1_W_2/ref_hw1_W_2_chr_1.txt -r hw1_W_2/reads_hw1_W_2_chr_1.txt -o test_output.txt -t hw1_W_2_chr_1
fi