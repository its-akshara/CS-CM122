#! /bin/bash

if [ "$1" = "1" ]
then
    python3 basic_assembly.py -r spectrum_A_1/reads_spectrum_A_1_chr_1.txt -o spectrum_A_1_output.txt -t spectrum_A_1_chr_1
elif [ "$1" = "2" ]
then
    python3 basic_assembly.py -r practice_A_2/reads_practice_A_2_chr_1.txt -o practice_A_2_output.txt -t practice_A_2_chr_1
else
    python3 basic_assembly.py -r reads_hw3all_A_3_chr_1.txt -o hw3all_A_3_chr_1.txt -t hw3all_A_3_chr_1
fi