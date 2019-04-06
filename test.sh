#!/bin/bash

## test run

# from stdin, CTCF DNA sequences
cat data/ctcf.meme | ./motto.py -s compact

# from inFile, lipoprotein, Amino acid sequences
./motto.py -i data/lipoAA.meme



