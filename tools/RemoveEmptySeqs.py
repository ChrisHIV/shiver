#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script removes from a fasta file any sequences
which contain nothing but "-", "?" and "N" characters. Output is printed to
stdout, suitable for redirection to a new fasta file.'''

import argparse
import os
import sys
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
args = parser.parse_args()

OutSeqs = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  empty = True
  for base in str(seq.seq):
    if not base in ["?", "-", "N"]:
      empty = False
      break
  if not empty:
    OutSeqs.append(seq)
    continue

SeqIO.write(OutSeqs, sys.stdout, "fasta")
