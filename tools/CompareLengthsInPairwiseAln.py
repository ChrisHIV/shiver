#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script takes a pairwise sequence alignment, and
reports the total number of bases (everything but the gap character "-") that
the first seq has before the second seq begins & after it finishes, and the
total number of bases the first seq has where the second seq has an internal
deletion. Both of these are signed integers, i.e. where the second seq has a
base and the first has a gap that counts as -1. We print these two values to
stdout in the order described.'''

import argparse
import os
import sys
import re
from Bio import AlignIO
import numpy as np
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
args = parser.parse_args()

# Read in the alignment
try:
  alignment = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise
AlignmentLength = alignment.get_alignment_length()

if len(alignment) != 2:
  print('Pairwise alignment expected. Quitting.', file=sys.stderr)
  exit(1)

# Our representation of each seq will be an array of bools: False for a gap char,
# True otherwise
seqs = np.array([np.array([base != '-' for base in seq.seq]) for seq in \
alignment])

def GetSeqStartAndEndPos(seq):
  '''Get the position of the first and last non-gap character in the seq.'''
  FirstBasePos = 0
  try:
    while not seq[FirstBasePos]:
      FirstBasePos += 1
  except IndexError:
    print('Encountered pure-gap sequence. Quitting', file=sys.stderr)
    quit(1)
  LastBasePos = len(seq) - 1
  while not seq[LastBasePos]:
    LastBasePos -= 1
  return FirstBasePos, LastBasePos

start1, end1 = GetSeqStartAndEndPos(seqs[0])
start2, end2 = GetSeqStartAndEndPos(seqs[1])

ExtraAtEndsFor1 = 0
ExtraAtEndsFor1 += sum(IsNotGap for IsNotGap in seqs[0][start1:start2] if IsNotGap)
ExtraAtEndsFor1 -= sum(IsNotGap for IsNotGap in seqs[0][start2:start1] if IsNotGap)
ExtraAtEndsFor1 += sum(IsNotGap for IsNotGap in seqs[0][end2:end1] if IsNotGap)
ExtraAtEndsFor1 -= sum(IsNotGap for IsNotGap in seqs[0][end1:end2] if IsNotGap)

ExtraInMiddleFor1 = 0
for pos in range(max(start1, start2), min(end1, end2)):
  Seq1isNotGap = seqs[0][pos]
  Seq2isNotGap = seqs[1][pos]
  ExtraInMiddleFor1 += int(Seq1isNotGap) - int(Seq2isNotGap)

print(ExtraAtEndsFor1, ExtraInMiddleFor1)

