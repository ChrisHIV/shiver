#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script analyses two input fasta files which should
contain the same set of sequences in the same order, except for differences in
(a) "?" and "-" characters, which may appear and disappear freely, and (b) the
case (lower or upper) of bases. For each sequence, we preserve the "?" and "-"
characters from the 'after' form of the sequence, but use the case from the
'before' form of the sequence. Output is written to stdout suitable for
redirection to a fasta-format file. The point of this script was the following
problem: manual editing of shiver consensuses using Geneious (where changing gap
placement also changed the number of "?" characters required at the ends) caused
loss of the meaningful upper/lower case distinction.'''

import argparse
import os
import sys
from Bio import SeqIO
from Bio import Seq
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('seqs_before', type=File)
parser.add_argument('seqs_after', type=File)
parser.add_argument('-1', '--first-seq-only', action='store_true', \
help='''Process only the first sequence in each of the two input files, leaving
the others as they are.''')

args = parser.parse_args()

seqs_before = list(SeqIO.parse(open(args.seqs_before), 'fasta'))
seqs_after  = list(SeqIO.parse(open(args.seqs_after),  'fasta'))

# Function to turn 1 into 1st, 2 into 2nd, etc.
ordinal = lambda n: "%d%s" % (n,"tsnrhtdd"[(n/10%10!=1)*(n%10<4)*n%10::4])

# Check same number of seqs
assert len(seqs_after) == len(seqs_before), \
'The two input files should have the same number of sequences'

# Iterate through the 1st seq in each file, then the 2nd etc.
for seq_num in range(len(seqs_before)):
  seq_before = str(seqs_before[seq_num].seq)
  seq_after = str(seqs_after[seq_num].seq)

  # Strip - and ? chars, then check they're identical
  seq_before_gapless = seq_before.replace("-", "").replace("?", "")
  seq_after_gapless  = seq_after.replace("-", "").replace("?", "")
  assert seq_before_gapless.upper() == seq_after_gapless.upper(), \
  "After stripping '-' and '?', and ignoring differences in case, the " \
  + str(ordinal(seq_num + 1)) + " seqs in the two input files differ"

  # Iterate through bases of the after seq. Use "-" or "?" where they occur;
  # otherwise use the base of the before seq.
  before_pos_gapless = 0
  seq_after_mod = ''
  for after_pos, after_base in enumerate(seq_after):
    if after_base == "-" or after_base == "?":
      seq_after_mod += after_base
    else:
      seq_after_mod += seq_before_gapless[before_pos_gapless]
      before_pos_gapless += 1

  assert seq_after_mod.upper() == seq_after.upper(), "Internal malfunction of "\
  + __file__ + ": modifying the seq changed more than its case. Please report "\
  "to Chris Wymant."

  seqs_after[seq_num].seq = Seq.Seq(seq_after_mod)

  if args.first_seq_only:
    break

SeqIO.write(seqs_after, sys.stdout, "fasta")

