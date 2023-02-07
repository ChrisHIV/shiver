#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script checks that all fasta files supplied as
arguments contain the same sequences, regardless of formatting. If any
difference is found we print "false" to stdout and exit with status 0; if no
difference is found we print "true" to stdout and exit with status 0. Any error
should result in different output printed to stdout/stderr and a non-zero exit
status.'''

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
parser.add_argument('FastaFile', type=File, nargs='+')
parser.add_argument('-C', '--case-sensitive', action='store_true', \
help='consider upper- and lower-case occurences of the same base to be '+\
'unequal, e.g. A != a. (By default, case is ignored.)')
parser.add_argument('-G', '--ignore-gaps', action='store_true', \
help='Ignore gaps when comparing sequences (by default, any difference in '+\
'gaps constitutes a difference).')
args = parser.parse_args()

# Check at least 2 arguments were supplied.
if len(args.FastaFile) < 2:
  print('At least two fasta files must be supplied as arguments. Quitting.', \
  file=sys.stderr)
  exit(1)

def FastaFileToSeqDict(FastaFile):
  '''Reads a fasta file into a dictionary.'''
  MyDict = {}
  for seq in SeqIO.parse(open(FastaFile),'fasta'):
    if seq.id in MyDict:
      print('Sequence', seq.id, 'occurs multiple times in', FastaFile+\
      '; sequence names ought to be unique. Quitting.', file=sys.stderr)
      exit(1)
    if args.ignore_gaps:
      seq.seq = seq.seq.ungap('-')
    SeqAsString = str(seq.seq)
    if not args.case_sensitive:
      SeqAsString = SeqAsString.upper()
    MyDict[seq.id] = SeqAsString
  return MyDict

# Check all files have the same sequences.
FirstFileDict = FastaFileToSeqDict(args.FastaFile[0])
for FastaFile in args.FastaFile[1:]:
  SeqDict = FastaFileToSeqDict(FastaFile)
  if SeqDict != FirstFileDict:
    print('false')
    exit(0)

print('true')
