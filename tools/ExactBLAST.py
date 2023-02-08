#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Search for exact matches of query sequences in a set of
sequences, and report the locations of the matches.'''

import argparse
import os
import sys
from Bio import SeqIO
import collections

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File, help="Containing the seqs in which"\
" we search for hits.")
parser.add_argument('-Q', '--query-seq', nargs='+',
help='One or more seqs to search for.')
parser.add_argument('-F', '--query-seq-file', type=File)
parser.add_argument('-RC', '--rev-comp', action='store_true',
help="Search for reverse complements of queries (not done by default).")
parser.add_argument('-G', '--gap-aware', action='store_true', help='''Search for
queries allowing for gaps inside the sequences. (By default we strip sequences
of their gaps before searching, discarding alignment information.) This is
slow.''')
args = parser.parse_args()

HaveSeqFile = args.query_seq_file != None
CommandLineSeqs = args.query_seq != None

if not CommandLineSeqs and not HaveSeqFile:
  print("You must specify some query seqs. Quitting.", file=sys.stderr)
  exit(1)
if CommandLineSeqs and HaveSeqFile:
  print("You must use either --query-seq or --query-seq-file, not both.",
  "Quitting.", file=sys.stderr)
  exit(1)

if CommandLineSeqs:
  QuerySeqs = args.query_seq
else:
  QuerySeqs = []
  for seq in SeqIO.parse(open(args.query_seq_file), 'fasta'):
    QuerySeqs.append(str(seq.seq))
    if args.rev_comp:
      QuerySeqs.append(str(seq.seq.reverse_complement()))


# Make a list of upper-case seqs and their ids. Strip gaps if desired.
if args.gap_aware:
  seqs = [(seq.id, str(seq.seq).upper()) for seq in \
  SeqIO.parse(open(args.FastaFile), 'fasta')]
else:
  seqs = [(seq.id, str(seq.seq).upper().replace("-", "")) for seq in \
  SeqIO.parse(open(args.FastaFile), 'fasta')]

# Thanks stack overflow:
def find_all(a_str, sub):
  '''Find all locations of a substring inside a string.'''
  start = 0
  while True:
    start = a_str.find(sub, start)
    if start == -1: return
    yield start
    start += len(sub)

for QuerySeq in QuerySeqs:

  hits = []
  QueryLen = len(QuerySeq)

  for ID, seq in seqs:
    SeqLen = len(seq)

    if args.gap_aware:
      for PositionMin1, base in enumerate(seq):

        # Skip gaps
        if base == "-":
          continue

        # For each query, work forward through both the query and the seq
        # while the two agree, stopping as soon as they don't, ignoring
        # gaps in the seq.
        # When it equals the primer length, we have found that query.
        matches = 0
        StepsForward = 0
        while seq[PositionMin1 + StepsForward] == QuerySeq[matches] \
        or seq[PositionMin1 + StepsForward] == "-":
          if not seq[PositionMin1+StepsForward] == "-":
            matches += 1
            if matches == QueryLen:
              # Found the query!
              hits.append(PositionMin1)
              break
          StepsForward += 1
          if PositionMin1 + StepsForward == SeqLen:
            break

    else:
      hits += find_all(seq, QuerySeq)

  print(QuerySeq + ':', ' '.join(map(str, hits)))

