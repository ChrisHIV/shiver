#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script analyses one or more fasta files,
printing the length of each sequence found therein.'''

import argparse
import os
import sys
import string
import re
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File, nargs="+")
parser.add_argument('-g', '--include-gaps', action='store_true', \
help='include gap characters, "-" and "?", ignored by default')
parser.add_argument('-C', '--ignore-lower-case', action='store_true', \
help="Doesn't count lower case letters.")
parser.add_argument('-1', '--first-seq-only', action='store_true', \
help='''Look at only the first sequence in each input fasta file.''')
parser.add_argument('-F', '--fragments', action='store_true', help='''Print the
length of the 'fragments' of the sequence, i.e. those parts that are separated
by one or more '?' characters (denoting missing sequence data). If a sequence
consists of N fragments, we will print N values followed by a zero (meaning that
the length of the (N+1)th fragment is zero).''')
parser.add_argument('--ignore-n', action='store_true', \
help='exclude the "N" and "n" characters, included by default')
parser.add_argument('-US', '--undetermined-start', action='store_true', help='''
Report the number of characters at the start of the sequence that are either
"N", "n" or "?". Except for --undetermined-end and --first-seq-only, other
options used in conjunction with this one will be ignored.''')
parser.add_argument('-UE', '--undetermined-end', action='store_true', help='''
Report the number of characters at the end of the sequence that are either
"N", "n" or "?". Except for --undetermined-start and --first-seq-only, other
options used in conjunction with this one will be ignored.''')
parser.add_argument('-LG', '--longest-gap', action='store_true', help='''Report
the length of the longest gap (the longest run of "-" characters). Except for
--first-seq-only, other options used in conjunction with this one will be
ignored.''')
args = parser.parse_args()

check_undetermined = args.undetermined_start or args.undetermined_end
check_both_undetermined = args.undetermined_start and args.undetermined_end

# Some options can't be used together.
if args.include_gaps and args.fragments:
  print('The --include-gaps and --fragments options cannot both be used at',
  'once. Quitting.', file=sys.stderr)
  exit(1)
if args.longest_gap and check_undetermined:
  print('The --longest-gap option cannot be used with the --undetermined-end',
  'or --undetermined-start options. Quitting.', file=sys.stderr)
  exit(1)

def get_max_match_length(string, pattern):
  "Of all matches of a pattern to a string, we report the largest length."
  match_lengths = [len(match) for match in re.findall(pattern, string)] 
  if len(match_lengths) == 0:
    return 0
  return max(match_lengths)


undetermined_start_regex = "^[nN\?]+"
undetermined_end_regex = "[nN\?]+$"

SeqLengths = []
for file_ in args.FastaFile:
  for seq in SeqIO.parse(open(file_),'fasta'):

    # Check for undetermined starts and ends if desired.
    seq_as_str = str(seq.seq)
    if check_undetermined:
      if args.undetermined_start:
        datum = [seq.id, get_max_match_length(seq_as_str,
        undetermined_start_regex)]
        if args.undetermined_end:
          datum.append(get_max_match_length(seq_as_str, undetermined_end_regex))
      else:
        datum = [seq.id, get_max_match_length(seq_as_str, undetermined_end_regex)]
      SeqLengths.append(datum)
      if args.first_seq_only:
        break
      continue

    # Check the longest gap if desired.
    if args.longest_gap:
      longest_gap_length = get_max_match_length(seq_as_str, "[-]+")
      SeqLengths.append([seq.id, longest_gap_length])
      if args.first_seq_only:
        break
      continue

    if not args.include_gaps:
      seq.seq = seq.seq.ungap("-")
      if not args.fragments:
        seq.seq = seq.seq.ungap("?")
    if args.ignore_n:
        seq.seq = seq.seq.ungap("n")
        seq.seq = seq.seq.ungap("N")
    if args.ignore_lower_case:
      seq.seq = ''.join(x for x in seq.seq if not x.islower())
    if args.fragments:
      FragLengths = [len(frag) for frag in seq.seq.split('?') if len(frag) > 0]
      FragLengths = sorted(FragLengths, reverse=True)
      FragLengths.append(0)
      SeqLengths.append([seq.id] + FragLengths)
    else:
      SeqLengths.append([seq.id, len(seq.seq)])
    if args.first_seq_only:
      break

for data in SeqLengths:
  print(' '.join(map(str,data)))

