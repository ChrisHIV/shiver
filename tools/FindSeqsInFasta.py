#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script retrieves searched-for sequences from a
fasta file. Output is printed to stdout in fasta format, with options to invert
the search, extract a window of alignment, and strip gaps (call with --help for
details). Use EITHER the --seq-names option to specify the names of the
sequences you're looking for at the command line, OR the --seq-name-file option
to specify a file that contains the names of the sequences you're looking for.
'''

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
# Define a function to convert from a comma-separated pair of positive integers
# as a string, to a list of two integers, as a type for the argparse.
def CoordPair(MyCoordPair):
  if MyCoordPair.count(',') != 1:
    raise argparse.ArgumentTypeError(MyCoordPair+\
    ' does not contain exactly 1 comma.')
  LeftCoord, RightCoord = MyCoordPair.split(',')
  try:
    LeftCoord, RightCoord = int(LeftCoord), int(RightCoord)
  except ValueError:
    raise argparse.ArgumentTypeError('Unable to understand the values in'+\
    MyCoordPair+' as integers.')
  if LeftCoord > RightCoord:
    raise argparse.ArgumentTypeError('The left value should not be greater '+\
    'than the right value in '+MyCoordPair)
  if LeftCoord < 1:
    raise argparse.ArgumentTypeError('The left value should be greater than'+\
    ' or equal to 1 in '+MyCoordPair)
  return [LeftCoord,RightCoord]

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('-N', '--seq-names', nargs='+', help='''Used to specify the
names of the sequences you're looking for, separated by whitespace.''')
parser.add_argument('-F', '--seq-name-file', type=File, help="""A file in which
each line contains the name of sequence you're looking for.""")
parser.add_argument('-v', '--invert-search', action='store_true', \
help='return all sequences except those searched for')
parser.add_argument('-W', '--window', type=CoordPair,\
help='A comma-separated pair of positive integers specifying the coordinates'+\
' (with respect to the alignment if the sequenced are aligned) of a window '+\
'outside of which the desired sequences will be truncated / trimmed / not '+\
'printed. e.g. specifying 2,10 means only the second to tenth positions in '+\
'desired sequence(s) will be printed.')
parser.add_argument('-g', '--gap-strip', action='store_true', \
help='Remove all gap characters ("-" and "?") before printing. NB if '+\
'multiple sequences are being retrieved from an alignment, this will probably'+\
' unalign them. NB if used in conjunction with -W, fewer bases than the '+\
'window width may be printed.')
parser.add_argument('-S', '--match-start', action='store_true', \
help='Sequences whose names begin with one of the strings-to-be-searched for '+\
'are returned.')
parser.add_argument('-M', '--ignore-missing', action='store_true', help='''
By default, if any named sequences are not found we stop with an error (unless
--match-start is used). With this option, we ignore such missing sequences and
simply print those that were found. Warning: if none of the desired sequences
were found, the ouput will be blank.''')
parser.add_argument('-L', '--min-length', type=int, help='''Exclude from the
results any sequence whose length (after removing any "-" or "?" characters) is
less than this value. Note that if this length criterion is the only search
criterion you require, i.e. you don't want to search by sequence name, you can
set the compulsory SequenceName argument to the empty string "" and use the
--match-start option.''')
parser.add_argument('-D', '--allow-duplicates', action='store_true', help='''
Used to specify that there may be duplicate names in the input sequences; for
each named searched for, return all matches.''')

args = parser.parse_args()

# Check sequence names were specified properly
if args.seq_name_file and args.seq_names:
  print('You must use either the --seq-names option or the --seq-name-file',
  'option, not both. Exiting.', file=sys.stderr)
  exit(1)
if (not args.seq_name_file) and (not args.seq_names):
  print('You must use exactly one the --seq-names option or the --seq-name-file',
  'option. Exiting.', file=sys.stderr)
  exit(1)

if args.seq_name_file:
  SeqNames = []
  with open(args.seq_name_file, 'r') as f:
    for line in f:
      SeqNames += line.split()
  if not SeqNames:
    print('Nothing but whitespace found in', args.seq_name_file + '.Exiting.',
    file=sys.stderr)
    exit(1)
else:
  SeqNames = args.seq_names

# Check all sequences to be searched for are unique
CounterObject = collections.Counter(SeqNames)
DuplicatedArgs = [i for i in CounterObject if CounterObject[i]>1]
if len(DuplicatedArgs) != 0:
  for DuplicatedArg in DuplicatedArgs:
    print('Sequence name', DuplicatedArg, 'was duplicated in the arguments.',\
    file=sys.stderr)
  print('All sequence names should be unique. Exiting.', file=sys.stderr)
  exit(1)

NumSeqsToSearchFor = len(SeqNames)

# Find the seqs
AllSeqNamesEncountered = []
SeqsWeWant = []
SeqsWeWant_names = []
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
  AllSeqNamesEncountered.append(seq.id)
  if args.match_start:
    ThisSeqWasSearchedFor = False
    for beginning in SeqNames:
      if seq.id[0:len(beginning)] == beginning:
        ThisSeqWasSearchedFor = True
        break
  else:
    ThisSeqWasSearchedFor = seq.id in SeqNames
  if ThisSeqWasSearchedFor and (not args.invert_search):
    if seq.id in SeqsWeWant_names and not args.allow_duplicates:
      print('Sequence', seq.id, 'occurs multiple times in', args.FastaFile+\
      '\nQuitting.', file=sys.stderr)
      exit(1)
    SeqsWeWant.append(seq)
    SeqsWeWant_names.append(seq.id)
    if not args.match_start and not args.allow_duplicates and \
    len(SeqsWeWant) == NumSeqsToSearchFor:
      break
  elif args.invert_search and (not ThisSeqWasSearchedFor):
    if seq.id in SeqsWeWant_names and not args.allow_duplicates:
      print('Sequence', seq.id, 'occurs multiple times in', args.FastaFile+\
      '\nQuitting.', file=sys.stderr)
      exit(1)
    SeqsWeWant.append(seq)
    SeqsWeWant_names.append(seq.id)

# Check we found some sequences for printing!
if (not args.ignore_missing) and SeqsWeWant == []:
  ErrorMsg = 'Searched in ' + args.FastaFile + ' for ' + \
  ' '.join(SeqNames)
  if args.invert_search:
    ErrorMsg += ' with the --invert-search option'
  if args.match_start:
    ErrorMsg += ' with the --match-start option'
  ErrorMsg += '; found nothing.'
  print(ErrorMsg, file=sys.stderr)
  exit(1)

# Check all specified seqs were encountered (unless only the beginnings of names
# were specified).
if not (args.match_start or args.ignore_missing):
  SeqsNotFound = [seq for seq in SeqNames \
  if not seq in AllSeqNamesEncountered]
  if len(SeqsNotFound) != 0:
    print('The following sequences were not found in', args.FastaFile+':', \
    ' '.join(SeqsNotFound) +'\nQuitting.', file=sys.stderr)
    exit(1)

# Trim to the specified window and/or gap strip, if desired
for seq in SeqsWeWant:
  if args.window != None:
    LeftCoord, RightCoord = args.window
    if RightCoord > len(seq.seq):
      print('A window', LeftCoord, '-', RightCoord, 'was specified but', \
      seq.id, 'is only', len(seq.seq), 'bases long. Quitting.', file=sys.stderr)
      exit(1)
    seq.seq = seq.seq[LeftCoord-1:RightCoord]
  if args.gap_strip:
    seq.seq = seq.seq.ungap("-").ungap("?")

# Skip too-short sequences if desired
if args.min_length:
  NewSeqsWeWant = []
  for seq in SeqsWeWant:
    if len(seq.seq.ungap("-").ungap("?")) >= args.min_length:
      NewSeqsWeWant.append(seq)
  SeqsWeWant = NewSeqsWeWant

SeqIO.write(SeqsWeWant, sys.stdout, "fasta")

