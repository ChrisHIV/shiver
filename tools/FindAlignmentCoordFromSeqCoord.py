#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This alignment translates positions with respect to a
named sequence to positions with respect to an alignment containing that
sequence. e.g. if your named seq started like this in the alignment: --a-g...
and you specified position 2, this script will tell you position 5.'''

import argparse
import os
import sys
from Bio import AlignIO
from ShiverFuncs import TranslateSeqCoordsToAlnCoords

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('SeqName')
parser.add_argument('PositionsInSeq', type=int, nargs='+',
help='One or more positive integers.')
args = parser.parse_args()

# Check all positions positive
if min(args.PositionsInSeq) < 1:
  print('All positions should be greater than zero. Quitting.', file=sys.stderr)
  exit(1)

alignment = AlignIO.read(args.alignment, "fasta")
ChosenSeq = None
for seq in alignment:
  if seq.id == args.SeqName:
    ChosenSeq = str(seq.seq)
    break
if ChosenSeq == None:
  print('Did not find', args.SeqName, 'in', args.alignment + '. Quitting.', 
  file=sys.stderr)
  exit(1)

# Count missing coverage as a gap
ChosenSeq = ChosenSeq.replace('?', '-')

# Check no positions are after the end of the seq
if max(args.PositionsInSeq) > len(ChosenSeq) - ChosenSeq.count('-'):
  print('You specified position', max(args.PositionsInSeq),
  'but', args.SeqName, 'is only', len(ChosenSeq) - ChosenSeq.count('-'),
  'bases long. Quitting.', file=sys.stderr)
  exit(1)

TranslatedCoords = TranslateSeqCoordsToAlnCoords(ChosenSeq, args.PositionsInSeq)
print(' '.join(map(str, TranslatedCoords)))
