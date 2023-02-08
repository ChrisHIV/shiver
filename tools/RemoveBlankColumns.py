#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script removes those columns/positions in alignment
that are wholly 'blank', i.e. consist solely of the gap character "-". Output is
printed to stdout.'''

import argparse
import os
import sys
from Bio import AlignIO
from ShiverFuncs import RemoveBlankColumns

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('-?', '--q-mark', action='store_true', help='''Also consider
the '?' character (assumed throughout Chris Wymant's code to mean missing
coverage) to be blank.''')
parser.add_argument('-U', '--uninformative', action='store_true', help='''Also
remove those columns in which all non-blank bases are all identical.''')
args = parser.parse_args()

# Are we checking just for "-", or also for "?"
BlankChars = '-'
if args.q_mark:
  BlankChars += '?'

alignment = AlignIO.read(args.FastaFile, "fasta")

alignment = RemoveBlankColumns(alignment, BlankChars, args.uninformative)

AlignIO.write(alignment, sys.stdout, 'fasta')
