#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script splits a fasta file up, printing each
sequence into a different file, named according to the sequence name (removing
any characters that are problematic for file names).'''

import argparse
import os
import sys
from Bio import SeqIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Define a function to check directories exist, as a type for the argparse.
def Dir(MyDir):
  if not os.path.isdir(MyDir):
    raise argparse.ArgumentTypeError(MyDir+\
    ' does not exist or is not a directory.')
  return MyDir

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('FastaFile', type=File)
parser.add_argument('OutputDir', type=Dir)
parser.add_argument('-G', '--gap-strip', action='store_true', \
help='Remove all gap characters ("-").')
args = parser.parse_args()

FilenameSafeChars = \
'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890_-.() '

def FilenameFromSeqID(FilenameSafeID):
  '''Join the output dir, the file name and its extension.'''
  return os.path.join(args.OutputDir, FilenameSafeID + '.fasta')

SeqDict = {}
for seq in SeqIO.parse(open(args.FastaFile),'fasta'):

  # Remove characters problematic for filenames from the seq id.
  # Check that something is left, that it's unique, and that the output file
  # does not already exist.
  FilenameSafeID = ''.join(char for char in seq.id if char in FilenameSafeChars)
  if FilenameSafeID == '':
    print('Sequence name', seq.id, 'in', args.FastaFile, 'containins no',
    'filename-safe characters. Quitting.', file=sys.stderr)
    exit(1)
  if FilenameSafeID in SeqDict:
    print('Sequence names', SeqDict[FilenameSafeID].id, 'and', seq.id, 'are',
    'identical after removing filename-unsafe characters. Please rename and', \
    'try again. Quitting.', file=sys.stderr)
    exit(1)
  if os.path.isfile(FilenameFromSeqID(FilenameSafeID)):
    print('The filename-safe version of', seq.id, 'is', FilenameSafeID + \
    ', but', FilenameFromSeqID(FilenameSafeID), 'already exists. We will not', \
    'overwrite it: please rename, move or delete it and try again. Quitting.', \
    file=sys.stderr)
    exit(1)

  # Gap strip if desired.
  if args.gap_strip:
    seq.seq = seq.seq.ungap("-")

  SeqDict[FilenameSafeID] = seq

# Check at least one seq was found.
if SeqDict == {}:
  print('No sequences found in', args.FastaFile+'. Quitting.', file=sys.stderr)
  exit(1)

# Write the files.
for FilenameSafeID, seq in list(SeqDict.items()):
  SeqIO.write(seq, FilenameFromSeqID(FilenameSafeID), "fasta")
