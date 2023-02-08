#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script removes reads with an 'identity' (the
fraction of bases which are mapped and agree with the reference) below a
specified value from a bam file.
'''

import os
import sys
import argparse
import pysam
from Bio import SeqIO
from ShiverFuncs import CalculateReadIdentity

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('RefFile', type=File)
parser.add_argument('InBamFile', type=File)
parser.add_argument('OutBamFile')
parser.add_argument('ReadIdentityThreshold', type=float)
parser.add_argument('-P', '--keep-pairs-only', action='store_true')
args = parser.parse_args()

# Check the identity threshold is between 0 and 1
if not 0 < args.ReadIdentityThreshold <= 1:
  print('The read identity threshold should be greater than zero, and less',\
  'than or equal to 1. Quitting.', file=sys.stderr)
  exit(1)

# Get the reference.
SeqList = list(SeqIO.parse(open(args.RefFile), 'fasta'))
if len(SeqList) != 1:
  print('There are', len(SeqList), 'sequences in', args.ref_file +\
  '. There should be exactly 1. Quitting.', file=sys.stderr)
  exit(1)
RefSeq = str(SeqList[0].seq)

InBam = pysam.AlignmentFile(args.InBamFile, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = InBam.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', args.InBamFile+'; found',\
  str(len(AllReferences))+'.Quitting.', file=sys.stderr)
  exit(1)
RefName = AllReferences[0]

OutBam = pysam.AlignmentFile(args.OutBamFile, "wb", template=InBam)

# Iterate through the reads
UnpairedReads = {}
for read in InBam.fetch(RefName):

  # Calculate the read's identity
  identity = CalculateReadIdentity(read, RefSeq)

  if identity >= args.ReadIdentityThreshold:

    # If we've seen the mate before, print this read and its mate; otherwise
    # record this read.
    if args.keep_pairs_only:
      if read.query_name in UnpairedReads:
        OutBam.write(UnpairedReads[read.query_name])
        OutBam.write(read)
        del UnpairedReads[read.query_name]
      else:
        UnpairedReads[read.query_name] = read

    else:
      OutBam.write(read)

