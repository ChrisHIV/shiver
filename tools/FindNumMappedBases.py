#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Prints the number of mapped bases for a bam file.
'''

import os
import collections
import sys
import argparse
import numpy
import pysam
from Bio import SeqIO
from ShiverFuncs import CalculateReadIdentity

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Define a comma-separated float pair object, as a type for the argparse.
def CommaSeparatedFloatPair(FloatPairAsString):
  try:
    values = FloatPairAsString.split(',')
    assert len(values) == 2
    values = [float(value) for value in values]
  except:
    raise argparse.ArgumentTypeError('Unable to understand ' +\
    FloatPairAsString + ' as a comma-separated pair of floats.')
  else:
    return values

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BamFile', type=File)
parser.add_argument('-I', '--identity-binning', type=CommaSeparatedFloatPair, \
help='''Breakdown the total count of mapped bases into contributions from reads
of different 'identity' (fractional agreement between a read and the reference
sequence). Use this option to specify two floats separated by a comma: the first
value being the minimum identity considered, the second being the identity step
width. e.g. specifying '0.5,0.05' means we'll report the total number of mapped
nucleotides from those reads with identity at least 0.5, then from those whose
identity is at least 0.55, etc. This option requires the reference sequence to
supplied using the --ref-file flag.''')
parser.add_argument('-R', '--ref-file', help='The file containing the '+\
'sequence of the reference (to which reads were mapped in the bam file).', \
type=File)
args = parser.parse_args()

BamFile = pysam.AlignmentFile(args.BamFile, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = BamFile.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', BamFileName+'; found',\
  str(len(AllReferences))+'.\nQuitting.', file=sys.stderr)
  exit(1)
RefName = AllReferences[0]

BinByIdentity = args.identity_binning != None

if BinByIdentity:

  if args.ref_file == None:
    print('The --identity-binning option requires the --ref-file option.', \
    'Quitting.', file=sys.stderr)
    exit(1)
    
  SeqList = list(SeqIO.parse(open(args.ref_file), 'fasta'))
  if len(SeqList) != 1:
    print('There are', len(SeqList), 'sequences in', args.ref_file +\
    '. There should be exactly 1. Quitting.', file=sys.stderr)
    exit(1)
  RefSeq = str(SeqList[0].seq)

  Min, BinWidth = args.identity_binning
  if not 0 <= Min < 1:
    print('The minimum for the --identity-binning option must be', \
    '0 <= Min < 1. Quitting.', file=sys.stderr)
    exit(1)
  Range = 1. - Min
  if not 0 < BinWidth < Range:
    print('The bin width for the --identity-binning option must be positive', \
    'and less than one minus the minimum (so that there is at least one bin).',\
    'Quitting.', file=sys.stderr)
    exit(1)
  NumBins = int(Range / BinWidth) + 1
  NumMappedBasesByReadIdentity = [0] * NumBins

  def Bin(x):
    '''returns 0 if x is in the first bin, 1 if it is in the second bin, ...'''
    if x >= 1:
      return NumBins-1
    if x <= Min:
      return 0
    return int(float(x - Min) / BinWidth)

else:
  NumMappedBases = 0

NumDone = 0
for read in BamFile.fetch(RefName):

  # Add the number of mapped bases to the total if that's all we're doing.
  if not BinByIdentity:
    NumMappedBases += read.query_alignment_length
    continue

  identity = CalculateReadIdentity(read, RefSeq)

  # Add the number of mapped bases to the bin for reads with this identity
  NumMappedBasesByReadIdentity[Bin(identity)] += \
  read.query_alignment_length

  
if not BinByIdentity:
  print(NumMappedBases)
  exit(0)

print('Minimum read identity, Number of mapped bases')
for i in range(NumBins):
  print(Min + i*BinWidth, sum(NumMappedBasesByReadIdentity[i:]), sep=',')


