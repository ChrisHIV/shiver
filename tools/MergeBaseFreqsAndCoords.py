#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script adds a different set of coordinates to a
file of base frequencies, in the format output by shiver (more specifically, by
the script shiver/tools/AnalysePileup.py which is run as part of shiver). Its
arguments are the base frequency file and either a file detailing the coordinate
transformation already performed by shiver (more specifically, by the script
shiver/tools/MergeAlignments.py), or a pairwise
alignment of the reference in the base frequency file (i.e. the reference to
which reads were mapped) with any other sequence. In the first case, the other
coordinate from shiver's transformation is added; in the second case, the
coordinate of the other sequence is added. Output is printed to stdout, suitable
for redirection to a csv file.'''

import sys
import os
import csv
import re
import argparse
from Bio import AlignIO

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BaseFreqFile', type=File)
parser.add_argument('--coords-file', '-C', type=File)
parser.add_argument('--pairwise-aln', '-P', type=File)
args = parser.parse_args()

# Check exactly one of --coords-file or --pairwise-aln has been used.
HavePairwiseAln = args.pairwise_aln != None
if HavePairwiseAln and args.coords_file != None:
  print('Specify either --coords-file or --pairwise-aln, not both. Quitting.',
  file=sys.stderr)
  exit(1)
if (not HavePairwiseAln) and args.coords_file == None:
  print('Specify exactly one of --coords-file or --pairwise-aln. Quitting.',
  file=sys.stderr)
  exit(1)

RefPosColumnInBaseFreqFile = 1
RefBaseColumnInBaseFreqFile = 2
AlnPosColumnInCoordsFile = 1
RefPosColumnInCoordsFile = 2


# Read in the BaseFreqFile.
BaseFreqData = []
with open(args.BaseFreqFile, 'r') as f:
  for line in f:
    BaseFreqData.append(line)
if len(BaseFreqData) < 2:
  print(args.coords_file, 'contains less than 2 lines; unexpected. Quitting.',
  file=sys.stderr)
  exit(1)
BaseFreqHeader = BaseFreqData[0]

RefPosToNewPosDict = {}

if HavePairwiseAln:

  # Find the reference name from the base freq file
  BaseFreqFirstColumnHeader = BaseFreqHeader.split(',',1)[0]
  RefRegex = '^position in (.*)'
  RefRegexCompiled = re.compile(RefRegex)
  match = RefRegexCompiled.match(BaseFreqFirstColumnHeader)
  if not match:
    print('The first field of the first line of', args.coords_file, "should",
    "begin with 'position in ', followed by the reference name. It doesn't.",
    "Quitting.", file=sys.stderr)
    exit(1)
  (RefName, ) = match.groups()

  # Read in the pairwise alignment.
  try:
    TwoSeqs = AlignIO.read(args.pairwise_aln, "fasta")
  except:
    print('Problem reading', args.pairwise_aln + ':',
    file=sys.stderr)
    raise
  if len(TwoSeqs) != 2:
    print('Found', len(TwoSeqs), 'sequences in', args.pairwise_aln + 
    "; there should be exactly 2. Quitting.", file=sys.stderr)
    exit(1)

  # Check one of the seqs has the right ref name; see which one.
  if TwoSeqs[0].id == RefName:
    if TwoSeqs[1].id == RefName:
      print('Both sequences in', args.pairwise_aln, 'are named', RefName, 
      "(the reference name mentioned in the base freq file); only one should"
      "be. Quitting.", file=sys.stderr)
      exit(1)
    RefSeq = TwoSeqs[0]
    OtherSeq = TwoSeqs[1]
  else:
    if TwoSeqs[1].id != RefName:
      print('Neither sequence in', args.pairwise_aln, 'is named', RefName, 
      "(the reference name mentioned in the base freq file); exactly one "
      "should be. Quitting.", file=sys.stderr)
      exit(1)
    RefSeq = TwoSeqs[1]
    OtherSeq = TwoSeqs[0]
  OutputString = 'Position in ' + OtherSeq.id + ',' + BaseFreqHeader
  RefSeq = str(RefSeq.seq)
  OtherSeq = str(OtherSeq.seq)
  GaplessRefSeq = RefSeq.replace('-','').upper()


  # These will be incremented to reflect the current position in each seq.
  # We start the OtherSeqPos at -1 because that's how we want to describe any
  # positions before the other seq starts.
  RefPos = 0
  OtherSeqPos = -1

  # Go through the alignment, and record the correspondence between coords in
  # the two seqs. We only need to record positions where neither has a gap.
  for AlnPosMin1 in range(TwoSeqs.get_alignment_length()):
    RefBase = RefSeq[AlnPosMin1]
    OtherSeqBase = OtherSeq[AlnPosMin1]
    if OtherSeqBase != '-':
      if OtherSeqPos == -1:
        OtherSeqPos += 2
      else:
        OtherSeqPos += 1
    if RefBase != '-':
      RefPos += 1
      if OtherSeqBase != '-':
        RefPosToNewPosDict[RefPos] = OtherSeqPos

else:

  OutputString = 'Position in alignment,' + BaseFreqHeader

  # Read in the coordinate file to get the correspondence between reference
  # coords and alignment coords.
  with open(args.coords_file, 'r') as f:
    reader = csv.reader(f, delimiter=',') #, quotechar='"')
    for LineNumberMin1, fields in enumerate(reader):
      if LineNumberMin1 == 0:
        continue
      try:
        RefPos = fields[RefPosColumnInCoordsFile-1]
        AlnPos = fields[AlnPosColumnInCoordsFile-1]
      except IndexError:
        print('Not enough fields on line', LineNumberMin1+1, 'in', \
        args.coords_file +'. Quitting.', file=sys.stderr)
        exit(1)
      if RefPos != '-':
        try:
          RefPos = int(RefPos)
        except IndexError:
          print('Failed to understand the reference position as an'
          ' integer on line', LineNumberMin1+1, 'in', \
          args.coords_file +'. Quitting.', file=sys.stderr)
          exit(1)
        RefPosToNewPosDict[RefPos] = AlnPos

# Iterate through the base freq data
LastRefPos = 0
for LineNumberMin2, line in enumerate(BaseFreqData[1:]):

  fields = line.split(',')

  # Get the reference position and base
  try:
    RefPos = fields[RefPosColumnInBaseFreqFile-1]
  except IndexError:
    print('Not enough fields on line', LineNumberMin2+2, 'in', \
    args.BaseFreqFile +'. Quitting.', file=sys.stderr)
    exit(1)
  try:
    RefBase = fields[RefBaseColumnInBaseFreqFile-1]
  except IndexError:
    print('Not enough fields on line', LineNumberMin2+2, 'in', \
    args.BaseFreqFile +'. Quitting.', file=sys.stderr)
    exit(1)

  # Skip positions where we have base freq data inside an insertion with respect
  # to the ref.
  if RefPos == 'NA':  
    continue

  # Check the ref position increments by 1.
  RefPos = int(RefPos)
  if RefPos != LastRefPos + 1:
    print('Error: the reference position jumped from', LastRefPos, 'to', RefPos,
    'on line', LineNumberMin2+2, 'in', args.BaseFreqFile + '. Quitting.',
    file=sys.stderr)
    quit(1)
  LastRefPos = RefPos

  # Check the reference base matches in the base freq file and in the alignment.
  if HavePairwiseAln and RefBase.upper() != GaplessRefSeq[RefPos-1]:
    print('Error: mismatch in reference base at position', str(RefPos) + \
    ':', RefBase, 'in', args.BaseFreqFile, 'but', GaplessRefSeq[RefPos-1],
    'in', args.pairwise_aln + '. Quitting.', file=sys.stderr)
    quit(1)

  try:  
    NewPos = RefPosToNewPosDict[RefPos]
  except KeyError:
    NewPos = '-'
  OutputString += str(NewPos) + ',' +line

print(OutputString.rstrip())

