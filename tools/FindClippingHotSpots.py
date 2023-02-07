#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script finds, at each position in a bam file, the
number and fraction of reads that are clipped from that position to their left
or right end. Having many such reads is a warning sign that the reference and
reads are so different that reads were not aligned correctly.'''

import os
import collections
import sys
import argparse
import pysam

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BamFile', type=File)
parser.add_argument('MinClipLength', type=int, help='''Don't count clipped ends
whose length is less than this. e.g. if you specify 3, we don't count cases
where 1 or 2 bases have been clipped. If you specify 1 or less, we count all
cases of clipping. A value a little larger than 1 is likely to be optimal,
depending on the mapper used, since a single SNP close to the end of the read
may result in the end being clipped, which is not a problem provided other reads
span this position.''')
parser.add_argument('-N', '--min-read-count', type=int, default=1, help='''Don't
report positions where the number of reads clipped is less than this value. (The
default value of 1 means any position with clipping is reported.)''')
args = parser.parse_args()

BamFile = pysam.AlignmentFile(args.BamFile, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = BamFile.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', BamFileName+'; found',\
  str(len(AllReferences))+'.\nQuitting.', file=sys.stderr)
  exit(1)
RefName = AllReferences[0]

# Get the length of the reference.
AllReferenceLengths = BamFile.lengths
if len(AllReferenceLengths) != 1:
  print('Pysam error: found one reference but', len(AllReferenceLengths), \
  'reference lengths.\nQuitting.', file=sys.stderr)
  exit(1)
RefLength = AllReferenceLengths[0]

ClipPositions = []
NumbersOfSpanningReads = [0 for i in range(RefLength)]
NumReads = 0

for read in BamFile.fetch(RefName):

  positions = read.get_reference_positions(full_length=True)

  # Shouldn't happen, but skip reads mapped to no position
  if not any(pos != None for pos in positions):
    continue

  NumReads += 1

  # If the left edge is clipped, find where.
  LeftMostMappedBase = 0
  while positions[LeftMostMappedBase] == None:
    LeftMostMappedBase += 1
  if LeftMostMappedBase > args.MinClipLength - 1:
    ClipPositions.append(positions[LeftMostMappedBase])

  # If the right edge is clipped, find where.
  RightMostMappedBase = len(positions)-1
  while positions[RightMostMappedBase] == None:
    RightMostMappedBase -= 1
  if RightMostMappedBase < len(positions) - args.MinClipLength:
    ClipPositions.append(positions[RightMostMappedBase]+1)

  # Update the counts of reads spanning each position
  for pos in range(positions[LeftMostMappedBase]+1, \
  positions[RightMostMappedBase]+1):
    NumbersOfSpanningReads[pos] += 1

if NumReads == 0:
  print('No mapped reads found in', args.BamFile, '. Quitting.', \
  file=sys.stderr)
  exit(1)

# Count how many times each position was recorded
ClipPositionCounts = collections.Counter(ClipPositions)

if args.min_read_count > 1:
  ClipPositionCounts = {key:value for key, value in list(ClipPositionCounts.items()) \
  if value >= args.min_read_count}

# Print the output
output = 'Reference position, Number of reads clipped, Percentage of spanning'+\
' reads clipped'

for pos, count in sorted(list(ClipPositionCounts.items()), key=lambda x:x[1], \
reverse=True):
  # 100% of reads overhanging the start or end of the reference are clipped.
  if pos == 0 or pos == RefLength:
    PercentageClipped = 100
  else:
    PercentageClipped = 100 * float(count) / \
    (count + NumbersOfSpanningReads[pos])
  output += '\n' + str(pos+1) + ',' + str(count) + ',%.3f' % PercentageClipped

print(output)
