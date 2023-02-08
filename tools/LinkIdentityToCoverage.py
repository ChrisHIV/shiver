#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script calculates the mean read identity
(fractional agreement with the mapping reference) over all positions in a bam
file, after grouping those positions with the same coverage (number of mapped
reads at that position). When there is a tendency for contaminant reads to make
up a bigger proportion of the reads at low coverage positions, there should also
be a tendency for mean read identity to decrease at low coverage. The value of
the coverage at which this happens is therefore informative of a possible
coverage threshold for calling a consensus base. The user is advised to run this
script separately on each bam file for which the relationship between coverage
and proportion of contaminant reads is expected to be roughly the same (for
example, over a data set sequenced in the same way), and then merge the results
with ~/shiver/tools/LinkIdentityToCoverage_CombineBams.py.'''

import os
import sys
import argparse
import pysam
import subprocess
from Bio import SeqIO
from ShiverFuncs import CalculateReadIdentity

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('BamFile', type=File)
parser.add_argument('RefFile', type=File)
parser.add_argument('-S', '--start', type=int, help='Used to specify a ' +\
'position in the reference before which we will ignore reads.')
parser.add_argument('-E', '--end',   type=int, help='Used to specify a ' +\
'position in the reference after which we will ignore reads.')
parser.add_argument('--double-count-overlaps', action="store_true", help="""For
paired reads, by default when we encounter the second in the pair we ignore
mapped positions that were contained in the first (i.e. counting the overlap
only once, which is what samtools mpileup does unless you specify
--ignore-overlaps and/or --min-BQ 0). With this option, we will count the
overlap twice.""")
args = parser.parse_args()

SeqList = list(SeqIO.parse(open(args.RefFile), 'fasta'))
if len(SeqList) != 1:
  print('There are', len(SeqList), 'sequences in', args.ref_file +\
  '. There should be exactly 1. Quitting.', file=sys.stderr)
  exit(1)
RefSeq = str(SeqList[0].seq)
RefLength = len(RefSeq)

# Try to make the bam index file.
if not os.path.isfile(args.BamFile+'.bai'):
  try:
    ExitStatus = subprocess.call(['samtools', 'index', args.BamFile])
    assert ExitStatus == 0
  except:
    print(args.BamFile, 'has no corresponding .bai file, and there was a',
    'problem calling samtools index to create this file. It is needed for',
    'pysam. Quitting.', file=sys.stderr)
    raise

BamFile = pysam.AlignmentFile(args.BamFile, "rb")

# Find the reference in the bam file; there should only be one.
AllReferences = BamFile.references
if len(AllReferences) != 1:
  print('Expected exactly one reference in', args.BamFile+'; found',\
  str(len(AllReferences))+'.Quitting.', file=sys.stderr)
  exit(1)
RefName = AllReferences[0]

# Make the start and end positions zero-based if we have them. If we only have
# one (i.e. start or end), set the other manually, to facilitate the comparison
# of positions. Check end >= start.
start = args.start
end = args.end
HaveStart = start != None
HaveEnd = end != None
HaveWindow = HaveStart or HaveEnd
if HaveStart:
  start -= 1
  if not HaveEnd:
    end == RefLength
  elif end < start:
    print('The end point should not be before the start point. Quitting.',
    file=sys.stderr)
    exit(1)
if HaveEnd:
  end -= 1
  if not HaveStart:
    start = 0

# At each reference position record (1) the number of reads mapped here and (2)
# the sum of the identity values of reads mapped here. Each read contributes to
# a count of 1, and its identity value, to every position to which it is mapped.
CoveragesByPos = [0 for pos in range(RefLength)]
IdentityTotalsByPos = [0 for pos in range(RefLength)]
PositionsByRead = {}
for read in BamFile.fetch(RefName):

  if read.is_unmapped or read.is_supplementary:
    continue

  MappedPositions = read.get_reference_positions(full_length=False)

  # If this read is paired and we've seen its mate already, ignore mapped
  # positions in this read that were contained in the mate. If we haven't seen
  # the mate already, record this read's positions.
  if (not args.double_count_overlaps) and read.is_paired:
    if read.query_name in PositionsByRead:
      MatePositions = PositionsByRead[read.query_name]
      MappedPositions = [pos for pos in MappedPositions if not pos in \
      MatePositions]
      if not MappedPositions:
        continue
    else:
      PositionsByRead[read.query_name] = MappedPositions

  identity = CalculateReadIdentity(read, RefSeq)

  CheckPositionsInWindow = False
  if HaveWindow:
    # Skip reads wholly outside the window of interest
    if MappedPositions[0] > end or MappedPositions[-1] < start:
      continue
    # Reads partially in, partially out the window need to be checked
    if MappedPositions[0] < start or MappedPositions[-1] > end:
      CheckPositionsInWindow = True
  if CheckPositionsInWindow:
    for pos in MappedPositions:
      if pos >= start:
        if pos > end:
          break
        CoveragesByPos[pos] += 1
        IdentityTotalsByPos[pos] += identity
  else:  
    for pos in MappedPositions:
      CoveragesByPos[pos] += 1
      IdentityTotalsByPos[pos] += identity  



# Combine identity totals for positions that have the same coverage, and count
# how many positions have that coverage.
IdentityTotalsByCoverage = {}
CoverageCounts = {}
for pos in range(RefLength):
  coverage = CoveragesByPos[pos]
  if coverage > 0:
    IdentityTotal = IdentityTotalsByPos[pos]
    if coverage in IdentityTotalsByCoverage:
      IdentityTotalsByCoverage[coverage] += IdentityTotal
      CoverageCounts[coverage] += 1
    else:
      IdentityTotalsByCoverage[coverage] = IdentityTotal
      CoverageCounts[coverage] = 1

# The total number of reads present in all positions with a given coverage
# (counting each read once per position-that-has-that-coverage, not once in
# total) is the value of the coverage multiplied by the number of
# positions-that-have-that-coverage. Dividing the total read identity associated
# with this coverage (again, with each read contributing once per 
# position-that-has-that-coverage) by that total number of reads gives the mean
# identity for that coverage.
outstring = 'Coverage,Number of positions with that coverage,Mean identity'
for coverage, count in sorted(list(CoverageCounts.items()), key=lambda x: x[0]):
  TotalNumReads = coverage * count
  MeanIdentity = float(IdentityTotalsByCoverage[coverage]) / TotalNumReads
  outstring += '\n' + str(coverage) + ',' + str(count) + ',' + str(MeanIdentity)

print(outstring)
