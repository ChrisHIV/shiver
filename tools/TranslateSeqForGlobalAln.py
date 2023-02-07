#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Called with a fasta file containing a consensus and its
reference for mapping (SEE IMPORTANT NOTE BELOW), and a coordinate file produced by fully processing a
sample with shiver, this script does the following: (1) excises insertions of
the consensus with respect to its reference, (2) excises insertions of the
reference (constructed out of the sample's contigs) with respect to the set of
existing references given as input to shiver, (3) adds gaps into the consensus
so that it is in the same coordinates system as the alignment of existing
references given as input to shiver, (4) prints the result to stdout, suitable
for redirection into a fasta file. Combining this output from multiple samples
into a single file gives an alignment.
IMPORTANT NOTE: this can only be used on shiver's first-round mapping result,
for which the mapping reference was constructed out of aligned contigs and real
references. It cannot be used with the second-round mapping result, for which
the mapping reference was the consensus from the first round of mapping. See
chapter 'The Global Alignment' in the shiver manual for more detail.'''
# NB, for brevity, 'pos' = position, 'aln' = alignment, 'seq' = sequence
GapChar = '-'

NoCoverageChar = '?'

import argparse
import os
import sys
import collections
from Bio import SeqIO
from Bio import Seq

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('ConsensusFile', type=File)
parser.add_argument('CoordsFile', type=File)
args = parser.parse_args()

# Exit if _remap is in the file name
if "_remap" in os.path.basename(args.ConsensusFile):
  print(__file__, 'was run on consensus file', args.ConsensusFile + \
  '; that file name contains the string "_remap". See the "IMPORTANT NOTE"',
  'part of the --help message for', __file__ + ". If", args.ConsensusFile,
  "really is the consensus produced by shiver's first-round mapping, rename it",
  "(just temporarily if desired) to bypass this saftey check. Quitting.",
  file=sys.stderr)
  exit(1)

# Read in the seq
SeqAsString = None
RefAsString = None
for InSeq in SeqIO.parse(open(args.ConsensusFile),'fasta'):
  if SeqAsString == None:
    seq = InSeq
    SeqAsString = str(InSeq.seq)
    continue
  RefAsString = str(InSeq.seq)
  break
if RefAsString == None:
  print('Did not find two sequences in', args.ConsensusFile + '. Quitting.', \
  file=sys.stderr)
  exit(1)
SeqLength = len(SeqAsString)
if len(RefAsString) != SeqLength:
  print(args.ConsensusFile, 'is not an alignment. Quitting.', file=sys.stderr)
  exit(1)

# Excise insertions of the seq relative to its ref
SeqNoInsertions = ''
for position,BaseInRef in enumerate(RefAsString):
  if BaseInRef != GapChar:
    SeqNoInsertions += SeqAsString[position]



SeqForGlobalAln = ''
PrevPosInAln = 0
PrevPosInRef = 0

# TODO: if first position encountered != 1, print a more helpful message

def CheckAndUpdatePos(pos, PrevPos, PosType):
  '''If pos is "-", returns PrevPos. Otherwise we check that
  int(pos) == PrevPos+1, and return this value.'''
  if pos == GapChar:
    return PrevPos
  try:
    pos = int(pos)
  except ValueError:
    # TODO: raise
    print('Error. Quitting')
    exit(1)
  if pos != PrevPos + 1:
    print('In', args.CoordsFile +', position', pos, 'in the', PosType, \
    'follows on from position', str(PrevPos) +'. Quitting.', file=sys.stderr)
    exit(1)
  if PosType == "reference" and pos > SeqLength:
    print('Encountered sequence position', pos, 'in', args.CoordsFile + \
    ', which is after the end of the sequence in', args.ConsensusFile, '('+ \
    str(SeqLength) +'bp long). Quitting.', file=sys.stderr)
    exit(1)
  return pos

# Iterate through the lines of the coords file
FirstLine = True
with open(args.CoordsFile, 'r') as f:
  for line in f:
    if FirstLine:
      FirstLine = False
      continue

    # Read in the fields
    fields = line.split(',')
    if len(fields) != 3:
      # TODO
      print('Error. Quitting')
      exit(1)
    PosInAln = fields[0]
    PosInRef = fields[1]

    # Check that, except for gaps, positions increment by 1.
    PrevPosInRef = CheckAndUpdatePos(PosInRef, PrevPosInRef, 'reference')
    PrevPosInAln = CheckAndUpdatePos(PosInAln, PrevPosInAln, 'alignment')

    # If this isn't a unique insertion in the sequence - which must be excised -
    # use the next base in the sequence.
    if PosInAln != GapChar:
      if PosInRef == GapChar:
        SeqForGlobalAln += GapChar
      else:
        SeqForGlobalAln += SeqNoInsertions[PrevPosInRef-1]


def PropagateNoCoverageChar(seq, LeftToRightDone=False):
  '''Replaces gaps that border "no coverage" by "no coverage".

  Where NoCoverageChars neighbour GapChars, propagate the former outwards until
  they touch bases on both sides (because insertions should only be called when
  the bases on either side are known). e.g.
  ACTG---?---ACTG
  becomes
  ACTG???????ACTG'''
  
  if LeftToRightDone:
    seq = seq[::-1]
  BaseToLeftIsNoCoverage = False
  ResultingSeq = ''
  for base in seq:
    if base == NoCoverageChar:
      BaseToLeftIsNoCoverage = True
      ResultingSeq += NoCoverageChar
    elif base == GapChar:
      if BaseToLeftIsNoCoverage:
        ResultingSeq += NoCoverageChar
      else:
        ResultingSeq += GapChar
    else:
      BaseToLeftIsNoCoverage = False
      ResultingSeq += base
  if LeftToRightDone:
    ResultingSeq = ResultingSeq[::-1]
  else:
    ResultingSeq = PropagateNoCoverageChar(ResultingSeq, True)
  return ResultingSeq

SeqForGlobalAln = PropagateNoCoverageChar(SeqForGlobalAln)

seq.seq = Seq.Seq(SeqForGlobalAln)
#SeqIO.write(seq, sys.stdout, "fasta")

FastaOut = SeqIO.FastaIO.FastaWriter(sys.stdout, wrap=50)
FastaOut.write_header()
FastaOut.write_record(seq)
FastaOut.write_footer()

