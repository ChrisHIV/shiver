#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script estimates which of the possible bases (A, C,
G or T) an IUPAC ambiguity code in a sequence alignment represents, by using the
most common of the possibilities at that position in the alignment. Upper/lower
case is preserved. Gaps are not changed, nor are "N" bases.'''

import argparse
import os
import sys
from Bio import AlignIO
from Bio import SeqIO
import collections
from AuxiliaryFunctions import IUPACdict

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('OutputFile')
parser.add_argument('-V', '--verbose', action='store_true',
help="Print a line for each ambiguity code estimated.")
args = parser.parse_args()

try:
  alignment = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise
AlignmentLength = alignment.get_alignment_length()

# Bases we'll leave untouched
OKbases = "ACGTNacgtn-?"

BaseCountsByPos = []
BaseCountTotalsByPos = []
for pos in range(AlignmentLength):

  # Sort the OKbases here by how common they are. (Check there are some!)
  BaseCounts = collections.Counter(alignment[:, pos])
  BasesObserved = list(BaseCounts.keys())
  for base in BasesObserved:
    if not base in OKbases:
      del BaseCounts[base]
  if len(BaseCounts) == 0:
    print('No "OK" bases, i.e.', OKbases + ", were observed at position", pos \
    + 1, "in", args.alignment + ". This code does not know how to estimate",
    "ambiguity codes in this case. Quitting.", file=sys.stderr)
    exit(1)
  BaseCountsByPos.append(BaseCounts.most_common())
  BaseCountTotalsByPos.append(sum(BaseCounts.values()))

# Make the alignment writable
MutableSeqList = [seq.seq.tomutable() for seq in alignment]
IDs = [seq.id for seq in alignment]

for row in range(len(MutableSeqList)):
  SeqAsStr = str(MutableSeqList[row])
  ID = IDs[row]
  for pos in range(AlignmentLength):

    # Check at this position whether ambiguity interpretation is needed.
    base = SeqAsStr[pos]
    if not base in OKbases:

      # Convert to upper case, but remember whether it was orinally lower
      IsLower = base == base.lower()
      UpperBase = base.upper()

      # Which bases does this ambiguity code mean?
      if not UpperBase in IUPACdict:
        print("Error: unexpected base", base, 'at position', pos + 1, 'for seq',
        ID, 'in', args.alignment + ". Quitting.", file=sys.stderr)
        exit(1)
      bases = IUPACdict[UpperBase]

      # Iterate through the OKbases from most common to least common, take the
      # first match. If there's no match, crash - we don't know what to do.
      MatchToBaseElsewhere = False
      for BaseElsewhere, count in BaseCountsByPos[pos]:
        if BaseElsewhere in bases:
          BaseToUse = BaseElsewhere
          CountToUse = count 
          MatchToBaseElsewhere = True
          break
      if not MatchToBaseElsewhere:
        BaseToUse = "N"
        print('Warning: position', pos + 1, 'for seq', ID, 'in', args.alignment,
        "is", base + ", which is the ambiguity code for", " or ".join(bases) + \
        ", however none of these bases were observed in other sequences at",
        "this position. Replacing this base by N.", file=sys.stderr)

      # Set the new base
      if IsLower:
        BaseToUse = BaseToUse.lower()
      MutableSeqList[row][pos] = BaseToUse
      if args.verbose:
        print('At position', pos + 1, 'for seq', ID, "base", base,
        "was changed to", BaseToUse, "(the latter appearing", CountToUse,
        "times amongst", BaseCountTotalsByPos[pos], "unambiguous bases here).")

SeqIO.write((SeqIO.SeqRecord(seq, id=IDs[i], description='') for i, seq in \
enumerate(MutableSeqList)), args.OutputFile, 'fasta')
