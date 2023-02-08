#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script aligns two base frequency files
of the format produced by shiver, each one resulting from mapping to a different
reference, figuring out the correspondance of positions based on alignment of
the two consensuses called from the two base frequency files. Note that
AlignBaseFreqFiles_ByConsensuses.py does the same thing but using an alignment
of the two mapping references, instead of the two consensuses; read the --help
for that script for a discussion of the differences between these two scripts.
Here, the consensuses must have been called from the two base frequency files
using a custom run of shiver/tools/CallConsensus.py with the argument
MinFracToCall set to any negative value, and with the options
--use-n-for-missing --keep-gaps-by-missing. The value of the MinCoverage
argument you used there must be specified here. Output is printed to stdout
suitable for redirection to a csv file.'''

import argparse
import os
import sys
from Bio import AlignIO
import itertools
from AuxiliaryFunctions import CallAmbigBaseIfNeeded

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File, help='An alignment that contains '+\
'both consensuses.')
parser.add_argument('seq1name')
parser.add_argument('seq2name')
parser.add_argument('BaseFreqsFile1', type=File)
parser.add_argument('BaseFreqsFile2', type=File)
parser.add_argument('MinCoverage', type=int, help='''The value of the
minimum coverage threshold used to call a lower case base when deriving the
consensuses from the base frequency files.''')
parser.add_argument('-C', '--coverage-only', action='store_true', help='''At
each position in the alignment, instead of printing the base frequencies for
each consensus, print their coverage (the total number of reads mapped there,
i.e. the sum of the base counts).''')
parser.add_argument('-CF', '--compare', action='store_true', help='''At each
position in the alignment
print a score in the range [0,1] to summarise how similar the two
sets of frequencies are to each other. That score is calculated as follows. The
count for "N" is ignored. For each base "A", "C", "G", "T", "-", we calculate
the absolute difference between the fraction of all reads each consensus has
with that base; we then sum these values, divide by 2, and subtract the
result from 1. This means a score of 1 is obtained if and only if the two sets
of frequencies agree perfectly, a score of zero is obtained if and only if there
are no bases in common, and all other situations give an intermediate value. A
score of NA is given at positions where one or both consensuses have no mapped
reads, or if one consensus has an insertion with respect to the other.''')
parser.add_argument('-CS', '--compare-simple', action='store_true', help='''Like
--compare, except that the following simple similarity metric is used: 0 for
disagreement on which of the five bases is most common, 1 for agreement.''')
parser.add_argument('--compare-snips-with-coverage', action='store_true',
help='''Amongst positions where the two sets of base frequenices disagree on
what base is most common (ignoring positions where either set thinks a gap is
most common), we count at how many positions BaseFreqsFile1 has higher coverage
and at how many positions BaseFreqsFile2 has equal or higher coverage. These two
values are printed to stdout in that order.''')
parser.add_argument('--start-pos-in-aln', type=int, default=1,
help='''For use with --compare-snips-with-coverage: specify a start
position in the alignment before which snips are not compared. (Default is
1.)''')
parser.add_argument('--end-pos-in-aln', type=int,
help='''For use with --compare-snips-with-coverage: specify an end
position in the alignment after which snips are not compared. (Default is
the end of the alignment.)''')

args = parser.parse_args()

# Check for unique seq names
if args.seq1name == args.seq2name:
  print('The two consensus names must be distinct. Quitting.', file=sys.stderr)
  exit(1)

# No commas in consensus names
if ',' in args.seq1name or ',' in args.seq2name:
  print('Consensus names may not contain commas. Rename the sequences in',
  args.alignment, 'and try again. Quitting.', file=sys.stderr)
  exit(1)

# Consensuses should not have "?" in them.
def CheckQmarks(seq, name):
  if "?" in seq:
    print('Question mark character in', name + \
    '; unexpected. To create the alignment of consensuses you should have',
    'aligned them in a way that preserves information on how many bases are',
    'in a stretch of missing coverage, e.g. replacing "?" by "N" and then',
    'aligning. Quitting.',
    file=sys.stderr)
    quit(1)

# Find the consensuses
alignment = AlignIO.read(args.alignment, "fasta")
seq1 = None
seq2 = None
for seq in alignment:
  if seq.id == args.seq1name:
    if seq1 != None:
      print('Found', args.seq1name, 'twice in', args.alignment + '. Quitting.',
      file=sys.stderr)
      quit(1)
    seq1 = str(seq.seq).upper()
    CheckQmarks(seq1, args.seq1name)
  if seq.id == args.seq2name:
    if seq2 != None:
      print('Found', args.seq2name, 'twice in', args.alignment + '. Quitting.',
      file=sys.stderr)
      quit(1)
    seq2 = str(seq.seq).upper()
    CheckQmarks(seq2, args.seq2name)

# If no end position was given, set it to be the end of the alignment.
if args.end_pos_in_aln == None:
  args.end_pos_in_aln = alignment.get_alignment_length() + 1

# Check the seqs were found
if seq1 == None:
  print(args.seq1name, 'not found in', args.alignment + \
  '. Quitting.', file=sys.stderr)
  exit(1)
if seq2 == None:
  print(args.seq2name, 'not found in', args.alignment + \
  '. Quitting.', file=sys.stderr)
  exit(1)

ExpectedBasesNoN = ['A', 'C', 'G', 'T', '-']

def GetFreqs(BaseFreqsFile, consensus):
  '''Converts freqs with respect to the consensus seq to freqs with respect to
  the alignment.'''

  GaplessConsensusUpper = consensus.replace('-','').upper()
  LenConsensus = len(GaplessConsensusUpper)

  # For each line of the base freqs file get the base freqs:
  FreqsInConsensus = []
  PosInConsensus = -1
  with open(BaseFreqsFile, 'r') as f:
    for LineNumMin1, line in enumerate(f):
      if LineNumMin1 == 0:
        continue
      fields = line.split(',')
      try:
        freqs = list(map(int, fields[2:]))
        assert len(freqs) == 6
      except (ValueError, AssertionError):
        print("Unexpected input format of", BaseFreqsFile + ". It appears that",
        "this is not a base frequency file of the format produced by shiver.",
        "Quitting.", file=sys.stderr)
        raise

      RefBase = fields[1]

      # Ignore the count for 'N'
      FreqsNoN = freqs[:-1]

      # Call the appropriate base from these freqs.
      coverage = sum(FreqsNoN)
      if coverage < args.MinCoverage:
        if RefBase == "-":
          continue
        BaseFromFreqs = 'N'
      else:
        MaxFreq = max(FreqsNoN)
        BasesWithMaxFreq = [ExpectedBasesNoN[i] for i, freq in \
        enumerate(FreqsNoN) if freq == MaxFreq]
        BaseFromFreqs = CallAmbigBaseIfNeeded(BasesWithMaxFreq, coverage,
        args.MinCoverage, BaseFreqsFile)
        # Skip to the next position if the consensus has a gap here.
        if BaseFromFreqs == "-":
          continue

      # Check we haven't got more base frequencies than positions in the
      # consensuses
      PosInConsensus += 1
      if PosInConsensus == LenConsensus:
        print('Line', LineNumMin1 + 1, 'in', BaseFreqsFile + ', which is\n' + \
        line + 'should correspond to position', PosInConsensus + 1, 'in the',
        'consensus, but the associated consensus in', args.alignment, 'is only',
        LenConsensus, 'bases long. Quitting.', file=sys.stderr)
        exit(1)

      # Check the base here agrees with the one in the consensus
      if BaseFromFreqs != GaplessConsensusUpper[PosInConsensus]:
        print('Line', LineNumMin1 + 1, 'in', BaseFreqsFile + ', which is\n' + \
        line + 'should correspond to position', PosInConsensus + 1, 'in the',
        'consensus, but these base frequencies correspond to a base call of "' \
        + BaseFromFreqs + '" whereas position', PosInConsensus + 1, 'in the',
        'associated consensus in', args.alignment, 'is "' + \
        GaplessConsensusUpper[PosInConsensus] + '". Quitting.', file=sys.stderr)
        exit(1)

      FreqsInConsensus.append(freqs)

  # Sanity check on the consensus freqs
  if len(FreqsInConsensus) != LenConsensus:
    print("Error: the consensus length inferred from", BaseFreqsFile, "is",
    len(FreqsInConsensus), "but the length of the associated consensus",
    "sequence in", args.alignment, "is", str(LenConsensus) + ". Quitting.",
    file=sys.stderr)
    exit(1)

  FreqsInAlignment = []
  AlnPosToConsensusPos = []
  LastPosFreqs = [0,0,0,0,0,0]
  ZeroBasedPosInConsensus = -1
  for base in consensus:
    if ZeroBasedPosInConsensus == LenConsensus - 1:
      # We're after the end of the seq in the alignment
      freqs = [0,0,0,0,0,0]
    elif base == '-':
      # If the seq hasn't started yet, all base counts equal zero; if this is a
      # gap inside the seq, reproduce the base counts from the last position.
      # This is just to facilitate the coverage calculation, for which we take 
      # the coverage (sum of counts) here to be equal to its last value before
      # the deletion, but we won't report the breakdown into different bases.
      freqs = LastPosFreqs
    else:
      ZeroBasedPosInConsensus += 1
      freqs = FreqsInConsensus[ZeroBasedPosInConsensus]
    FreqsInAlignment.append(freqs)
    LastPosFreqs = freqs
    if base == '-':
      AlnPosToConsensusPos.append('-')
    else:
      AlnPosToConsensusPos.append(ZeroBasedPosInConsensus + 1)
  return FreqsInAlignment, AlnPosToConsensusPos, GaplessConsensusUpper

AllSeq2freqs, seq2PosConversions, GaplessSeq2upper = \
GetFreqs(args.BaseFreqsFile2, seq2)
AllSeq1freqs, seq1PosConversions, GaplessSeq1upper = \
GetFreqs(args.BaseFreqsFile1, seq1)

# Set up the csv file header
outstring = 'Alignment position,Position in ' + args.seq1name + \
',Position in ' + args.seq2name
if args.coverage_only:
  outstring += ',coverage for ' + args.seq1name + ',coverage for ' + \
  args.seq2name
else:
  outstring += ',A count for ' + args.seq1name + \
  ',C count for ' + args.seq1name + ',G count for ' + args.seq1name + \
  ',T count for ' + args.seq1name + ',gap count for ' + args.seq1name + \
  ',N count for ' + args.seq1name + ',A count for ' + args.seq2name + \
  ',C count for ' + args.seq2name + ',G count for ' + args.seq2name + \
  ',T count for ' + args.seq2name + ',gap count for ' + args.seq2name + \
  ',N count for ' + args.seq2name
if args.compare_simple:
  outstring += ',agreement on the most common base?'
if args.compare:
  outstring += ',base frequency similarity score'



NumPosWithHigherCovIn1 = 0
NumPosWithHigherCovIn2 = 0

# Record each row of the csv file
for PosMin1, (seq1freqs, seq2freqs) in enumerate(zip(AllSeq1freqs,
AllSeq2freqs)):
  PosInSeq1 = seq1PosConversions[PosMin1]
  PosInSeq2 = seq2PosConversions[PosMin1]
  outstring += '\n' + str(PosMin1+1) + ',' + str(PosInSeq1) + ',' + \
  str(PosInSeq2) 

  seq1cov = sum(seq1freqs[:5])
  seq2cov = sum(seq2freqs[:5])

  if args.coverage_only:
    outstring += ',' + str(seq1cov) + ',' + str(seq2cov)
  else:
    if PosInSeq1 == '-':
      seq1freqs = ['-'] * 6
    if PosInSeq2 == '-':
      seq2freqs = ['-'] * 6
    outstring += ',' + ','.join(map(str,seq1freqs)) + ',' + \
    ','.join(map(str,seq2freqs))

  if args.compare or args.compare_simple or args.compare_snips_with_coverage:
    if PosInSeq1 == '-' or PosInSeq2 == '-':
      SimScoreBin = 'NA'
      SimScoreCont = 'NA'
    else:
      seq1freqs = seq1freqs[:5]
      seq2freqs = seq2freqs[:5]  
      if seq1cov == 0 or seq2cov == 0:
        SimScoreBin = 'NA'
        SimScoreCont = 'NA'
      else:
        if args.compare_simple or args.compare_snips_with_coverage:
          BaseInSeq1 = seq1[PosMin1]
          BaseInSeq2 = seq2[PosMin1]
          if BaseInSeq1 == BaseInSeq2:
            SimScoreBin = 1
          else:
            SimScoreBin = 0
          # Compare snps with coverage, if desired, if both consensuses have
          # enough coverage:
          if SimScoreBin == 0 and args.compare_snips_with_coverage and \
          PosMin1+1 >= args.start_pos_in_aln and \
          PosMin1+1 <= args.end_pos_in_aln and \
          BaseInSeq1 != "N" != BaseInSeq2:
            #seq1cov >= args.MinCoverage and \
            #seq2cov >= args.MinCoverage:
            #print("disagreement at aln pos", PosMin1+1, BaseInSeq1, BaseInSeq2)
            if seq1cov > seq2cov:
              NumPosWithHigherCovIn1 += 1
            else:
              NumPosWithHigherCovIn2 += 1
        if args.compare:
          SimScoreCont = 0
          for i in range(5):
            SimScoreCont += \
            abs(float(seq1freqs[i])/seq1cov - float(seq2freqs[i])/seq2cov)
          SimScoreCont = 1 - SimScoreCont/2
    if args.compare_simple:
      outstring += ',' + str(SimScoreBin)
    if args.compare:
      outstring += ',' + str(SimScoreCont)

# Print output
if args.compare_snips_with_coverage:
  simple_total_diffs = 0
  for PosMin1, (base1, base2) in enumerate(zip(seq1, seq2)):
    if base1 != base2 and \
    PosMin1+1 >= args.start_pos_in_aln and \
    PosMin1+1 <= args.end_pos_in_aln and \
    base1 != "N" != base2 and \
    base1 != "-" != base2:
      simple_total_diffs += 1
  if simple_total_diffs != NumPosWithHigherCovIn1 + NumPosWithHigherCovIn2:
    print("Internal malfunction of", __file__ + ": found", simple_total_diffs,
    "positions disagreeing",
    "between", args.seq1name, "and", args.seq2name, "in", args.alignment,
    "but found", NumPosWithHigherCovIn1 + NumPosWithHigherCovIn2,
    "positions disagreeing from our site by site comparison of base",
    "frequencies. Quitting.", file=sys.stderr)
    exit(1)
    #print("WARNING: found", simple_total_diffs, "positions disagreeing",
    #"between", args.seq1name, "and", args.seq2name, "in", args.alignment,
    #"in a simple comparison excluding any position where either had an 'N'",
    #"(or a gap); but found", NumPosWithHigherCovIn1 + NumPosWithHigherCovIn2,
    #"positions disagreeing from our site by site comparison of base",
    #"frequencies, which excludes those 'N' positions due to insufficient",
    #"coverage but not due to an exact tie of a gap with another base.",
    #"Another possible reason for this discrepancy is an error in this code.",
    #"Investigate!", file=sys.stderr)
  print(NumPosWithHigherCovIn1, NumPosWithHigherCovIn2)
else:
  print(outstring)
