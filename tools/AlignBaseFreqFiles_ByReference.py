#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script aligns two base frequency files
of the format produced by shiver, each one resulting from mapping to a different
reference, figuring out the correspondance of positions based on alignment of
the two references. Note that AlignBaseFreqFiles_ByConsensuses.py does the same
thing but using an alignment of the two consensuses called from the two base
frequency files, instead of the two mapping references. Note that the accuracy
of the output of either of these scripts relies on the accuracy of the input
alignment. If using this script with mapping references that have indels with
respect to each other, alignment ambiguity at the site of the indels can cause a
slight shift/stagger in the correspondence between positions close to the site,
which will lead to an overestimate of the number of positions where the base
frequencies disagree on what base is most common. Therefore if your two base
frequency files were derived by mapping the same reads to two different
references and your aim is to count sites where base frequencies differ, you
should use AlignBaseFreqFiles_ByConsensuses.py instead (because the consensuses
should be nearly identical, and so their alignment should give a more accurate
correspondence of positions). If your two base frequency files were derived by
mapping the reads from two different samples and part of the genome is missing
for one or both samples, this script is likely to be better, because an
alignment of the references instead of the consensuses does not suffer from the
ambiguity due to missing sequence. Output of this script is printed to stdout
suitable for redirection to a csv file.'''

import argparse
import os
import sys
from Bio import AlignIO
import itertools

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File, help='An alignment that contains '+\
'both references.')
parser.add_argument('ref1name')
parser.add_argument('ref2name')
parser.add_argument('BaseFreqsFile1', type=File)
parser.add_argument('BaseFreqsFile2', type=File)
parser.add_argument('-C', '--coverage-only', action='store_true', help='''At
each position in the alignment, instead of printing the base frequencies for
each reference, print their coverage (the total number of reads mapped there,
i.e. the sum of the base counts).''')
parser.add_argument('-CF', '--compare', action='store_true', help='''At each
position in the alignment
print a score in the range [0,1] to summarise how similar the two
sets of frequencies are to each other. That score is calculated as follows. The
count for "N" is ignored. For each base "A", "C", "G", "T", "-", we calculate
the absolute difference between the fraction of all reads each reference has
with that base; we then sum these values, divide by 2, and subtract the
result from 1. This means a score of 1 is obtained if and only if the two sets
of frequencies agree perfectly, a score of zero is obtained if and only if there
are no bases in common, and all other situations give an intermediate value. A
score of NA is given at positions where one or both references have no mapped
reads, or if one reference has an insertion with respect to the other.''')
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
parser.add_argument('--min-coverage', type=int, default=1,
help='''For use with --compare-snips-with-coverage: specify a minimum coverage
below which SNPs are not compared.''')

args = parser.parse_args()

# Check for unique ref names
if args.ref1name == args.ref2name:
  print('You cannot specify the same reference twice. Quitting.',
  file=sys.stderr)
  exit(1)

# No commas in reference names
if ',' in args.ref1name or ',' in args.ref2name:
  print('Reference names may not contain commas. Rename the sequences in',
  args.alignment, 'and try again. Quitting.', file=sys.stderr)
  exit(1)

# Find the consensus and its ref
alignment = AlignIO.read(args.alignment, "fasta")
ref1seq = None
ref2seq = None
for seq in alignment:
  if seq.id == args.ref1name:
    if ref1seq != None:
      print('Found', args.ref1name, 'twice in', args.alignment + '. Quitting.',
      file=sys.stderr)
      quit(1)
    ref1seq = str(seq.seq)
  if seq.id == args.ref2name:
    if ref2seq != None:
      print('Found', args.ref2name, 'twice in', args.alignment + '. Quitting.',
      file=sys.stderr)
      quit(1)
    ref2seq = str(seq.seq)

# If no end position was given, set it to be the end of the alignment.
if args.end_pos_in_aln == None:
  args.end_pos_in_aln = alignment.get_alignment_length() + 1

# Check the refs were found
if ref1seq == None:
  print(args.ref1name, 'not found in', args.alignment + \
  '. Quitting.', file=sys.stderr)
  exit(1)
if ref2seq == None:
  print(args.ref2name, 'not found in', args.alignment + \
  '. Quitting.', file=sys.stderr)
  exit(1)

def GetFreqs(BaseFreqsFile, RefSeq):
  '''Converts freqs with respect to the ref to freqs with respect to the
  alignment.'''

  GaplessRefSeq = RefSeq.replace('-','').upper()
  FreqsInRef = []
  with open(BaseFreqsFile, 'r') as f:
    LastRefPos = 0
    for LineNumMin1, line in enumerate(f):
      if LineNumMin1 == 0:
        continue
      fields = line.split(',')
      RefPos = fields[0]
      if RefPos == 'NA':  
        continue
      RefPos = int(RefPos)
      if RefPos != LastRefPos + 1:
        print('Error: the reference position jumped to', RefPos, 'on line',
        LineNumMin1+1, 'in', BaseFreqsFile + '. Quitting', file=sys.stderr)
        quit(1)
      RefBase = fields[1]
      if RefBase.upper() != GaplessRefSeq[RefPos-1]:
        print('Error: mismatch in reference base at position', str(RefPos) + \
        ':', RefBase, 'in', BaseFreqsFile, 'but', GaplessRefSeq[RefPos-1],
        'in', args.alignment + '. Quitting.', file=sys.stderr)
        quit(1)
      freqs = list(map(int, fields[2:]))
      assert len(freqs) == 6
      FreqsInRef.append(freqs)
      LastRefPos += 1
  RefLength = len(GaplessRefSeq)
  assert len(FreqsInRef) == RefLength
  FreqsInAlignment = []
  AlnPosToRefPos = []
  LastPosFreqs = [0,0,0,0,0,0]
  ZeroBasedPosInRef = -1
  for base in RefSeq:
    if ZeroBasedPosInRef == RefLength - 1:
      # We're after the end of the ref in the alignment
      freqs = [0,0,0,0,0,0]
    elif base == '-':
      # If the ref hasn't started yet, all base counts equal zero; if this is a
      # gap inside the ref, reproduce the base counts from the last position.
      # This is just to facilitate the coverage calculation, for which we take 
      # the coverage (sum of counts) here to be equal to its last value before
      # the deletion, but we won't report the breakdown into different bases.
      freqs = LastPosFreqs
    else:
      ZeroBasedPosInRef += 1
      freqs = FreqsInRef[ZeroBasedPosInRef]
    FreqsInAlignment.append(freqs)
    LastPosFreqs = freqs
    if base == '-':
      AlnPosToRefPos.append('-')
    else:
      AlnPosToRefPos.append(ZeroBasedPosInRef+1)
  return FreqsInAlignment, AlnPosToRefPos

ref1freqs, ref1PosConversions = GetFreqs(args.BaseFreqsFile1, ref1seq)
ref2freqs, ref2PosConversions = GetFreqs(args.BaseFreqsFile2, ref2seq)

# Set up the csv file header
outstring = 'Alignment position,Position in ' + args.ref1name + \
',Position in ' + args.ref2name
if args.coverage_only:
  outstring += ',coverage for ' + args.ref1name + ',coverage for ' + \
  args.ref2name
else:
  outstring += ',A count for ' + args.ref1name + \
  ',C count for ' + args.ref1name + ',G count for ' + args.ref1name + \
  ',T count for ' + args.ref1name + ',gap count for ' + args.ref1name + \
  ',N count for ' + args.ref1name + ',A count for ' + args.ref2name + \
  ',C count for ' + args.ref2name + ',G count for ' + args.ref2name + \
  ',T count for ' + args.ref2name + ',gap count for ' + args.ref2name + \
  ',N count for ' + args.ref2name
if args.compare_simple:
  outstring += ',agreement on the most common base?'
if args.compare:
  outstring += ',base frequency similarity score'



NumPosWithHigherCovIn1 = 0
NumPosWithHigherCovIn2 = 0

# Record each row of the csv file
for PosMin1, (ref1freqs, ref2freqs) in enumerate(zip(ref1freqs,
ref2freqs)):
  PosInRef1 = ref1PosConversions[PosMin1]
  PosInRef2 = ref2PosConversions[PosMin1]
  outstring += '\n' + str(PosMin1+1) + ',' + str(PosInRef1) + ',' + \
  str(PosInRef2) 

  if args.coverage_only:
    ref1cov = sum(ref1freqs)
    ref2cov = sum(ref2freqs)
    outstring += ',' + str(ref1cov) + ',' + str(ref2cov)
  else:
    if PosInRef1 == '-':
      ref1freqs = ['-'] * 6
    if PosInRef2 == '-':
      ref2freqs = ['-'] * 6
    outstring += ',' + ','.join(map(str,ref1freqs)) + ',' + \
    ','.join(map(str,ref2freqs))

  if args.compare or args.compare_simple or args.compare_snips_with_coverage:
    if PosInRef1 == '-' or PosInRef2 == '-':
      SimScoreBin = 'NA'
      SimScoreCont = 'NA'
    else:
      ref1freqs = ref1freqs[:5]
      ref2freqs = ref2freqs[:5]  
      ref1cov = sum(ref1freqs)
      ref2cov = sum(ref2freqs)
      if ref1cov == 0 or ref2cov == 0:
        SimScoreBin = 'NA'
        SimScoreCont = 'NA'
      else:
        if args.compare_simple or args.compare_snips_with_coverage:
          MaxFreq1 = max(ref1freqs)
          BasesWithMaxCount1 = [i for i,count in enumerate(ref1freqs) \
          if count == MaxFreq1]
          MaxFreq2 = max(ref2freqs)
          BasesWithMaxCount2 = [i for i,count in enumerate(ref2freqs) \
          if count == MaxFreq2]
          if BasesWithMaxCount1 == BasesWithMaxCount2:
            SimScoreBin = 1
          else:
            SimScoreBin = 0
          # Compare snps with coverage, if desired, in the specified window of
          # the alignment, when both references have enough coverage, and when
          # neither of them calls a gap as the most common base.
          if SimScoreBin == 0 and args.compare_snips_with_coverage and \
          PosMin1+1 >= args.start_pos_in_aln and \
          PosMin1+1 <= args.end_pos_in_aln and \
          ref1cov >= args.min_coverage and \
          ref2cov >= args.min_coverage and \
          BasesWithMaxCount1 != [4] and \
          BasesWithMaxCount2 != [4]:
            if ref1cov > ref2cov:
              NumPosWithHigherCovIn1 += 1
            else:
              NumPosWithHigherCovIn2 += 1
        if args.compare:
          SimScoreCont = 0
          for i in range(5):
            SimScoreCont += \
            abs(float(ref1freqs[i])/ref1cov - float(ref2freqs[i])/ref2cov)
          SimScoreCont = 1 - SimScoreCont/2
    if args.compare_simple:
      outstring += ',' + str(SimScoreBin)
    if args.compare:
      outstring += ',' + str(SimScoreCont)

# Print output
if args.compare_snips_with_coverage:  
  print(NumPosWithHigherCovIn1, NumPosWithHigherCovIn2)
else:
  print(outstring)
