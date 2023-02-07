#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Aligns more sequences to a pairwise alignment in which
the first sequence contains missing coverage, i.e. the "?" character. Output is
printed to stdout suitable for redirection to a fasta-format file. The pairwise 
alignment is nominally the consensus (called from parsing mapped reads) and the
reference used for mapping, in that order. Alignment is performed using mafft;
mafft does not know what missing coverage is, hence the need for this program.
How it works: we replace missing coverage by gaps, realign, match the consensus
fragments after (which in general contain new gaps) with those before, then
replace the appropriate gaps by missing coverage.'''

import argparse
import os
import sys
import re
import copy
from Bio import SeqIO
from Bio import Seq
import itertools
import subprocess
import collections
from AuxiliaryFunctions import PropagateNoCoverageChar

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('OtherSeqsToBeAdded', type=File)
parser.add_argument('SeqPairWithMissingCov', type=File)
parser.add_argument('-F', '--addfragments', action='store_true', \
help='Call mafft with --addfragments instead of --add.')
parser.add_argument('-S', '--swap-seqs', action='store_true', \
help='Swap the consensus and the reference in the output (so that the order ' +\
'becomes consensus then reference then the other sequences added).')
parser.add_argument('-T1', '--temp-file-1', \
default='temp_AlnToMissingCovPair1.fasta')
parser.add_argument('-T2', '--temp-file-2', \
default='temp_AlnToMissingCovPair2.fasta')
parser.add_argument('--x-mafft', default='mafft', help=\
'''The command required to call mafft. You may include extra mafft options in
this command (except --add or --addfragments, since one of these must be used
and that is addressed with the -F argument to this script; and --preservecase is
always included here), which need to separated by white space as usual and then
the whole thing needs to be surrounded with one pair of quotation marks (so that
the mafft binary and its options are kept together as one option for this
script). If you include a path to your mafft binary (necessary if it is not in
the $PATH variable of your terminal), it may not include whitespace, since
whitespace is interpreted as separating raxml options. (The default value is
mafft).''')
args = parser.parse_args()

# Find the consensus and its ref
ConsensusFound = False
RefFound = False
for seq in SeqIO.parse(open(args.SeqPairWithMissingCov),'fasta'):
  if not ConsensusFound:
    consensus = seq
    ConsensusFound = True
    continue
  if not RefFound:
    ref = seq
    RefFound = True
    continue
  print('Found three sequences in', args.SeqPairWithMissingCov+\
  '; expected only two. Quitting.', file=sys.stderr)
  exit(1)
if not RefFound:
  print('Less than two sequences found in', args.SeqPairWithMissingCov+\
  '; expected two. Quitting.', file=sys.stderr)
  exit(1)

# Check consensus.id != ref.id
if consensus.id == ref.id:
  print("The consensus and the ref should have different names in",
  args.SeqPairWithMissingCov + ". Quitting.", file=sys.stderr)
  exit(1)

# Check the consensus and its ref are aligned. Skip positions with no bases.
ConsensusAsString = str(consensus.seq)
RefAsString = str(ref.seq)
if len(ConsensusAsString) != len(RefAsString):
  print(args.SeqPairWithMissingCov, 'is not an alignment - seq lengths', \
  'differ. Quitting.', file=sys.stderr)
  exit(1)
NewConsensusAsString = ''
NewRefAsString = ''
for ConsensusBase, RefBase in zip(ConsensusAsString, RefAsString):
  if RefBase == '-' and (ConsensusBase == '-' or ConsensusBase == '?'):
    continue
  else:
    NewConsensusAsString += ConsensusBase
    NewRefAsString += RefBase
ConsensusAsString = NewConsensusAsString
RefAsString = NewRefAsString
ref.seq = Seq.Seq(RefAsString)

# Replaces gaps that border "no coverage" by "no coverage".
ConsensusAsString = PropagateNoCoverageChar(ConsensusAsString)

# Check all seq IDs are unique.
IDsOfSeqsToBeAdded = [seq.id for seq in \
SeqIO.parse(open(args.OtherSeqsToBeAdded),'fasta')]
if consensus.id in IDsOfSeqsToBeAdded:
  print('A sequence in', args.OtherSeqsToBeAdded, 'is called', consensus.id +
  ', like the consensus in', args.SeqPairWithMissingCov + 
  '. Rename to avoid such a clash. Quitting.', file=sys.stderr)
  exit(1)

# Align
ConsensusNoMissingCov = copy.copy(consensus)
ConsensusNoMissingCovStr = ConsensusAsString.replace('?', '-')
ConsensusNoMissingCov.seq = Seq.Seq(ConsensusNoMissingCovStr)
if args.swap_seqs:
  seqs = [ref, ConsensusNoMissingCov]
else:
  seqs = [ConsensusNoMissingCov, ref]
SeqIO.write(seqs, args.temp_file_1, "fasta")
if args.addfragments:
  AddOption = '--addfragments'
else:
  AddOption = '--add'
with open(args.temp_file_2, 'w') as f:
  try:
    ExitStatus = subprocess.call(args.x_mafft.split() + ['--preservecase',
    AddOption, args.OtherSeqsToBeAdded, args.temp_file_1], stdout=f)
    assert ExitStatus == 0
  except:
    print('Problem calling mafft. Quitting.', file=sys.stderr)
    raise

# Read in the aligned seqs. Note which one is the consensus. Check all the 
# expected seqs are recovered.
AlignedSeqs = []
for i, seq in enumerate(SeqIO.parse(open(args.temp_file_2),'fasta')):
  AlignedSeqs.append(seq)
  if seq.id == consensus.id:
    ConsensusPosition = i
if sorted([seq.id for seq in AlignedSeqs]) != \
sorted([consensus.id, ref.id] + IDsOfSeqsToBeAdded):
  print('Error: different sequences found in', args.temp_file_2, \
  'compared to', args.SeqPairWithMissingCov, 'and', args.OtherSeqsToBeAdded + \
  '. Quitting.', file=sys.stderr)
  exit(1)

# Check the consensus only has changes in gaps.
PostAlignmentConsensus = AlignedSeqs[ConsensusPosition]
PostAlignmentConsensusAsStr = str(PostAlignmentConsensus.seq)
if PostAlignmentConsensusAsStr.replace('-','') != \
ConsensusNoMissingCovStr.replace('-',''):
  print('Error:', consensus.id, 'contains different bases before and after', \
  'alignment. Quitting.', file=sys.stderr)
  exit(1)

# To be used shortly
def CharToInt(char):
  if char == '?':
    return 0
  if char == '-':
    return 1
  return 2

# To be used shortly
def SplitSeqByGapsAndMissingCov(seq):
  '''Split up a sequence into runs of missing coverage, runs of gaps, and runs
  of bases.'''
  for i, base in enumerate(seq):
    BaseType = CharToInt(base)
    if i == 0:
      ListOfBitsOfSeq = [base]
      ListOfBitTypes = [BaseType]
      continue
    if BaseType == ListOfBitTypes[-1]:
      ListOfBitsOfSeq[-1] += base
    else:
      ListOfBitsOfSeq.append(base)
      ListOfBitTypes.append(BaseType)
  return ListOfBitsOfSeq, ListOfBitTypes

# Split up the consensus, pre- and post-alignment, into 'bits', namely runs of
# bases ('BitType' 2), runs of gaps ('BitType' 1) and runs of missing coverage
# ('BitType' 0).
PreAlnConsensusBits, PreAlnConsensusBitTypes = \
SplitSeqByGapsAndMissingCov(ConsensusAsString)
PostAlnConsensusBits, PostAlnConsensusBitTypes = \
SplitSeqByGapsAndMissingCov(PostAlignmentConsensusAsStr)

NewConsensus = ''
ProgThroughPostAln = 0
InsidePreAlnBaseRun = False
for ProgThroughPreAln, PreAlnConsensusBit in enumerate(PreAlnConsensusBits):
  PreAlnConsensusBitType  = PreAlnConsensusBitTypes[ProgThroughPreAln]
  PostAlnConsensusBit     = PostAlnConsensusBits[ProgThroughPostAln]
  PostAlnConsensusBitType = PostAlnConsensusBitTypes[ProgThroughPostAln]

  # A PreAln BitType 0 or 1 should become a 1 after alignment. We want the 
  # PostAln length of the bit (it could be longer due to accommodating a new 
  # insertion), but the PreAln BitType.
  if PreAlnConsensusBitType in [0,1]:
    if PostAlnConsensusBitType != 1:
      print('Error running', sys.argv[0] + ': gap or missing coverage became', \
      'something other than a gap after alignment. Please report to Chris', \
      'Wymant (google for current email address).', file=sys.stderr)
      exit(1)
    NewConsensus += PreAlnConsensusBit[0] * len(PostAlnConsensusBit)
    ProgThroughPostAln += 1
    continue

  # A PreAln BitType 2 can either stay the same, or be chopped into alternating
  # 1s and 2s i.e. bits of sequence separated by gaps. We want to check this is
  # the case, then just use the new form (chopped into gap-separated-bits
  # appropriately). The new form should start with a 2, not a 1, unless we're
  # currently dealing with the very first bit, i.e. pre-aln the consensus began
  # with sequence but post-aln it begins with some gaps. We want that number og
  # gaps, but replaced by '?' characters.
  NewForm = ''
  if PostAlnConsensusBitType == 1:
    if ProgThroughPreAln != 0:
      print('Error running', sys.argv[0] + ': a consensus sequence fragment', \
      "before alignment turns into something starting with a gap, but its not",\
      'the first sequence fragment. Please report to Chris', \
      'Wymant (google for current email address).', file=sys.stderr)
      exit(1)
    NewForm += '?' * len(PostAlnConsensusBit)
    ProgThroughPostAln += 1
  PreAlnConsensusBitLength = len(PreAlnConsensusBit)
  while len(NewForm.replace('-', '').replace('?', '')) < \
  PreAlnConsensusBitLength:
    NewForm += PostAlnConsensusBits[ProgThroughPostAln]
    ProgThroughPostAln += 1
  if NewForm.replace('-', '').replace('?', '') != PreAlnConsensusBit:
    print('Error running', sys.argv[0] + ': unable to match a split fragment', \
    'of consensus to the same fragment pre-alignment. Please report to Chris', \
    'Wymant (google for current email address).', file=sys.stderr)
    exit(1)
  NewConsensus += NewForm

# If the PreAlnConsensus ended with sequence but the PostAlnConsensus ends with
# gaps, we won't have met those gaps in the iteration above. Append as many '?'
# as there are such gaps to the consensus.
NumPostAlnConsensusBits = len(PostAlnConsensusBits)
if ProgThroughPostAln == NumPostAlnConsensusBits - 1 and \
PostAlnConsensusBitTypes[-1] == 1:
  NewConsensus += '?' * len(PostAlnConsensusBits[-1])
elif ProgThroughPostAln != NumPostAlnConsensusBits:
  print('Error running', sys.argv[0] + ': something went wrong matching up', \
  'the consensus before and after alignment. Please report to Chris', \
  'Wymant (google for current email address).', file=sys.stderr)
  exit(1)

# Output sanity checks
if NewConsensus.replace('?', '-') != PostAlignmentConsensusAsStr:
  print('Error running', sys.argv[0] + ': something went wrong replacing "-"', \
  'characters by of "?" characters. Please report to Chris', \
  'Wymant (google for current email address).', file=sys.stderr)
  exit(1)
if '?-' in NewConsensus or '-?' in NewConsensus:
  print('Error running', sys.argv[0] + ': found "-"', \
  'character next to "?" character in output. Please report to Chris', \
  'Wymant (google for current email address).', file=sys.stderr)
  exit(1)

AlignedSeqs[ConsensusPosition].seq = Seq.Seq(NewConsensus)

SeqIO.write(AlignedSeqs, sys.stdout, "fasta")

