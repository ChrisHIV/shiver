#!/usr/bin/env python


## Author: Tanya Golubchik and Chris Wymant chris.wymant@bdi.ox.ac.uk
## Acknowledgement: written with funding from the ERC Advanced Grant PBDR-339251 
##
## Overview
ExplanatoryMessage = '''This script takes as input a sequence alignment and the
names of some sequences in it to consider for correction. Correction consists of
splitting the sequence into two at each place where there is a large enough gap,
then applying a minimum length threshold for sequences (or parts of sequence) to
keep. The intended usage is for a set contigs aligned to a set of references;
the process of alignment may have introduced large gaps into the contigs that
should not be interpreted as genuine deletions.'''


import os.path
import sys
import argparse
import collections
import re
import itertools
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from ShiverFuncs import GetSeqStartAndEndPos, RemoveBlankColumns

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment', type=File)
parser.add_argument('NameOfSeqToCorrect', nargs='+')
parser.add_argument( '-G', '--split-gap-size', default=100, type=int, help='''
The gap size at which a contig will be split into two. (Default: 100.)''')
parser.add_argument( '-M', '--min-contig-size', default=80, type=int, help='''
The minimum length of a contig (or a piece of a contig, after contig splitting) 
for it to be kept in the alignment. (Default: 80.)''')
parser.add_argument('-O', '--trim-overhangs', action="store_true", help='''With
this option we trim (by replacing with gaps) bases in one of the sequences to
correct that starts before or ends after any of the other sequences.''')
args = parser.parse_args()

ContigNames = args.NameOfSeqToCorrect

# Check size args are positive
if args.min_contig_size <= 0:
  print('The minimum contig size must be greater than 0. Quitting.',
  file=sys.stderr)
  exit(1)
if args.split_gap_size <= 0:
  print('The split gap size must be greater than 0. Quitting.', file=sys.stderr)
  exit(1)

# Check all contig names are unique
CounterObject = collections.Counter(ContigNames)
DuplicatedArgs = [i for i in CounterObject if CounterObject[i]>1]
if len(DuplicatedArgs) != 0:
  for DuplicatedArg in DuplicatedArgs:
    print('Contig name', DuplicatedArg, 'was duplicated in the arguments.',\
    file=sys.stderr)
  print('All contig names should be unique. Exiting.', file=sys.stderr)
  exit(1)

# Read in the alignment
try:
  AlignedSeqs = AlignIO.read(args.alignment, "fasta")
except:
  print('Problem reading', args.alignment + ':', file=sys.stderr)
  raise

# Tanya's function :)
def split_parts(c, min_contig_size, split_gap_size, prefix=''):
  '''Split sequence at gaps >= split_gap_size, retaining subsequences >=
  min_contig_size. Rename where we keep >1 subsequence.'''
  if type(c) is SeqRecord:
    prefix = c.name
    c = str(c.seq)
  aln_length = len(c)
  pat = '(-{%s,})' % split_gap_size
  contig_bits = re.compile(pat).split(c)
  num_contig_bits_to_keep = sum(1 for bit in contig_bits if bit and \
  len(bit) - bit.count('-') >= min_contig_size)
  rename = num_contig_bits_to_keep > 1
  pos = 0
  i = 0
  for s in contig_bits:
    if s and len(s) - s.count('-') >= min_contig_size:
      i += 1
      if rename:
        label = '{0}_{1}'.format(prefix, i)
      else:
        label = prefix
      yield SeqRecord(seq=Seq('-' * pos + s + \
      '-' * (aln_length - (pos + len(s)))), id=label, name=label,
      description='')
    pos = pos + len(s)

# Check we find all contigs. Find the start of the first and end of the last
# non-contig sequence.
ContigsFound = {contig:False for contig in ContigNames}
NumRefSeqs = 0
AlignmentLength = AlignedSeqs.get_alignment_length()
FirstRefStart = AlignmentLength
LastRefEnd = 0
for seq in AlignedSeqs:
  if seq.id in ContigsFound:
    ContigsFound[seq.id] = True
  else:
    if args.trim_overhangs:
      RefStart, RefEnd = GetSeqStartAndEndPos(str(seq.seq))
      FirstRefStart = min(FirstRefStart, RefStart)
      LastRefEnd    = max(LastRefEnd, RefEnd)
    NumRefSeqs += 1

# Check if any of the named contigs were not found.
MissingContigs = {contig for contig, found in list(ContigsFound.items()) if not found}
if MissingContigs:
  print("The following contigs were not found in " + args.alignment + ":",
  " ".join(MissingContigs) + '\nQuitting.', file=sys.stderr)
  exit(1)


# Iterate through all seqs in the alignment, splitting contigs where necessary.
SeqsForOutput = []
for seq in AlignedSeqs:
  if seq.id in ContigsFound:
    if args.trim_overhangs:
      TrimmedSeq = "-" * FirstRefStart + str(seq.seq)[FirstRefStart:LastRefEnd \
      + 1] + "-" * (AlignmentLength - LastRefEnd - 1)
      assert len(TrimmedSeq) == AlignmentLength, \
      "Internal malfunction of overhang trimming"
      seq.seq = Seq(TrimmedSeq)
    SeqsForOutput += list(split_parts(seq, args.min_contig_size,
    args.split_gap_size))
  else:
    SeqsForOutput.append(seq)

# It's possible that after splitting contigs and imposing a minimum length
# threshold, there are no contigs left. Exit with status 3 - shiver's reserved
# non-zero exit status to indicate a lack of HIV data.
NumContigs = len(SeqsForOutput) - NumRefSeqs
if NumContigs == 0:
  print("After splitting contigs at gaps of length at least",
  args.split_gap_size, "and discarding contigs of length less than " + \
  str(args.min_contig_size) + ", no contigs were left. Quitting.",
  file=sys.stderr)
  exit(3)

# Remove pure-gap columns and print the output.
OutputAlignment = AlignIO.MultipleSeqAlignment(SeqsForOutput)
OutputAlignment = RemoveBlankColumns(OutputAlignment)
AlignIO.write(OutputAlignment, sys.stdout, 'fasta')

