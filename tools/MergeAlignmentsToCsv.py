#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script takes one or more sequence alignment files as arguments,
and a chosen reference. All other sequences are represented by what they have at each position where the reference has a base; in general this is a kmer, because sequences can have insertions with respect to the reference. e.g. for sequence = ACGT, ref = A--T, we report that the sequence has ACG at the position where the reference has A. Output is printed to stdout, suitable for redirection to a .csv file.'''

import argparse
import os
import sys
from Bio import AlignIO
from collections import OrderedDict
import numpy as np

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('RefName')
parser.add_argument('alignment', type=File, nargs='+')
parser.add_argument('-O', '--offset-for-positions', type=int, help='''Use this
to specify an integer to add the reference position when reporting the results,
i.e. if the first base of the reference present in the alignment is not really
the reference's first base, because some sequence to the left has been chopped
off. If you are including multiple alignments, the same offset will be used for
all of them.''', default=0)
parser.add_argument('-S', '--skip-string', help='''Use this to specify a string;
we will skip any sequence we find whose ID contains that string.''')

args = parser.parse_args()

# The skip string shouldn't be in the ref name.
HaveSkipString = args.skip_string != None
if HaveSkipString and args.skip_string in args.RefName:
  print('The string specified with --skip-string should not be found in',
  RefName + '. Quitting.', file=sys.stderr)
  exit(1)

def GetSeqToRefComparison(File):
  '''TODO'''

  alignment = AlignIO.read(File, 'fasta')
  AlnLength = alignment.get_alignment_length()

  if HaveSkipString:
    alignment = [seq for seq in alignment if not args.skip_string in seq.id]

  # Check that there are no duplicate seq names
  AllSeqNames = []
  for seq in alignment:
    if seq.id in AllSeqNames:
      print('Sequence', seq.id, 'occurs multiple times in', File + \
      '. Sequence names should be unique. Quitting.', file=sys.stderr)
      exit(1)
    AllSeqNames.append(seq.id)

  # Get the ref seq
  if not args.RefName in AllSeqNames:
    print(args.RefName, 'does not appear in', File + '. Quitting.', 
    file=sys.stderr)
    exit(1)
  RefPosition = AllSeqNames.index(args.RefName)
  RefSeq = str(alignment[RefPosition].seq)

  # Initialise a matrix with a row for each seq and a column for each non-gap
  # character of the ref. Set the data type of each element as a generic python
  # object. It will be a string, but of unknown length (at most the alignment
  # length but perhaps better not to reserve that much memory).
  # TODO: find the longest deletion in the ref (after lstrip("-")) and use that
  # +1 as the string length. Check that, using this with different files,
  # mergin the results doesn't truncate longer strings.
  RefLength = AlnLength - RefSeq.count("-")
  if RefLength == 0:
    print(args.RefName, 'in', File, 'contains no bases. Quitting.', 
    file=sys.stderr)
    exit(1)
  NumSeqs = len(AllSeqNames)
  matrix = np.empty(shape=(NumSeqs, RefLength), dtype='object')

  # Now populate the matrix: each element will be what that sequence has at that
  # non-gap char of the ref, may be a base or a gap, or more generally a kmer,
  # because that non-gap char of the ref may be followed by gaps in the ref.
  AlnAsArray = np.array([str(seq.seq) for seq in alignment])
  RefPos0based = -1
  for AlnPos0based, RefBaseHere in enumerate(RefSeq):

    if RefBaseHere == '-':
      continue
    RefPos0based += 1

    # If the reference has one or more gap chars after this position, we want to
    # consider all those positions together with this one. 
    SizeOfRefDeletionAfterThisPos = 0
    while AlnPos0based + SizeOfRefDeletionAfterThisPos < AlnLength - 1 and \
    RefSeq[AlnPos0based + SizeOfRefDeletionAfterThisPos + 1] == '-':
      SizeOfRefDeletionAfterThisPos += 1

    for i in range(NumSeqs):
      KmerHere = AlnAsArray[i][AlnPos0based : AlnPos0based + \
      SizeOfRefDeletionAfterThisPos + 1].replace("-", "")
      if i == RefPosition and KmerHere != RefBaseHere:
        print("Malfunction of ", sys.argv[0], ': inconsistent determination of',
        ' what the reference has at alignment position ', AlnPos0based + 1,
        ' for file ', File, '. Quitting.', sep='', file=sys.stderr)
        exit(1)
      if not KmerHere:
        KmerHere = "-"
      matrix[i, RefPos0based] = KmerHere

  return matrix, AllSeqNames


for FileNum, FastaFile in enumerate(args.alignment):

  matrix, KeyForRowNames = GetSeqToRefComparison(FastaFile)

  if FileNum == 0:
    AllResults = matrix
    KeyForAllRowNames = KeyForRowNames

  else:

    # Check that, for any sequence we're seeing in this file but have already
    # seen (which will always be the case for the reference, but possibly others
    # if present in multiple files), we get the same result. If so, remove them 
    # from the new matrix and its key before updating all the results.
    PreviouslyEncounteredSeqNames = \
    set(KeyForAllRowNames).intersection(KeyForRowNames)
    PreviouslyEncounteredSeqNameIndices = sorted([KeyForRowNames.index(name) \
    for name in PreviouslyEncounteredSeqNames], reverse=True)
    for index in PreviouslyEncounteredSeqNameIndices:
      name = KeyForRowNames[index]
      PositionInPreviousResults = KeyForAllRowNames.index(name)
      OldResult = AllResults[PositionInPreviousResults]
      NewResult = matrix[index]
      if np.array_equal(NewResult, OldResult):
        #matrix = matrix[:index] + matrix[index + 1:]
        matrix = np.delete(matrix, index, 0)
        KeyForRowNames = KeyForRowNames[:index] + \
        KeyForRowNames[index + 1:]
      else:
        if name == args.RefName:
          print('In ', FastaFile, ', the reference sequence ', name,
          ' is different from a in the previous alignment arguments (',
          ' '.join(args.alignment[:FileNum]), '). The reference sequence ',
          'should be identical in all files. Quitting.', sep='',
          file=sys.stderr)
        else:
          print('In ', FastaFile, ', obtained a different result for sequence ',
          name, ' than in one or more of the previous alignment arguments (',
          ' '.join(args.alignment[:FileNum]), '). If the same sequence appears',
          ' in more than one input file, its alignment relative to the ',
          'reference should be identical each time. Quitting.', sep='',
          file=sys.stderr)
        exit(1)

    # Add new results (if there are any) to those so far.
    if matrix.size:
      AllResults = np.append(AllResults, matrix, axis=0)
      KeyForAllRowNames += KeyForRowNames

# Add a row of reference positions, and swap the reference sequence to be the
# first row after that. The temporary assignment is ugly but necessary.
FinalRefPosition = KeyForAllRowNames.index(args.RefName)
RefLength = AllResults.shape[1]
RefPositions = np.array([np.arange(args.offset_for_positions + 1, args.offset_for_positions + RefLength + 1)])
temp = np.append(RefPositions, np.array([AllResults[FinalRefPosition]]), axis=0)
temp = np.append(temp, AllResults[:FinalRefPosition], axis=0)
AllResults = np.append(temp, AllResults[FinalRefPosition + 1:], axis=0)

AllResults = AllResults.transpose()

np.savetxt(sys.stdout, AllResults, delimiter=",", fmt='%s', comments='',
header='Position in ' + args.RefName + ',Base in ' + args.RefName + ',' + \
','.join('kmer in ' + name for name in KeyForAllRowNames if \
name != args.RefName))
    


