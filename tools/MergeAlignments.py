#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251 
## TODO:
## Overview: this script excises insertions and introduces gaps into  takes as its arguments one 'main' alignment and one
## 'paried' alignment (containing two sequences), both in fasta format. The
## paired alignment contains one sequence which is present in the main
## alignment (the reference) and one which is not.
## 1) In the paired alignment, insertions of the other sequence relative to the
## reference are removed, so that the reference there contains no gaps.
## 2) In the main alignment, unique insertions of the reference relative to all
## other sequences are removed; such positions are also removed from both
## sequences in the paired alignment       Then this gapless  The two versions of the sequence 
## which is present in both are compared; gaps which are present in the main
## alignment but not the paired alignment are inserted into the latter, so
## that the sequence which is not in the main alignment can be added to it.
##
## Output is printed to stdout, suitable for redirection to a target fasta file.
##
################################################################################
## USER INPUT
## Setting the bool below to True means that any gap found in a reference in one
## of the paired alignments will cause an exit with an error. Setting it to
## False means any such gaps are removed, along with the corresponding insertion
## in the paired sequence.
GapInRefIsError = False

GapChar = '-'
# Wrap the sequence in the output fasta file to this length per line
FastaSeqLineLength=50
################################################################################

import os.path, sys
from inspect import currentframe, getframeinfo
from itertools import groupby
from operator import itemgetter
from AuxiliaryFunctions import ReadSequencesFromFile, PropagateNoCoverageChar
import argparse

parser = argparse.ArgumentParser(description='A script for aligning a '+\
'sequence in a pair-wise alignment file to another alignment, by comparing '+\
'the positions of gaps of a second sequence - the reference - which is in '+\
'both alignments.', epilog = 'Either the -d option or the -e option (but not'+\
' both) must be specified.')
parser.add_argument('MainAlignmentFile', help='The fasta file containing an '+\
'alignment of the reference and (many?) other sequences.')
parser.add_argument('PairedAlignmentFile', help='The fasta file containing an'+\
' alignment of the reference and one other sequence - the one you want to add'+\
' to the main alignment file.')
parser.add_argument('-e', '--excise', action='store_true', \
help='To be used when the alignment you want to add to is the main alignment'+\
" *without* the reference (and so the reference's unique insertions in that "+\
"alignment are excised from both alignments).")
parser.add_argument('-d', '--dont-excise', action='store_true', \
help='To be used when the alignment you want to add to '+\
"is the main alignment as is, i.e. including the reference.")
parser.add_argument('-L', '--log-file', help='Used to specify a log file '
'describing the coordinate transformation.')

args = parser.parse_args()

# Check whether we're excising or not.
if (args.excise and args.dont_excise) or \
(not args.excise and not args.dont_excise):
  print('Either the -d option or the -e option (but not both) must be '+\
  'specified.', file=sys.stderr)
  exit(1)

# Rename arguments for brevity / clarity. 
MainAlnFile = args.MainAlignmentFile
PairedAlnFile = args.PairedAlignmentFile
ExciseUniqueInsertionsOfRefInMainAlignment = args.excise

# Check that the arguments exist and are files
for InputFile in [MainAlnFile, PairedAlnFile]:
  if not os.path.isfile(InputFile):
    print(InputFile, 'does not exist or is not a file.', file=sys.stderr)
    exit(1)

# Read in the sequences from the main alignment file (into a dictionary)
MainAlnSeqDict, MainAlnSeqLength = ReadSequencesFromFile(MainAlnFile)
MainAlnSeqNames = list(MainAlnSeqDict.keys())
MainAlnSeqs     = list(MainAlnSeqDict.values())


# Read in the sequences from the paired alignment file
PairedAlnSeqDict, PairedAlnSeqLength = ReadSequencesFromFile(PairedAlnFile)

# Check it has two sequences
if len(PairedAlnSeqDict) != 2:
  print('File', PairedAlnFile, 'contains', len(PairedAlnSeqDict),\
  'sequences; two were expected.\nQuitting.', file=sys.stderr)
  exit(1)
Seq1name, Seq2name = list(PairedAlnSeqDict.keys())

# Check that one of the sequences is in the main alignment file (the 'Ref')
# and one is not (the 'SeqToAdd').
if Seq1name in MainAlnSeqNames:
  if Seq2name in MainAlnSeqNames:
    print('Both sequences from', PairedAlnFile, 'are present in', MainAlnFile+\
    '. Only one should be.\nQuitting.', file=sys.stderr)
    exit(1)
  else:
    RefSeqName, SeqToAddName = Seq1name, Seq2name
else:
  if not Seq2name in MainAlnSeqNames:
    print('Neither sequence from', PairedAlnFile, 'is present in', MainAlnFile+\
    '. Exactly one should be.\nQuitting.', file=sys.stderr)
    exit(1)
  else:
    SeqToAddName, RefSeqName = Seq1name, Seq2name

SeqToAdd = PairedAlnSeqDict[SeqToAddName]
RefSeqFromPair = PairedAlnSeqDict[RefSeqName]
RefSeqFromMain = MainAlnSeqDict[RefSeqName]

# Compare the duplicated sequence to its duplicate in the main alignment.
# Does the reference sequence in the paired alignment contain gaps? If so we
# we either exit, or cut out that insertion, according to the user-specifed
# bool.
# We keep track of how coordinates are changing as we excise those insertions,
# in NewToOldDict.
if GapChar in RefSeqFromPair:
  if GapInRefIsError:
    print('Gap in the version of', RefSeqName, 'in the paired alignment!',\
    'Unexpected.\nQuitting.', file=sys.stderr)
    exit(1)
  else:
    # Loop through every position in RefSeqFromPair. If it's a gap, we ignore 
    # that position in SeqToAdd.
    SeqToAdd_NoInsertions = ''
    RefSeqFromPair_NoGaps = ''
    NewPosition = 0
    NewToOldDict = {}
    for position,BaseInRef in enumerate(RefSeqFromPair):
      if BaseInRef != GapChar:
        RefSeqFromPair_NoGaps += BaseInRef
        SeqToAdd_NoInsertions += SeqToAdd[position]
        NewToOldDict[NewPosition] = position
        NewPosition += 1
    RefSeqFromPair = RefSeqFromPair_NoGaps
    SeqToAdd = SeqToAdd_NoInsertions
else:
  NewToOldDict = {pos:pos for pos in range(len(RefSeqFromPair))}

# TODO: how to keep track of the positions where RefSeqFromPair has a gap, and
# insert them, in order, into TranslationRecord?
# Don't bother. The 'from' coordinates then become hard to interpret. Simply
# ignoring such positions means the 'from' coords are with respect to the
# refseq (and the 'to' coords are with respect to the final alignment, i.e. the
# output.)
 
# Check that the two versions of the ref seq differ only with regards to gaps
# and upper/lower case:
if RefSeqFromMain.replace(GapChar,'').upper() != RefSeqFromPair.upper():
  print('Removing the gaps from the version of sequence', RefSeqName, 'in',\
  MainAlnFile, 'gives a different sequence to the one in the paired',\
  'alignment file. Unexpected.\nQuitting.', file=sys.stderr)
  #print(RefSeqFromMain.replace(GapChar,''), RefSeqFromPair)
  exit(1)

# Loop through every base in RefSeqFromMain. a) If it's a GapChar, we need to 
# add a gap at that position in SeqToAdd, unless every reference in the main
# alignment has a gap there in which case we just skip. b) If it's not a
# GapChar, and is a  unique insertion (i.e. at that position in the main 
# alignment only that  sequence has a base), we want to excise that base from 
# SeqToAdd. c) If it's neither a GapChar nor a unique insertion, we use that 
# base from SeqToAdd.
SeqToAdd_WithGaps = ''
ReferenceWithoutInsertions = ''
PositionInSeqToAdd = 0
PositionInFinalAln = 0
TranslationRecord = '"position w.r.t. final alignment","position w.r.t. ref","base"'
for MainPosition,BaseInMainRef in enumerate(RefSeqFromMain):
  if BaseInMainRef == GapChar:
    EveryBaseHereIsAGap = False
    if ExciseUniqueInsertionsOfRefInMainAlignment:
      EveryBaseHereIsAGap = True
      for SeqName, seq in list(MainAlnSeqDict.items()):
        if SeqName == RefSeqName:
          continue
        elif seq[MainPosition] != GapChar:
          EveryBaseHereIsAGap = False
          break
    if EveryBaseHereIsAGap:
      continue
    SeqToAdd_WithGaps += GapChar
    ReferenceWithoutInsertions += GapChar
    PositionInFinalAln += 1
    TranslationRecord += '\n'+str(PositionInFinalAln)+',-,-'
  else:
    SeqBase = SeqToAdd[PositionInSeqToAdd]
    RefInMainHasUniqueInsertion = False
    if ExciseUniqueInsertionsOfRefInMainAlignment:
      RefInMainHasUniqueInsertion = True
      for SeqName, seq in list(MainAlnSeqDict.items()):
        if SeqName == RefSeqName:
          continue
        elif seq[MainPosition] != GapChar:
          RefInMainHasUniqueInsertion = False
          break
    OldPosition = NewToOldDict[PositionInSeqToAdd]
    RefPosition = PositionInSeqToAdd
    if not RefInMainHasUniqueInsertion:
      SeqToAdd_WithGaps += SeqBase
      ReferenceWithoutInsertions += RefSeqFromPair[PositionInSeqToAdd]
      PositionInFinalAln += 1
      TranslationRecord += '\n' +str(PositionInFinalAln) +',' +\
      str(RefPosition+1) +',' +SeqBase
    else:
      TranslationRecord += '\n-,' +str(RefPosition+1)+','+SeqBase 
    PositionInSeqToAdd += 1

if args.log_file != None:
  with open(args.log_file, 'w') as f:
    f.write(TranslationRecord)

FinalSeqToAdd = PropagateNoCoverageChar(SeqToAdd_WithGaps)

# Thanks Stackoverflow:
def insert_newlines(string, every=FastaSeqLineLength):
    lines = []
    for i in range(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

print('>'+SeqToAddName)
print(insert_newlines(FinalSeqToAdd))
#print('>'+RefSeqName+'_UniqueInsertionsExcised')
#print(insert_newlines(ReferenceWithoutInsertions))



