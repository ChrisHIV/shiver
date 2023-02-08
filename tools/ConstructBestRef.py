#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251 
##
## Overview
ExplanatoryMessage = '''We construct the optimal reference for a given sample by
flattening its (de novo assembled) contigs with a set of references.
How it works:
(1) We flatten the contigs, taking the base (or gap) of the longest one where
they disagree (expecting that the de novo assembly makes the longest contig
first with the dominant i.e. most numerous set of short reads).
(2) We compare every reference to the flattened contigs and count the number of
agreeing bases. The comparison is not done in the gaps between contigs nor
before/after the reference/contigs. When we are inside a contig and inside the
reference, however, gaps count - here they are taken to represent a genuine
deletion, so the gap character holds equal footing with a base.
(3) Starting with the reference with the highest score, we elongate it both
directions using references or progressively lower scores. We stop when both
edges reach the edges of the alignment, or when a reference is reached which has
a score of zero. This defines the elongated best reference.
(4) We fill in gaps before, between and after (but not within) the contigs using
the elongated best reference. This defines the constructed best reference.
'''

################################################################################
## USER INPUT
# The character that indicates a gap (missing base). If there are multiple such
# characters, put the one to be used in the output first.
GapChar = '-'
# Excise unique insertions present in the contigs but not any of the references?
ExciseUniqueInsertions = False
################################################################################

FloatComparisonTolerance = 1e-6

import os.path
import sys
import collections
import argparse
from AuxiliaryFunctions import ReadSequencesFromFile, PropagateNoCoverageChar

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# For using new lines in the argument help
class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()  
        return argparse.HelpFormatter._split_lines(self, text, width)

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage, \
formatter_class=SmartFormatter)
parser.add_argument('AlignmentOfContigsToRefs', type=File)
parser.add_argument('OutputFile')
parser.add_argument('ContigName', nargs='+')
parser.add_argument('-AS', '--always-use-sequence', action='store_true',
help='''By default, at each position we always use whatever the longest contig
has there, even if it's a deletion (i.e. an internal gap) and a shorter contig
has a base there. With this option we'll always use a base from the contigs if
possible, i.e. if a shorter contig has an insertion relative to a longer contig
we will use the insertion. This protects against cases where a long contig has
an erroneous internal gap due to misalignment, but may introduce artefactual
insertions in the flattened contigs due to misalignment. e.g. By default the
contigs GGGGA-CC- and --TG-ACCT would be flattened to GGGGA-CCT (using the
longer first contig wherever possible); with this option they would be flattened
to GGGGAACCT, giving a double-A that is seen in neither contig.''')
parser.add_argument('-L', '--contigs-length', action='store_true', \
help='Simply print the length of the flattened contigs (i.e. the number of '+\
'positions in the alignment that are inside at least one contig, excluding '+\
'positions where all contigs have a gap), and exit.')
parser.add_argument('-P', '--print-best-score', action='store_true', \
help='Simply print the fractional identity and the name of the reference ' +\
'with the highest fractional identity, then exit.')
parser.add_argument('-S1', '--summarise-contigs-1', action='store_true', \
help='Simply print the length and gap fraction of each contig, then exit.')
parser.add_argument('-S2', '--summarise-contigs-2', action='store_true', \
help='Simply print the length and gap fraction of each contig in an alignment'+\
' in which only the best reference is kept with the contigs (i.e. all other ' +\
'references are discarded, then pure-gap columns are removed), then exit.')
parser.add_argument('-C', '--compare-contigs-to-consensus', \
help='''R|Use this option to specify the name of the consensus
sequence; use it when the AlignmentOfContigsToRefs argument is the
alignment of contigs, consensus and mapping reference.
We put each position in the alignment in one of the
following four categories, then print the counts of 
each. (1) At least one of the flattened contigs agrees
with the consensus. (2) All contigs disagree with the
consensus. (3) At least one contig has a base and the
reference has "?". (4) There is no contig coverage but
the consensus has a base. We do not count positions
where the consensus has "?" and all contigs have gaps or
no coverage, nor positions where the consensus has a gap
and [at least one contig has a gap or no contig has
coverage], since such positions constitute trivial
agreement. We also do not count any position where the
consensus has an 'N'.The table below may help clarify
the categorisation of positions:

                          Consensus:
                       |  base  | gap |  ?  |  N
          -------------|--------|-----|-----|-----
         just bases    | 1 or 2 |  2  |  3  | n/a
Contigs: bases + gaps  | 1 or 2 | n/a |  3  | n/a
         just gaps     |   2    | n/a | n/a | n/a
         no coverage   |   4    | n/a | n/a | n/a
''')
parser.add_argument('-C2', '--check-contig-snps', help='''Use this option to
specify the name of the consensus sequence; use it when the AlignmentOfContigsToRefs argument
is an alignment containing the contigs and the consensus. Amongst all positions
at which more than one base (not including gaps) is represented, we count the
number where the longest contig's base agrees with the consensus (the first
value we'll print to stdout), and the number where it does not (the second value
we'll print to stdout).''')

args = parser.parse_args()

# Shorthand
AlignmentFile = args.AlignmentOfContigsToRefs
ContigNames = args.ContigName
CompareContigsToConsensus = args.compare_contigs_to_consensus != None
CheckContigSNPs = args.check_contig_snps != None
if CompareContigsToConsensus:
  ConsensusName = args.compare_contigs_to_consensus
elif CheckContigSNPs:
  ConsensusName = args.check_contig_snps
else:
  ConsensusName = None

# You would think stripping leading and trailing whitespace off the contig names
# would do nothing, because they are parsed taking whitespace as a delimiter.
# However, if the Mac \r character is amongst them that messes things up. This
# solves it.
ContigNames = [name.strip() for name in ContigNames]

# Can't use these two options together
if CompareContigsToConsensus and CheckContigSNPs:
  print('You cannot use both --compare-contigs-to-consensus and',
  '--check-contig-snps options at once. Quitting.', file=sys.stderr)
  exit(1)

# Check all contig names are unique
CounterObject = collections.Counter(ContigNames)
DuplicatedContigNames = [i for i in CounterObject if CounterObject[i]>1]
if len(DuplicatedContigNames) != 0:
  for ContigName in DuplicatedContigNames:
    print('Contig name', ContigName, 'was duplicated in the arguments.', \
    file=sys.stderr)
  print('All contig names should be unique. Exiting.', file=sys.stderr)
  exit(1)

# Check the consensus name does not match one fo the contig names.
if ConsensusName != None and ConsensusName in ContigNames:
  print('The consensus name should not be the same as one of the contig', \
  'names. Quitting.', file=sys.stderr)
  exit(1)

# Read in the sequences from the alignment file (into a dictionary)
AllSeqsDict, AlignmentLength = ReadSequencesFromFile(AlignmentFile)

# Check the consensus is found
if ConsensusName != None:
  if not ConsensusName in AllSeqsDict:
    print(ConsensusName, 'not found in', AlignmentFile + '. Quitting.', \
    file=sys.stderr)
    exit(1)
  ConsensusSeq = AllSeqsDict[ConsensusName]

# Separate sequences into references and contigs
RefDict = {}
ContigDict = {}
for SeqName in AllSeqsDict:
  if SeqName in ContigNames:
    ContigDict[SeqName] = AllSeqsDict[SeqName]
  else:
    RefDict[SeqName] = AllSeqsDict[SeqName]

# Check we found all the contigs
MissingContigs = [contig for contig in ContigNames if not contig in ContigDict]
if len(MissingContigs) != 0:
  for ContigName in MissingContigs:
    print('Contig "'+ContigName+'" was not found in the alignment file.', \
    file=sys.stderr)
  print('Exiting.', file=sys.stderr)
  exit(1)

# A function we'll need more than once:
def FindSeqStartAndEnd(SeqName, seq, AlignmentLength, FileName):
  '''Find the 0-based positions of the start and end of the sequence.'''
  StartOfSeq = 0
  try:
    while seq[StartOfSeq] == GapChar:
      StartOfSeq += 1
  except IndexError:
    print(SeqName, "in", FileName, \
    "has no bases - it's just one big gap. Quitting.", file=sys.stderr)
    exit(1)
  EndOfSeq = AlignmentLength-1
  while seq[EndOfSeq] == GapChar:
    EndOfSeq -= 1
  return StartOfSeq,EndOfSeq

# Find the start, end, length and gap fraction of each contig.
# Find the start of the left-most ('first') contig and the end of the right-most
# ('last') contig.
ContigStartsAndEnds = {}
ContigLengths = {}
ContigGapFractions = {}
TotalGapsInContigs = 0
for ContigName,ContigSeq in list(ContigDict.items()):
  StartOfContig, EndOfContig = FindSeqStartAndEnd(ContigName, ContigSeq, \
  AlignmentLength, AlignmentFile)
  ContigStartsAndEnds[ContigName] = [StartOfContig,EndOfContig]
  NumGapsInAlignedContig = ContigSeq.count(GapChar)
  NumInternalGaps = ContigSeq[StartOfContig:EndOfContig+1].count(GapChar)
  TotalGapsInContigs += NumInternalGaps
  ContigLengths[ContigName] = AlignmentLength -NumGapsInAlignedContig
  ContigGapFractions[ContigName] = \
  float(NumInternalGaps)/(EndOfContig+1 - StartOfContig)

if args.summarise_contigs_1:
  for ContigName, length in list(ContigLengths.items()):
    sys.stdout.write(str(length) + ' ' + str(ContigGapFractions[ContigName]))
  exit(0)

#print(' '.join(map(str, sorted(ContigGapFractions.values(), reverse=True) )))
#print(float(TotalGapsInContigs)/sum(ContigLengths.values()))
#exit(0)

AllContigStarts = []
AllContigEnds = []
for [ContigStart,ContigEnd] in list(ContigStartsAndEnds.values()):
  AllContigStarts.append(ContigStart)
  AllContigEnds.append(ContigEnd)
StartOfFirstContig = min(AllContigStarts)
EndOfLastContig = max(AllContigEnds)

#GappiestContigName, GapFraction = \
#sorted(ContigGapFractions.items(), key=lambda x:x[1], reverse=True)[0]

# Flatten the contigs.
# Repeat the gap character until the start of the first contig.
# Then at each position: if all contigs agree on what's there, use that;
# otherwise, use either [whatever is present in the longest contig there] OR
# [whatever is present in the longest contig there provided it's not a gap], as
# desired. Then repeat the gap character until the end.
FlattenedContigsSeq = GapChar * StartOfFirstContig
for position in range(StartOfFirstContig,EndOfLastContig+1):
  DictOfBasesHere = {}
  for ContigName,ContigSeq in list(ContigDict.items()):
    DictOfBasesHere[ContigName] = ContigSeq[position]
  BasesHere = set(DictOfBasesHere.values())
  if len(BasesHere) == 1:
    BaseHere = list(DictOfBasesHere.values())[0]
  else:
    LengthOfLongestDesiredContig = 0
    for ContigName in DictOfBasesHere:
      StartOfContig, EndOfContig = ContigStartsAndEnds[ContigName]
      if StartOfContig <= position <= EndOfContig and \
      ContigLengths[ContigName] > LengthOfLongestDesiredContig:
        ThisContigsBase = DictOfBasesHere[ContigName]
        if not args.always_use_sequence or ThisContigsBase != GapChar:
          BaseHere = ThisContigsBase
          LengthOfLongestDesiredContig = ContigLengths[ContigName]
    if LengthOfLongestDesiredContig == 0:
      print('Malfunction of', sys.argv[0] + ": we failed to figure out what",
      "base to use from the contigs at ", position+1, 'in', AlignmentFile + \
      '. Please report to Chris Wymant. Quitting.', file=sys.stderr)
      exit(1)
  FlattenedContigsSeq += BaseHere
FlattenedContigsSeq += GapChar * (AlignmentLength - EndOfLastContig -1)

if args.contigs_length:
  print(len(FlattenedContigsSeq) - FlattenedContigsSeq.count(GapChar))
  exit(0)

# Make a list, of the same length of the alignment, of integers: each one
# counting the number of contigs with coverage there. Gaps inside contigs get
# counted as coverage; gaps between contigs get a count of 0.
ContigCoverageByPosition = [0 for n in range(0,AlignmentLength)]
for [start,end] in list(ContigStartsAndEnds.values()):
  for position in range(start,end+1):
    ContigCoverageByPosition[position] += 1

# Compare contigs to consensus to count positions in the four categories
# ('cats') described in the help above. Convert all bases to upper case, and
# replace any gap char that neighbours a '?' char in the consensus by a '?'.
if ConsensusName != None:
  categories = []
  ConsensusSeq = ConsensusSeq.upper()
  ConsensusSeq = PropagateNoCoverageChar(ConsensusSeq)
  FlattenedContigsSeq = FlattenedContigsSeq.upper()
  for ContigName in ContigDict:
    ContigDict[ContigName] = ContigDict[ContigName].upper()

  if CheckContigSNPs:
    NumPosLongestContigCorrect = 0
    NumPosLongestContigIncorrect = 0
    for pos in range(0,AlignmentLength):
      if ContigCoverageByPosition[pos] < 2:
        continue
      ContigBases = []
      for ContigName, ContigSeq in list(ContigDict.items()):
        StartOfContig, EndOfContig = ContigStartsAndEnds[ContigName]
        if StartOfContig <= pos <= EndOfContig:
          ContigBase = ContigSeq[pos]
          if ContigBase != '-':
            ContigBases.append(ContigSeq[pos])
      if len(set(ContigBases)) >= 2:
        # There are at least two different bases in the contigs here.
        #print('pos', pos+1, 'consensus:', ConsensusSeq[pos], 'longest contig:', 
        #FlattenedContigsSeq[pos], 'contig bases here:', ContigBases) 
        if FlattenedContigsSeq[pos] == ConsensusSeq[pos]:
          NumPosLongestContigCorrect += 1
        else:
          NumPosLongestContigIncorrect += 1
    print(NumPosLongestContigCorrect, NumPosLongestContigIncorrect)
    exit(0)


  if CompareContigsToConsensus:
    for pos in range(0,AlignmentLength):
      ConsensusBase = ConsensusSeq[pos]
      if ConsensusBase == 'N': 
        cat = None
      elif ContigCoverageByPosition[pos] == 0:
        if ConsensusBase == '?' or ConsensusBase == GapChar:
          cat = None
        else:
          cat  = 4
      else:

        # Find all bases (or gaps) inside a contig here.
        ContigBases = []
        for ContigName, ContigSeq in list(ContigDict.items()):
          StartOfContig, EndOfContig = ContigStartsAndEnds[ContigName]
          if StartOfContig <= pos <= EndOfContig:
            ContigBases.append(ContigSeq[pos])
        if ContigBases == []:
          print('Malfunction of', sys.argv[0] + ": we've lost track of what", \
          "contigs are present at position", pos+1, 'in', AlignmentFile + \
          '. Please report to Chris Wymant. Quitting.', file=sys.stderr)
          exit(1)

        if ConsensusBase == '?':
          if all(ContigBase == GapChar for ContigBase in ContigBases):
            cat = None
          else:
            cat = 3
        elif ConsensusBase == GapChar:
          if all(ContigBase != GapChar for ContigBase in ContigBases):  
            cat = 2
          else:
            cat = None
        else:
          if any(ContigBase == ConsensusBase for ContigBase in ContigBases):  
            cat = 1
          else:
            cat = 2
          
      categories.append(cat)
    CatCounts = [0,0,0,0]
    for cat in categories:
      if cat != None:
        CatCounts[cat-1] += 1
    print(' '.join(map(str, CatCounts)))
    exit(0)
    

#TotalContigCoverage = sum([1 for base in FlattenedContigsSeq if base != GapChar])
#sys.stdout.write(str(len(ContigDict)) +' '+ str(TotalContigCoverage) +' ')
#exit(0)

# Count the number of positions where exactly one contig has a base and no
# reference has a base - 'unique insertions'. Replace such positions by gaps if
# desired.
# Count how many deletions there are inside contigs.
FlattenedContigsSeq_NoNewInsertions = ''
LengthOfUniqueInsertions = 0
LengthOfDeletions = 0
for position in range(0,AlignmentLength):
  NumContigsHere = ContigCoverageByPosition[position]
  ContigBase = FlattenedContigsSeq[position]
  if NumContigsHere > 0 and ContigBase == GapChar:
    LengthOfDeletions += 1
  if NumContigsHere == 1:
    NoRefHasBaseHere = True
    for RefName,RefSeq in list(RefDict.items()):
      if RefSeq[position] != GapChar:
        NoRefHasBaseHere = False
        break
    if NoRefHasBaseHere:
      LengthOfUniqueInsertions += 1
    if ExciseUniqueInsertions and NoRefHasBaseHere:
      FlattenedContigsSeq_NoNewInsertions += GapChar
    else:
      FlattenedContigsSeq_NoNewInsertions += ContigBase
  else:
    FlattenedContigsSeq_NoNewInsertions += ContigBase
if ExciseUniqueInsertions:
  FlattenedContigsSeq = FlattenedContigsSeq_NoNewInsertions

# For each reference, find the fraction of positions with both reference and
# contig coverage where the two are in agreement. Record the reference start and
# end.
ListOfRefsAndScores = []
for RefName,RefSeq in list(RefDict.items()):
  NumBasesAgreeing = 0
  OverlapLength = 0
  StartOfRef, EndOfRef = FindSeqStartAndEnd(RefName, RefSeq, AlignmentLength,
  AlignmentFile)
  for position in range(StartOfRef, EndOfRef+1):
    if ContigCoverageByPosition[position] > 0:
      OverlapLength += 1
      if RefSeq[position] == FlattenedContigsSeq[position]:
        NumBasesAgreeing += 1
  if OverlapLength == 0:
    FractionalAgreement = 0
  else:
    FractionalAgreement = float(NumBasesAgreeing)/OverlapLength
  ListOfRefsAndScores.append([RefName, StartOfRef, EndOfRef, \
  FractionalAgreement])
  #NumBasesAgreeing])

# Check at least one reference matches at least one position!
if all(item[3] < FloatComparisonTolerance for item in ListOfRefsAndScores):
  print('No reference matches the contigs at any position! We assume this is',\
  'is an error.\nQuitting.', file=sys.stderr)
  exit(1)

# Sort the references - the ones closest to the contigs first.
ListOfRefsAndScores = sorted(ListOfRefsAndScores, key=lambda x:x[3], \
reverse=True)

BestRefName, BestRefStart, BestRefEnd, BestRefScore = ListOfRefsAndScores[0]

if args.print_best_score:
  print(BestRefScore, BestRefName)
  exit(0)

if args.summarise_contigs_2:
  ContigsWithBestRef = [RefDict[BestRefName]] + list(ContigDict.values())
  NumSeqs = len(ContigsWithBestRef)
  for column in range(AlignmentLength-1,-1,-1):
    if all(seq[column] == '-' for seq in ContigsWithBestRef):
      for i in range(NumSeqs):
        ContigsWithBestRef[i] = ContigsWithBestRef[i][:column] + \
        ContigsWithBestRef[i][column+1:]
  ThisAlignmentLength = len(ContigsWithBestRef[0])
  assert all(len(seq) == ThisAlignmentLength for seq in ContigsWithBestRef)
  #for i in range(NumSeqs):
  #  print('>'+str(i+1))
  #  print(ContigsWithBestRef[i])
  GapFracsAndLengths = []
  for contig in ContigsWithBestRef[1:]:
    start, end = FindSeqStartAndEnd('contig', contig, ThisAlignmentLength,
    AlignmentFile)
    NumGaps = contig[start:end+1].count('-')
    GapFrac = float(NumGaps)/(end - start + 1)
    Length = end - start + 1 - NumGaps
    GapFracsAndLengths.append(str(Length) +',' + str(GapFrac))
  print('\n'.join(GapFracsAndLengths))
  exit(0)

# Start with the best ref, and iteratively extend it using longer references
# with lower scores. Break if we reach a reference with zero score, or if our
# construct becomes as long as the alignment.
ElongatedRef = RefDict[BestRefName]
ElongatedRefStart, ElongatedRefEnd = BestRefStart, BestRefEnd
for RefName, StartOfRef, EndOfRef, NumBasesAgreeing in ListOfRefsAndScores[1:]:
  if NumBasesAgreeing == 0 or \
  (ElongatedRefStart == 0 and ElongatedRefEnd == AlignmentLength-1):
    break
  ThisRef = RefDict[RefName]
  if StartOfRef < ElongatedRefStart:
    ElongatedRef = \
    ThisRef[:ElongatedRefStart] + ElongatedRef[ElongatedRefStart:]
    ElongatedRefStart = StartOfRef
  if EndOfRef > ElongatedRefEnd:
    ElongatedRef = \
    ElongatedRef[:ElongatedRefEnd+1] + ThisRef[ElongatedRefEnd+1:]
    ElongatedRefEnd = EndOfRef

# Fill in gaps in contig coverage using the elongated best reference from above.
ConstructedRef = ''
for position in range(0,AlignmentLength):
  if ContigCoverageByPosition[position] > 0:
    ConstructedRef += FlattenedContigsSeq[position]
  else:
    ConstructedRef += ElongatedRef[position]

# Inserting line breaks: thanks Stackoverflow:
def insert_newlines(string, every=50):
  lines = []
  for i in range(0, len(string), every):
    lines.append(string[i:i+every])
  return '\n'.join(lines)

LengthOfFlattenedContigs = len(FlattenedContigsSeq) - FlattenedContigsSeq.count('-')
LengthOfConstructedRef = len(ConstructedRef) - ConstructedRef.count('-')
NumBasesOfElongatedRefUsed = LengthOfConstructedRef - LengthOfFlattenedContigs
print('Info: using', BestRefName, '(elongated with other longer references if',
'needed) to provide', NumBasesOfElongatedRefUsed, 'bases to fill in gaps',
'before/between/after contigs. This reference was the best match for the',
'contigs:', str(100*BestRefScore) + '% of positions in agreement.')

# Print output.

with open(args.OutputFile, 'w') as f:
  f.write('>ContigsFlattenedWith_'+BestRefName + '\n' + \
  insert_newlines(ConstructedRef) + '\n>' + BestRefName + '_elongated' + \
  insert_newlines(ElongatedRef) + '\n')

