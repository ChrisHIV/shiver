#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: this script is called from the command line with a pileup file and
## a directory of reference genome files as its two arguments. We iterate
## through the pileup file and calculate the percentages of bases found at each
## position. The reference is read in from the directory, and compared to the
## information in the pileup file as a sanity check, and also to fill in any
## positions missing from the pileup file. The output columns are:
## Position with respect to the reference; the reference's base; number of short
## reads mapping here; a string listing the bases found here, in order from the
## most common to the least common; then the frequencies of those bases, in that
## same order.
## Extra positions are added in between two reference positions if the most
## common insertion size (with 0 always a possibility) is greater than zero;
## only bases inside an insertion of exactly the right size will count towards
## base frequencies there.
## Quality information is ignored.

################################################################################
# USER INPUT
# The filename extension for reference genome files
RefGenomeExtension = '.fasta'
################################################################################

import os.path, sys, collections, operator
from AuxiliaryFunctions import ReadSequencesFromFile

# Check this file is called from the command line with two arguments
if len(sys.argv) != 3:
  print('This script,', os.path.basename(sys.argv[0])+', requires two',\
  'arguments: firstly a pileup file, then EITHER a fasta file containing the',\
  'reference OR a directory of reference fasta files (one of which we want).', \
  file=sys.stderr)
  exit(1)
PileupFile = sys.argv[1]
ReferenceFileOrDir = sys.argv[2]

# Check that the PileupFile exists and is a file
if not os.path.isfile(PileupFile):
  print(PileupFile, 'does not exist or is not a file.', file=sys.stderr)
  exit(1)

def ReadReferenceFromFile(File):
  '''Read in all sequences in the reference file; check there is only one.'''
  AllSequences, ReferenceLength = ReadSequencesFromFile(File,False)
  if len(AllSequences) != 1:
    print('Found', len(AllSequences), 'sequences in', ReferenceFile+\
    '; expected 1.\nQuitting.', file=sys.stderr)
    exit(1)
  return list(AllSequences.items()), ReferenceLength

# If a reference file was given, read it in; otherwise just check a directory
# was specified.
RefMustBeFoundInDir=False
if os.path.isdir(ReferenceFileOrDir):
  RefMustBeFoundInDir=True
else:
  if not os.path.isfile(ReferenceFileOrDir):
    print(ReferenceFileOrDir, 'is neither a file nor a directory.\nQuitting.', \
    file=sys.stderr)
    exit(1)
  ReferenceFile = ReferenceFileOrDir
  [(RefNameInFasta, RefSeq)], ReferenceLength = \
  ReadReferenceFromFile(ReferenceFileOrDir)
  RefSeq = RefSeq.upper()
  
ExpectedBases = ['A', 'C', 'G', 'T', '-', 'N']
MissingCoverageBaseCounts = [0] * len(ExpectedBases)
MissingCoverageBaseCountsAsStr = ','.join(map(str,MissingCoverageBaseCounts))

def ProcessBaseCounts(BaseCounts, ReferencePosition, ReferenceBase, \
IsInsertion, ReferenceLength):
  '''Counts each kind of base present, with sanity checks. Returns
  a string of the bases present in order of their counts, and the counts.'''

  # If this isn't an insertion, we may have gaps: rename them from * to -
  if (not IsInsertion) and '*' in BaseCounts:
    BaseCounts['-'] = BaseCounts['*']
    del BaseCounts['*']

  # If this isn't an insertion, '.' and ',' both count as the reference base.
  if not IsInsertion:
    ForwardReadsMatchReference  = '.' in BaseCounts
    BackwardReadsMatchReference = ',' in BaseCounts
    if ForwardReadsMatchReference or BackwardReadsMatchReference:
      if ReferenceBase in BaseCounts:
        print('Unexpected behavifour: reference base', ReferenceBase, \
        'explicitly mentioned in the pileup string, instead of "," or ".".'+\
        '\nQuitting.', file=sys.stderr)
        exit(1)
      NumBasesMatchingReference = 0
      if ForwardReadsMatchReference:
        NumBasesMatchingReference += BaseCounts['.']
        del BaseCounts['.']
      if BackwardReadsMatchReference:
        NumBasesMatchingReference += BaseCounts[',']
        del BaseCounts[',']
      BaseCounts[ReferenceBase] = NumBasesMatchingReference

  # Warn about unexpected bases.
  for base in list(BaseCounts.keys()):
    if not base in ExpectedBases:
      warning = 'WARNING: unexpected base '+base+' occurs '+\
      str(BaseCounts[base])+' times '
      if IsInsertion:
        warning += 'in an insertion '
      warning += 'at reference position '+str(ReferencePosition)+'.'
      print(warning, file=sys.stderr)

  NumberOfReads = sum(BaseCounts.values())

  # A list of things we want to record for this position.
  # Firstly, what position is it with respect to the reference.
  # Then, what base does the reference have.
  # Then, how many short reads do we have here.
  # Then, a string listing the bases found here, in order from the most common
  # to the least common.
  # Then, the frequencies of those bases, in that same order.
  if IsInsertion:
    ReferencePosition = 'NA'
  SummaryList = [ReferencePosition, ReferenceBase]
  for base in ExpectedBases:
    try:
      count = BaseCounts[base]
    except KeyError:
      count = 0
    SummaryList.append(count)

  return SummaryList
  # TODO: Change references to
  # specific positions in SummaryList.




# Analyse the pileup file, line by line.
InfoFromAllPositions = []
with open(PileupFile, 'r') as f:
  for LineNumberMin1,line in enumerate(f):

    # Separate the line into fields based on whitespace
    fields = line.split()
    NumFields = len(fields)

    # Check the number of fields is four (if coverage is zero and the samtools
    # version was < 1.4), five, or six (final quality field is optional):
    if not NumFields in [4,5,6]:
      print('Expected 4, 5 or 6 fields; encountered', NumFields, 'on line', \
      str(LineNumberMin1+1)+'.\nQuitting.', file=sys.stderr)
      exit(1)

    # If there's no coverage, check the fields are consistent with that:
    # depending on the version of samtools, there may or may not be a fifth
    # 'placeholder' field.
    try:
      NumReads = int(fields[3])
    except ValueError:
      print('On line', str(LineNumberMin1+1) + ', could not understand the',
      'fourth field,', fields[3] + ', as an integer. Quitting.',
      file=sys.stderr)
      exit(1)
    assert NumReads >= 0, 'Number of mapped reads must be positive'
    NoCoverage = NumReads == 0
    if NumFields > 4:
      PileupString = fields[4]
      if NoCoverage and PileupString != '*':
        print('On line ', LineNumberMin1+1, ', unexpected fifth field "',
        PileupString, '" given that the number of mapped reads is zero (',
        "expected either nothing or the placeholder '*'). Quitting.", sep='',
        file=sys.stderr)
        exit(1)

    # On line 1, read the reference name. Check the file exists then read it in.
    if LineNumberMin1 == 0:
      RefNameInPileup = fields[0]
      if RefMustBeFoundInDir:
        ReferenceFile = \
        os.path.join(ReferenceFileOrDir, RefNameInPileup+RefGenomeExtension)
        if not os.path.isfile(ReferenceFile):
          print('Line 1 quotes a reference', RefNameInPileup+', but', \
          ReferenceFile, 'does not exist or is not a file.\nQuitting.', \
          file=sys.stderr)
          exit(1)
        [(RefNameInFasta, RefSeq)], ReferenceLength = \
        ReadReferenceFromFile(ReferenceFile)
        RefSeq = RefSeq.upper()

    # On lines after line 1, check the reference is the same.
    elif fields[0] != RefNameInPileup:
      print('ERROR: a different reference is reported on the line\n'+\
      line.rstrip()+'\nQuitting.', file=sys.stderr)
      exit(1)

    # Check the reference base here reported by the pileup file matches the one
    # in the reference file.
    BasePosition = int(fields[1])
    ReferenceBase = fields[2].upper()
    if ReferenceBase != RefSeq[BasePosition-1]:
      print('The pileup file', PileupFile, 'reports a base "'+ReferenceBase+\
      '" at position', BasePosition, 'of the reference, but', ReferenceFile, \
      'has a base', RefSeq[BasePosition-1], 'here.\nQuitting.', file=sys.stderr)
      exit(1)

    # If no coverage, we'll record '?' for the string describing the bases
    # appearing here and 'NA' for the frequencies (explained later), and skip.
    if NoCoverage:
      SummaryList = [BasePosition, ReferenceBase] + MissingCoverageBaseCounts
      InfoFromAllPositions.append(SummaryList)
      continue

    # We will only consider adding a new column(s) in between two reference
    # positions if more than half of the reads have an insertion here. If not,
    # we'll save time by not keeping track of the insertion information.
    NumReadsWithInsertion = PileupString.count('+')
    NumReadsWithoutInsertion = NumReads - NumReadsWithInsertion
    MostReadsHaveInsertion = NumReadsWithInsertion > NumReadsWithoutInsertion

    # Iterate through PileupString, modifying the 
    # base counts by interpreting the pileup format appropriately.
    # The character following a ^ should be ignored. $ should be ignored.
    # After a '+' or '-' then a number, that number of characters should be
    # ignored.
    PileupString_OnlyBases = ''
    SkipThisManyBases = 0
    insertions = {}
    for position,char in enumerate(PileupString):
      if SkipThisManyBases > 0:
        SkipThisManyBases -= 1
        continue
      if char == '$':
        continue
      if char == '^':
        SkipThisManyBases = 1
        continue
      if char == '-' or char == '+':
        # Next in the string will be an unknown number of digits, together
        # comprising a number specifying the indel size.
        IndelSizeNumDigits = 0
        while PileupString[position+1+IndelSizeNumDigits] in \
        ['0','1','2','3','4','5','6','7','8','9']:
          IndelSizeNumDigits += 1
        IndelSize = int(PileupString[position+1:position+1+IndelSizeNumDigits])
        SkipThisManyBases = IndelSizeNumDigits+IndelSize

        # If most reads have insertions, we record the insertions: as a
        # dictionary of lists, indexed by the length of the insertion. e.g.
        # insertions[1] = ['A','A','C']; insertions[2] =  ['GG','CT'] etc.
        if char == '+' and MostReadsHaveInsertion:
          insertion = PileupString[position+1+IndelSizeNumDigits:\
          position+1+IndelSizeNumDigits+IndelSize]
          if IndelSize in insertions:
            insertions[IndelSize].append(insertion)
          else:
            insertions[IndelSize] = [insertion]
        continue

      # If we get to here, the character is just a regular base.
      PileupString_OnlyBases += char

    # Start by counting all the different characters at this position.
    # Ignore case.
    PileupString_OnlyBases = PileupString_OnlyBases.upper()
    BaseCounts = collections.Counter(PileupString_OnlyBases)

    # Process the base counts. ('False' means this position isn't an insertion.)
    SummaryList = ProcessBaseCounts(BaseCounts, BasePosition, ReferenceBase, \
    False, ReferenceLength)

    # Check that our total number of bases matches what the pileup file says.
    OurNumReads = sum(SummaryList[2:])
    if OurNumReads != NumReads:
      print('Error in the interpretation of the pileup string at line',\
      str(LineNumberMin1+1)+': counted', OurNumReads, \
      ' bases whereas the pileup file itself says', str(NumReads)+\
      '.\nQuitting.', file=sys.stderr)
      exit(1)

    InfoFromAllPositions.append(SummaryList)

    # Find the most common insertion size here.
    MostCommonInsertionSize = 0
    NumberOfReadsWithMostCommonInsertionSize = NumReadsWithoutInsertion
    for InsertionSize,AllInsertionsThatSize in list(insertions.items()):
      if len(AllInsertionsThatSize) > NumberOfReadsWithMostCommonInsertionSize:
        MostCommonInsertionSize = InsertionSize
        NumberOfReadsWithMostCommonInsertionSize = len(AllInsertionsThatSize)

    # If the most common insertion size here is greater than zero, add in that
    # many columns between this reference position and the next. For these
    # columns, we consider only those reads with an insertion of exactly that
    # size, ignoring all other reads (and not counting them in the 'number of 
    # reads mapping here' field).
    if MostCommonInsertionSize > 0:
      InsertionsConsidered = insertions[MostCommonInsertionSize]
      for i in range(0,len(InsertionsConsidered)):
        InsertionsConsidered[i] = InsertionsConsidered[i].upper()
      for PositionInInsertion in range(0,MostCommonInsertionSize):
        BasesHere = [insertion[PositionInInsertion] for insertion in \
        InsertionsConsidered]
        BaseCounts = collections.Counter(BasesHere)
        SummaryList = ProcessBaseCounts(BaseCounts, BasePosition, \
        '-', True, ReferenceLength)
        InfoFromAllPositions.append(SummaryList)

# Check we've got information from at least one position
if len(InfoFromAllPositions) == 0:
  print('Found no pileup information. Quitting.', file=sys.stderr)
  exit(1)

print('position in ', RefNameInPileup, ',base in ', RefNameInPileup, \
',A count,C count,G count,T count,gap count,N count', sep='')

RightMostReferencePositionSoFar = 0
for PositionWithPileup in InfoFromAllPositions:
  ReferencePosition = PositionWithPileup[0]

  if ReferencePosition == 'NA':
    # This is an insertion with respect to the reference  
    print(','.join(map(str,PositionWithPileup)))
    continue

  if ReferencePosition != RightMostReferencePositionSoFar+1:
    # We have skipped some positions with respect to the reference, due to
    # having no pileup information there. Let's fill in those blanks.  
    for SkippedPosition in range(RightMostReferencePositionSoFar+1,\
    ReferencePosition):
      print(SkippedPosition, RefSeq[SkippedPosition-1], \
      MissingCoverageBaseCountsAsStr, sep=',')

  print(','.join(map(str,PositionWithPileup)))
  RightMostReferencePositionSoFar = ReferencePosition

# Include any missing lines after the data finishes
if RightMostReferencePositionSoFar != ReferenceLength:
  for SkippedPosition in range(RightMostReferencePositionSoFar+1,\
    ReferenceLength+1):
      print(SkippedPosition, RefSeq[SkippedPosition-1], \
      MissingCoverageBaseCountsAsStr, sep=',')


'''
# Now print the information from all positions with respect to the reference.
# Positions that have no coverage may not appear in the pileup file:
# we print '?' for the bases-appearing-here string and 'NA' for the frequencies.
NumPostionsWithPileup = len(InfoFromAllPositions)
NumPostionsWithPileupSoFar = 0
NextPositionWithPileup = InfoFromAllPositions[NumPostionsWithPileupSoFar][0]
for PositionMin1,RefBase in enumerate(ReferenceSeq):

  if NextPositionWithPileup = 'NA'
  if (PositionMin1+1) == NextPositionWithPileup:
    # We've got pileup information for this position.
    print ' '.join(map(str,InfoFromAllPositions[NumPostionsWithPileupSoFar]))
    NumPostionsWithPileupSoFar += 1
    if NumPostionsWithPileupSoFar == NumPostionsWithPileup:
      break
    NextPositionWithPileup = InfoFromAllPositions[NumPostionsWithPileupSoFar][0]

  else:
    print PositionMin1+1, RefBase, 0, '? NA'
'''
