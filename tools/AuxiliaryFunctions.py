

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: here we define some functions: for reading in all sequences from
## a file into a dictionary, for reading in patient details from a file
## into a dictionary, and for checking if two bases match while allowing for one
## or both to be ambiguous.

import os, sys, collections

# A dictionary of what the IUPAC ambiguous letters mean
IUPACdict = {}
IUPACdict['M'] = ['A','C']
IUPACdict['R'] = ['A','G']
IUPACdict['W'] = ['A','T']
IUPACdict['S'] = ['C','G']
IUPACdict['Y'] = ['C','T']
IUPACdict['K'] = ['G','T']
IUPACdict['V'] = ['A','C','G']
IUPACdict['H'] = ['A','C','T']
IUPACdict['D'] = ['A','G','T']
IUPACdict['B'] = ['C','G','T']
IUPACdict['N'] = ['A','C','G','T']
#for AmbigLetter in IUPACdict:
#  IUPACdict[AmbigLetter] = sorted(IUPACdict[AmbigLetter])


# The reverse dictionary: the set of all letters which could mean 'a', the set
# of all letters which could mean 'c', etc.
acgt = ['A','C','G','T']
ReverseIUPACdict = {}
for UnambigLetter in acgt:
  ReverseIUPACdict[UnambigLetter] = [UnambigLetter]
  for AmbigLetter,TargetLetters in list(IUPACdict.items()):
    if UnambigLetter in TargetLetters:
      ReverseIUPACdict[UnambigLetter].append(AmbigLetter)


# Another kind of reverse dictionary, giving the ambiguity code for a set of
# unambiguous bases.
ReverseIUPACdict2 = {}
ReverseIUPACdict2['AG']  = 'R'
ReverseIUPACdict2['CT']  = 'Y'
ReverseIUPACdict2['AC']  = 'M'
ReverseIUPACdict2['GT']  = 'K'
ReverseIUPACdict2['AT']  = 'W'
ReverseIUPACdict2['CG']  = 'S'
ReverseIUPACdict2['CGT'] = 'B'
ReverseIUPACdict2['AGT'] = 'D'
ReverseIUPACdict2['ACT'] = 'H'
ReverseIUPACdict2['ACG'] = 'V'
ReverseIUPACdict2['ACGT'] = 'N'


def InterpretIUPAC(MyDict):
  '''The argument should be a dictionary whose keys are (some subset of) the DNA
  letters and IUPAC ambiguity codes, with the values being numeric. The value
  for an ambiguity code key is divided equally between the letters involved in
  the ambiguity code.'''

  keys = list(MyDict.keys())
  UpdatedDict = {}
  for UnambigLetter in acgt:
    if UnambigLetter in keys:
      UpdatedDict[UnambigLetter] = MyDict[UnambigLetter]

  for AmbigLetter,TargetLetters in list(IUPACdict.items()):
    if AmbigLetter in keys:
      WeightPerTargetLetter = float(MyDict[AmbigLetter])/len(TargetLetters)
      for TargetLetter in TargetLetters:
        if TargetLetter in UpdatedDict:
          UpdatedDict[TargetLetter] += WeightPerTargetLetter
        else:
          UpdatedDict[TargetLetter] = WeightPerTargetLetter
  return UpdatedDict




# Two bases match if: they are equal, or if one is an ambiguity code coding for
# the other, or if they are both ambiguity codes coding for at least one shared
# base. 
def BaseMatch(base1,base2):
  '''Checks if base1 = base2, allowing for one or both to be ambiguous.'''
  if base1 == base2:
    return True
  if base1 in IUPACdict:
    if base2 in IUPACdict:
      return any(i in IUPACdict[base1] for i in IUPACdict[base2])
    return base2 in IUPACdict[base1]
  if base2 in IUPACdict:
    return base1 in IUPACdict[base2]
  return False


def CallAmbigBaseIfNeeded(bases, coverage, MinCovForUpper, BaseFreqFile):
  '''If several bases are supplied, calls an ambiguity code. Uses upper/lower
  case appropriately.'''

  bases = ''.join(sorted(bases))
  NumBases = len(bases)
  assert NumBases > 0, 'CallAmbigBaseIfNeeded function called with no bases.'
  if len(bases) == 1:
    BaseHere = bases
  else:
    # If a gap is one of the things most common at this position, call an 'N';
    # otherwise, find the ambiguity code for this set of bases.
    if "-" in bases:
      BaseHere = 'N'
    else:  
      try:
        BaseHere = ReverseIUPACdict2[bases]
      except KeyError:
        print('Unexpected set of bases', bases, 'found in', BaseFreqFile, \
        ', not found amonst those for which we have ambiguity codes, namely:', \
        ' '.join(list(ReverseIUPACdict2.keys())) + '. Quitting.', file=sys.stderr)
        raise
  if coverage < MinCovForUpper - 0.5:
    return BaseHere.lower()
  else:
    return BaseHere.upper()


def PropagateNoCoverageChar(seq, LeftToRightDone=False):
  '''Replaces gaps that border "no coverage" by "no coverage".

  Where NoCoverageChars neighbour GapChars, propagate the former outwards until
  they touch bases on both sides (because deletions should only be called when
  the bases on either side are known). e.g.
  ACTG---?---ACTG
  becomes
  ACTG???????ACTG'''
  
  if LeftToRightDone:
    seq = seq[::-1]
  BaseToLeftIsNoCoverage = False
  ResultingSeq = ''
  for base in seq:
    if base == '?':
      BaseToLeftIsNoCoverage = True
      ResultingSeq += '?'
    elif base == '-':
      if BaseToLeftIsNoCoverage:
        ResultingSeq += '?'
      else:
        ResultingSeq += '-'
    else:
      BaseToLeftIsNoCoverage = False
      ResultingSeq += base
  if LeftToRightDone:
    ResultingSeq = ResultingSeq[::-1]
  else:
    ResultingSeq = PropagateNoCoverageChar(ResultingSeq, True)
  return ResultingSeq


def ReadSequencesFromFile(DataFile,IsAlignment=True):
  '''Reads in all sequences from a file into a dictionary.'''

  # Check that the first argument exists and is a file
  if not os.path.isfile(DataFile):
    print(DataFile, 'does not exist or is not a file.', file=sys.stderr)
    exit(1)

  # Check that the second argument is a bool
  if type(IsAlignment) != type(True):
    print('Function ReadSequencesFromFile called with a second argument that',\
    'is not a bool.\nQuitting.', file=sys.stderr)
    exit(1)

  # Read in all sequences
  AllSequences = {}
  HaveReachedFirstSequence = False
  with open(DataFile, 'r') as f:
    for line in f:

      # Strip whitespace
      ThisLine = line.strip()

      # Ignore blank lines
      if ThisLine == '':
        continue

      # If we're at the start of a new sequence, check its name hasn't already
      # been encountered, then start an empty string to which we'll append the
      # sequence contained in subsequent lines of the file.
      if ThisLine[0] == '>':
        NameOfCurrentSequence = ThisLine[1:].split()[0]
        if NameOfCurrentSequence in AllSequences:
          print('Found a second sequence titled', NameOfCurrentSequence+\
          '; sequence names should be unique.\nQuitting.', file=sys.stderr)
          exit(1)
        AllSequences[NameOfCurrentSequence] = ''
        HaveReachedFirstSequence = True
        continue

      # If we're not at the start of a new sequence and we haven't read any
      # sequences yet, there's nothing to do (i.e. we ignore everything that comes
      # before the name of the first sequence).
      elif not HaveReachedFirstSequence:
        continue

      # If we're here, we should read in sequence data.
      AllSequences[NameOfCurrentSequence] += ThisLine


  # Check we have at least one sequence
  if len(AllSequences) == 0:
    print('No sequences found in', DataFile+'.\nQuitting.', file=sys.stderr)
    exit(1)


  # Check all sequences have the same length, if they're supposed to
  FirstSequenceName, FirstSequence = list(AllSequences.items())[0]
  SequenceLength = len(FirstSequence)
  if IsAlignment:
    for SequenceName, Sequence in list(AllSequences.items()):
      if len(Sequence) != SequenceLength:
        print(SequenceName, 'has length', len(Sequence), 'whereas', \
        FirstSequenceName, 'has length', str(SequenceLength)+\
        '. Aligned sequences were expected.\nQuitting.', file=sys.stderr)
        exit(1)

  return AllSequences, SequenceLength




def ReadSequencesFromFile_ordered(DataFile,IsAlignment=True):
  '''Reads in all sequences from a file into a list of items [name,sequence].'''

  # Check that the first argument exists and is a file
  if not os.path.isfile(DataFile):
    print(DataFile, 'does not exist or is not a file.', file=sys.stderr)
    exit(1)

  # Check that the second argument is a bool
  if type(IsAlignment) != type(True):
    print('Function ReadSequencesFromFile called with a second argument that',\
    'is not a bool.\nQuitting.', file=sys.stderr)
    exit(1)

  # Read in all sequences
  AllSequences = []
  HaveReachedFirstSequence = False
  with open(DataFile, 'r') as f:
    for line in f:

      # Strip whitespace
      ThisLine = line.strip()

      # Ignore blank lines
      if ThisLine == '':
        continue

      # If we're at the start of a new sequence, check its name hasn't already
      # been encountered, then start an empty string to which we'll append the
      # sequence contained in subsequent lines of the file.
      if ThisLine[0] == '>':
        NameOfCurrentSequence = ThisLine[1:].split()[0]
        AllSequences.append([NameOfCurrentSequence,''])
        HaveReachedFirstSequence = True
        continue

      # If we're not at the start of a new sequence and we haven't read any
      # sequences yet, there's nothing to do (i.e. we ignore everything that comes
      # before the name of the first sequence).
      elif not HaveReachedFirstSequence:
        continue

      # If we're here, we should read in sequence data.
      AllSequences[-1][1] += ThisLine


  # Check we have at least one sequence
  if len(AllSequences) == 0:
    print('No sequences found in', DataFile+'.\nQuitting.', file=sys.stderr)
    exit(1)

  # Check all sequences have the same length, if they're supposed to
  FirstSeqName =AllSequences[0][0]
  FirstSeqLength = len(AllSequences[0][1])
  if IsAlignment:
    OtherSeqs = [item[1] for item in AllSequences[1:]]
    if not all(len(OtherSeq) == FirstSeqLength for OtherSeq in OtherSeqs):
      for [SeqName,seq] in AllSequences[1:]:
        if len(item[1]) != FirstSeqLength:
          print(SeqName, 'has length', len(seq), 'whereas', \
          FirstSeqName, 'has length', FirstSeqLength, \
      'Aligned sequences were expected.\nQuitting.', file=sys.stderr)
      exit(1)

  return AllSequences, FirstSeqLength





def ReadPatientFile(OneLinePerPatientOnly, filename):
  '''Read in patient data from a csv file.

  If the bool OneLinePerPatientOnly is True, and the same patient is found on
  two or more lines in this file, we throw an error and quit.
  '''

  AllPatientsDict = {}
  with open(filename, 'r') as f:
    CurrentLineNumber=0
    for line in f:
      CurrentLineNumber += 1

      if CurrentLineNumber == 1:

        # Let the file itself choose the names of the fields. Strip whitespace.
        fields = [field.strip() for field in line.split(',')]
        NumFields = len(fields)

        # Complain if any field except the first one is called 'ID'.
        if 'ID' in fields[1:]:
          print('Column', fields[1:].index('ID')+2, 'of',filename,\
          'is titled "ID" - we reserve this title for the first column (the',\
          'patient ID). Exiting.', file=sys.stderr)
          exit(1)
        fields[0] = 'ID'

        # Check all field names are unique
        CounterObject = collections.Counter(fields)
        DuplicatedFieldNames = [i for i in CounterObject if CounterObject[i]>1]
        if len(DuplicatedFieldNames) != 0:
          for DuplicatedFieldName in DuplicatedFieldNames:
            print('The field name "'+ DuplicatedFieldName+\
            '" is duplicated in file', filename, file=sys.stderr)
          print('Each field should have a unique name. Exiting.', \
          file=sys.stderr)
          exit(1)

      else:
        # After reading the first line (with the field names) we hit data (csv).
        data = line.split(',')

        # Strip leading and trailing whitespace from every field
        for i in range(0,len(data)):
          data[i] = data[i].strip()

        # Ignore empty lines
        if all(datum == '' for datum in data):
          continue

        # Trim single or double quotes from around the patient ID
        ID = data[0].strip('''"' ''')

        # Check for duplicated patients
        if OneLinePerPatientOnly and ID in AllPatientsDict:
          print('ERROR: on line number', CurrentLineNumber, 'of', filename+ \
          ', patient', ID, 'was encountered a second time. Quitting.', \
          file=sys.stderr)
          exit(1)

        # Check this line has the expected number of items of data
        if len(data) != NumFields:
          print('WARNING: Expected',NumFields, 'fields; encountered',len(data),\
          "in the following line, which we're therefore skipping:")
          print(line)
          continue

        # If each patient might have multiple lines of data in this file,
        # create a list for each field and append to it. If this line
        # corresponds to a patient we haven't seen yet, the list needs 
        # initialising.
        if not OneLinePerPatientOnly:
          if not ID in AllPatientsDict:
            AllPatientsDict[ID] = {}
            for i,field in enumerate(fields):
              if field == 'ID':
                continue
              AllPatientsDict[ID][field] = [data[i].strip()]
          else:
            # We've seen this patient in this file already, so append to the
            # existing list instead of creating it.
            for i,field in enumerate(fields):
              if field == 'ID':
                continue
              AllPatientsDict[ID][field].append(data[i].strip())

        else:
          # OneLinePerPatientOnly is True
          ThisPatient = {}
          for i,field in enumerate(fields):
            if field == 'ID':
              continue
            ThisPatient[field] = data[i].strip()
          AllPatientsDict[ID] = ThisPatient

  if len(AllPatientsDict) == 0:
    print('No patients found in file', filename+'.\nQuitting.', file=sys.stderr)
    exit(1)

  return AllPatientsDict

