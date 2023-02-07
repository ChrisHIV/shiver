#!/usr/bin/env python

#
## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview: this script reads in one blast file of hits for forward reads and
## one for backward reads (from paired-read data). It finds read pairs for which
## either a) both reads blast to something other than the named sequence, or b)
## one read blasts to the named sequence and its mate blasts to something else
## but that second blast has a better evalue. The names of the reads in these
## pairs are written, without preserving their original order, to two files: one
## for forward reads and one for reverse reads.
## The intended use is for processing the result of blasting paired-read data
## to a database consisting of something that looks like your sample plus other
## sequences that should be recognised as contamination. (For example, if de
## novo assembly has been done with these reads, some of the resulting contigs
## may be identified as contamination -- we want to find the reads that
## correspond to those contigs, in order to remove them.)


################################################################################
# USER INPUT
# What columns are the fields we want?
column_qacc = 1
column_sacc = 2
column_evalue = 4
################################################################################

# Import some modules we'll need.
import os.path
import sys

# Check that this script was called from the command line with two arguments.
if len(sys.argv) != 5:
  sys.stderr.write('Incorrect number of arguments given. Correct usage:\n'+\
  sys.argv[0] +' BlastFileFor1reads BlastFileFor2reads RefName '+\
  'OutFileBasename\nQuitting\n')
  exit(1)
BlastFileFor1reads = sys.argv[1]
BlastFileFor2reads = sys.argv[2]
RefName = sys.argv[3]
OutFileBasename = sys.argv[4]

# Check that the arguments exist and are files
for DataFile in [BlastFileFor1reads,BlastFileFor2reads]:
  if not os.path.isfile(DataFile):
    sys.stderr.write(DataFile +' does not exist or is not a file. Quitting.\n')
    exit(1)

RightmostColumn = max([column_qacc,column_sacc,column_evalue])

def ReadBlastFile(BlastFile):
  '''Reads in the qacc, sacc, and evalue from a blast output file.
  The qacc is expected to end in /1 or /2, and this is removed.'''
  DictOfHits = {}
  with open(BlastFile, 'r') as f:
    for line in f:
      fields = line.split(',')
      try:
        ReadName = fields[column_qacc-1]
        DesiredFields = [fields[column_sacc-1],fields[column_evalue-1]]
      except IndexError:
        sys.stderr.write('The following line in ' +BlastFile +' has too few '+\
        'fields (less than ' +str(RightmostColumn)+'):\n' +line+'Quitting.\n')
        exit(1)
      if len(ReadName) < 2 or (not ReadName[-2:] in ['/1','/2']):
        sys.stderr.write('Unexpected format for read ' +ReadName +'; read '+\
        'names are expected to end in /1 or /2.\nQuitting.\n')
        exit(1)
      ReadName = ReadName[:-2]
      if ReadName in DictOfHits:
        sys.stderr.write('Encountered read ' +ReadName +' (having trimmed '+\
        '"/1" or "/2" from the name) a second time in ' +BlastFile +\
        '.\nQuitting.\n')
        exit(1)
      DictOfHits[ReadName] = DesiredFields
  return DictOfHits

# Read in the blast files.
HitsFor1reads = ReadBlastFile(BlastFileFor1reads)
HitsFor2reads = ReadBlastFile(BlastFileFor2reads)

# Find contaminant read pairs
ContaminantReadPairs = []
for ReadName, [read1Hit, read1Evalue] in list(HitsFor1reads.items()):

  try:
    [read2Hit,read2Evalue] = HitsFor2reads[ReadName]
  except KeyError:
    continue
  read1HitsRef = read1Hit == RefName
  read2HitsRef = read2Hit == RefName
  if read1HitsRef and read2HitsRef:
    continue
  elif read1HitsRef or read2HitsRef:
    try:
      read1Evalue = float(read1Evalue)
      read2Evalue = float(read2Evalue)
    except ValueError:
      sys.stderr.write('Failed to understand one of the evalues for read ' +\
      ReadName +', namely ' +str(read1Evalue).strip() +' and ' +\
      str(read2Evalue).strip()+', as a float.\nQuitting.\n')
      exit(1)
    if read1HitsRef:
      RefEvalue         = read1Evalue
      ContaminantEvalue = read2Evalue
    else:
      ContaminantEvalue = read1Evalue
      RefEvalue         = read2Evalue    
    if ContaminantEvalue < RefEvalue:
      ContaminantReadPairs.append(ReadName)
  else:
    ContaminantReadPairs.append(ReadName)

# Write the output
OutFile_1reads = OutFileBasename+'_1.txt'
OutFile_2reads = OutFileBasename+'_2.txt'
HaveSomeReads = len(ContaminantReadPairs) > 0
with open(OutFile_1reads, 'w') as f:
  if HaveSomeReads:
    f.write('\n'.join([ReadName+'/1' for ReadName in ContaminantReadPairs]) + "\n")
with open(OutFile_2reads, 'w') as f:
  if HaveSomeReads:
    f.write('\n'.join([ReadName+'/2' for ReadName in ContaminantReadPairs]) + "\n")

# Antiquated code, from when non-contaminant read pairs were printed:
'''
## NB we discard pairs where one read blasted to the reference and the other
## did not blast to anything. Our reasons are two-fold. 1) Keeping such 'pairs'
## would involve making the assumption that the read missing from the blast 
## table really exists; it might not. 2) Though we do keep reads where
## one read hits the reference and the other hits a contaminant (if the
## contaminant hit is not as good as the reference hit), it's not inconsistent
## to discard pairs where one hits the reference and the other doesn't blast at
## all. In the former case, the contaminant hit could be masking a hit to the
## reference that's almost as good; in the latter case, one can be sure that the
## missing read (if it exists) failed to blast to the reference at all.

# Find the read pairs that blast best to the reference.
GoodReadPairs = []
for ReadName, [read1Hit, read1Evalue] in HitsFor1reads.items():

  try:
    [read2Hit,read2Evalue] = HitsFor2reads[ReadName]
  except KeyError:
    if read1Hit == RefName:
      GoodReadPairs.append(ReadName)
    continue
  read1HitsRef = read1Hit == RefName
  read2HitsRef = read2Hit == RefName
  if read1HitsRef and read2HitsRef:
    GoodReadPairs.append(ReadName)
  elif read1HitsRef or read2HitsRef:
    try:
      read1Evalue = float(read1Evalue)
      read2Evalue = float(read2Evalue)
    except ValueError:
      sys.stderr.write('Failed to understand one of the evalues for read ' +\
      ReadName +', namely ' +str(read1Evalue).strip() +' and ' +\
      str(read2Evalue).strip()+', as a float.\nQuitting.\n')
      exit(1)
    if read1HitsRef:
      RefEvalue         = read1Evalue
      ContaminantEvalue = read2Evalue
    else:
      ContaminantEvalue = read1Evalue
      RefEvalue         = read2Evalue    
    if RefEvalue <= ContaminantEvalue:
      GoodReadPairs.append(ReadName)

for ReadName, [read2Hit, read2Evalue] in HitsFor2reads.items():
  if (not ReadName in HitsFor1reads) and read2Hit == RefName:
    GoodReadPairs.append(ReadName)
'''
