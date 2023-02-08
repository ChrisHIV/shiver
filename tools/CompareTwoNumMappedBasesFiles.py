#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Finds the percentage increase, from file1 to file 2, in
corresponding numbers of mapped bases (as a function of identity threshold) in
files produced by
$ ./SummariseBam.py bam1 -R ref1.fasta -I foo,bar > file1
$ ./SummariseBam.py bam2 -R ref2.fasta -I foo,bar > file2
Output is printed to stdout suitable for redirection to a csv file.
'''

import os
import collections
import sys
import argparse
import numpy

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('file1', type=File)
parser.add_argument('file2', type=File)

args = parser.parse_args()

def ReadFromFile(MyFile):
  with open(MyFile, 'r') as f:
    for LineNumberMin1, line in enumerate(f):
      values = line.split(',')
      if LineNumberMin1 == 0:
        AllValues = [values]
      else:
        values[1] = int(values[1])
        AllValues.append(values)
  return AllValues

file1values = ReadFromFile(args.file1)
file2values = ReadFromFile(args.file2)

if file1values[0] != file2values[0]:
  print('Differing header lines found in', args.file1, 'and', args.file2 + \
  '. Quitting.', file=sys.stderr)
  exit(1)

if [value[0] for value in file1values] != [value[0] for value in file2values]:
  print('Differing first columns in', args.file1, 'and', args.file2 + \
  '. Quitting.', file=sys.stderr)
  exit(1)

print(','.join(file1values[0]).rstrip())
for i in range(1,len(file1values)):
  IdentityThreshold, count1 = file1values[i]
  IdentityThreshold, count2 = file2values[i]
  if count1 == 0:
    if count2 == 0:
      PercentageIncrease = 'NA'
    else:
      PercentageIncrease = 'inf'
  else:
    PercentageIncrease = (count2 - count1) * 100. / count1
  print(IdentityThreshold, PercentageIncrease, sep=',')


