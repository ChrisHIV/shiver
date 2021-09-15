from __future__ import print_function

import argparse
import os

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This trivial script replaces all whitespace (spaces,
tabs, an \\r character and/or a \\n character) at the end of each line in the
input file with a single \\n character. Output is printed to stdout suitable for
redirection to another file.'''


# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('InputFile', type=File)
args = parser.parse_args()

with open(args.InputFile, 'r') as f:
  for line in f:
    print(line.rstrip())
