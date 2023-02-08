#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''Takes in a pairwise alignment of seqs, the first one
(nominally the consensus) possibly containing the missing coverage character '?'
and ambiguity codes, and one more seq (nominally a reference used for mapping to
derive the consensus). This script uses the reference to fill in gaps in the
consensus: where the consensus has '?' or 'N' we use the reference base; where
the consensus has an ambiguity code we take the first (in alphabetical order) of
the bases A, C, G or T that that code represents; otherwise, we use the base (or
gap) of the consensus. Output is printed to stdout suitable for redirection to a
fasta file.''' 

import argparse
import os
import sys
from Bio import SeqIO
from Bio import Seq
import itertools
from AuxiliaryFunctions import PropagateNoCoverageChar, IUPACdict

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('ConsensusWithRef', type=File)
parser.add_argument('-N', '--output-seq-name')
args = parser.parse_args()

# Find the consensus and its ref
ConsensusFound = False
RefFound = False
for seq in SeqIO.parse(open(args.ConsensusWithRef),'fasta'):
  if not ConsensusFound:
    consensus = seq
    ConsensusFound = True
    continue
  if not RefFound:
    ref = seq
    RefFound = True
    continue
  print('Found three sequences in', args.ConsensusWithRef+\
  '; expected only two. Quitting.', file=sys.stderr)
  exit(1)
if not RefFound:
  print('Less than two sequences found in', args.ConsensusWithRef+\
  '; expected two. Quitting.', file=sys.stderr)
  exit(1)

# Check the consensus and its ref are aligned with no pure-gap columns.
ConsensusAsString = str(consensus.seq).upper()
RefAsString = str(ref.seq).upper()
if len(ConsensusAsString) != len(RefAsString):
  print(args.ConsensusWithRef, 'is not an alignment - seq lengths', \
  'differ. Quitting.', file=sys.stderr)
  exit(1)

# The main bit.
NewConsensus = ''
ExpectedBases = ['A', 'C', 'G', 'T', '-']
for ConsensusBase, RefBase in zip(ConsensusAsString, RefAsString):
  if ConsensusBase == '?' or ConsensusBase == 'N':
    NewConsensus += RefBase
  elif ConsensusBase in ExpectedBases:
    NewConsensus += ConsensusBase
  elif ConsensusBase in IUPACdict:
    NewConsensus += IUPACdict[ConsensusBase][0]
  else:
    print('Encountered unexpected base', ConsensusBase, 'in', \
    args.ConsensusWithRef + '. Quitting.', file=sys.stderr)
    exit(1)

# Ungap, and rename if desired.
consensus.seq = Seq.Seq(NewConsensus).ungap('-')
if args.output_seq_name != None:
  consensus.id = args.output_seq_name
consensus.description = ''

SeqIO.write(consensus, sys.stdout, "fasta")

