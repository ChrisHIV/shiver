#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script finds the consensus sequence, i.e. the
sequence defined by taking what is most common at each position, from a set of
input sequences that are aligned in the csv format produced by
shiver/tools/MergeAlignmentsToCsv.py (i.e. each column after the first two
corresponds to one of the input sequences; each row corresponds to one base in a
reference sequence, and what each input sequence has at the position of that
base can be a base, a gap or a kmer).'''

import argparse
import os
import sys
from collections import Counter
from re import sub
from Bio import Seq  
from Bio import SeqIO  

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('alignment_csv', type=File)
parser.add_argument('--fractional-missingness-threshold', type=float, help='''
The fractional missingness (the fraction of sequences that have only "N" and/or
"-" at this position) above which we will call nothing, skipping onto the next
position. The default is 0.5.''', default=0.5)
parser.add_argument('--output_seq_name', default='GlobalConsensus',
help="The default is 'GlobalConsensus'.")
args = parser.parse_args()

consensus = ""

with open(args.alignment_csv, 'r') as f:
  for lin_num_min_1, line in enumerate(f):

    fields = line.strip().split(",")

    # Record the number of fields on line 1. Check it's >=3.
    if lin_num_min_1 == 0:
      num_fields = len(fields)
      num_seqs = num_fields - 2
      if num_fields < 3:
        print("Only", num_fields, "fields on line number", lin_num_min_1 + 1,
        "in", args.alignment_csv + "; expected at least 3. Quitting.",
        file=sys.stderr)
        exit(1)
      continue

    # Check for consistent number of fields.
    if len(fields) != num_fields:
      print("Error:", len(fields), "fields on line number", lin_num_min_1 + 1,
      "in", args.alignment_csv + ", c.f.", num_fields,  "on line 1. Quitting.",
      file=sys.stderr)
      exit(1)

    # Check the first two fields are what we expect.
    ref_pos = int(fields[0])
    ref_base = fields[1]
    assert ref_base in "ACGT", "Unexpected base" + ref_base + "on line" + \
    str(lin_num_min_1 + 1)

    # Count all of the different kmers at this position. Remove "N" and "-" from
    # every kmer, skip empty ones, merge ones that are now identical. 
    unprocessed_kmers = fields[2:]
    unprocessed_kmer_counts = Counter(unprocessed_kmers)
    kmer_counts = Counter()
    kmer_len_counts = Counter()
    for unprocessed_kmer, count in list(unprocessed_kmer_counts.items()):
      kmer = sub("[N-]+", "", unprocessed_kmer)
      len_kmer = len(kmer)
      if len_kmer > 0:
        kmer_counts[kmer] += count
      if any(base != "N" for base in unprocessed_kmer):
        kmer_len_counts[len_kmer] += count

    # Skip positions where too many kmers were wholly undetermined (nothing but
    # "N" and "-").
    num_known_kmers = sum(kmer_counts.values())
    missingness = 1 - float(num_known_kmers) / num_seqs
    if missingness > args.fractional_missingness_threshold:
      continue

    # Append the most common kmer to our consensus seq so far.
    most_common_kmer = kmer_counts.most_common(1)[0][0]
    len_of_most_common_kmer = len(most_common_kmer)
    true_most_common_len = kmer_len_counts.most_common(1)[0][0]
    if len(most_common_kmer) != true_most_common_len:
      print("Warning: at position", ref_pos, "the most common kmer is",
      most_common_kmer, "which has length", str(len_of_most_common_kmer) + \
      ", but the most common length amongst all kmers here is",
      str(true_most_common_len) + ". This is a bit irritating but an",
      "unavoidable result of simply using the most common kmer at each",
      "position. Continuing.", file=sys.stderr)
    consensus += most_common_kmer

out_seq = SeqIO.SeqRecord(seq=Seq.Seq(consensus), id=args.output_seq_name,
description='')
SeqIO.write([out_seq], sys.stdout, "fasta")
