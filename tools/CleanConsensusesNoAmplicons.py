#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script 'cleans' a sequence alignment. The "?"
character is replaced by "N". Multiple sequences for the same patient & date
(determined by assuming sequence IDs start by the patient ID, then "_", then the
date, then optionally "_" and anything else) are flattened into one sequence by
using the base of the longest sequence (the most bases excluding N or gaps) if
is not an N, otherwise the base of the second-longest sequence if it is not an
N, etc. Then any gap character neighbouring an N is iteratively replaced by an
N. Any wholly undetermined (purely N) sequences are removed.'''

import argparse
import os
import sys
from Bio import AlignIO
from Bio import Seq  
from Bio import SeqIO  
import collections
import pandas
import datetime
from re import sub
from ShiverFuncs import TranslateSeqCoordsToAlnCoords

# Define a function to check files exist, as a type for the argparse.
def File(_file):
  if not os.path.isfile(_file):
    raise argparse.ArgumentTypeError(_file +' does not exist or is not a file.')
  return _file

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('global_aln', type=File, help='''An alignment of all
consensuses to be cleaned. It should contain no ambiguity codes; Ambiguity codes
can be estimated using shiver/tools/EstimateAmbiguousBases.py. If sequence names
contain the string "_consensus" and/or "_MinCov", the name is truncated by
removing the string and every thing after it. After that, each sequence name
should be either the ID of a patient, or the ID of a patient followed by an
underscore then a unique string, to distinguish multiple sequences associated
with the same patient.''')
parser.add_argument('seq_based_blacklist', type=File, help='''A plain-text file
listing sequences to be blacklisted.''')
parser.add_argument('output_file')
parser.add_argument('--date_sampled_csv', type=File, help='''A csv-format file
containing a column "PATIENT" and a column "Date.sampled". For patients with
sequences from multiple dates, we'll choose the date closest to
"Date.sampled". If this file is not provided, we include all dates.''')
parser.add_argument('-V', '--verbose', action='store_true')
args = parser.parse_args()

have_dates = args.date_sampled_csv != None

def preprocess_seq(seq):
  '''(1) replaces '?' by N; (2) converts to upper case; (3) checks for ambiguous
  bases; (4) checks whether it's wholly undetermined.'''
  seq_as_str = str(seq.seq)
  initial_length = len(seq_as_str)
  seq_as_str = sub("\?", "N", seq_as_str)
  assert len(seq_as_str) == initial_length
  seq_as_str = seq_as_str.upper()
  if any(not base in "ACGTN-" for base in seq_as_str):
    unexpected_bases = set(base for base in seq_as_str if not base in "ACGTN-")
    print('Seq', seq.id, 'contains a base other than A, C, G, T, N or -',
    '(specifically: ', ', '.join(unexpected_bases) + '). This is unexpected.',
    'Have you run shiver/tools/EstimateAmbiguousBases.py on your',
    'alignment first? Quitting.', file=sys.stderr)
    exit(1)
  wholly_undetermined = all(base == "N" or base == "-" for base in seq_as_str)
  return seq_as_str, wholly_undetermined

def preprocess_seq_id(seq_id):
  "Removes things we don't want from the end of seq ids."
  seq_id = seq_id.split('_consensus', 1)[0]
  #seq_id = seq_id.split('_MinCov', 1)[0]
  #seq_id = seq_id.split('_remap', 1)[0]
  #seq_id = seq_id.split('_BaseFreqs_', 1)[0]
  #seq_id = seq_id.split('.csv', 1)[0]
  return seq_id

def get_beehive_id(seq_id):
  "Return whatever's before the first underscore if there is one."
  return seq_id.split("_", 1)[0]

def PropagateNoCoverageChar(seq, LeftToRightDone=False):
  '''Iteratively replace gaps that border N by N.

  Where N neighbours a gap, propagate N outwards until it touches bases on both
  sides (because deletions should only be called when the bases on either side
  are known). e.g.
  ACTG---N---ACTG
  becomes
  ACTGNNNNNNNACTG'''
  
  assert type(seq) == type("a string")
  if LeftToRightDone:
    seq = seq[::-1]
  BaseToLeftIsNoCoverage = False
  ResultingSeq = ''
  for base in seq:
    if base == 'N':
      BaseToLeftIsNoCoverage = True
      ResultingSeq += 'N'
    elif base == '-':
      if BaseToLeftIsNoCoverage:
        ResultingSeq += 'N'
      else:
        ResultingSeq += '-'
    else:
      BaseToLeftIsNoCoverage = False
      ResultingSeq += base
  if LeftToRightDone:
    ResultingSeq = ResultingSeq[::-1]
  else:
    ResultingSeq = PropagateNoCoverageChar(ResultingSeq, True)
  assert not "N-" in ResultingSeq and not "-N" in ResultingSeq
  return ResultingSeq

def num_known_bases(seq):
  return len(seq) - seq.count("N") - seq.count("-")


# Read the alignment
try:
  alignment = AlignIO.read(args.global_aln, "fasta")
except:
  print('Problem reading', args.global_aln + ':', file=sys.stderr)
  raise
alignment_length = alignment.get_alignment_length()

# Read the seq-based blacklist.
with open(args.seq_based_blacklist, 'r') as f:
  seq_blacklist = set(line.strip() for line in f)

# Iterate through all seqs.
seq_dict = collections.defaultdict(dict)
blacklisted_seqs_found = set([])
for seq in alignment:

  seq.id = preprocess_seq_id(seq.id)

  # Check expected seq id format
  if seq.id.count("_") == 0:
    print('Error: expected seq ids to be the patient id, then "_", then the',
    'date, then optionally another "_" plus another string; seq id', seq.id,
    'does not have this format. Quitting.', file=sys.stderr)
    exit(1)

  # Skip seqs that are blacklisted or wholly undetermined. 
  if seq.id in seq_blacklist:
    blacklisted_seqs_found.add(seq.id)
    if args.verbose:
      print("Skipping sequence", seq.id, "which was blacklisted.")
    continue
  seq_as_str, wholly_undetermined = preprocess_seq(seq)
  if wholly_undetermined:
    if args.verbose:
      print("Skipping sequence", seq.id, "which is wholly undetermined.")
    continue

  # Store all seqs for the same patient in one dict within the larger dict for
  # all patients. Check seq ids are unique.
  beehive_id = seq.id.split("_", 1)[0]
  if seq.id in seq_dict[beehive_id]:
    print('Error: encountered seq id', seq.id, "a second time. Quitting.",
    file=sys.stderr)
    exit(1)
  seq_dict[beehive_id][seq.id] = seq_as_str

# Warn about blacklisted seqs that were not found.
blacklisted_seqs_not_found = seq_blacklist - blacklisted_seqs_found
if len(blacklisted_seqs_not_found) != 0:
  print('Warning: of', len(seq_blacklist), 'blacklisted seqs specified in',
  args.seq_based_blacklist + ', the following', len(blacklisted_seqs_not_found),
  'were not encountered:', ' '.join(blacklisted_seqs_not_found),
  file=sys.stderr)

if have_dates:
  date_sampled_df = pandas.read_csv(args.date_sampled_csv, na_filter=False)
  try:
    pats_in_date_sampled_df = date_sampled_df["PATIENT"]
  except KeyError:
    print('Error: no column called "PATIENT" found in', args.date_sampled_csv + \
    ". Quitting.", file=sys.stderr)
    exit(1)
  try:
    dates_in_date_sampled_df = date_sampled_df["Date.sampled"]
  except KeyError:
    print('Error: no column called "Date.sampled" found in',
    args.date_sampled_csv + ". Quitting.", file=sys.stderr)
    exit(1)
  #pats_in_date_sampled_df = list(pats_in_date_sampled_df)
  #dates_in_date_sampled_df = list(dates_in_date_sampled_df)

# Iterate through patients and the set of seqs for each patient.
final_seqs_dict = {}
for beehive_id, sub_seq_dict in seq_dict.items():

  # If this patient only has one sequence, nothing to do except restrict the 
  # seq id to the patient id and the date.
  if len(sub_seq_dict) == 1:
    seq_id, seq_as_str = list(sub_seq_dict.items())[0]
    seq_id_pieces = seq_id.split("_")
    seq_id_to_use = beehive_id + "_" + seq_id_pieces[1]
    assert not seq_id_to_use in final_seqs_dict
    final_seqs_dict[seq_id_to_use] = seq_as_str
    continue

  # Now we have multiple sequences for this patient. First, group them by date.
  seqs_by_date = collections.defaultdict(list)
  for seq_id, seq_as_str in sub_seq_dict.items():
    date = seq_id.split("_", 2)[1]
    seqs_by_date[date].append(seq_as_str)

  # Merge multiple seqs from the same date: rank them by length, and at each
  # position use the base of the best seq that doesn't have an N there.
  merged_seqs_dict = {}
  for date, seqs in seqs_by_date.items():
    num_seqs = len(seqs)
    if num_seqs == 1:
      best_seq = seqs[0]
    else:
      if args.verbose:
        print("Merging", num_seqs, "seqs for", beehive_id, "for date", date)
      seqs_by_length = sorted(seqs, key=lambda seq : num_known_bases(seq),
      reverse=True)
      best_seq = ""
      for pos in range(alignment_length):
        best_base = "N"
        for i in range(num_seqs):
          base = seqs_by_length[i][pos]
          if base != "N":
            best_base = base
            break
        best_seq += best_base
    merged_seqs_dict[date] = best_seq

  # If we're including all dates for each pat, do so now...
  if not have_dates:
    for date, seq in list(merged_seqs_dict.items()):
      seq_id_to_use = beehive_id + "_" + date
      assert not seq_id_to_use in final_seqs_dict
      final_seqs_dict[seq_id_to_use] = seq
    continue
  # ... otherwise we go on to choose one date per pat.

  # If there is only one date, no choice:
  if len(merged_seqs_dict) == 1:
    date, seq = list(merged_seqs_dict.items())[0]
    seq_id_to_use = beehive_id + "_" + date
    assert not seq_id_to_use in final_seqs_dict
    final_seqs_dict[seq_id_to_use] = seq
    continue

  # Now we have multiple dates. If NA is one of them, remove it and check again
  # whether there's only one left.
  try:
    del merged_seqs_dict["NA"]
  except KeyError:
    pass
  else:
    if len(merged_seqs_dict) == 1:
      date, seq = list(merged_seqs_dict.items())[0]
      seq_id_to_use = beehive_id + "_" + date
      assert not seq_id_to_use in final_seqs_dict
      final_seqs_dict[seq_id_to_use] = seq
      if args.verbose:
        print("For", beehive_id, "using date", date,
        "as the only other date is NA")
      continue

  # Now we have multiple non-NA dates. Must choose one.

  # Get the date sampled from the csv.
  which_rows_this_pat = [row_num for row_num, pat in \
  enumerate(pats_in_date_sampled_df) if pat == beehive_id]
  if len(which_rows_this_pat) != 1:
    print('Error: expected exactly 1 row in', args.date_sampled_csv,
    "with a PATIENT value of", beehive_id + "; found " + \
    str(len(which_rows_this_pat)) + ". Quitting.", file=sys.stderr)
    exit(1)
  this_pat_row_num = which_rows_this_pat[0]
  date_sampled = dates_in_date_sampled_df[this_pat_row_num]

  # We may or may not use this, but try to convert it to a date object. If we
  # can't, warn and skip the patient.
  try:
    date_sampled_date_object = datetime.datetime.strptime(date_sampled, '%Y-%m-%d')
  except ValueError:
    print('Warning: could not understand the Date.sampled - "' + date_sampled\
    + '" - for patient', beehive_id, "in", args.date_sampled_csv + \
    ". Skipping patient.", file=sys.stderr)
    continue

  # If the date sample is a perfect match to one of the seq dates, use that one;
  # otherwise use the seq date closest to the date sampled. 
  if date_sampled in merged_seqs_dict:
    date_to_use = date_sampled
    if args.verbose:
      print("For", beehive_id, "with seq dates",
      " and ".join(list(merged_seqs_dict.keys())) + ", using date", date_sampled,
      "as this perfectly matches the Date.sampled in", args.date_sampled_csv)
  else:
    date_to_use = None
    for seq_date in list(merged_seqs_dict.keys()):
      if date_to_use == None:
        seq_date_object = datetime.datetime.strptime(seq_date, '%Y-%m-%d')
        min_date_diff = abs((date_sampled_date_object - seq_date_object).days)
        date_to_use = date
      else:
        seq_date_object = datetime.datetime.strptime(seq_date, '%Y-%m-%d')
        date_diff = abs((date_sampled_date_object - seq_date_object).days)
        if date_diff < min_date_diff:
          min_date_diff = date_diff
          date_to_use = date
    if args.verbose:
      print("For", beehive_id, "with seq dates",
      " and ".join(list(merged_seqs_dict.keys())) + ", using date", date_to_use,
      "as this is the closest match to the Date.sampled -", date_sampled,
      "- in", args.date_sampled_csv)

  # Record the seq to use.
  seq_to_use = merged_seqs_dict[date_to_use]
  seq_id_to_use = beehive_id + "_" + date_to_use
  assert not seq_id_to_use in final_seqs_dict
  final_seqs_dict[seq_id_to_use] = seq_to_use

# Iteratively replace gaps that border N by N.
for seq_id, seq in final_seqs_dict.items():
  final_seqs_dict[seq_id] = PropagateNoCoverageChar(seq)

# Output the final seqs, sorted by name.
sorted_seq_objects = [SeqIO.SeqRecord(Seq.Seq(seq), id=seq_id, description='') \
for seq_id, seq in sorted(list(final_seqs_dict.items()), key=lambda x:x[0])]
SeqIO.write(sorted_seq_objects, args.output_file, 'fasta')

