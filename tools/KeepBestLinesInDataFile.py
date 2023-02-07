#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script focusses on two fields in an input data
file: an ID field and a field used for sorting (containing a float). For any ID
that appears in multiple rows, we keep only one of these rows: by default the
one with the smallest value of the sorting field. (In the event of an exact tie
we use the row encountered first.) The specific purpose we have in mind is
keeping the blast hit with the smallest evalue.
This script emulates terminal commands like "sort MyData.csv -t, -k1,1 -k2,2 |
sort -t, -k1,1 -u --merge" (here keeping one occurence of each unique value in
field one, based on the value of field 2), which in my experience work on linux
but fail on MacOS (see https://www.biostars.org/p/144569/#325003).'''

import argparse
import os
import sys
from collections import OrderedDict

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

# Set up the arguments for this script
parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('in_file', type=File)
parser.add_argument('out_file')
parser.add_argument('--separator', default=',', help='''Used to specify
the character that separates fields in the data file. (By default this is a
comma.)''')
parser.add_argument('-N', '--num_fields', type=int, default=9, help='''Used to
specify the expected number of fields in the data file. (By default this is
%(default)s.)''')
parser.add_argument('-I', '--id_field', type=int, default=1, help='''Used to
specify which field in the file is the ID. (By default this is 1. i.e. the first
field.)''')
parser.add_argument('-S', '--sort_field', type=int, default=4, help='''Used
to specify the field in the data file whose value we will sort by. (By default
this is %(default)s.)''')
parser.add_argument('-L', '--keep_largest', action='store_true', help='''Used
to specify that we should keep the largest value in the sort field for each id.
(By default we keep the smallest.)''')
parser.add_argument('-O', '--order_by_id', action='store_true', help='''Used to
specify that the output file should be sorted by id, rather than left in the
original order.''')
parser.add_argument('-H', '--header', action='store_true', help='''Used to
specify that the first line is a header line (which we therefore always include
in the output, excluding it from the sort-based choice of lines). By default we
assume there is no header line (as is the case in blast output).''')
args = parser.parse_args()

if args.sort_field > args.num_fields:
  print("Error: sort_field cannot be larger than num_fields. Quitting.",
  file=sys.stderr) 
  exit(1)
if args.id_field > args.num_fields:
  print("Error: id_field cannot be larger than num_fields. Quitting.",
  file=sys.stderr) 
  exit(1)

rows_to_keep = OrderedDict()

with open(args.in_file, 'r') as f:
  for lin_num_min_1, line in enumerate(f):

    if lin_num_min_1 == 0 and args.header:
      header = line
      continue

    # Split into fields.
    fields = line.split(args.separator)
    assert len(fields) == args.num_fields, \
    "Error: expected {n} fields in data file; found {m}.".format(\
    n=args.num_fields, m=len(fields))

    # Get the sort value as a float.
    try:
      sort_value = float(fields[args.sort_field - 1])
    except ValueError:
      print("Error: could not understand value", fields[args.sort_field - 1],
      "in field", args.sort_field, "on line", lin_num_min_1 + 1, "in",
      args.in_file, "as a float. Quitting.", file=sys.stderr) 
      exit(1)

    id_ = fields[args.id_field - 1]

    # Record this line unless we've seen it before with a preferable sort value.
    if id_ in rows_to_keep:
      prev_sort_value = rows_to_keep[id_][1]
      if args.keep_largest:
        if sort_value > prev_sort_value:
          rows_to_keep[id_] = (line, sort_value)
      else:
        if sort_value < prev_sort_value:
          rows_to_keep[id_] = (line, sort_value)
    else:
      rows_to_keep[id_] = (line, sort_value)

with open(args.out_file, "w") as f:
  if args.header:
    f.write(header)
  if args.order_by_id:
    for value in sorted(list(rows_to_keep.values()),
    key=lambda x: x[0].split(",", 1)[0]):
      f.write(value[0])
  else:
    for value in list(rows_to_keep.values()):
      f.write(value[0])

