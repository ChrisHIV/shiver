#!/usr/bin/env python


## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251
##
## Overview:
ExplanatoryMessage = '''This script merges the output from running
shiver/tools/LinkIdentityToCoverage.py on several (usually many) bam files,
producing a csv file and a plot of how mean read identity varies with coverage.
It also reports the mean number of positions per bam file that have a coverage
equalling or exceeding each coverage value considered, showing the loss in the
amount of data as one increases a coverage threshold for calling a consensus.'''

################################################################################
# USER INPUT
XaxisLabel='coverage'
YaxisLabel='mean read identity'
colour = 'blue'
################################################################################

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import argparse
import os
import collections
import sys
import numpy as np
try:
  import pandas as pd
except ImportError:
  print("This script requires your python installation to have the pandas",
  "module installed. Search for 'python pandas installation' or similar.")
  raise

# Define a function to check files exist, as a type for the argparse.
def File(MyFile):
  if not os.path.isfile(MyFile):
    raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
  return MyFile

parser = argparse.ArgumentParser(description=ExplanatoryMessage)
parser.add_argument('OutFileStem', help="We'll produce .pdf and .csv "
"extensions from this file stem")
parser.add_argument('DataFile', type=File, nargs='*', help='''One or more csv
files produced by running shiver/tools/LinkIdentityToCoverage.py on one or more
bam files.''')
parser.add_argument('--csv-for-replotting', help="""If you have already run this
script and created the csv file output, and now just want to replot the same
data with different plotting parameters, use this option to specify the existing
csv file. In this case, you shouldn't specify anything for the DataFile
argument.""")

parser.add_argument('--title', help='For the plot.',
default='''The read identity averaged over all genome positions in all samples,
grouping all positions with the same coverage value.
Point annotations show the mean number of genome positions
per sample with a coverage of at least that value.''')
parser.add_argument('-AS', '--axis-font-size', type=int,
help='For the plot. The default is 11.', default=11)
parser.add_argument('-TS', '--title-font-size', type=int,
help='For the plot. The default is 11.', default=11)
parser.add_argument('--point-size', type=int,
help='For the plot. The default is 15.', default=15)
parser.add_argument('--point-thickness', type=float,
help='For the plot. The default is 0.5.', default=0.5)
parser.add_argument('--point-marker', help='For the plot. The default is "+".',
default="+")
parser.add_argument('--point-alpha', type=float,
help='For the plot. The default is 0.6.', default=0.6)
parser.add_argument('--point-font-size', type=int,
help='For the plot. The default is 3.', default=3)
parser.add_argument('-XM', '--x-min-max', help='The minimum and maximum for '\
'the x axis in the plot, specified together as a comma-separated pair of '\
'numbers.')
parser.add_argument('-YM', '--y-min-max', help='The minimum and maximum for '\
'the y axis in the plot, specified together as a comma-separated pair of '\
'numbers.')

args = parser.parse_args()

if len(args.DataFile) == 0 and args.csv_for_replotting == None:
  print('You must either specify something for the DataFile argument or for',
  'the --csv-for-replotting argument. Quitting.', file=sys.stderr)
  exit(1)

def ParseNumberPair(arg, ArgName):
  MinMax = arg.split(',')
  if len(MinMax) != 2:
    print(ArgName, 'should be used to specify a comma-separated pair of',
    'numbers. Quitting.', file=sys.stderr)
    exit(1)
  try:
    Min = float(MinMax[0])
    Max = float(MinMax[1])
  except ValueError:
    print(ArgName, 'should be used to specify a comma-separated pair of',
    'numbers. Quitting.', file=sys.stderr)
    exit(1)
  return min(Min, Max), max(Min, Max)

# Get plot limits
if args.x_min_max:
  Xmin, Xmax = ParseNumberPair(args.x_min_max, '--x-min-max')
if args.y_min_max:
  Ymin, Ymax = ParseNumberPair(args.y_min_max, '--y-min-max')

OutFilePdf = args.OutFileStem + ".pdf"

if args.csv_for_replotting == None:

  # Check the csv doesn't exist already.
  OutFileCsv = args.OutFileStem + ".csv"
  if os.path.isfile(OutFileCsv):
    print(OutFileCsv, 'exists already; rename, move or delete it. Quitting to',
    'prevent overwriting.', file=sys.stderr)
    exit(1)

  # Record how many times we encounter each coverage
  CoverageCounts = {}
  # Keep a running count of the total read identity associated with each coverage,
  # 'scaled' in that it's missing a multiplication by coverage
  ScaledTotalIdentitiesByCoverage = {}
  for DataFile in args.DataFile:
    with open(DataFile, 'r') as f:
      for LineNumMin1, line in enumerate(f):

        if LineNumMin1 == 0:
          continue
        fields = line.split(',')
        coverage = int(fields[0])
        NumPosWithThatCov = int(fields[1])
        MeanIdentity = float(fields[2])
        # 'scaled' in that it's missing a multiplication by coverage
        ScaledTotalIdentity = NumPosWithThatCov * MeanIdentity

        if coverage in CoverageCounts:
          CoverageCounts[coverage] += NumPosWithThatCov
          ScaledTotalIdentitiesByCoverage[coverage] += ScaledTotalIdentity
        else:
          CoverageCounts[coverage] = NumPosWithThatCov
          ScaledTotalIdentitiesByCoverage[coverage] = ScaledTotalIdentity

  # Constructing vectors in a loop.
  # Firstly, all the different coverage values encountered, from largest
  # to smallest.
  # Secondly, in the same order, the mean identity for the corresponding coverage.
  # Thirdly, in the same order, the total number of positions with coverage
  # exceeding (or equalling) that coverage.
  x = np.empty(len(CoverageCounts), dtype=int)
  y = np.empty(len(CoverageCounts))
  NumbersOfPositionsExceedingCoverages = np.empty(len(CoverageCounts))
  RunningTotalNumberOfPositions = 0
  for i, (coverage, count) in enumerate(sorted(list(CoverageCounts.items()),
  key = lambda x : x[0], reverse = True)):
    x[i] = coverage
    y[i] = ScaledTotalIdentitiesByCoverage[coverage] / count
    RunningTotalNumberOfPositions += count
    NumbersOfPositionsExceedingCoverages[i] = RunningTotalNumberOfPositions

  x = np.flipud(x)
  y = np.flipud(y)
  NumbersOfPositionsExceedingCoverages = \
  list(reversed(NumbersOfPositionsExceedingCoverages))

  # Divide NumbersOfPositionsExceedingCoverages by the number of samples, to get
  # a mean value per sample.
  MeanNumPosExceedingCoveragesPerSample = [float(NumPos) / len(args.DataFile) \
  for NumPos in NumbersOfPositionsExceedingCoverages]

  # Write the csv output
  df = pd.DataFrame({'coverage' : x})
  df['mean identity (over all positions in all samples)'] = y
  df['mean (over all samples) number of positions with at least that coverage'] = \
  MeanNumPosExceedingCoveragesPerSample
  df.to_csv(OutFileCsv, index=False)

else:

  # Read in the csv input
  df = pd.read_csv(args.csv_for_replotting)
  x = df['coverage']
  y = df['mean identity (over all positions in all samples)']
  MeanNumPosExceedingCoveragesPerSample = \
  df['mean (over all samples) number of positions with at least that coverage']

ax = plt.figure().add_subplot(111)

plt.scatter(x, y, marker=args.point_marker, s=args.point_size,
linewidth=args.point_thickness, alpha=args.point_alpha)

for i, MeanNumPos in enumerate(MeanNumPosExceedingCoveragesPerSample):
  ax.annotate("{:.0f}".format(MeanNumPos), (x[i], y[i]),
  fontsize = args.point_font_size)


if args.x_min_max:
  ax.set_xlim(xmin=Xmin, xmax=Xmax)
if args.y_min_max:
  ax.set_ylim(ymin=Ymin, ymax=Ymax)


ax.set_xscale('log')
plt.xlabel(XaxisLabel, fontsize=args.axis_font_size)
plt.ylabel(YaxisLabel, fontsize=args.axis_font_size)
ax.tick_params(axis='both', which='major', labelsize=args.axis_font_size)
plt.title(args.title, fontsize=args.title_font_size)
plt.tight_layout()
plt.savefig(OutFilePdf)
