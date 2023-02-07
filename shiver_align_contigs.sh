#!/usr/bin/env bash

set -u
set -o pipefail

UsageInstructions=$(echo '
Arguments for this script:
(1) the initialisation directory you created using the shiver_init.sh command;
(2) the configuration file, containing all your parameter choices etc.;
(3) a fasta file of contigs (output from processing the short reads with an
assembly program);
(4) A sample ID ("SID") used for naming the output from this script (a sensible
choice might be the contig file name minus its path and extension).
If this script completes successfully it will produce a .blast file (detailing
the contigs that have blast hits). If the .blast file is not empty it will
produce a fasta file of those contigs with hits and another fasta file of
these contigs aligned to the refereces; if cutting and/or reversing of contigs
is necessary, two more fasta files are produced - the cut/reversed contigs on
their own and also aligned to references (i.e. there will be two files of the
contigs on their own and two files of the contigs aligned to references).
')

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=4
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
InitDir="$1"
ConfigFile="$2"
ContigFile="$3"
SID="$4"

# Check InitDir exists. Remove a trailing slash, if present.
if [ ! -d "$InitDir" ]; then
  echo "$InitDir does not exist or is not a directory. Quitting." >&2
  exit 1
fi
InitDir=$(cd "$InitDir"; pwd)

BlastDatabase="$InitDir/ExistingRefsBlastDatabase"
RefAlignment="$InitDir/ExistingRefAlignment.fasta"

# Source the shiver funcs, check files exist, source the config file, check it.
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.sh'
CheckFilesExist "$ContigFile" "$RefAlignment"
CheckConfig "$ConfigFile" false true false || \
{ echo "Problem with $ConfigFile. Quitting." >&2 ; exit 1 ; }

# Print how this script was called, and what the config file parameter values
# were.
echo '###############################################'
echo "Info:" $(basename "$0") "was called thus:"
echo $0 $@
echo "With these config file parameter values:"
awk '{if (substr($0, 1, 23) == "# Suffixes we'\''ll append") {exit};
if (substr($0, 1, 1) == "#") {next} else if (NF == 0) {next} else print}' \
"$ConfigFile"
echo '###############################################'
echo

# Out files we'll make
BlastFile="$SID$BlastSuffix"
MergedBlastFile="$SID$MergedBlastSuffix"
RawContigAlignment="$SID"'_raw_wRefs.fasta'
CutContigAlignment="$SID"'_cut_wRefs.fasta'

# Extract just the HIV contigs (those that blast to the refs) and put them in
# $RawContigFile1
GetHIVcontigs "$ContigFile" "$ContigsNoShortOnes" "$BlastFile" "$RawContigFile1"
GetHIVcontigsStatus=$?
if [[ $GetHIVcontigsStatus == 3 ]]; then
  echo "No HIV contigs to analyse for $SID. Quitting." >&2 
  exit 3
elif [[ $GetHIVcontigsStatus != 0 ]]; then
  echo "Problem encountered while checking the contigs in $ContigFile and"\
  "extracting those thought to be HIV. Quitting." >&2
  exit 1
fi

# Align the HIV contigs.
OldMafft=false
SwapContigsToTop=true
AlignContigsToRefs "$mafft" '--quiet' "$RawContigFile1" "$RefAlignment" \
"$RawContigAlignment" "$SwapContigsToTop" "$OldMafft" || \
{ echo 'Problem aligning the raw contigs to refs. Quitting.' >&2 ; exit 1 ; }

PrintAlnLengthIncrease "$RefAlignment" "$RawContigAlignment" || \
{ echo "Problem checking the alignment length increase after adding the raw" \
"contigs to the existing references. Quitting." >&2 ; exit 1 ; }

# Run the contig cutting & flipping code
"$python" "$Code_CorrectContigs" "$BlastFile" "$ContigMinBlastOverlapToMerge" \
-C "$ContigFile" -O "$CutContigFile" -B "$MergedBlastFile" || \
{ echo "Problem encountered running $Code_CorrectContigs. Quitting." >&2 ; \
exit 1; }

# Set up the arguments for when we cut the aligned contigs.
if [[ "$TrimToKnownGenome" == "true" ]]; then
  CutAlignedContigsArgs="--split-gap-size $MinGapSizeToSplitGontig \
  --min-contig-size $MinContigFragmentLength --trim-overhangs"
else
  CutAlignedContigsArgs="--split-gap-size $MinGapSizeToSplitGontig \
  --min-contig-size $MinContigFragmentLength"
fi

# If the contigs needed correcting, align the corrected contigs. 
if [[ -f "$CutContigFile" ]]; then

  # Align
  AlignContigsToRefs "$mafft" '--quiet' "$CutContigFile" "$RefAlignment" \
  "$TempContigAlignment3" "$SwapContigsToTop" "$OldMafft" || \
  { echo 'Problem aligning the cut/modified contigs to refs. Quitting.' >&2 ;
  exit 1 ; }

  # Split gappy contigs after alignment.
  CutContigNames=$(awk '/^>/ {print substr($1,2)}' "$CutContigFile")
  "$python" "$Code_CutAlignedContigs" "$TempContigAlignment3" $CutContigNames \
  $CutAlignedContigsArgs > "$CutContigAlignment"
  CutAlignedContigsStatus=$?
  if [[ $CutAlignedContigsStatus == 3 ]]; then
    echo "After contig correction, all (pieces of) contigs for $SID were below"\
    "the length threshold. Quitting." >&2
    exit 3
  elif [[ $CutAlignedContigsStatus != 0 ]]; then
    echo 'Problem splitting gappy contigs after alignment. Quitting.' >&2
    exit 1
  fi

  PrintAlnLengthIncrease "$RefAlignment" "$CutContigAlignment" || \
  { echo "Problem checking the alignment length increase after adding the cut"\
  "contigs to the existing references. Quitting." >&2 ; exit 1 ; }

else
  # If we don't have a cut contig file, the split-gappy-contigs code might still
  # want to split the aligned raw contigs. If so, name it the same as aligned
  # cut contigs for consistency.
  RawContigNames=$(awk '/^>/ {print substr($1,2)}' "$RawContigFile1")
  "$python" "$Code_CutAlignedContigs" "$RawContigAlignment" $RawContigNames \
  $CutAlignedContigsArgs > "$CutContigAlignment"
  CutAlignedContigsStatus=$?
  if [[ $CutAlignedContigsStatus == 3 ]]; then
    echo "After contig correction, all (pieces of) contigs for $SID were below"\
    "the length threshold. Quitting." >&2
    exit 3
  elif [[ $CutAlignedContigsStatus != 0 ]]; then
    echo 'Problem splitting gappy contigs after alignment. Quitting.' >&2
    exit 1
  fi

  equal=$("$python" "$Code_CheckFastaFileEquality" "$RawContigAlignment" \
  "$CutContigAlignment") || { echo "Problem running"\
  "$Code_CheckFastaFileEquality. Quitting." >&2 ; exit 1 ; }
  if [[ "$equal" == "true" ]]; then
    rm "$CutContigAlignment"
  elif [[ "$equal" == "false" ]]; then
    PrintAlnLengthIncrease "$RefAlignment" "$CutContigAlignment" || \
    { echo "Problem checking the alignment length increase after adding the"\
    "cut contigs to the existing references. Quitting." >&2 ; exit 1 ; }
  else
    echo "Unexpected output \"$equal\" from $Code_CheckFastaFileEquality."\
    "Quitting." >&2
    exit 1
  fi

fi

