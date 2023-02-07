#!/usr/bin/env bash

set -u
set -o pipefail

UsageInstructions=$(echo '
Arguments for this script:
(1) an output directory for the initialisation files.
(2) the configuration file, containing all your parameter choices etc.;
(3) your chosen alignment of references;
(4) a fasta file of the adapters used in sequencing;
(5) a fasta file of the primers used in sequencing.
')

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=5
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
OutDir="$1"
ConfigFile="$2"
RefAlignment="$3"
adapters="$4"
primers="$5"

# Source the shiver funcs, check files exist, source the config file, check it.
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.sh'
CheckFilesExist "$ConfigFile" "$RefAlignment"
CheckConfig "$ConfigFile" true false false || \
{ echo "Problem with $ConfigFile. Quitting." >&2 ; exit 1 ; }

# If OutDir does not exist, try to create it.
if [ ! -d "$OutDir" ]; then
  mkdir "$OutDir"
fi || { echo 'Unable to create the specified output directory. (NB it should' \
'not exist already.) Quitting.' >&2 ; exit 1; }

# Remove a trailing slash, if present.
OutDir=$(cd "$OutDir"; pwd)

# Check that OutDir does not have whitespace in it
if [[ "$OutDir" =~ ( |\') ]]; then
  echo "Your specified directory $OutDir contains whitespace; we need to be"\
  "able to build a blast database in there, and unfortunately, blast cannot"\
  "handle whitespace in paths (stupid, I know). Try again with a different"\
  "directory. Quitting." >&2
  exit 1;
fi

# Check OutDir is empty
if ! find "$OutDir"/ -maxdepth 0 -empty | read v; then
  echo "$OutDir is not empty. Delete or move its contents. Quitting." >&2
  exit 1;
fi

# Some files we'll create
NewRefAlignment="$OutDir"/'ExistingRefAlignment.fasta'
RefList="$OutDir"/'ExistingRefNamesSorted.txt'
UngappedRefs="$OutDir"/'ExistingRefsUngapped.fasta'
database="$OutDir"/'ExistingRefsBlastDatabase'

# Copy the three input fasta files into the initialisation directory, removing
# pure-gap columns from RefAlignment.
"$python" "$Code_RemoveBlankCols" "$RefAlignment" > "$NewRefAlignment" || \
{ echo "Problem removing pure-gap columns from $RefAlignment. Quitting." >&2 ; \
exit 1; }
cp "$adapters" "$OutDir"/'adapters.fasta'
cp "$primers" "$OutDir"/'primers.fasta'

if [[ "$TrimPrimerWithOneSNP" == "true" ]]; then
  "$python" "$Code_AddSNPsToSeqs" "$OutDir/primers.fasta" \
  "$OutDir/PrimersWithSNPs.fasta" || { echo "Problem generating all possible"\
  "variants of the sequences in $primers differing by a single base mutation."\
  "Quitting." >&2; exit 1; }
fi

# List all names in the reference alignment
awk '/^>/ {print substr($1,2)}' "$NewRefAlignment" | sort > "$RefList"

# Check that RefAlignment has some sequences, that their IDs are unique, and
# that their IDs don't contain commas.
NumRefs=$(wc -l "$RefList" | awk '{print $1}')
if [[ $NumRefs -eq 0 ]]; then
  echo "$RefAlignment contains no sequences. Quitting." >&2
  exit 1;
fi
NumUniqueIDs=$(uniq "$RefList" | wc -l)
if [[ $NumUniqueIDs -ne $NumRefs ]]; then
  echo "$RefAlignment contains some identically named sequences. Rename"\
  'these and try again. Quitting.' >&2
  exit 1;
fi
RefNames=$(cat "$RefList")
if [[ "$RefNames" == *","* ]]; then
  echo "Reference names must not contain commas. Quitting." >&2
  exit 1
fi

# Ungap RefAlignment
"$python" "$Code_UngapFasta" "$NewRefAlignment" > "$UngappedRefs" || \
{ echo "Problem ungapping $RefAlignment. Quitting." >&2 ; exit 1; }

# Create the blast database
"$BlastDBcommand" -dbtype nucl -in "$UngappedRefs" -input_type fasta -out \
"$database" || \
{ echo 'Problem creating a blast database out of' \
"$OutDir/ExistingRefsUngapped.fasta. Quitting." >&2 ; exit 1; }

# Make a set of files each containing one (ungapped) sequence from the reference
# alignment.
IndividualRefDir="$OutDir"/'IndividualRefs'
mkdir -p "$IndividualRefDir" || \
{ echo "Problem making the directory $IndividualRefDir. Quitting." >&2 ; \
exit 1; }
"$python" "$Code_SplitFasta" -G "$RefAlignment" "$IndividualRefDir" || { echo "Problem" \
"splitting $RefAlignment into one file per sequence. Quitting." >&2 ; exit 1; }
