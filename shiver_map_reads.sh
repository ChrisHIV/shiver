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
choice might be the contig file name minus its path and extension);
(5) the blast file created by the shiver_align_contigs.sh command;
(6) either the alignment of contigs to refs produced by the
shiver_align_contigs.sh command, or a fasta file containing a single reference
to be used for mapping;
(7) the forward reads;
(8) the reverse reads.
')

################################################################################
# PRELIMINARIES

# Check for the right number of arguments. Assign them to variables.
NumArgsExpected=8
if [ "$#" -ne "$NumArgsExpected" ]; then
  echo $UsageInstructions
  echo "$#" 'arguments specified;' "$NumArgsExpected" 'expected. Quitting' >&2
  exit 1
fi
InitDir="$1"
ConfigFile="$2"
RawContigsFile="$3"
SID="$4"
ContigBlastFile="$5"
FastaFile="$6"
reads1="$7"
reads2="$8"

# Check InitDir exists. Remove a trailing slash, if present.
if [ ! -d "$InitDir" ]; then
  echo "$InitDir does not exist or is not a directory. Quitting." >&2
  exit 1
fi
InitDir=$(cd "$InitDir"; pwd)

RefList="$InitDir"/'ExistingRefNamesSorted.txt'
ExistingRefAlignment="$InitDir"/'ExistingRefAlignment.fasta'
adapters="$InitDir"/'adapters.fasta'
primers="$InitDir"/'primers.fasta'

# Source the shiver funcs, check files exist, source the config file, check it.
ThisDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
source "$ThisDir"/'shiver_funcs.sh'
CheckFilesExist "$ConfigFile" "$reads1" "$reads2" "$RawContigsFile" \
"$ContigBlastFile" "$FastaFile" "$RefList" "$ExistingRefAlignment" "$adapters" \
"$primers"
CheckConfig "$ConfigFile" false false true || \
{ echo "Problem with $ConfigFile. Quitting." >&2 ; exit 1 ; }

# Check the primer + SNPs file exists, if it is needed.
if [[ "$TrimPrimerWithOneSNP" == "true" ]]; then
  PrimersToUse="$InitDir/PrimersWithSNPs.fasta"
  if [[ ! -f "$PrimersToUse" ]]; then
    echo "The file $PrimersToUse does not exist; it should have been created"\
    "when you ran shiver_init.sh. Perhaps the variable 'TrimPrimerWithOneSNP'"\
    "in the config file was set to false when you ran shiver_init.sh, and you"\
    "changed it to true afterwards? You can generate this file from your"\
    "primers by running" >&2
    echo "$python" "$Code_AddSNPsToSeqs $InitDir/primers.fasta"\
    "$InitDir/PrimersWithSNPs.fasta" >&2
    echo "(Though be careful copying and pasting that if any file paths"\
    "contain whitespace.) Quitting." >&2
    exit 1
  fi
else
  PrimersToUse="$primers"
fi

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

# Some files we'll create
TheRef="$SID$OutputRefSuffix"
consensus="$SID"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2.fasta"
ConsensusForGlobalAln="$SID"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2$GlobalAlnSuffix"
consensusWcontigs="$SID"'_consensus_MinCov_'"$MinCov1"'_'"$MinCov2"'_wContigs.fasta'
MappedContaminantReads="$SID$MappedContaminantReadsSuffix"
cleaned1reads="$SID$ReadsPreMapping1Suffix"
cleaned2reads="$SID$ReadsPreMapping2Suffix"
CoordsDict="$SID$CoordsDictSuffix"
BaseFreqs="$SID$BaseFreqsSuffix"
BaseFreqsWGlobal="$SID$BaseFreqsWGlobalSuffix"
################################################################################

# Check we have some reads
CheckNonEmptyReads "$reads1" || { echo 'Quitting.' >&2; exit 3; }

################################################################################
# CONSTRUCT A REFERENCE, OR USE THE ONE SUPPLIED

# List the HIV contigs
awk -F, '{print $1}' "$ContigBlastFile" | sort | uniq > "$HIVContigsListOrig"

# FastaFile should be either a single seq, which we map to as is, or else an
# alignment of contigs to real refs.
RefIsInAlignment=true
NumSeqsInFastaFile=$(grep -e '^>' "$FastaFile" | wc -l)
if [ "$NumSeqsInFastaFile" -eq 0 ]; then
  echo 'Error: there are no sequences in' "$FastaFile". 'Quitting.' >&2
  exit 1

elif [ "$NumSeqsInFastaFile" -eq 1 ]; then

  # Try to find the sequence in FastaFile in ExistingRefAlignment.
  RefName=$(awk '/^>/ {print substr($1,2)}' "$FastaFile")
  "$python" "$Code_FindSeqsInFasta" "$ExistingRefAlignment" -g -N "$RefName" > \
  "$RefFromAlignment" || \
  { echo "Could not find seq $RefName in $ExistingRefAlignment; that's OK," \
  'but after mapping we will not be able to produce a version of the' \
  'consensus seq with gaps suitable for a global alignment. Continuing.' ; \
  RefIsInAlignment=false ; }

  # Compare the sequence in FastaFile to the one in ExistingRefAlignment.
  if $RefIsInAlignment; then
    equal=$("$python" "$Code_CheckFastaFileEquality" "$RefFromAlignment" "$FastaFile") ||
    { echo 'Problem running' "$Code_CheckFastaFileEquality"'. Quitting.' >&2 ; \
    exit 1 ; }
    if [[ "$equal" == "false" ]]; then
      echo 'Seq' "$RefName" 'differs between' "$FastaFile" 'and' \
      "$ExistingRefAlignment; that's OK," \
      'but after mapping we will not be able to produce a version of the' \
      'consensus seq with gaps suitable for a global alignment. Continuing.' ;
      RefIsInAlignment=false
    elif [[ "$equal" != "true" ]]; then
      echo "Unexpected output \"$equal\" from $Code_CheckFastaFileEquality."\
      "Quitting." >&2
      exit 1
    fi
  fi

  # Set the flag appropriate for use of a real ref, when it comes to
  # coordinate translation for a global alignment.
  GlobalAlignExcisionFlag='-d'

  cp "$FastaFile" "$TheRef"
  cp "$ExistingRefAlignment" "$TempRefAlignment"

  # Extract those contigs that have a blast hit.
  NumHIVContigsOrig=$(wc -w "$HIVContigsListOrig" | awk '{print $1}')
  if [[ $NumHIVContigsOrig -gt 0 ]]; then
    "$python" "$Code_FindSeqsInFasta" "$RawContigsFile" -F "$HIVContigsListOrig" > \
    "$RawContigFile2" || \
    { echo 'Problem extracting the HIV contigs. Quitting.' >&2 ; exit 1 ; }
  fi

else

  ContigToRefAlignment="$FastaFile"

  # ContigToRefAlignment should contain the same set of sequences in the input
  # existing reference alignment, plus the contigs.
  awk '/^>/ {print substr($1,2)}' "$ContigToRefAlignment" | sort \
  | sed $'s/\r//g' > "$AllSeqsInAln"
  MissingRefs=$(comm -1 -3 "$AllSeqsInAln" "$RefList")
  NumMissingRefs=$(echo $MissingRefs | wc -w)
  if [ $NumMissingRefs -gt 0 ]; then
    echo "Error: the following references from $ExistingRefAlignment are"\
    "missing from $ContigToRefAlignment: $MissingRefs. Quitting." >&2
    exit 1
  fi
  comm -2 -3 "$AllSeqsInAln" "$RefList" > "$HIVContigsListUser"
  NumHIVContigsUser=$(wc -w "$HIVContigsListUser" | awk '{print $1}')
  if [ $NumHIVContigsUser -eq 0 ]; then
    echo "Error: no contigs found in $ContigToRefAlignment. Quitting" >&2
    exit 1
  fi

  # Extract just the existing references (i.e. everything but the contigs) from
  # ContigToRefAlignment, and check that they are the same as in
  # ExistingRefAlignment. Also extract just the contigs, stripping gaps, ready
  # for later.
  "$python" "$Code_FindSeqsInFasta" "$ContigToRefAlignment" -F "$HIVContigsListUser" -v > \
  "$TempRefAlignment" &&
  "$python" "$Code_FindSeqsInFasta" "$ContigToRefAlignment" -F "$HIVContigsListUser" -g > \
  "$RawContigFile2" || { echo 'Problem separating the contigs and existing'\
  "refs in $ContigToRefAlignment. Quitting." >&2 ; exit 1 ; }
  "$python" "$Code_RemoveBlankCols" "$TempRefAlignment" > "$AlignmentForTesting" || \
  { echo "Problem removing pure-gap columns from $TempRefAlignment (which was"\
  "created by removing the contigs from $ContigToRefAlignment - that's"\
  "probably the problematic file). Quitting." >&2 ; exit 1; }
  equal=$("$python" "$Code_CheckFastaFileEquality" "$AlignmentForTesting" \
  "$ExistingRefAlignment") || { echo "Problem running"\
  "$Code_CheckFastaFileEquality. Quitting." >&2 ; exit 1 ; }
  if [[ "$equal" == "false" ]]; then
    echo "The reference sequences in $ContigToRefAlignment are different from"\
    "those in $ExistingRefAlignment. If you modify your alignment of contigs"\
    "to existing reference sequences, you should only modify the contig"\
    "sequences, not the references. Quitting." >&2
    exit 1
  elif [[ "$equal" != "true" ]]; then
    echo "Unexpected output \"$equal\" from $Code_CheckFastaFileEquality."\
    "Quitting." >&2
    exit 1
  fi

  # Construct the tailored ref
  HIVcontigNames=$(cat "$HIVContigsListUser")
  "$python" "$Code_ConstructRef" "$ContigToRefAlignment" "$GappyRefWithExtraSeq" \
  $HIVcontigNames || \
  { echo 'Failed to construct a ref from the alignment. Quitting.' >&2 ; \
  exit 1 ; }

  # Extract just the constructed ref (the first sequence)
  awk '/^>/{if(N)exit;++N;} {print;}' "$GappyRefWithExtraSeq" > "$RefWithGaps"

  # Remove any gaps from the reference
  "$python" "$Code_UngapFasta" "$RefWithGaps" > "$TheRef" || \
  { echo 'Gap stripping code failed. Quitting.' >&2 ; exit 1 ; }

  RefName=$(awk '/^>/ {print substr($1,2)}' "$TheRef")

  # Set the flag appropriate for use of a constructed ref, when it comes to
  # coordinate translation for a global alignment.
  GlobalAlignExcisionFlag='-e'

  # Create a version of the alignment of contigs to real refs, with the contigs
  # replaced by the constructed ref, ready for coordinate translation later.
  cat "$RefWithGaps" >> "$TempRefAlignment"

fi

################################################################################

################################################################################
# TRIM & CLEAN READS

# Copy the reads to the working directory.
if [[ $(dirname "$reads1") != "." ]]; then
  NewReads1=$(basename "$reads1")
  if [[ -f "$NewReads1" ]]; then
    echo "A file called $NewReads1 exists in the working directory, $PWD,"\
    "already; we will not overwrite it with $reads1. Quitting." >&2; exit 1;
  fi
  cp -i "$reads1" . || { echo "Failed to copy $reads1 to the working"\
  "directory. Quitting." >&2; exit 1; }
  reads1="$NewReads1"
fi
if [[ $(dirname "$reads2") != "." ]]; then
  NewReads2=$(basename "$reads2")
  if [[ -f "$NewReads2" ]]; then
    echo "A file called $NewReads2 exists in the working directory, $PWD,"\
    "already; we will not overwrite it with $reads2. Quitting." >&2; exit 1;
  fi
  cp -i "$reads2" . || { echo "Failed to copy $reads2 to the working"\
  "directory. Quitting." >&2; exit 1; }
  reads2="$NewReads2"
fi

# Unzip the reads if they end in .gz.
if [[ "$reads1" == *.gz ]]; then
  gunzip -f "$reads1" &&
  NewReads1="${reads1%.gz}" &&
  ls $NewReads1 > /dev/null || { echo "Problem unzipping $reads1. Quitting.">&2;
  exit 1; }
  reads1="$NewReads1"
fi
if [[ "$reads2" == *.gz ]]; then
  gunzip -f "$reads2" &&
  NewReads2="${reads2%.gz}" &&
  ls $NewReads2 > /dev/null || { echo "Problem unzipping $reads2. Quitting.">&2;
  exit 1; }
  reads2="$NewReads2"
fi

# Check all 1 read seq ids end in /1, and 2 reads in /2. Check there are
# no tabs in the seq id lines.
CheckReadNames "$reads1" 1 && CheckReadNames "$reads2" 2 || \
{ echo 'Problem with read names. Quitting.' >&2 ; exit 1 ; }

HaveModifiedReads=false

# Read trimming:
if [[ "$TrimReadsForAdaptersAndQual" == "true" ]]; then

  # Trim adapters and low-quality bases
  echo 'Now trimming reads - typically a slow step.'
  $trimmomatic PE -quiet -threads $NumThreadsTrimmomatic \
  "$reads1" "$reads2" "$reads1trim1" "$reads1trimmings" "$reads2trim1" \
  "$reads2trimmings" ILLUMINACLIP:"$adapters":"$IlluminaClipParams" \
  $BaseQualityParams || \
  $trimmomatic PE -threads $NumThreadsTrimmomatic \
  "$reads1" "$reads2" "$reads1trim1" "$reads1trimmings" "$reads2trim1" \
  "$reads2trimmings" ILLUMINACLIP:"$adapters":"$IlluminaClipParams" \
  $BaseQualityParams || { echo 'Problem running trimmomatic. Quitting.' >&2 ; \
  exit 1 ; }

  HaveModifiedReads=true
  reads1="$reads1trim1"
  reads2="$reads2trim1"

fi
if [[ "$TrimReadsForPrimers" == "true" ]]; then

  # Trim primers
  "$fastaq" 'sequence_trim' --revcomp "$reads1" "$reads2" "$reads1trim2" \
  "$reads2trim2" "$PrimersToUse" || \
  { echo 'Problem running fastaq. Quitting.' >&2 ; exit 1 ; }
  echo "fastaq completed successfully."

  HaveModifiedReads=true
  reads1="$reads1trim2"
  reads2="$reads2trim2"
fi

# Check some reads are left after trimming
if $HaveModifiedReads; then
  CheckNonEmptyReads "$reads1" ||
  { echo 'No reads left after trimming. Quitting.' >&2; exit 3; }
fi

# If RawContigsFile is empty we cannot do read cleaning, so switch it from true
# to false if necessary, and warn.
if [[ "$CleanReads" == "true" ]]; then

  # List all the contigs and the HIV ones.
  awk '/^>/ {print substr($1,2)}' "$RawContigsFile" | sort > "$AllContigsList"

  # Check there are some contigs
  NumContigs=$(wc -l "$AllContigsList" | awk '{print $1}')
  if [ "$NumContigs" -eq 0 ]; then
    echo "WARNING: the config variable 'CleanReads' was set to 'true', but"\
    "cleaning is not possible as $RawContigsFile is empty. Setting"\
    "'CleanReads' to 'false' and continuing." >&2
    CleanReads=false
  fi
fi

# If we're not cleaning reads:
if [[ "$CleanReads" != "true" ]]; then

  # If we have trimmed the reads, change their name to the final output name for
  # the reads, to indicate that something has been done to them (we don't want
  # their name to continue beginning 'temp').
  # If we haven't trimmed, then just make the cleaned reads variables point to
  # the unprocessed reads: we're not doing anything to the read files provided
  # as input so there's no need to rename to indicate that they're shiver
  # output worth saving separately.
  if $HaveModifiedReads; then
    mv "$reads1" "$cleaned1reads"
    mv "$reads2" "$cleaned2reads"
  else
    cleaned1reads="$reads1"
    cleaned2reads="$reads2"
  fi
else

  # Check that there aren't any contigs appearing in the blast file & missing
  # from the file of contigs.
  NumUnknownContigsInBlastHits=$(comm -1 -3 "$AllContigsList" \
  "$HIVContigsListOrig" | wc -l | awk '{print $1}')
  if [ "$NumUnknownContigsInBlastHits" -ne 0 ]; then
    echo 'Error: the following contigs are named in' "$ContigBlastFile"\
    'but are not in' "$RawContigsFile"':' >&2
    comm -1 -3 "$AllContigsList" "$HIVContigsListOrig" >&2
    echo 'Quitting.' >&2
    exit 1
  fi

  # Find contigs without a blast hit: includes both contaminants and contigs
  # that are too short...
  comm -3 "$AllContigsList" "$HIVContigsListOrig" > "$DiscardedContigNames"
  NumContaminantContigs=$(wc -l $DiscardedContigNames | awk '{print $1}')

  # ...and now remove from these contigs those that are too short, leaving only
  # contaminants.
  if [ "$NumContaminantContigs" -gt 0 ]; then
    "$python" "$Code_FindSeqsInFasta" "$RawContigsFile" -F "$DiscardedContigNames" \
    --min-length "$MinContigLength" > "$RefAndContaminantContigs" ||
    { echo "Problem extracting contaminant contigs from $RawContigsFile." \
    "Quitting." >&2; exit 1; }
    ContaminantContigNames=$(awk '/^>/ {print substr($1,2)}' \
    "$RefAndContaminantContigs")
    NumContaminantContigs=$(echo $ContaminantContigNames | wc -w)
  fi

  # If there are no contaminant contigs, we don't need to clean.
  # We create a blank mapping file to more easily keep track of the fact that
  # there are no contaminant reads in this case.
  if [ "$NumContaminantContigs" -eq 0 ]; then
    echo 'There are no contaminant contigs: read cleaning unnecessary.'
    echo -n > "$MappedContaminantReads"
    if $HaveModifiedReads; then
      mv "$reads1" "$cleaned1reads"
      mv "$reads2" "$cleaned2reads"
    else
      cleaned1reads="$reads1"
      cleaned2reads="$reads2"
    fi

  # We enter this scope if there are some contaminant contigs:
  else

    # Make a blast database out of the contaminant contigs and the ref.
    cat "$TheRef" >> "$RefAndContaminantContigs"
    "$BlastDBcommand" -dbtype nucl -in "$RefAndContaminantContigs" \
    -input_type fasta -out "$BlastDB" || \
    { echo 'Problem creating a blast database. Quitting.' >&2 ; exit 1 ; }

    # Convert fastq to fasta.
    "$fastaq" to_fasta "$reads1" "$reads1asFasta" &&
    "$fastaq" to_fasta "$reads2" "$reads2asFasta" || \
    { echo 'Problem converting the reads from fastq to fasta. Quitting.' >&2 ; \
    exit 1 ; }

    # Blast the reads.
    echo 'Now blasting the reads - typically a slow step.'
    "$BlastNcommand" -query "$reads1asFasta" -db "$BlastDB" -out \
    "$reads1blast1" -max_target_seqs 1 -outfmt \
    '10 qacc sacc sseqid evalue pident qstart qend sstart send' &&
    "$BlastNcommand" -query "$reads2asFasta" -db "$BlastDB" -out \
    "$reads2blast1" -max_target_seqs 1 -outfmt \
    '10 qacc sacc sseqid evalue pident qstart qend sstart send' || \
    { echo 'Problem blasting' "$ContigFile"'. Quitting.' >&2 ; exit 1 ; }

    # For multiple blast hits, keep the one with the highest evalue
    # TODO: test what blast does with fasta headers that have comments in them -
    # does it include them too?
    "$python" "$Code_KeepBestLinesInDataFile" "$reads1blast1" "$reads1blast2" &&
    "$python" "$Code_KeepBestLinesInDataFile" "$reads2blast1" "$reads2blast2" || 
    { echo "Problem extracting the best blast hits using"\
    "$Code_KeepBestLinesInDataFile. Quitting." >&2 ; exit 1 ; }

    # Find the read pairs that blast best to something other than the reference.
    "$python" "$Code_FindContaminantReadPairs" "$reads1blast2" "$reads2blast2" \
    "$RefName" "$BadReadsBaseName" && ls "$BadReadsBaseName"_1.txt \
    "$BadReadsBaseName"_2.txt > /dev/null 2>&1 || \
    { echo 'Problem finding contaminant read pairs using' \
    "$Code_FindContaminantReadPairs. Quitting." >&2 ; exit 1 ; }

    # If none of the read pairs blast better to contaminant contigs than the
    # reference, we just duplicate the original short read files.
    NumContaminantReadPairs=$(wc -l "$BadReadsBaseName"_1.txt | \
    awk '{print $1}')
    if [ "$NumContaminantReadPairs" -eq 0 ]; then
      echo 'There are no contaminant read pairs.'
      echo -n > "$MappedContaminantReads"
      if $HaveModifiedReads; then
        mv "$reads1" "$cleaned1reads"
        mv "$reads2" "$cleaned2reads"
      else
        cleaned1reads="$reads1"
        cleaned2reads="$reads2"
      fi

    # We enter this scope if there are some read pairs that blast better to
    # contaminant contigs than the reference.
    else

      # Check every read has a mate.
      # TODO: move the 'unpaired' check right to the beginning?
      if ! cmp <(awk '{if ((NR-1)%4==0) print substr($1,2,length($1)-3)}' \
      "$reads1" | sort) \
      <(awk '{if ((NR-1)%4==0) print substr($1,2,length($1)-3)}' \
      "$reads2" | sort); then
        echo 'At least one read in' "$reads1" 'or' "$reads2" 'is unpaired.' \
        'Quitting.' >&2 ; exit 1 ;
      fi

      # Extract the non-contaminant read pairs
      "$python" "$Code_FindReadsInFastq" -v -s "$reads1" "$BadReadsBaseName"_1.txt > \
      "$cleaned1reads" &&
      "$python" "$Code_FindReadsInFastq" -v -s "$reads2" "$BadReadsBaseName"_2.txt > \
      "$cleaned2reads" || \
      { echo 'Problem extracting the non-contaminant reads using' \
      "$Code_FindReadsInFastq"'. Quitting.' >&2 ; exit 1 ; }

      # Check some reads are left after cleaning
      CheckNonEmptyReads "$cleaned1reads" ||
      { echo 'No reads left after cleaning. Quitting.' >&2; exit 3; }

      # Map the contaminant reads to the reference, to measure how useful the
      # cleaning procedure was.
      if [[ "$MapContaminantReads" == "true" ]]; then
        "$python" "$Code_FindReadsInFastq" -s "$reads1" "$BadReadsBaseName"_1.txt > \
        "$BadReadsBaseName"_1.fastq &&
        "$python" "$Code_FindReadsInFastq" -s "$reads2" "$BadReadsBaseName"_2.txt > \
        "$BadReadsBaseName"_2.fastq || \
        { echo 'Problem extracting the contaminant reads using' \
        "$Code_FindReadsInFastq. Quitting." >&2 ; exit 1 ; }
        BamOnly=true
        map "$BadReadsBaseName"_1.fastq "$BadReadsBaseName"_2.fastq "$TheRef" \
        "$MappedContaminantReads" "$BamOnly" || { echo "Problem mapping the" \
        "contaminant reads to $RefName using smalt. Quitting." >&2 ; exit 1 ; }
      fi

      HaveModifiedReads=true

    fi
  fi
fi

# Prepend modified read filenames by temp_ if desired.
if [[ "$KeepPreMappingReads" == "false" ]] && $HaveModifiedReads; then
  mv "$cleaned1reads" "temp_$cleaned1reads"
  mv "$cleaned2reads" "temp_$cleaned2reads"
  cleaned1reads="temp_$cleaned1reads"
  cleaned2reads="temp_$cleaned2reads"
fi

################################################################################
# MAP

OldMafft=false

# Do the mapping
BamOnly=false
map "$cleaned1reads" "$cleaned2reads" "$TheRef" "$SID" "$BamOnly"
MapStatus=$?
if [[ $MapStatus == 3 ]]; then
  echo "Quitting." >&2
  exit 3
elif [[ $MapStatus != 0 ]]; then
  echo "Problem mapping to the reference in $TheRef. Quitting." >&2
  exit 1
fi

# Add gaps and excise unique insertions, to allow this consensus to be added to
# a global alignment with others.
if $RefIsInAlignment; then
  "$python" "$Code_MergeAlignments" "$GlobalAlignExcisionFlag" -L "$CoordsDict" \
  "$TempRefAlignment" "$consensus" > "$ConsensusForGlobalAln" ||
  { echo 'Problem translating the coordinates of the consensus for the'\
  'global alignment. Quitting.' >&2 ; exit 1 ; }

  # Add the global alignment coordinates to the base frequencies file.
  "$python" "$Code_MergeBaseFreqsAndCoords" "$BaseFreqs" -C "$CoordsDict" > \
  "$BaseFreqsWGlobal" || { echo 'Problem adding the global alignment'\
  'coordinates to the base frequencies file. Quitting.' >&2 ; exit 1 ; }

fi

if [[ "$remap" == "true" ]]; then

  echo 'Beginning remapping to the consensus.'

  # New file names
  NewSID="$SID"'_remap'
  NewRef="$NewSID$OutputRefSuffix"
  NewRefName="$SID"'_ConsensusRound1_GapsFilled'

  # Fill in any gaps in the consensus with the corresponding part of the orginal
  # reference for mapping.
  "$python" "$Code_FillConsensusGaps" "$consensus" '--output-seq-name' \
  "$NewRefName" > "$NewRef" || { echo 'Problem'\
  'filling in gaps in the consensus with the corresponding part of the orginal'\
  'reference for mapping. Quitting.' >&2 ; exit 1 ; }

  # Map!
  BamOnly=false
  map "$cleaned1reads" "$cleaned2reads" "$NewRef" "$NewSID" "$BamOnly" ||
  { echo 'Problem remapping to the consensus from the first round of mapping.'\
  'Quitting.' >&2 ; exit 1 ; }

fi

