#!/bin/sh

# This requires a separate gff file for each protein that contains all the exons of that protein.
# These gff files are in the directory defined by $GFFDIR
# This uses the strand (+ vs - in column 7 of the gff) to determine what is upstream or downstream.
# The AGAT script agat_sp_extract_sequences.pl is then used to extract the sequence upstream and downstream
# of the splice. No sequence is taken from the upstream side of the first exon and the downstream side of the
# of the last exon since this is not a splice site.
# Seqtk is used to remoformat the line widths in the fasta files. 
#
# Matthew Kweskin, kweskinm@si.edu

INFASTA=./Abisporus_varbisporusH97.v2_AssembledScaffolds.fasta
OUTDIRUP=./up-fasta      # Output directory for sequence upstream of the splice
OUTDIRDO=./do-fasta      # Output directory for sequence downstream of the splice
GFFDIR=./gff-subsets     # Directory with GFF files for each protein (created by get-gff-subset.sh)

mkdir -p $OUTDIRUP $OUTDIRDO

# Repeat for each gff in $GFFDIR

for INGFF in $GFFDIR/*.gff; do
  # Get the output prefix from the current gff file
  uce=$(basename $INGFF .gff)

  # Get the strand from the 7th column of the gff
  strand=`head -n 1 $INGFF | awk '{print $7}'`

  # Use the strand to determine what is upstream vs downstream.
  # The .up14 and .do6 files are gff files used with agat_sp_extract_sequences.pl
  # 
  # The sed command removes either the first or last line (depending on upstream or downstream)
  # so that data sequence is not included from the start or end of the protein because that's not
  # an intron/exon boundry

  if [[ $strand == "+" ]]; then
    sort -n $INGFF | sed '1d;' $INGFF > $INGFF.up14
    sort -n $INGFF | sed '$d;' $INGFF > $INGFF.do6
  else
    sort -n $INGFF | sed '$d;' $INGFF > $INGFF.up14
    sort -n $INGFF | sed '1d;' $INGFF > $INGFF.do6
  fi

  # Use agat_sp_extract_sequences.pl to get the fasta sequence based on the .do6 and do14 gff files
  agat_sp_extract_sequences.pl -g $INGFF.do6 -f $INFASTA -t cds --split --do 6 --output $OUTDIRDO/$uce-cds-do6.fasta >/dev/null
  agat_sp_extract_sequences.pl -g $INGFF.up14 -f $INFASTA -t cds --split --up 14 --output $OUTDIRUP/$uce-cds-up14.fasta >/dev/null

  #convert to single line fasta using seqtk
  seqtk seq -l0 $OUTDIRDO/$uce-cds-do6.fasta > $OUTDIRDO/$uce-cds-do6.fasta.tmp
  mv $OUTDIRDO/$uce-cds-do6.fasta.tmp $OUTDIRDO/$uce-cds-do6.fasta

  seqtk seq -l0 $OUTDIRUP/$uce-cds-up14.fasta > $OUTDIRUP/$uce-cds-up14.fasta.tmp
  mv $OUTDIRUP/$uce-cds-up14.fasta.tmp $OUTDIRUP/$uce-cds-up14.fasta
done
