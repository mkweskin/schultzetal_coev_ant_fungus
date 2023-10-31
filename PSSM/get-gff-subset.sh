#!/bin/bash
OUTDIR=./gff-subsets
INGFF=./Abisporus_varbisporusH97.v2.FilteredModels3.gff
INFASTA=./Abisporus_varbisporusH97.v2_AssembledScaffolds.fasta
LOCIPROT=./loci_protein

mkdir -p $OUTDIR

#Loop through the file that contains uce id and protein id
while read x; do
  uce=`echo $x |  awk '{ print $1 }'`
  proteinid=`echo $x |  awk '{ print $2 }'`
  grep "proteinId ${proteinid};" $INGFF >$OUTDIR/$uce.gff || echo "Error for $uce, proteinID $proteinid not found"
done < $LOCIPROT
