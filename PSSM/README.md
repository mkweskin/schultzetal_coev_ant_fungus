# Scripts used to create a custom PSSM file for exonerate

These scripts were used to calculate base frequences on the 3' and 5' sides of an intron/extron boundry.

**These scripts require gnu sed and grep. If you have a Mac, install gnu sed with: conda install -c conda-forge --override-channels sed grep**

## Get separate GFFs for each protein part of the UCE dataset

`get-gff-subset.sh`: uses a list of UCE loci and their associated protein code as input. Creates a separate gff file for each protein with a separate line for each exon of that protein.

## Get fasta sequence from splice sites

`pos-strand.sh`: uses the gff files from `get-gff-subset.sh` to produce fasta sequences upstream and downstream of the exon/intron splices.

This requires:

- AGAT: <https://github.com/NBISweden/AGAT>  Dainat J. AGAT: Another Gff Analysis Toolkit to handle annotations in any GTF/GFF format.  
(Version v0.7.0). Zenodo. <https://www.doi.org/10.5281/zenodo.3552717>
- seqtk: <https://github.com/lh3/seqtk>  Shen W, Le S, Li Y, Hu F (2016) SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE 11(10): e0163962. <https://doi.org/10.1371/journal.pone.0163962>

## Get upstream/downstream bases

These `egrep` regex expresseions get the last 9 nucleotides of the downstream fastas and the first 15 of the upstream fastas

```bash
egrep --no-filename -o "[ACGT]{9}$" do-fasta/*.fasta > down-stream-9bases
egrep --no-filename -o "^[ACGT]{15}" up-fasta/*.fasta > up-stream-15bases
```

The percentage of each nucleotide at each of the positions was manually calculated to create the PSSM file for exonerate.
