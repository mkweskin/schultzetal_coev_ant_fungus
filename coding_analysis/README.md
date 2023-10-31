# Creating alignments of coding sequence

After we used `phyluce_probe_annotations_from_genomes` to get protein IDs of UCE that mapped to the agabi_varbur_2 genome, we used the agabi_varbur_2 protein sequences to get the coding region of our fungal samples.

`exonerate_hits.py` from HybPiper <https://github.com/mossmatters/HybPiper> (version 1.3.1) was used to run exonerate (version 2.4.0).

For each UCE that mapped to a coding region of the  genome, `exonerate_hits.py` was run separately for each sample that had a UCE sequence for that locus.

- The amino acid sequence of the protein was used as the 'proteinfile' (in `exonerate_hits.py` terminology)
- The unaligned UCE sequence was used as the 'assemblyfile' (in `exonerate_hits.py` terminology)
  - A single nucleotide sequence is given in each `exonerate_hits.py` run.
  - HybPiper expects the nucleotide fasta headers to be in a specific format. To allow the program to run, this header was used for each UCE input file: `>NODE_1_length_2770_cov_7.873250`. With the `exonerate_hits.py` options used (see below), this header information is not used.
- The additional options were given to the `exonerate_hits.py` run: `-t 65 --depth_multiplier 0 --length_pct 90`
  - `--depth_multiplier 0` this uses the coverage depth criterion which isn't available for this analysis.

After exonerate was run for each sample in a locus, the results (as amino acids) were concantenated and converted into single line fasta files.

Any stop codons (`*`) at the end of a sequence were removed from this concantenated fasta: `sed -ri 's/\*$//' UCE-XX.fasta`

Any sequences with remaining (interior) stop codons were removed from the analysis.

The first and last three amino acids were removed from the remaining sequences: `sed -r 's/^[^>]..(.*)...$/\1/' UCE-XX-nostop-trim3.fasta`

This fasta was aligned with mafft using the `--auto` setting and then ends were trimmed with `phyluce_align_get_trimmed_alignments_from_untrimmed` (trim settings:   `--proportion 0.60 --threshold 0.60 --max_divergence 0.3 --min-length 90 --window 20`)

The trimming process can leave columns in the alignment that are all gaps. These were removed with Trimal: `trimal -in <inputfile> -out <outputfile> --noallgaps`

For nucleotide analysis of the coding data, the AA alignment was preserved, but translated to the original nucleotide sequence using the script `back_trans.py`.

In pseudocode this process is:

```text
FOR each UCE locus with protein sequence
    FOR each sample with UCE data
        run exonerate_hits.py with Agabi_varbisH97_2 protein sequence as 'proteinfile' and sample's nucleotide sequence as 'assemblyfile'
    END FOR
    concatenate protein sequences of samples
    Remove stop codons from end of AA sequence
    Remove any samples that have interior stop codons
    Remove first and last three amino acids
    Align amino acids with MAFFT
    Trim ends with `phyluce_align_get_trimmed_alignments_from_untrimmed`
    Remove columns that are all gaps with `trimal`
    IF number of remaining samples < 70% of total
        Remove from further analysis
    ELSE
        Use for phylogenetic analyses
    END IF
    Back translate aligned and trimmed AA sequences to nucleotides with back_trans.py
END FOR
Use amino acid alignments for phylogenomic analyses
Use nucleotide alignments for phylogenomic anlayses
```
