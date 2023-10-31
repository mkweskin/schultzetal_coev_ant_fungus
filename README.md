# schultzetal_coev_ant_fungus

This contains scripts used in the paper "The coevolution of ant-fungus agriculture" (working title), Schultz, et al., submitted.

Please cite this paper if you use these scripts.

```text
.
├── PSSM
│   ├── README.md
│   ├── get-gff-subset.sh
│   └── pos-strand.sh
├── README.md
├── build_protein_annot_table.R
├── coding_analysis
│   ├── README.md
│   └── back_trans.py
├── plot_dists_rev.R
└── probe_annotations_from_genomes
    ├── README.md
    ├── get_uniques.sh
    ├── phyluce_probe_annotations_from_genomes
    └── phyluce_probe_annotations_from_genomes_1.7.3
```

- `PSSM/`: Contains scripts used for calculating a custom Position Specific Score Matrix (PSSM) for use with exonerate. See the [README.md](PSSM/README.md) in that directory for more details.
- `README.md`: This file.
- `build_protein_annot_table.R`: Used to build a table of UCE and metadata for the associated metadata.
- `coding_analysis`: Process of using exonerate to protein coding sequence from nucleotide UCE data. Includes steps to create coding alignments. See the [README.md](coding_analysis/README.md) in that directory for more details.
- `plot_dists_rev.R`: Compares alignments where the mafft feature `--adjustdirection` was enabled and disabled. It produces histograms of pairwise genetic distances of the aligned sequences highlighting sequences reverse complemented by mafft. Requires the R package [seqinr](https://seqinr.r-forge.r-project.org). R 4.0.3 and seqinr 4.2_5 were used in the publication.
- `probe_annotations_from_genomes/`: Contains `phyluce_probe_annotations_from_genomes` a modified version of `phyluce_probe_slice_sequence_from_genomes` from Phyluce 1.6.7 that, given a gff annotation file, returns annotations where UCE loci were mapped. In addition to phyluce 1.6.7, requires gffutils (version 0.9 for Python 2.7 used). See the [README.md](probe_annotations_from_genomes/README.md) in that directory for more details. The file `phyluce_probe_annotations_from_genomes_1.7.3` is the same script for phyluce 1.7.3.
