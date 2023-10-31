# `phyluce_probe_annotations_from_genomes`

## Installing (for phyluce 1.6.7)

### Create a phyluce 1.6.7 conda environment

```bash
(base)$ conda create -c conda-forge -c bioconda --override-channels -n phyluce167_annot phyluce=1.6.7 gffutils=0.9
(base)$ conda activate phyluce167_annot
```

### Download and install `phyluce_probe_annotations_from_genomes`

```basj
(phyluce167_annot)$ curl -O ... #TODO: add url to the script
(phyluce167_annot)$ cp phyluce_probe_annotations_from_genomes $CONDA_PREFIX/bin/ 
(phyluce167_annot)$ chmod a+x $CONDA_PREFIX/bin/phyluce_probe_annotations_from_genomes
```

### Test `phyluce_probe_annotations_from_genomes`

```bash
(phyluce167_annot)$ phyluce_probe_annotations_from_genomes
usage: phyluce_probe_annotations_from_genomes [-h] --conf CONF --lastz LASTZ
                                              --output OUTPUT
                                              [--name-pattern PATTERN]
                                              [--probe-prefix PROBE_PREFIX]
                                              [--probe-regex PROBE_REGEX]
                                              [--exclude EXCLUDE [EXCLUDE ...]]
                                              [--verbosity {INFO,WARN,CRITICAL}]
                                              [--contig_orient]
                                              (--flank FLANK | --probes PROBES)
phyluce_probe_annotations_from_genomes: error: argument --conf is required
```

## Running `phyluce_probe_annotations_from_genomes`

### Convert reference to `.2bit`

See the [Phyluce tutorial-4](https://phyluce.readthedocs.io/en/latest/tutorials/tutorial-4.html#convert-genomes-to-2bit-format) for this step

### Map probes to genome

```bash
(phyluce167_annot)$ phyluce_probe_run_multiple_lastzs_sqlite \
    --db agabi_varbish97_2.sqlite \
    --output matched_fungal_genome-lastz \
    --scaffoldlist agabi_varbish97_2 \
    --genome-base-path $PWD/genomes \
    --probefile ./agaricales-v1-attine-cultivar-SPECIFIC-probe-list-FINAL.fasta \
    --cores 1 \
    --coverage 67 \
    --identity 65
```

### Retreive annotations

Requires a conf file that has the path to the scaffolds in 2bit format (`[scaffolds]` section) and the gff annotation file (`[gffs]`), see `genomes.conf`

```bash
(phyluce167_annot)$ phyluce_probe_annotations_from_genomes \
    --lastz matched_fungal_genome-lastz \
    --conf genomes.conf \
    --flank 0 \
    --name-pattern "agaricales-v1-attine-cultivar-SPECIFIC-probe-list-FINAL.fasta_v_{}.lastz.clean" \
    --output features

(phyluce167_annot)$ head -n 1 features/agabi_varbish97_2.gff
scaffold_10	JGI	exon	526047	526419	.	+	.	transcriptId "153259"; name "fgenesh2_pm.10_#_144"; uce uce-500
```

The resulting file is a gff that only includes annotations from the regions where probes matched. The UCE locus is appended to column 9, the attributes (`; uce uce-500`)

*All annotations are including regardless of the feature type (gff column 3)*

### Generate list of UCE and Protein

The output of this section is a table of UCE loci with the JGI protein ID. This is used to retreive the peptide sequence of the CDS in downstreams steps.

Get the `CDS` annotations from the gff file

```bash
(phyluce167_annot)$ cd features
(phyluce167_annot)$ grep CDS agabi_varbish97_2.gff > agabi_varbur_2_harvest.features-0bp.CDS
```

**The following steps rely on specifics of the `Abisporus_varbisporusH97.v2.FilteredModels3.gff` file and will likely not work with other gff files**

In the CDS annotations of `Abisporus_varbisporusH97.v2.FilteredModels3.gff`, e.g.

```bash
head -n 1 agabi_varbur_2_harvest.features-0bp.CDS
scaffold_10	JGI	CDS	526047	526419	.	+	0	exonNumber "1"; name "fgenesh2_pm.10_#_144"; proteinId "153259"; uce uce-500
```

(Note: in the following steps gnu sed is required. If you are using a Mac, install gnu sed with: `conda install -c conda-forge --override-channels sed`)

For further processing we are changing the CDS file so:

1. The uce ID is moved to the first column:

```bash
(phyluce167_annot)$ sed -r -i 's/(.*) (uce-[0-9]+)$/\2\t\1/' agabi_varbur_2_harvest.features-0bp.CDS
```

1. The protein ID (without quotes is now the final column)

```bash
(phyluce167_annot)$ sed -r -i 's/proteinId "([0-9]+)".*/proteinId "\1"\t\1/' agabi_varbur_2_harvest.features-0bp.CDS
```

The result:

```bash
(phyluce167_annot)$ head -n 1 agabi_varbur_2_harvest.features-0bp.CDS
uce-500	scaffold_10	JGI	CDS	526047	526419	.	+	0	exonNumber "1"; name "fgenesh2_pm.10_#_144"; proteinId "153259"	153259
```

#### Generate list of UCE with their protein ID

In this step, a sqlite database, `CDS.sqlite`, is created that contains the CDS data.

This database is queried for UCE that have more than one proteinid that were found. These UCE are not included in downstream analyses of coding regions.

```bash
(phyluce167_annot)$ ../get_uniques.sh agabi_varbur_2_harvest.features-0bp.CDS
Creating sqlite db with CDS info
Exporting all the loci into the text file all_loci_CDS
Finding loci where >1 proteinid was found
uce-281:        2 ProteinIDs found, removing this locus from the database
Exporting UNIQUE loci into the text file loci_protein
Cleaning up intermediate files (import.sql, export-all.sql, all_loci_CDS, export-unique.sql)
Output: loci_protein (Lists each UCE locus with one protein ID found, and its protein ID)
```

### Get protein sequence from JGI protein ID

The JGI assembly file `Abisporus_varbisporus_H97.v2.FilteredModels3.proteins.fasta` gives the protein sequences of the *A. bisporus var. bisporus* genome.

The fasta headers include the protein ID. E.g. `>jgi|Agabi_varbisH97_2|189137|estExt_fgenesh2_kg.C_10001`

The proteins fasta was converted to single-line fasta format so that grep could be used to output the protein sequence:

```bash
mkdir -p baits
while read x; do
  uce=$(echo $x | awk '{print $1}')
  protid=$(echo $x | awk '{print $2}')
  grep -A 1 "|$protid|" Abisporus_varbisporusH97.v2.FilteredModels3.proteins.singleline.fasta >baits/$uce.fasta
done <loci_protein
```
