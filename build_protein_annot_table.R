# Takes a table of UCE loci and the JGI protein ID and returns NCBI protein information (protein ID, interproscan annotations)
# This was run in R 4.3.1 on macOS and these version from tidyverse:
# ── Attaching core tidyverse packages ───────────────────────────────────────────────────────────────────────────── tidyverse 2.0.0 ──
# ✔ dplyr     1.1.3     ✔ readr     2.1.4
# ✔ forcats   1.0.0     ✔ stringr   1.5.0
# ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
# ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
# ✔ purrr     1.0.2     
#
# Matthew Kweskin, kweskinm@si.edu


library(tidyverse)

# loci_protein is from my pipeline that matches probes to the reference and pullls out protein IDs from the gff file (from JGI)
#  locus:           UCE locus ID
#  gff_protein_id:  Proetin from the gff
loci_protein <- read.delim("loci_protein", header = FALSE, col.names = c("locus", "gff_protein_id"))

# Locus.tag is for merging with NCBI data
loci_protein <- loci_protein %>% mutate(Locus.tag = paste("AGABI2DRAFT_", gff_protein_id, sep = ""))

# proteins_858_31528.csv is from NCBI's genome DB from Agaricus bisporus bispours H97
# https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/858/31528%7CAgaricus%20bisporus%20var.%20bisporus%20H97/Un/
ncbi_protein_df <- read.csv("proteins_858_31528.csv")

# Merge the ncbi data with the protein code form the gff
loci_protein_ncbi <- loci_protein %>% left_join(y = ncbi_protein_df, by = c("Locus.tag" = "Locus.tag"))

# The number of lines with N/A from the NCBI protein info
cat(paste("This is the number of protein IDs from UCEs that don't match a protein ID in NCBI: ", (loci_protein_ncbi %>% filter(is.na(Protein.Name)) %>% count())))
# There are 7

# And print those lines
cat("The protein IDs not in ncbi are:")
loci_protein_ncbi %>% filter(is.na(Protein.Name)) %>% select(locus, gff_protein_id)

cat(paste("Lines where the NCBI protein name *isn\'t* the boring \"hypothetical protein...\"", loci_protein_ncbi %>% filter(str_starts(Protein.Name, "hypothetical protein", negate = TRUE)) %>% count()))
# There are only 233 like this

# Export list of NCBI protein IDs to a text doc, protein_list.csv
cat("writing \'protein_list.csv\' which has a list of NCBI protein IDs...")
loci_protein_ncbi %>% select(Protein.product) %>%  drop_na() %>% write_csv(file = "protein_list.csv", col_names = FALSE)

######
# Get protein annotations from NCBI
######

# Upload list to Batch Entrez: https://www.ncbi.nlm.nih.gov/sites/batchentrez
# Search in the 'protein' database.
# Then download all the records as "GenPept" format (proteins.gb)
#
# 2023-10-31: The batch entrez search returned two errors, XP_006457855.1 and XP_006454128.1
#             search for these manually and also download them in genpept format. 
#
# Import .gb files into Geneious and then export the documents as .gff3 (choose "Don't export sequences" and "Export all qualifiers as attributes...")
# Run this on the command line to capture the annotation information that's important:
#   egrep -v "Protein|Site|source|^#" "2 documents from annots.gff" | cut -d'  ' -f 1,3,9 | sort -u >proteins-annotations.tsv
#   Note: this is a tab character in the cut command ------------------------/\ 
#   Note: "2 documents from annots.gff" is my Geneious export file name

#####
# Extract IPR codes from CDS
#####

# Read in the annotations that came from NCBI->Geneious->egrep filtering->proteins-annotations.tsv
prot_annot <- read_delim(file = "proteins-annotations.tsv", col_names = c("Protein.product","Type","Notes"), delim = "\t")

# Split the Notes column by the text "InterPro:" which proceeds each IPR code
# With separate(), I give up to eight IPR_* column names. If there are more than 8 "InterPro:" annotations, you would need to list more of these IPR_ column names 
CDS_IPR <- prot_annot %>% filter(Type == "CDS", grepl("InterPro", Notes)) %>% separate(Notes, into = c("ignore", "IPR_1", "IPR_2", "IPR_3", "IPR_4", "IPR_5", "IPR_6", "IPR_7", "IPR_8", "rest"), sep="InterPro:", extra = "merge")

# Check that all the ipr codes were found (`rest` should be all NA)
if(CDS_IPR %>%  filter(! is.na(rest)) %>% count() != 0) 
  stop("Error: There were additional IPR codes found that have not been caught by `separate`")
  cat("All the IPR codes were found... proceed")

# cleanup the CDS_IPR table...
# Drop "ignore" and "rest" columns
CDS_IPR <- CDS_IPR %>% select(-c(ignore, rest))
# Merge back all the CDS_IPR parts
CDS_IPR <- CDS_IPR %>% unite('IPR_merged', IPR_1, IPR_2, IPR_3, IPR_4, IPR_5, IPR_6, IPR_7, IPR_8)
# Remove everything after the last ";"
CDS_IPR <- CDS_IPR %>% mutate(IPR_merged = str_replace(IPR_merged, "^(.*);.*$", "\\1"))
# Replace ";db_xref=_IPR..." with the IPR code and a comma (this creates the comma separated list of IPR codes)
CDS_IPR <- CDS_IPR %>% mutate(IPR_merged = str_replace_all(IPR_merged, ";db_xref=.(IPR[0-9]+)", ", \\1"))
# There's sometimes some extra codes after the last ";", this removes that
CDS_IPR <- CDS_IPR %>% mutate(IPR_merged = str_replace(IPR_merged, ";.*$", ""))

# Instead of a comma separated listed, create a table (CDS_IPR_each) with each IPR as a separate row
CDS_IPR_each <- data.frame(Protein.product = vector(mode="character"), Type = vector(mode="character"), IPR = vector(mode="character"))
for (i in 1:nrow(CDS_IPR)){
  if(! grepl(",", CDS_IPR[i, ]$IPR_merged) ){
    CDS_IPR_each[nrow(CDS_IPR_each) + 1,] <- CDS_IPR[i, ]
  } else{
    IPRs <- unlist(strsplit(CDS_IPR[i, ]$IPR_merged, ", ", fixed = TRUE))
    for (j in 1:length(IPRs)){
      CDS_IPR_each[nrow(CDS_IPR_each) + 1,] <- c(CDS_IPR[i, ]$Protein.product, CDS_IPR[i, ]$Type, IPRs[j])
    }
  }
}

# Write out CDS_IPR_each for use with downstream analyses
write_delim(file = "CDS_IPR_each.tsv", CDS_IPR_each, delim = "\t")

#####
# Extract note="..." info from NCBI annotations
#####

# Example: "...note=similar to hypothetical protein CC1G_03415 [Coprinopsis cinerea okayama7#130]..."
CDS_Notes <- prot_annot %>% filter(Type == "CDS", grepl("note=", Notes)) %>% separate(Notes, into = c("ignore", "note_1", "rest"), sep="note=", extra = "merge")

# Check that all the ipr codes were found (`rest` should be all NA)
if(CDS_Notes %>%  filter(! is.na(rest)) %>% count() != 0) 
  stop("Error: There were additional note tags found that have not been caught by `separate`")
cat("All notes were found, proceed...")

# Clean up CDS_Notes
# Drop "ignore" and "rest" columns
CDS_Notes <- CDS_Notes %>% select(-c(ignore, rest))
# Remove everything after the last ";"
CDS_Notes <- CDS_Notes %>% mutate(note_1 = str_replace(note_1, ";.*$", ""))


#####
# Add Interpro IPR info (type and name)
#####

# Download interpro tsv table of codes and definitions
entry_list <- read_delim(file="https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list", delim = "\t")
# Merge the IPR data for each gene product with the interproscan info for that code
CDS_IPR_description_each <- left_join( CDS_IPR_each, entry_list, by = c("IPR" = "ENTRY_AC") )

# Add a column that has the IPR and description combined
CDS_IPR_description_each <- CDS_IPR_description_each %>%
  mutate(descriptionCombined = paste(IPR, " (", ENTRY_TYPE, ": ", ENTRY_NAME, ")", sep = "")) %>% 
  mutate(descriptionCombined = ifelse(is.na(ENTRY_TYPE), IPR, descriptionCombined))

# Make a new df that has each protein, and a formatted list of the IPR descriptions
CDS_IPR_description_summary <- CDS_IPR_description_each %>%
  group_by(Protein.product) %>%
  mutate(allDescriptions = paste0(descriptionCombined, collapse = ", ")) %>%
  distinct(Protein.product, allDescriptions)

# Add this to the table of loci (loci_protein_ncbi)
loci_protein_ncbi <- loci_protein_ncbi %>%
       # remove the ".1" from the end of the protein product ID to merge with CDS_IPR_description_summary
       mutate(Protein.product.merge = str_replace(Protein.product, "\\.\\d$", "")) %>%
       left_join(CDS_IPR_description_summary, by = c("Protein.product.merge" = "Protein.product"))

# Output this as the final table
write_delim(loci_protein_ncbi %>% mutate(locus_id = str_replace(locus, "^(uce|UCE)-", "")) %>% select(locus, locus_id, gff_protein_id, Protein.product, Protein.Name, allDescriptions), file = "loci_protein_IPR.tsv", delim = "\t")