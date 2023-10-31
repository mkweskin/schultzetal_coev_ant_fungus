#!/usr/bin/env Rscript

# plot_dsts_rev.R
# Creates a histogram of the mean genetic distance for each taxon in an alignment.
# Highlights taxa that have been reversed by mafft (have "_R_" at the taxon name) in the histogram.
# This operates on a directory of alignments in fasta format, producing PNG output for each alignment and outputs summary stats for each locus to the screen
# Matthew Kweskin, kweskinm@si.edu

#*******************************************************************************************************************************
# The code for the computation of p_distance is used with permission from remove_misaligned.R:
# https://github.com/marekborowiec/dorylinae_phylogenomics/blob/master/remove_misaligned.R by Marek Borowiec
# Marek L Borowiec Convergent Evolution of the Army Ant Syndrome and Congruence in Big-Data Phylogenetics Systematic Biology,
#   Volume 68, Issue 4, July 2019, Pages 642–656 https://doi.org/10.1093/sysbio/syy088
#*******************************************************************************************************************************

# Usage:
# Rscript plot_dsts_rev.R <input_dir_1> <input_dir_2> <output_dir>
#
# <input_dir_1> directory with fasta alignment files, some of which may have been reversed by mafft with the --adjustdirection or --adhjustaccurate optoins (taxa names start with "_R_")
# <input_dir_2> directory with fasta alignment files, none of which ARE NOT reversed (i.e., --adjustdirection or --adhjustaccurate were not specified)
# <output_dir>  directory where PNGs of the histograms and a summary .tsv file are output.
#   If the output directory doesn't exist, it will be created.
#
# WARNINGS: - The two input directories should have the same names and the taxa in the alignments must match.
#           - The taxa names should match, except for the addtion of "_R_" for those that were reverse complemented by mafft.
#           - The alignment names MUST end in .fa or .fasta

# Report file description:
# file1_name:   Name of alignment file from directory1 (potentially with reversals)
# file1_taxon:  Reversed taxon in the alignment (there is one row per reversed taxon per file)
# file1_mean:   Mean pairwise genetic distance for all taxa in alignment from file1
# file1_sd:     Standard deviation of all pairwise distances for file1
# file2_name:   Alignment with NO reversals
# file2_taxon:  Name of taxon that was reversed in file1
# file2_mean:   Mean pairwise genetic distance for all taxa in alignment from file2
# file2_sd:     Standard deviation of all pairwise distances for file2
# better_file:  In which file is the reversed sequence closer to the mean (file1 means it's 'better' reversed, file2 means it's 'better' not-reversed) 
# diff_sd:      Differences in SDs between file1 and file2 for this reversed sequence. Positive means reversed is better, negative means not-reversed is better.
#               The farther this value is from 0, the greater the difference between the two files


# Required R package:
#  - seqinr

# disable scientific notation
## FROM remove_misaligned.R, Borowiec 2019
options(scipen=999)

# Load libraries
#   Use install.packages("seqinr") if seqinr is not already installed
suppressPackageStartupMessages(library("seqinr"))


## Get row means (taxon means)
## from a matrix of uncorrected p-distances
## FROM remove_misaligned.R, Borowiec 2019
get_p_distance <- function(alignment) {
  # matrix with pairwise identity
  mat <- dist.alignment(alignment, matrix="identity")
  # matrix with uncorrected p-distances
  p_mat <- mat*mat
  # means for each row
  mns <- rowMeans(as.matrix(p_mat), na.rm=T)
  
  return(mns)
}

has_revs <- function(alignment){
  # Checks the given alignment for taxa that start with "_R_"
  # Returns True if at least one taxa starts with "_R_", otherwise False
  for (name in unlist(alignment["nam"])) {
    if (startsWith(name, "_R_")){
      return(TRUE)
    }
  }
  return(FALSE)
}

plot_hist <- function(means, name){
  # Output the basic histogram
  hist(means, breaks = 50, main = name, xlab="Mean genentic distance to other taxa. Dashed blue line is mean for all taxa, Solid red line(s) are for taxa starting with _R_. Further to the right is more distant.")
}

plot_compare <- function(alignment1, alignment2, outpath, fasta_name, outfile) {
  
  # First use alignment1, this has the reversed reads prefixed with "_R_"
  # We will track the names of reversed sequences (rev_names), as well as
  # the mean (p_matrix_means_1_mean) and SD (p_matrix_means_1_sd) of alignment1

  rev_names <- c()
 
  # Prepare output to a png
  png( paste(outpath, "/", fasta_name, ".png", sep = ""), width = 1920, height = 2160)
  par(mfrow=c(2,1))
  
  # get row means of p-distance matrix
  p_matrix_means_1 <- get_p_distance(alignment1)

  # Plot histogram for alignment1
  plot_hist(p_matrix_means_1, fasta_name)

  # get the mean and standard deviation of the means for all taxa
  p_matrix_means_1_mean <- mean(p_matrix_means_1)
  p_matrix_means_1_sd <- sd(p_matrix_means_1)

  # Add vertical line for mean of all taxa
  abline(v=p_matrix_means_1_mean, col="blue", lty="dashed")

  # Add vertical line for each taxon starting with "_R_"
  for (name in unlist(alignment1["nam"])) {
    if (startsWith(name, "_R_")){
      rev_names <- append(rev_names, name)
      abline(v=p_matrix_means_1[name], col="red", lty="solid")
    }
  }
  
  # Get info from alignment2  
  if (! is.null(alignment2) ){
    # get row means of p-distance matrix
    p_matrix_means_2 <- get_p_distance(alignment2)
    plot_hist(p_matrix_means_2, outpath)

    # get the mean and standard deviation of the means
    p_matrix_means_2_mean <- mean(p_matrix_means_2)
    p_matrix_means_2_sd <- sd(p_matrix_means_2)
    
    # Add vertical line for mean
    abline(v=p_matrix_means_1_mean, col="blue", lty="dashed")
    
    # Go through the revsered taxa found in alignment 1 and output info from alignment2
    names_2 <- unlist(alignment2["nam"])
  }else{
    names_2 <- NULL
  }
  
  for (rev_name in rev_names){
    sd_rev_1 <- (p_matrix_means_1[rev_name] - p_matrix_means_1_mean) / p_matrix_means_1_sd
    rev_name_noprefix <- substring(rev_name, first = 4)
    if ( is.element(rev_name, names_2) | is.element(rev_name_noprefix, names_2)){
      if (! is.element(rev_name, names_2) ){
        name_in_fasta2 <- rev_name_noprefix
      }else{
        name_in_fasta2 <- rev_name
      }
      # Add red vertical line for the sequences that were reversed in alignment1
      abline(v=p_matrix_means_2[name_in_fasta2], col="red", lty="solid")
      
      # Number of sd away from the mean for this taxon in fasta1 and fasta2
      sd_rev_2 <- (p_matrix_means_2[name_in_fasta2] - p_matrix_means_2_mean) / p_matrix_means_2_sd
      
      # The amount in differnce in standard deviations of this taxon in the two alignments
      # Negative: file2 is "better"
      # Positive: file1 is "better"
      diff_sd <- sd_rev_2 - sd_rev_1
      if (diff_sd < 0){
        better <- "file2"
      }else{
        better <- "file1"
      }
      
      write(paste(fasta_name, rev_name, p_matrix_means_1[rev_name], sd_rev_1, name_in_fasta2, p_matrix_means_2[name_in_fasta2], sd_rev_2, better, diff_sd, sep = "\t"), outfile, append = TRUE)
    }else{
      # name not found in alignment2, print what we have
      write(paste(fasta_name, rev_name, p_matrix_means_1[rev_name], sd_rev_1, "NA", "NA", "NA", "NA", "NA", sep = "\t"), outfile, append = TRUE)
    }
      
  }
  
  ignorewarningmessage <- dev.off()
}


# Arguments are:
#   Input directory 1
#   Input directory 2
#   Output directory

args = commandArgs(trailingOnly=TRUE)
if(length(args) < 3){
  stop("Expecting three arguments: path1 path2 outpath")
  quit(save = "no")
}
path1 <- args[1]
path2 <- args[2]
outpath <- args[3]

# read in file names from dir1 (MUST end in .fasta or .fa)
fastas1 <- dir(path=path1, pattern="\\.fasta|\\.fa$")

if (length(fastas1) == 0){
  stop(paste("ERROR: No fasta files (*.fa or *.fasta) found in", path1))
  quit(save = "no")
}

# Create output dir if it doesn't exist
if (! dir.exists(outpath)){
  dir.create(outpath)
}

# Prepare report file
outfile <-paste(outpath,"/", "report.tsv", sep="")
write(paste("file1_name", "file1_taxon", "file1_mean", "file1_sd", "file2_name", "file2_taxon", "file2_mean", "file2_sd", "better_file", "diff_sd", sep = "\t"), outfile, append = TRUE)

# Loop through fasta files in path1
for (fasta in fastas1){
  alignment_1 <- read.alignment(file=paste(path1, "/", fasta, sep = ""), format="fasta")
  # Only anlayze alignments where there's a taxon starting with "_R_"
  if (has_revs(alignment_1)){
    message(paste("Found reversed sequences:", fasta))
    # Read in the alignments to an object
    if (file.exists(paste(path2,"/", fasta, sep = ""))){
      alignment_2 <- read.alignment(file=paste(path2, "/", fasta, sep = ""), format="fasta")
    }else{
      warning(paste("Warning expected file is missing from path2: ", path2, "/", fasta, sep = ""))
      warning("Output from the file from path1 with be generated")
      alignment_2 <- NULL
    }
    plot_compare(alignment_1, alignment_2, outpath, fasta, outfile)
  }else{
    message(paste("No reversed sequences found, skipping:", fasta))
  }
}
