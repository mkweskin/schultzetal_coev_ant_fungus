#!/usr/bin/env python3
"""
File: back_trans.py
Matthew Kweskin, kweskinm@si.edu

Description:
Takes an amino acid alignment and returns the nucleotide sequence (given in a separate fasta file).
Preserves the alignment of the original alignment

Requires:
Python3, bioconda
"""

import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils.CheckSum import seguid

def get_args():
    parser = argparse.ArgumentParser(description="""Back translates an AA alignment given a fasta of unalignmed nucleotide sequences.""")
    parser.add_argument(
        'inalignment',
        type=str,
        help="""Name of the AA alignment file"""
    )
    parser.add_argument(
        'innucleotide',
        type=str,
        help="""Name of the fasta file with the nucleotide sequence of sequences found in the alignment"""
    )
    parser.add_argument(
        '--informat',
        type=str,
        choices=["fasta","fastq","phylip","nexus","clustalw"],
        default="fasta",
        help="""Input format of the alignment (fasta, fastq, phylip, nexus, clustalw). Default: fasta"""
    )
    return parser.parse_args()

def main():
    args = get_args()

    align = AlignIO.read(args.inalignment,args.informat)
    align_length = align.get_alignment_length()

    unaligned_dict = SeqIO.to_dict(SeqIO.parse(args.innucleotide,args.informat))

    for aa_record in align:
        if aa_record.id in unaligned_dict.keys():
            nt_record = unaligned_dict[aa_record.id]
        else:
            print("ERROR: Sequence name not found in unaligned sequences")

        aa_pos=0
        nt_pos=0
        left_gaps=0
        aa_seq = aa_record.seq
        nt_seq = nt_record.seq
        nt_aligned = ""
        first_AA_found = False
        first_AA_match = False

        # Count how many gaps before the first AA
        while not first_AA_found:
            if aa_seq[aa_pos] == "-" or aa_seq[aa_pos] == "?":
                left_gaps += 1
                aa_pos += 1
            else:
                first_AA_found = True

        while aa_pos < align_length:
            if aa_seq[aa_pos] == "-" or aa_seq[aa_pos] == "?":
                nt_aligned += "---"
                aa_pos += 1
            else:
                codon = Seq(nt_seq[nt_pos:nt_pos+3])
                codon_trans = codon.translate()
                if codon_trans == aa_seq[aa_pos]:
                    nt_aligned += codon
                    if first_AA_match == False:
                        first_AA_match = True
                        aa_pos_match_start = aa_pos
                        nt_pos_natch_start = nt_pos
                    aa_pos += 1
                    nt_pos += 3
                elif first_AA_match == False:
                    nt_pos += 3
                else:
                    # Reset, this was a false start, return to previous match
                    aa_pos = aa_pos_match_start
                    nt_pos = nt_pos_natch_start + 3
                    nt_aligned = ""
                    first_AA_match = False
            #print (aa_pos, nt_pos, aa_pos_match_start, nt_pos_natch_start)
        nt_aligned = ("---" * left_gaps) + nt_aligned
        print(">"+aa_record.id)
        print(nt_aligned)

if __name__ == '__main__':
    main()
