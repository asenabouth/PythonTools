#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 10:43:35 2018

Converts FASTA sequences into GTF records 
for use with the Cell Ranger pipeline.

@author: a.senabouth
"""
# Import argument-parsing modules
import sys, getopt

# Import Biopython modules
from Bio import SeqIO
        
def parseInput(input_file):
    # Read in fasta file
    fasta_records = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))
    fasta_lines = []
    gtf_lines = []
    
    # Template strings
    fasta_header_template = ">{chromosome} dna:chromosome chromosome:GRCh38:{chromosome}:1:{length}:1 REF"
    gtf_template = """{chromosome}\thavana\tgene\t1\t{length}\t.\t+\t.\tgene_id "{chromosome}"; gene_name "{chromosome}"; gene_source "ensembl_havana"; gene_biotype "protein_coding";\n{chromosome}\thavana\ttranscript\t1\t{length}\t.\t+\t.\tgene_id "{chromosome}"; transcript_id "{chromosome}"; gene_name "{chromosome}"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "{chromosome}"; transcript_source "havana";\n{chromosome}\thavana\texon\t1\t{length}\t.\t+\t.\tgene_id "{chromosome}"; transcript_id "{chromosome}"; exon_number "1"; gene_name "{chromosome}"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; transcript_name "{chromosome}"; transcript_source "havana"; exon_id "{chromosome}_exon";"""

    # Iterate through record, create new record
    for header, record in fasta_records.items():
        seq = str(record.seq.lower())
        seq_length = len(seq)
        new_header = fasta_header_template.format(chromosome = header, length = seq_length)
        new_gtf = gtf_template.format(chromosome = header, length = seq_length)
        
        # Add to fasta_lines
        fasta_lines.append(new_header)
        fasta_lines.append(seq)
        
        # Add to gtf_lines
        gtf_lines.append(new_gtf)
    
    return(fasta_lines, gtf_lines)
    
def main(argv):
    # Put the variables here
    input_file = None
    output_filename = None
    
    # Start checks
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=", "ofile="])
    except getopt.GetoptError:
        print 'ConvertFASTAtoGTF.py -i <inputfile> -o <outputfilename>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'ConvertFASTAtoGTF.py -i <inputfile> -o <outputfilename>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            input_file = arg
        elif opt in ("-o", "--ofile"):
            output_filename = arg

    if input_file and output_filename:
        return (input_file, output_filename)
    else:
        return None, None

if __name__ == "__main__":
    # Call function to sort out arguments
    input_file, output_filename = main(sys.argv[1:])
    
    if not any([input_file, output_filename]):
        print ("Please specify an input,fasta output and gtf output files.")
        sys.exit(2)
    
    # Call function to parse the input file
    fasta_lines, gtf_lines = parseInput(input_file)
    
    # Write processed data out
    fasta_output = output_filename + ".fa"
    gtf_output = output_filename + ".gtf"

    # Add lines
    fasta_lines = [line + "\n" for line in fasta_lines]
    with open(fasta_output, "w") as fasta_handle:
        fasta_handle.writelines(fasta_lines)
    
    gtf_lines = [line + "\n" for line in gtf_lines]
    with open(gtf_output, "w") as gtf_handle:
        gtf_handle.writelines(gtf_lines)
    
    print "Conversion complete!"