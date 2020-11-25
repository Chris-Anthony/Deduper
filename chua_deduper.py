#!/usr/bin/env python

print()
print("Initializing...")

import argparse
import re

def get_args():
    '''Argparse code to allow for in-line command interface.'''

    parser = argparse.ArgumentParser(description="""Reads in a sorted SAM file. Remove PCR duplicates keeping the first PCR record seen.""")
    parser.add_argument('-f','--file',  action='store', nargs='?', type=str, 
                    required=True, help='Name of sorted SAM file')
    parser.add_argument('-u','--umi',  action='store', nargs='?', type=str, 
                required=False, default=False, 
                help='Enter file containing the list of UMIs (unset if randomers)')
    parser.add_argument('-p','--paired',  action='store_true', 
                required=False, help='Include if reads are paired', dest="paired")
    parser.add_argument('-s','--sort',  action='store_true', 
                required=False, help='Include if SAM file needs to be sorted', dest="sort")

    return parser.parse_args()

def checkCIGAR(CIGAR:str, POS:int, FLAG:str)->int:
    '''Checks the CIGAR string for soft clipping. Updates leftmost position as needed for the forward read
    or returns the backcalculated rightmost postion for the reverse read.'''

    # forward strand or unstranded
    soft_clip = re.search("^([0-9]+)S", CIGAR)
    
    if soft_clip:
        POS -= int(soft_clip.group(1))

        if FLAG == "reverse": # if paired is specified
            length = 0
            for x in re.findall("([0-9]+)", CIGAR):
                length += int(x)

            remove_length = 0 # removes insertions (and others) because they do not "consume" reference
            for y in re.findall("([0-9]+)[IHP=X]", CIGAR):
                remove_length += int(y)

            POS = POS + length - remove_length
        
    return POS

def main(file1:str, paired:bool, umi:str, sort:bool):
    '''Reads in a sorted SAM file. Remove PCR duplicates keeping the first PCR record seen.'''

    # if a UMI file was listed, create a dictionary with UMI barcodes
    if umi_file != False:
        barcodes = dict()

        with open(umi_file, "r") as umi_fh:
            for line in umi_fh:
                barcodes[line.rstrip("\n")] = []

    # Creates a new file with the output name the input name but with "_deduped" appended
    file_name = re.search(r"(.*)(\.sam)", file1)
    # deduped file
    output_deduped = ''.join([file_name.group(1), "_deduped", file_name.group(2)])
    out_deduped = open(output_deduped, "w")
    # duplicate file
    output_duplicate = ''.join([file_name.group(1), "_duplicate", file_name.group(2)])
    out_duplicate = open(output_duplicate, "w")
    # low quality file
    output_lowqual = ''.join([file_name.group(1), "_lowqual", file_name.group(2)])
    out_lowqual = open(output_lowqual, "w")

    # Create an empty dictionary to hold record identifiers
    uniq_dict = dict()

    # Regex to pull barcode information
    barcode = re.compile("[ATCGN]{4,}")

    # sort SAM file
    if sort == True:
        import pysam
        input_sort = ''.join([file_name.group(1), "_sorted", file_name.group(2)])
        pysam.sort("-o", input_sort, file1)
        file1 = input_sort

    # empty variables
    count_uniq, count_dup, count_lowqual = 0, 0, 0

    with open(file1, "r") as fh:
        for line in fh:
            if line[0] == "@":
                # print header lines to output file(s)
                out_deduped.writelines(line)
                out_duplicate.writelines(line)
                out_lowqual.writelines(line)
            else:
                record = line.split("\t")

                # Find UMI barcode
                QNAME = record[0]
                UMI = re.findall(barcode, QNAME)[0] # UMI barcode

                FLAG = record[1] # strandedness info
                if ((int(FLAG) & 32) == 32):
                    FLAG = "forward"
                elif ((int(FLAG) & 16) == 16):
                    FLAG = "reverse"
                    
                RNAME = record[2] # chromosome/scaffold/group

                # Find leftmost position, update as necessary
                CIGAR = record[5]
                POS = record[3] 
                POS = checkCIGAR(CIGAR, int(POS), FLAG) # updated leftmost position

                # Create a key unique to the read
                key = tuple([UMI, RNAME, POS, FLAG])

                # Write to output files
                if umi_file != False: # for UMIs
                    if UMI in barcodes.keys():
                        if key in uniq_dict.keys():
                            uniq_dict[key] += 1
                            count_dup += 1
                            out_duplicate.writelines(line)
                        else: 
                            uniq_dict[key] = 1
                            count_uniq += 1
                            out_deduped.writelines(line)
                    else:
                        count_lowqual += 1
                        out_lowqual.writelines(line)
                else: # for randomers
                    if re.search("N", UMI):
                        count_lowqual += 1
                        out_lowqual.writelines(line)
                    else:
                        if key in uniq_dict.keys():
                            uniq_dict[key] += 1
                            count_dup += 1
                            out_duplicate.writelines(line)
                        else: 
                            uniq_dict[key] = 1
                            count_uniq += 1
                            out_deduped.writelines(line)

    # Close files
    out_deduped.close()
    out_duplicate.close()
    out_lowqual.close()

    print("Number of unique reads:", count_uniq)
    print("Number of duplicate reads:", count_dup)
    print("Number of misindexed reads:", count_lowqual)

    print("Complete.")

if __name__ =="__main__":
    args = get_args()
    file1, paired, umi_file, sort = args.file, args.paired, args.umi, args.sort

    main(file1, paired, umi_file, sort)