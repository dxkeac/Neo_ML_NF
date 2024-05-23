#!/usr/bin/env python

import argparse
import sys
import re

def combine_hla_ii(in_file, out_file):
    hla_list = []
    comb = []
    comb_ab = []
    with open(in_file, "r") as infile:
        for line in infile:
            hla = line.strip()
            hla_list.append(hla)
                #print(hla1)

    hla_list_copy = hla_list
    for hla in hla_list:
        for hla_copy in hla_list_copy:
        #print(hla2)
            comb.append(hla + "-" + hla_copy)

    for comb_hla in comb:
        if re.search( r'DPA.*?DPB', comb_hla, re.I):
        #print(comb_hla)
        #if "DPA1" in comb_hla and "DPB1" in comb_hla:
            comb_ab.append(comb_hla)
        #elif "DQA1" in comb_hla and "DQB1" in comb_hla:
        elif re.search( r'DQA.*?DQB', comb_hla, re.I):
            comb_ab.append(comb_hla)
        elif comb_hla.count("HLA") == 1 and "D" not in comb_hla:
            comb_ab.append(comb_hla)
        else:
            pass
    with open(out_file, "w") as outfile:
        #hla_ii_no_dups = list(dict.fromkeys(hla_ii))
        outfile.write("\n".join(comb_ab))



if __name__ == "__main__":

    #usage = __doc__.split("\n\n\n")
    parser = argparse.ArgumentParser(description="Combine alpha and beta chain for HLA-II alleles")
    parser.add_argument("--hla_ii_list", required=True, help="HLA II list one per line")
    parser.add_argument("--sample_name", required=True, help="Sample name")
    parser.add_argument("--output_dir", required=True, help="Path to the output directory")
    #parser.add_argument("--supported_list", required=True, help="Supported list")
    #parser.add_argument("--model_list", required=True, help="Secondary list")
    args = parser.parse_args()
    # Parse arguments
    in_file = args.hla_ii_list
    out_file = args.output_dir + args.sample_name + "_hla_ii_combine.txt"

    combine_hla_ii(in_file, out_file)
