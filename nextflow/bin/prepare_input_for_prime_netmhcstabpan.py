#!/usr/bin/env python
# coding: utf-8

import os
import argparse
import sys
import glob
import json,csv
import pandas as pd
#import numpy as np

def prepare_input(path_pvacseq_out, sample_name, path_out):
    print(sample_name)
    for name in ['/deepimmuno','/prime','/netMHCstabpan']:
        if not os.path.exists(path_out+ name):
            os.makedirs(path_out + name)

    # read pvacseq file
    file = glob.glob(path_pvacseq_out + '/' + '*all_epitopes.tsv')
    print(file)
    All = pd.read_csv(file[0],sep='\t',low_memory=False)
    #All[All['HLA Allele'].str.contains('HLA-A|HLA-B|HLA-C')]
    #All_c = All.rename(columns=str.lower)
    All_c = All[All['HLA Allele'].str.contains('HLA-A|HLA-B|HLA-C')].rename(columns=str.lower)

    # copy into All_c
    #All_c=All
    All_c['hla allele'] = All_c['hla allele'].str.replace(':','')
    All_c.rename(columns={'epitope seq':'mt epitope seq'}, inplace = True)
    if 'peptide length' in All_c.columns:
        print("peptide length")
    else:
        All_c['peptide length'] = All_c['mt epitope seq'].str.len()

    # check binding_stability by prime - prepare input
    hla = All_c['hla allele'].astype('category').cat.categories.values
    hla = pd.DataFrame(hla)
    All_c['peptide length'] = All_c['peptide length'].astype(int)
    #sub = All_c[All_c['hla allele']==All_c['hla allele'].astype('category').cat.categories[0]][All_c['peptide length']<15][All_c['peptide length']>7]['mt epitope seq']
    sub = All_c[(All_c['peptide length']<15) & (All_c['peptide length']>7)]['mt epitope seq'].astype('category').cat.categories.values
    sub = pd.DataFrame(sub)
    sub.to_csv(path_out + '/prime/' + sample_name +'_pep.txt', header=False, index=False)
    hla.to_csv(path_out + '/prime/' + sample_name +'_hla.txt', header=False, index=False)

    # check binding_stability by deepimmuno - prepare input
    deep = All_c[All_c['peptide length']==9|10].iloc[:,[18,14]]
    deep.to_csv(path_out + '/deepimmuno/' + sample_name +'_pep_hla.csv', header=False, index=False)

    # check netMHCstabpan - prepare input
    #sub = All_c[All_c['hla allele']==All_c['hla allele'].astype('category').cat.categories[0]]
    #for l in sub['peptide length'].astype('category').cat.categories.values:
    for l in All_c['peptide length'].astype('category').cat.categories.values:
        sub_pep = All_c[All_c['peptide length']==l]['mt epitope seq'].astype('category').cat.categories.values
        sub_pep = pd.DataFrame(sub_pep)
        sub_pep.to_csv(path_out + '/netMHCstabpan/' + sample_name +'_len_'+ str(l) + '.txt', header=False, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare input files for prime and netmhcstabpan")
    parser.add_argument("-i", "--in_file_dir", required=True, help="Path including *all_epitopes.tsv file from pvacseq")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    #parser.add_argument("--median_mt_score", type = float, default = 500, help = "Median ic50 binding affinity of the mutant epitope across all prediction algorithms used. From pvacseq [default: %(default)s]")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    prepare_input(args.in_file_dir, args.sample_name, args.output_dir)
    print('Preparing input files done')