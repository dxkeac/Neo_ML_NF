#!/usr/bin/env python

import pandas as pd
import re
import argparse
import sys
import os
#import math

def merge_hla_i_ii_peptide(hla_i_pep, hla_ii_pep, sample_name, pep_length, suffix):

    pd_hla_i_pep = pd.read_csv(hla_i_pep, sep='\t', low_memory=False)
    pd_hla_ii_pep = pd.read_csv(hla_ii_pep, sep='\t', low_memory=False)

    #pd_epitopes.drop(pd_epitopes.loc[:,'median wt percentile':'cterm_7mer_gravy_score'].columns[1:-2],axis=1,inplace=True)
    #pd_epitopes_pep= pd.merge(pd_epitopes,pd_fasta_id_pep_mt_wt, left_on='index', right_on='id', how='left')
    #pd_epitopes_pep = pd_epitopes_pep.rename(columns=str.lower)

    pd_hla_i_ii_pep = pd.merge(pd_hla_i_pep, pd_hla_ii_pep, on=['chromosome', 'start', 'stop', 'reference', 'variant', 'transcript'], how='left')
    #pd_file1_2= pd.merge(pd_file1, pd_file2, on=['chromosome', 'start', 'stop', 'reference', 'variant', 'transcript'], how='inner')

    pd_hla_i_pep['hla_ii_pep'] = "F"
    i = 0
    for row in range(pd_hla_i_ii_pep.shape[0]):
        i += 1
        print(i)
        tsv_line = pd_hla_i_ii_pep.iloc[row]
        if 'D' in str(tsv_line['hla allele_y']):
            pd_hla_i_pep.loc[row,'hla_ii_pep'] = 'T'

    pd_hla_i_ii_pep_v = pd.merge(pd_hla_i_pep, pd_hla_ii_pep, on=['chromosome', 'start', 'stop', 'reference', 'variant', 'transcript'], how='left')

    pd_hla_i_pep.to_csv('pvacseq_'+ sample_name + '_' + str(pep_length) + 'aa_peptide_' + str(suffix) + '.txt', header=True, index=False, sep='\t')
    pd_hla_i_pep.to_excel('pvacseq_'+ sample_name + '_' + str(pep_length) + 'aa_peptide_' + str(suffix) + '.xlsx', header=True, index=False, sheet_name=sample_name)
    pd_hla_i_ii_pep_v.to_csv('pvacseq_'+ sample_name + '_' + str(pep_length) + 'aa_peptide_' + str(suffix) + '_verbose.txt', header=True, index=False, sep='\t')
    pd_hla_i_ii_pep_v.to_excel('pvacseq_'+ sample_name + '_' + str(pep_length) + 'aa_peptide_' + str(suffix) + '_verbose.xlsx', header=True, index=False, sheet_name=sample_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge hla i and ii peptides")
    parser.add_argument("-i", "--pvacseq_hla_i", required=True, help="hla i peptide file")
    parser.add_argument("-t", "--pvacseq_hla_ii", required=True, help="hla ii peptide file")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    parser.add_argument("-l", "--peptide_length", required=True, help="specific length")
    parser.add_argument("-x", "--output_suffix", required=True, help="suffix for output file")
    #parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    pep_length = int(args.peptide_length)

    merge_hla_i_ii_peptide(args.pvacseq_hla_i, args.pvacseq_hla_ii, args.sample_name, pep_length, args.output_suffix)
    print('merging hla i and ii peptides done')