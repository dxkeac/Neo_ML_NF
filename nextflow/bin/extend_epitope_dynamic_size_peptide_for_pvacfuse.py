#!/usr/bin/env python

import pandas as pd
import re
import argparse
import sys
import os
#import math


def check_peptide(mt_epitope, mt_peptide, wt_aa):
    if mt_epitope in mt_peptide and mt_peptide in wt_aa:
        return 'T'
    else:
        return 'F'

def get_specific_length_peptide(sample_name, pep_length, pd_epitopes_pep2, suffix):
    i=0
    for row in range(pd_epitopes_pep2.shape[0]):
        i += 1
        print(i)
        tsv_line = pd_epitopes_pep2.iloc[row]
        if tsv_line['variant type']=='missense' or tsv_line['variant type']=='inframe_ins' or tsv_line['variant type']=='inframe_del':
            tsv_line_inframe = tsv_line
            pd_epitopes_pep2.loc[row,'mt_peptide_val'] = check_peptide(tsv_line_inframe['mt epitope seq'], tsv_line_inframe['mt_peptide_sequence'], tsv_line_inframe['wildtype_amino_acid_sequence'])
        elif tsv_line['variant type']=='FS':
            tsv_line_fs = tsv_line
            if tsv_line_fs['frameshift_amino_acid_sequence'].count(tsv_line_fs['mt epitope seq']) == 1:
                ind = tsv_line_fs['frameshift_amino_acid_sequence'].find(tsv_line_fs['mt epitope seq'])
                fs_epitope_center_ind = ind + int(len(tsv_line_fs['mt epitope seq'])/2)
                pd_epitopes_pep2.loc[row,'mt_peptide_sequence'] = tsv_line_fs['frameshift_amino_acid_sequence'][fs_epitope_center_ind-pep_length:fs_epitope_center_ind+pep_length+1]
                pd_epitopes_pep2.loc[row,'mt_peptide_val'] = check_peptide(tsv_line_fs['mt epitope seq'], tsv_line_fs['mt_peptide_sequence'], tsv_line_fs['frameshift_amino_acid_sequence'])
            else:
                pd_epitopes_pep2.loc[row,'mt_peptide_val'] = check_peptide(tsv_line_fs['mt epitope seq'], tsv_line_fs['mt_peptide_sequence'], tsv_line_fs['frameshift_amino_acid_sequence'])
        elif tsv_line['variant type']=='inframe_fusion':
            tsv_line_inframe = tsv_line
            pd_epitopes_pep2.loc[row,'mt_peptide_val'] = check_peptide(tsv_line_inframe['mt epitope seq'], tsv_line_inframe['mt_peptide_sequence'], tsv_line_inframe['fusion_amino_acid_sequence'])
        elif tsv_line['variant type']=='frameshift_fusion':
            tsv_line_fus = tsv_line
            if tsv_line_fus['fusion_amino_acid_sequence'].count(tsv_line_fus['mt epitope seq']) == 1:
                ind = tsv_line_fus['fusion_amino_acid_sequence'].find(tsv_line_fus['mt epitope seq'])
                fs_epitope_center_ind = ind + int(len(tsv_line_fus['mt epitope seq'])/2)
                left_pos = 0 if 0 > (fs_epitope_center_ind-pep_length) else (fs_epitope_center_ind-pep_length)
                pd_epitopes_pep2.loc[row,'mt_peptide_sequence'] = tsv_line_fus['fusion_amino_acid_sequence'][left_pos:fs_epitope_center_ind+pep_length+1]
                pd_epitopes_pep2.loc[row,'mt_peptide_val'] = check_peptide(tsv_line_fus['mt epitope seq'], tsv_line_fus['mt_peptide_sequence'], tsv_line_fus['fusion_amino_acid_sequence'])
            else:
                pd_epitopes_pep2.loc[row,'mt_peptide_val'] = check_peptide(tsv_line_fus['mt epitope seq'], tsv_line_fus['mt_peptide_sequence'], tsv_line_fus['fusion_amino_acid_sequence'])
        else:
            pd_epitopes_pep2.loc[row,'mt_peptide_val'] = 'Novel variant type'

    pd_epitopes_pep2['mt_peptide_sequence_len'] = pd_epitopes_pep2['mt_peptide_sequence'].str.len()


    if 'mhcflurry wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['mhcflurry wt score', 'mhcflurry mt score', 'mhcflurry wt percentile', 'mhcflurry mt percentile'], axis=1,inplace=True)

    if 'mhcnuggetsi wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['mhcnuggetsi wt score', 'mhcnuggetsi mt score'], axis=1,inplace=True)

    if 'netmhccons wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['netmhccons wt score', 'netmhccons mt score', 'netmhccons wt percentile', 'netmhccons mt percentile'], axis=1,inplace=True)

    if 'netmhcpan wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['netmhcpan wt score', 'netmhcpan mt score', 'netmhcpan wt percentile', 'netmhcpan mt percentile'], axis=1,inplace=True)

    if 'pickpocket wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['pickpocket wt score', 'pickpocket mt score', 'pickpocket wt percentile', 'pickpocket mt percentile'], axis=1,inplace=True)

    if 'netmhc wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['netmhc wt score', 'netmhc mt score', 'netmhc wt percentile', 'netmhc mt percentile'], axis=1,inplace=True)

    if 'smm wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['smm wt score', 'smm mt score', 'smm wt percentile', 'smm mt percentile'], axis=1,inplace=True)

    if 'smmpmbec wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['smmpmbec wt score', 'smmpmbec mt score', 'smmpmbec wt percentile', 'smmpmbec mt percentile'], axis=1,inplace=True)

    #if 'mt_peptide_val' in pd_epitopes_pep2.columns:
        #pd_epitopes_pep2.drop(['mt_peptide_val'], axis=1,inplace=True)
    if 'mhcnuggetsii wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['mhcnuggetsii wt score', 'mhcnuggetsii mt score', 'mhcnuggetsii wt percentile', 'mhcnuggetsii mt percentile'], axis=1,inplace=True)

    if 'netmhciipan wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['netmhciipan wt score', 'netmhciipan mt score', 'netmhciipan wt percentile', 'netmhciipan mt percentile'], axis=1,inplace=True)

    if 'nnalign wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['nnalign wt score', 'nnalign mt score', 'nnalign wt percentile', 'nnalign mt percentile'], axis=1,inplace=True)

    if 'smmalign wt score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['smmalign wt score', 'smmalign mt score', 'smmalign wt percentile', 'smmalign mt percentile'], axis=1,inplace=True)

    if 'cterm_7mer_gravy_score' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['cterm_7mer_gravy_score', 'max_7mer_gravy_score', 'difficult_n_terminal_residue', 'c_terminal_cysteine', 'c_terminal_proline', 'cysteine_count','n_terminal_asparagine','asparagine_proline_bond_count'], axis=1,inplace=True)

    if 'index' in pd_epitopes_pep2.columns:
        pd_epitopes_pep2.drop(['index', 'id', 'wildtype_amino_acid_sequence', 'frameshift_amino_acid_sequence', 'fusion_amino_acid_sequence'], axis=1,inplace=True)

    pd_epitopes_pep2.to_csv('pvacfuse_'+ sample_name + '_flank' + str(pep_length) + 'aa_peptide_' + str(suffix) + '.txt', header=True, index=False, sep='\t')
    pd_epitopes_pep2.to_excel('pvacfuse_'+ sample_name + '_flank' + str(pep_length) + 'aa_peptide_' + str(suffix) + '.xlsx', header=True, index=False, sheet_name=sample_name)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get peptides of a specific length")
    parser.add_argument("-e", "--pvacfuse_all_epitopes_tsv", required=True, help="all_epitopes.tsv from pvacfuse run")
    parser.add_argument("-f", "--pvacfuse_fasta_tsv", required=True, help="manufacturability.tsv from pvacfuse generate_protein_fasta")
    parser.add_argument("-t", "--pvacfuse_parse_tsv", required=True, help="sample_name.tsv from pvacfuse run")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    parser.add_argument("-l", "--peptide_length", required=True, help="specific length")
    parser.add_argument("-x", "--output_suffix", required=True, help="suffix for output file")
    #parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    #all_epitopes='xxx_tumor_DNA_all_epitopes_HLA_I.tsv'
    all_epitopes = args.pvacfuse_all_epitopes_tsv
    pd_epitopes = pd.read_csv(all_epitopes, sep='\t', low_memory=False)
    pd_epitopes.rename(columns={'Mutation':'Index'}, inplace = True)
    pd_epitopes = pd_epitopes.rename(columns=str.lower)

    #fasta = 'xxx_flank_13_phased.fasta.manufacturability.tsv'
    fasta = args.pvacfuse_fasta_tsv
    pd_fasta = pd.read_csv(fasta, sep='\t', low_memory=False)
    pd_fasta_id_pep = pd_fasta.iloc[:,[0,1]]
    #pd_fasta_id_pep_mt_wt = pd.merge(pd_fasta_id_pep_mt,pd_fasta_id_pep_wt,on='id')

    pvac_tsv = args.pvacfuse_parse_tsv
    pd_tsv = pd.read_csv(pvac_tsv, sep='\t', low_memory=False)
    pd_tsv_id_pep = pd_tsv.loc[:,['wildtype_amino_acid_sequence','frameshift_amino_acid_sequence','fusion_amino_acid_sequence','index']]

    #pd_epitopes.drop(pd_epitopes.loc[:,'median wt percentile':'cterm_7mer_gravy_score'].columns[1:-2],axis=1,inplace=True)
    pd_epitopes_pep = pd.merge(pd_epitopes,pd_fasta_id_pep, left_on='index', right_on='id', how='left')
    pd_epitopes_pep = pd_epitopes_pep.rename(columns=str.lower)
    pd_epitopes_pep2 = pd.merge(pd_epitopes_pep,pd_tsv_id_pep, on='index', how='left')
    pd_epitopes_pep2.rename(columns={'peptide_sequence':'mt_peptide_sequence','epitope seq':'mt epitope seq'}, inplace = True)

    pd_epitopes_pep2['mt_peptide_val'] = "F"
    #pd_epitopes_pep2['wt_peptide_val'] = "F"
    pep_length = int(args.peptide_length)

    get_specific_length_peptide(args.sample_name, pep_length, pd_epitopes_pep2, args.output_suffix)
    print('specific length peptide processing done')