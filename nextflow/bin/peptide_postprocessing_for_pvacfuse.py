#!/usr/bin/env python

import argparse
import sys
import os
#import re
import pandas as pd
#pip install openpyxl
#import openpyxl

#if len(sys.argv)<2:
#   print("You need to specify file!")
#   sys.exit()

def peptide_postprocessing(in_file, sample_name, output_dir, median_mt_score, rank_stab, rank_bestallele, thalf, score_bestallele):
    #in_pd = pd.read_csv( infile, sep = '\t', header = 0)
    in_pd = pd.read_csv( in_file, sep = '\t', header = 0, low_memory=False)
    #in_pd.head(0)
    print(sample_name)
    print(in_file)
    in_pd.rename(columns={'median score':'median mt score'}, inplace = True)
    in_pd_f = in_pd[
    (in_pd['median mt score'] < median_mt_score) & 
    (in_pd['%rank_stab'] <= rank_stab) & 
    (in_pd['%rank_bestallele'] <= rank_bestallele) & 
    (in_pd['thalf(h)'] >= thalf) &
    (in_pd['score_bestallele'] >= score_bestallele)
    ].sort_values(by = ['score_bestallele','pred','%rank_stab'], ascending=[False,False,True])
    #(in_pd['median fold change'] >= median_fold_change) & 
    #(in_pd['score_bestallele'] >= score_bestallele) & 
    #(in_pd['gene expression'] >= gene_expression) & 
    #(in_pd['transcript expression'] >= transcript_expression) & 
    #(in_pd['tumor dna vaf'] >= tumor_dna_vaf) & 
    #(in_pd['normal vaf'] < normal_dna_vaf) & 
    #(in_pd['tumor rna vaf'] >= tumor_rna_vaf)
    #].sort_values(by = ['score_bestallele','pred','%rank_stab','gene expression'], ascending=[False,False,True,False])
    #os.getcwd()
    os.chdir(output_dir)
    #in_pd_f.drop(in_pd_f.loc[:,'mhcflurry score':'smmpmbec percentile'].columns,axis=1,inplace=True)
    in_pd_f.drop(in_pd_f.loc[:,'best percentile method':'cterm_7mer_gravy_score'].columns[1:-2],axis=1,inplace=True)
    del in_pd_f['identity']

    #in_pd_f.drop_duplicates(subset=['hla allele','mt epitope seq'], keep='first', inplace=True)
    in_pd_f.to_csv(sample_name + '_peptide_postprocessing' + '.txt', header=True, index=False, sep='\t')
    in_pd_f.to_excel(sample_name + '_peptide_postprocessing' + '.xlsx', header=True, index=False, sheet_name=sample_name)
    #in_pd_f.to_csv(output_dir + '/' + sample_name + '_peptide_postprocessing' + '.txt', header=True, index=False, sep='\t')
    #in_pd_f.to_excel(output_dir + '/' + sample_name + '_peptide_postprocessing' + '.xlsx', header=True, index=False, sheet_name=sample_name)

    #path_out = output_dir
    #patient= sample_name
    #if patient == patients[0]:
    #    in_pd_f.to_excel("all_peptide_postprocessing.xlsx", header=True, index=False, sheet_name=patient)
    #else:
    #    writer = pd.ExcelWriter(path_out + '/' + 'all_peptide_postprocessing.xlsx', mode='a', engine='openpyxl')
    #    in_pd_f.to_excel(writer,sheet_name=patient, header=True, index=False)
    #    writer.save()
    #    writer.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Peptide postprocessing based on merged results from pvacseq, prime, and netmhcstabpan")
    parser.add_argument("-i", "--in_file", required=True, help="Merged results from pvacseq, prime, and netmhcstabpan")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    parser.add_argument("--median_mt_score", type = float, default = 500, help = "Median ic50 binding affinity of the mutant epitope across all prediction algorithms used. From pvacseq [default: %(default)s]")
    parser.add_argument("--rank_stab", type = float, default = 5, help = "percent Random - percent Rank of predicted stability score to a set of 200.000 random natural 9mer peptides. From netmhcstabpan [default: %(default)s]")
    parser.add_argument("--rank_bestallele", type = float, default = 3.633, help = "Lowest percent Rank for PRIME score across the different alleles. [default: %(default)s]")
    parser.add_argument("--thalf", type = float, default = 0.39, help = "Predicted half life of the pMHC complex (in hours). From netmhcstabpan [default: %(default)s]")
    #parser.add_argument("--median_fold_change", type = float, default = 0.759, help = "Median WT Score / Median MT Score. NA if there is no WT Epitope Seq. From pvacseq [default: %(default)s]")
    parser.add_argument("--score_bestallele", type = float, default = 0.011951, help = "PRIME Score corresponding to allele with the lowest percent Rank. [default: %(default)s]")
    #parser.add_argument("--gene_expression", type = float, default = 0.878, help = "Gene expression value for the annotated gene containing the variant. From pvacseq [default: %(default)s]")
    #parser.add_argument("--transcript_expression", type = float, default = 0.365, help = "Transcript expression value for the annotated transcript containing the variant. From pvacseq [default: %(default)s]")
    #parser.add_argument("--tumor_dna_vaf", type = float, default = 0.1, help = "Tumor DNA variant allele frequency (VAF) at this position. From pvacseq [default: %(default)s]")
    #parser.add_argument("--normal_dna_vaf", type = float, default = 0.01, help = "Normal DNA variant allele frequency (VAF) at this position. From pvacseq [default: %(default)s]")
    #parser.add_argument("--tumor_rna_vaf", type = float, default = 0.01, help = "Tumor RNA variant allele frequency (VAF) at this position. From pvacseq [default: %(default)s]")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    peptide_postprocessing(args.in_file, args.sample_name, args.output_dir, args.median_mt_score, args.rank_stab, args.rank_bestallele, args.thalf, args.score_bestallele)
    print('peptide postprocessing done')
