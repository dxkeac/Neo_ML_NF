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

def peptide_postprocessing(in_file, sample_name, output_dir, median_mt_score, rank_stab, rank_bestallele, thalf, median_fold_change, score_bestallele, gene_expression, transcript_expression, tumor_dna_vaf, normal_dna_vaf, tumor_rna_vaf):
    in_pd = pd.read_csv( in_file, sep = '\t', header = 0, low_memory=False)
    #in_pd.head(0)
    print(sample_name)
    print(in_file)

    #median_8m
    in_pd_fc = in_pd.loc[(1-in_pd['median wt score'].isnull()).astype(bool),:]
    #bool_series = pd.notnull(data["Gender"])
    in_pd_fc_f = in_pd_fc[
    (in_pd_fc['tumor dna vaf'] >= tumor_dna_vaf) &
    (in_pd_fc['normal vaf']  < 0.02) &
    (in_pd_fc['tumor rna vaf'] >= 0.01) &
    (in_pd_fc['gene expression'] >= 0.6) &
    (in_pd_fc['transcript expression'] >= 0.06) &
    (in_pd_fc['median mt score'] <= 500) &
    (in_pd_fc['median mt percentile'] <= 2) &
    (in_pd_fc['median wt score'] >= 6) &
    (in_pd_fc['median fold change'] >= 0.5) &
    (in_pd_fc['pred'] >= 0.5) &
    (in_pd_fc['thalf(h)'] >= 0.1) &
    (in_pd_fc['%rank_stab'] <= 5.5) &
    (in_pd_fc['%rank_bestallele'] <= 2.1) &
    (in_pd_fc['score_bestallele'] >= 0.01) &
    (in_pd_fc['bigmhc'] >= 0.01)
    ]

    #in_pd_no_mis = in_pd[in_pd['variant type'].str.contains('inframe_ins|inframe_del|FS')]
    in_pd_no_fc = in_pd.loc[in_pd['median wt score'].isnull(),:]
    in_pd_no_fc_f = in_pd_no_fc[
    (in_pd_no_fc['tumor dna vaf'] >= tumor_dna_vaf) &
    (in_pd_no_fc['normal vaf']  < 0.02) &
    (in_pd_no_fc['tumor rna vaf'] >= 0.01) &
    (in_pd_no_fc['gene expression'] >= 0.6) &
    (in_pd_no_fc['transcript expression'] >= 0.06) &
    (in_pd_no_fc['median mt score'] <= 500) &
    (in_pd_no_fc['median mt percentile'] <= 2) &
    (in_pd_no_fc['pred'] >= 0.5) &
    (in_pd_no_fc['thalf(h)'] >= 0.1) &
    (in_pd_no_fc['%rank_stab'] <= 5.5) &
    (in_pd_no_fc['%rank_bestallele'] <= 2.1) &
    (in_pd_no_fc['score_bestallele'] >= 0.01) &
    (in_pd_no_fc['bigmhc'] >= 0.01)
    ]

    in_pd_f = pd.concat([in_pd_fc_f,in_pd_no_fc_f])
    in_pd_f['rank_prime'] = in_pd_f['score_bestallele'].rank(method='first',ascending=False)
    in_pd_f['rank_bigmhc'] = in_pd_f['bigmhc'].rank(method='first',ascending=False)
    in_pd_f['rank_median_fold_change'] = in_pd_f['median fold change'].rank(method='first',ascending=False)

    in_pd_f['min_rank_final'] = in_pd_f[['rank_bigmhc','rank_prime','rank_median_fold_change']].min(axis=1)
    #in_pd_f['min_rank_final'] = in_pd_f[['rank_bigmhc','rank_prime','rank_median_fold_change']].median(axis=1)
    #in_pd_fs = in_pd_f.sort_values(by = ['%rank_bestallele','pred','%rank_stab','gene expression'], ascending=[True,False,True,False])
    in_pd_fs = in_pd_f.sort_values(by = ['min_rank_final','rank_prime','rank_bigmhc','rank_median_fold_change'], ascending=[True,True,True,True])

    os.chdir(output_dir)
    #in_pd_f.drop(in_pd_f.loc[:,'mhcflurry wt score':'smmpmbec mt percentile'].columns,axis=1,inplace=True)
    #in_pd_f.drop(in_pd_f.loc[:,'median wt percentile':'cterm_7mer_gravy_score'].columns[1:-2],axis=1,inplace=True)
    del in_pd_fs['identity']

    in_pd_fs.to_csv(sample_name + '_peptide_postprocessing' + '_AF' + str(tumor_dna_vaf) + '.txt', header=True, index=False, sep='\t')
    in_pd_fs.to_excel(sample_name + '_peptide_postprocessing' + '_AF' + str(tumor_dna_vaf) + '.xlsx', header=True, index=False, sheet_name=sample_name)
    #in_pd_f.to_csv(output_dir + '/' + sample_name + '_peptide_postprocessing' + '.txt', header=True, index=False, sep='\t')
    #in_pd_f.to_excel(output_dir + '/' + sample_name + '_peptide_postprocessing' + '.xlsx', header=True, index=False, sheet_name=sample_name)

    #path_out = output_dir
    #patient= sample_name
    #if patient == patients[0]:
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
    parser.add_argument("--median_fold_change", type = float, default = 0.759, help = "Median WT Score / Median MT Score. NA if there is no WT Epitope Seq. From pvacseq [default: %(default)s]")
    parser.add_argument("--score_bestallele", type = float, default = 0.011951, help = "PRIME Score corresponding to allele with the lowest percent Rank. [default: %(default)s]")
    parser.add_argument("--gene_expression", type = float, default = 0.878, help = "Gene expression value for the annotated gene containing the variant. From pvacseq [default: %(default)s]")
    parser.add_argument("--transcript_expression", type = float, default = 0.365, help = "Transcript expression value for the annotated transcript containing the variant. From pvacseq [default: %(default)s]")
    parser.add_argument("--tumor_dna_vaf", type = float, default = 0.1, help = "Tumor DNA variant allele frequency (VAF) at this position. From pvacseq [default: %(default)s]")
    parser.add_argument("--normal_dna_vaf", type = float, default = 0.01, help = "Normal DNA variant allele frequency (VAF) at this position. From pvacseq [default: %(default)s]")
    parser.add_argument("--tumor_rna_vaf", type = float, default = 0.01, help = "Tumor RNA variant allele frequency (VAF) at this position. From pvacseq [default: %(default)s]")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    peptide_postprocessing(args.in_file, args.sample_name, args.output_dir, args.median_mt_score, args.rank_stab, args.rank_bestallele, args.thalf, args.median_fold_change, args.score_bestallele, args.gene_expression, args.transcript_expression, args.tumor_dna_vaf, args.normal_dna_vaf, args.tumor_rna_vaf)
    print('peptide postprocessing done')
