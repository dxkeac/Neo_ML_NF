#!/usr/bin/env python
# coding: utf-8

### merge results from different softwares output

import argparse
import os
import sys
import glob
import json,csv
import pandas as pd
import numpy as np
from functools import reduce


def merge_result(sample_name, path_pvacseq_out, path_prime_netmhc_out, path_bigmhc_out, path_merge_result_out):
    print(sample_name)
    # read file
    file = glob.glob(path_pvacseq_out + '/' + '*all_epitopes*.tsv')
    All = pd.read_csv(file[0],sep='\t', low_memory=False)
    #All= All.rename(columns=str.lower)
    All_c = All[All['HLA Allele'].str.contains('HLA-A|HLA-B|HLA-C')].rename(columns=str.lower)
    #All = All.values.tolist ()
    All_c.rename(columns={'epitope seq':'mt epitope seq'}, inplace = True)

    # read prime results
    prime = pd.read_csv(path_prime_netmhc_out + '/prime/' + sample_name +'_prime_out_all.txt',sep='\t', low_memory=False)
    prime['BestAllele'] = prime['BestAllele'].apply(lambda i: ('HLA-'+ i[0] + '*' + i[1:3] + ':' + i[3:5]))
    prime = prime.rename(columns=str.lower)
    prime.columns=['mt epitope seq', '%rank_bestallele', 'score_bestallele','%rankbinding_bestallele', 'hla allele']
    prime = prime.drop_duplicates(ignore_index=True)
    #prime = prime.values.tolist ()

    # read netMHCstabpan results
    netmhc = pd.read_csv(path_prime_netmhc_out + '/netMHCstabpan/' + sample_name +'_netmhcstabpan_out_all.txt', comment='#', sep='\t',low_memory=False)
    if netmhc.shape[1] < 2 :
        for i in range(10):
            netmhc[netmhc.keys()[0]] = netmhc[netmhc.keys()[0]].str.replace('  ',' ')
        netmhc2 = pd.DataFrame(netmhc[netmhc.keys()[0]].str.split(' ').tolist())
        netmhc = netmhc2.iloc[:,range(2,8)]
    #netmhc.columns=['pos','hla allele','mt epitope seq','identity', 'pred', 'thalf(h)', '%rank_stab', 'bindlevel']
    netmhc.columns=['hla allele','mt epitope seq','identity', 'pred', 'thalf(h)', '%rank_stab']
    netmhc = netmhc.drop_duplicates(ignore_index=True)
    #netmhc = netmhc.values.tolist ()

    # read bigmhc results
    bigmhc = pd.read_csv(path_bigmhc_out,sep=',', low_memory=False)
    #bigmhc = bigmhc[["mhc","pep","bigmhc_im"]]
    bigmhc = bigmhc.iloc[:,[0,1,4]]
    bigmhc.columns=['hla allele', 'mt epitope seq', 'bigmhc']
    bigmhc = bigmhc.drop_duplicates(ignore_index=True)
    #bigmhc = bigmhc.values.tolist ()

    dfs=[All_c, netmhc, prime, bigmhc]
    final = reduce(lambda left,right: pd.merge(left,right,how='left',on=['mt epitope seq','hla allele']), dfs)
    nan_value = float("NaN")
    final.replace("", nan_value, inplace=True)
    final.dropna(how='all', axis=1, inplace=True)
    final.to_csv(path_merge_result_out + '/' + sample_name + '_pvacseq_prime_netmhcstabpan_bigmhc' + '.txt', index=False,sep='\t')

    #DFs.append(final.values.tolist())

#s = json.dumps(DFs)
#with open(path_+'tesla_results_all_patients.json', "w") as text_file:
#    text_file.write(s)  

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge results from pvacseq, prime, and netmhcstabpan")
    #parser.add_argument("-i", "--in_file", required=True, help="Merge results from pvacseq, prime, and netmhcstabpan")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    #parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    parser.add_argument("--path_pvacseq_out", required=True, help="Path including *all_epitopes.tsv file from pvacseq")
    parser.add_argument("--path_prime_netmhcstabpan_out", required=True, help="Path including output file from prime and netmhcstabpan")
    parser.add_argument("--bigmhc_out", required=True, help="The output file from bigmhc")
    parser.add_argument("--path_merge_result_out", required=True, help="Path to the output directory")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    #args.in_file: '/home/dxk/data/peptide_analysis/'+ patient + '_pvacseq_prime_netmhcstabpan.txt'
    merge_result(args.sample_name, args.path_pvacseq_out, args.path_prime_netmhcstabpan_out, args.bigmhc_out, args.path_merge_result_out)
    print('Merging results done')