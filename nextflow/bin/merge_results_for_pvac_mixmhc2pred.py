#!/usr/bin/env python
# coding: utf-8

import argparse
import os
#import sys
#import glob
import pandas as pd

#hla_l = []
def convert_hla(hla):
    hla_l = hla.split("_")
    if hla.count("__") == 1:
        hla_star = hla_l[0] + '*' + hla_l[1] + ':' + hla_l[2] + '-' + hla_l[4] + '*' + hla_l[5] + ':' + hla_l[6]
    elif hla.count("__") == 0:
        hla_star = hla_l[0] + '*' + hla_l[1] + ':' + hla_l[2]
    else:
        print("bad hla")
    return hla_star

file = 'patient_mixmhc2pred_output.txt'
mixmhc2pred_df = pd.read_csv(file,sep='\t',comment='#')
#df=df.iloc[:,np.r_[0:5]]
mixmhc2pred_df = mixmhc2pred_df.iloc[:,0:5]

#mixmhc2pred_df['BestAllele'] = mixmhc2pred_df['BestAllele'].apply(lambda i: ('HLA-'+ i[0] + '*' + i[1:3] + ':' + i[3:5]))
mixmhc2pred_df['BestAllele'] = mixmhc2pred_df['BestAllele'].apply(lambda i:convert_hla(i))
mixmhc2pred_df = mixmhc2pred_df.rename(columns=str.lower)
mixmhc2pred_df.columns=['mt epitope seq', 'context', 'bestallele', '%rank_best', 'core_best']
mixmhc2pred_df = mixmhc2pred_df.drop_duplicates(ignore_index=True)



dfs=[All_c, netmhc, prime]

final = reduce(lambda left,right: pd.merge(left,right,how='left',on=['mt epitope seq','hla allele']), dfs)

nan_value = float("NaN")
final.replace("", nan_value, inplace=True)
final.dropna(how='all', axis=1, inplace=True)
final.to_csv(path_merge_result_out + '/' + sample_name + '_pvacseq_prime_netmhcstabpan' + '.txt', index=False,sep='\t')
