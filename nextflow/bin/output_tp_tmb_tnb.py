#!/usr/bin/env python

import argparse
import sys
import os
import pandas as pd


def TP_TMB_TNB(in_file_TP, in_file_TMB, in_file_NEO, sample_name, output_dir):
    pd_in_tp = pd.read_csv( in_file_TP, sep = ',', header = 0, low_memory=False)
    pd_in_tmb = pd.read_csv( in_file_TMB, sep = ',', header = 0, low_memory=False)
    pd_in_neo = pd.read_csv( in_file_NEO, sep = '\t', header = 0, low_memory=False)
    #print(sample_name)
    tp = round(pd_in_tp['Purity'].tolist()[0],2)
    tmb = round(pd_in_tmb['somatic.rate.ontarget'].tolist()[0],2)
    callable_bases = pd_in_tmb['callable.bases.ontarget'].tolist()[0]
    tnb = round(pd_in_neo.shape[0]/callable_bases*1000000,2)
    dic = {"Tumor purity":[tp],
           "TMB":[tmb],
           "TNB":[tnb]}
    pd_out = pd.DataFrame(dic)
    #pd_out = pd_out.T
    os.chdir(output_dir)
    pd_out.to_csv(sample_name + '_TP_TMB_TNB' + '.txt', header=True, index=False, sep='\t')
    pd_out.to_excel(sample_name + '_TP_TMB_TNB' + '.xlsx', header=True, index=False, sheet_name=sample_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Outputting TP (tumor purity), TMB (tumor mutation burden), and TNB (tumor neoantigen burden)")
    parser.add_argument("-p", "--in_file_TP", required=True, help="Input file from purecn")
    parser.add_argument("-m", "--in_file_TMB", required=True, help="Input file from purecn")
    parser.add_argument("-n", "--in_file_NEO", required=True, help="Input file from peptide_postprocessing")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    #args.in_file: '/home/dxk/data/peptide_analysis/'+ patient + '_pvacseq_prime_netmhcstabpan.txt'
    TP_TMB_TNB(args.in_file_TP, args.in_file_TMB, args.in_file_NEO, args.sample_name, args.output_dir)
    print('Outputting TP, TMB, and TNB done')
