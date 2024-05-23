#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys
import glob
import pandas as pd
import numpy as np

def merge_results(sample_name, input_path, prime_path, netmhcstabpan_path, mixmhcpred_path, prime, netmhcstabpan):

    pep = glob.glob(input_path + '/prime/*_pep.txt')
    hla = glob.glob(input_path + '/prime/*_hla.txt')
    hla_= pd.read_csv(hla[0],sep='\t',header=None)

    # Prime test (only for MHC-I)

    if prime ==True:        
        print('running prime')
        for index, row in hla_.iterrows():
            iname = "".join(row[0].replace('HLA-','').split('*'))
            run = prime_path + '/PRIME -i '+ input_path + '/prime/' + sample_name  +'_pep.txt -o '+ input_path + '/prime/' + sample_name + '_' + iname +'_out.txt -a '+ iname + ' -mix ' + mixmhcpred_path + '/MixMHCpred'
            os.system(run)
            size = os.path.getsize(input_path + '/prime/' + sample_name + '_' + iname +'_out.txt')
            #10byte
            if size < 10:
                rm_empty = 'rm -rf ' + input_path + '/prime/' + sample_name + '_' + iname + '_out.txt'
                os.system(rm_empty)
        # Compare data

        print('merging prime')
        files = glob.glob(input_path + '/prime/*_out.txt')
        for file in files:
            df = pd.read_csv(file,sep='\t',comment='#')
            #df=df.iloc[:,np.r_[0:5]]
            df = df.iloc[:,0:5]
            if file==files[0]:
                DF=df 
            else:
                #DF=DF.append(df)
                DF=pd.concat([DF,df])
        DF.to_csv(input_path +'/prime/'+ sample_name + '_prime_out_all.txt', header=True, index=False,sep='\t')

    # netMHCstabpan
    if netmhcstabpan ==True:    
        print('runing netmhcstabpan')
        for index, row in hla_.iterrows():
            iname = "".join(row[0].split('*'))
            hla = iname[0:7]+':'+iname[7:9]

            files = glob.glob(input_path+'/netMHCstabpan/'+'*.txt')
            for i in files:
                l = i.replace('.txt','').replace(input_path,'').replace('/netMHCstabpan/','').split('_')[2]
                run = netmhcstabpan_path + '/netMHCstabpan -p -a ' + hla + ' -l ' + l + ' '+ i +' > '+ input_path + '/netMHCstabpan/' + sample_name + '_len_' + l + '_' + iname + '_out'
                os.system(run)


        # merge netMHCstabpan

        print('merging netmhcstabpan')
        files = glob.glob(input_path + '/netMHCstabpan/' + '*_out')
        for i in files:
            #df = pd.read_csv(i, comment='#', skipfooter=5, engine='python')
            df = pd.read_csv(i, comment='#', header = None, skipfooter=5, engine='python')
            df = df[4:]
            df=df.drop_duplicates()
            #df=df.iloc[2:]
            #df=df.iloc[:,0].str.split(" ", n=100, expand=True)
            #df=df.iloc[:,[4,6,11,17,23,30,39,44]]
            df.columns=['pos HLA peptide Identity Pred Thalf(h) %Rank_Stab BindLevel']
            if i==files[0]:
                DF=df
            else:
                #DF=DF.append(df)
                DF=pd.concat([DF,df])
        DF.to_csv(input_path + '/netMHCstabpan/' + sample_name + '_netmhcstabpan_out_all.txt', header=True, index=False,sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Running prime and netmhcstabpan")
    parser.add_argument("-i", "--in_file_dir", required=True, help="Output path of prepare_input_args.py")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample name")
    #parser.add_argument("-o", "--output_dir", required=False, help="Path to the output directory")
    parser.add_argument("--prime_path", required=True, help = "Path to PRIME-2.0. (xxx/software/PRIME-2.0)")
    parser.add_argument("--mixmhcpred_path", required=True, help = "Path to MixMHCpred-2.2. (xxx/software/MixMHCpred-2.2)")
    parser.add_argument("--netmhcstabpan_path", required=True, help = "Path to netMHCstabpan-1.0. (xxx/software/netMHCstabpan-1.0)")
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    merge_results(args.sample_name, args.in_file_dir, args.prime_path, args.netmhcstabpan_path, args.mixmhcpred_path, True, True)
    print('Running prime and netmhcstabpan done')
