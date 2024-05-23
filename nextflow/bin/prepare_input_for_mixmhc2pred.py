#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv) < 4:
    print("Need 3 parameters: sample name, all epitopes tsv from pvacseq, and parsed tsv from pvacseq")
    sys.exit(1)

sample_name = sys.argv[1]

all_epitopes = sys.argv[2]
pd_epitopes = pd.read_csv(all_epitopes, sep='\t', low_memory=False)
pd_epitopes = pd_epitopes.rename(columns=str.lower)


pvac_tsv = sys.argv[3]
pd_tsv = pd.read_csv(pvac_tsv, sep='\t', low_memory=False)
pd_tsv_id_pep = pd_tsv.loc[:,['wildtype_amino_acid_sequence','frameshift_amino_acid_sequence','fusion_amino_acid_sequence','index']]

pd_epitopes_pep = pd.merge(pd_epitopes, pd_tsv_id_pep, on='index', how='left')
#pd_epitopes_pep_mis = pd_epitopes_pep[pd_epitopes_pep['variant type']=='missense']
pd_epitopes_pep_cont_copy = pd_epitopes_pep[['variant type','protein position','mutation position', 'mt epitope seq','wt epitope seq','wildtype_amino_acid_sequence','frameshift_amino_acid_sequence','fusion_amino_acid_sequence']]
pd_epitopes_pep_cont = pd_epitopes_pep_cont_copy.copy()
pd_epitopes_pep_cont['context'] = "No"

#pd_epitopes_pep_cont.loc[:,'context'] = 'No'

#i=0
for row in range(pd_epitopes_pep_cont.shape[0]):
    #i += 1
    print(row)
    tsv_line = pd_epitopes_pep_cont.iloc[row]
    pp = str(tsv_line['protein position'])
    ind = int(pp.split('-')[0]) - int(tsv_line['mutation position'])
    if tsv_line['variant type']=='missense' or tsv_line['variant type']=='inframe_ins' or tsv_line['variant type']=='inframe_del':
        epi_length = len(tsv_line['wt epitope seq'])
        if tsv_line['wildtype_amino_acid_sequence'].count(tsv_line['wt epitope seq']) == 1 or "True":
            #pp = str(tsv_line['protein position'])
            #pp_limit = 1 if (int(pp.split('-')[0]) - 30) < 0 else int(pp.split('-')[0]) - 30
            #ind = tsv_line['wildtype_amino_acid_sequence'].find(tsv_line['wt epitope seq'], pp_limit)
            if (ind-3) < 0 and (ind+epi_length+3) <= len(tsv_line['wildtype_amino_acid_sequence']):
                pre = (3-ind)*'-' + tsv_line['wildtype_amino_acid_sequence'][0:ind] + tsv_line['mt epitope seq'][:3]
                suf = tsv_line['mt epitope seq'][-3:] + tsv_line['wildtype_amino_acid_sequence'][ind+epi_length:ind+epi_length+3]
                context = pre + suf
                pd_epitopes_pep_cont.loc[row,'context'] = context
                #print(row)
            elif (ind-3) >= 0 and (ind+epi_length+3) > len(tsv_line['wildtype_amino_acid_sequence']):
                pre = tsv_line['wildtype_amino_acid_sequence'][ind-3:ind] + tsv_line['mt epitope seq'][:3]
                suf = tsv_line['mt epitope seq'][-3:] + tsv_line['wildtype_amino_acid_sequence'][ind+epi_length:] + (3-len(tsv_line['wildtype_amino_acid_sequence'][ind+epi_length:]))*'-'
                context = pre + suf
                pd_epitopes_pep_cont.loc[row,'context'] = context
                #print(row)
            elif (ind-3) >= 0 and (ind+epi_length+3) <= len(tsv_line['wildtype_amino_acid_sequence']):
                pre = tsv_line['wildtype_amino_acid_sequence'][ind-3:ind] + tsv_line['mt epitope seq'][:3]
                suf = tsv_line['mt epitope seq'][-3:] + tsv_line['wildtype_amino_acid_sequence'][ind+epi_length:ind+epi_length+3]
                context = pre + suf
                pd_epitopes_pep_cont.loc[row,'context'] = context
                #print(row)
            else:
                continue
        else:
            continue
    elif tsv_line['variant type']=='FS':
        epi_length = len(tsv_line['mt epitope seq'])
        if tsv_line['frameshift_amino_acid_sequence'].count(tsv_line['mt epitope seq']) == 1 or "True":
            #ind = tsv_line['frameshift_amino_acid_sequence'].find(tsv_line['mt epitope seq'])
            if (ind-3) < 0 and (ind+epi_length+3) <= len(tsv_line['frameshift_amino_acid_sequence']):
                pre = (3-ind)*'-' + tsv_line['frameshift_amino_acid_sequence'][0:ind] + tsv_line['mt epitope seq'][:3]
                suf = tsv_line['mt epitope seq'][-3:] + tsv_line['frameshift_amino_acid_sequence'][ind+epi_length:ind+epi_length+3]
                context = pre + suf
                pd_epitopes_pep_cont.loc[row,'context'] = context
                #print(row)
            elif (ind-3) >= 0 and (ind+epi_length+3) > len(tsv_line['frameshift_amino_acid_sequence']):
                pre = tsv_line['frameshift_amino_acid_sequence'][ind-3:ind] + tsv_line['mt epitope seq'][:3]
                suf = tsv_line['mt epitope seq'][-3:] + tsv_line['frameshift_amino_acid_sequence'][ind+epi_length:] + (3-len(tsv_line['frameshift_amino_acid_sequence'][ind+epi_length:]))*'-'
                context = pre + suf
                pd_epitopes_pep_cont.loc[row,'context'] = context
                #print(row)
            elif (ind-3) >= 0 and (ind+epi_length+3) <= len(tsv_line['frameshift_amino_acid_sequence']):
                pre = tsv_line['frameshift_amino_acid_sequence'][ind-3:ind] + tsv_line['mt epitope seq'][:3]
                suf = tsv_line['mt epitope seq'][-3:] + tsv_line['frameshift_amino_acid_sequence'][ind+epi_length:ind+epi_length+3]
                context = pre + suf
                pd_epitopes_pep_cont.loc[row,'context'] = context
                #print(row)
            else:
                continue
        else:
            continue
    else:
        continue

pd_epitopes_pep_cont_2col = pd_epitopes_pep_cont[['mt epitope seq','context']]
pd_epitopes_pep_cont_2col = pd_epitopes_pep_cont_2col.drop_duplicates(ignore_index=True)
pd_epitopes_pep_cont_2col.to_csv(sample_name + '_mixmhc2pred_input_peptide_context' + '.txt', header=False, index=False, sep='\t')
