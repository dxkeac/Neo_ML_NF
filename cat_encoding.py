import numpy as np
import pandas as pd
import argparse
import sys
import os
from numpy import unique
from collections import Counter

from sklearn.preprocessing import FunctionTransformer
import warnings

#from Utils.GlobalParameters import *

#warnings.filterwarnings(action='ignore', category=UserWarning)


class CatEncoder:
    def __init__(self, feature):
        self.feature = feature
        self.float_encoding = {}
        self.int_encoding = {}
        self.unknown_float_code = -1
        self.unknown_int_code = -1
        self.nr_classes = 0
        return
    def fit(self, x, y):
        pos_cnt = Counter([x_sel for x_sel, y_sel in zip(x, y) if y_sel == 1])
        tot_cnt = Counter(x)
        self.float_encoding = {}
        if len(tot_cnt.keys()) == 2 and (True in tot_cnt.keys()) and (False in tot_cnt.keys()):
            self.float_encoding[True] = 1.0
            self.float_encoding[False] = 0.0
            self.unknown_float_code = 0.0
        else:
            s = 0.0
            for l in tot_cnt:
                if l in pos_cnt:
                    self.float_encoding[l] = pos_cnt[l] / tot_cnt[l]
                    s += self.float_encoding[l]
                else:
                    self.float_encoding[l] = 0.0
            for l in tot_cnt:
                self.float_encoding[l] /= s
            self.unknown_float_code = 0.0
        self.float_encoding = dict(sorted(self.float_encoding.items(), key=lambda item: item[1]))
        self.nr_classes = len(self.float_encoding.keys())
        self.encode_int_from_float()
        return self
    def encode_int_from_float(self):
        self.int_encoding = {}
        idx = 1
        for label in self.float_encoding:
            self.int_encoding[label] = idx
            idx += 1
        self.unknown_int_code = 0
    def encode_float(self, label):
        if label in self.float_encoding:
            return self.float_encoding[label]
        else:
            return self.unknown_float_code
    def encode_int(self, label):
        if label in self.int_encoding:
            return self.int_encoding[label]
        else:
            return self.unknown_int_code
    def transform(self, values, cat_type='float'):
        if cat_type == 'float':
            return list(map(self.encode_float, values))
        elif cat_type == 'int':
            return list(map(self.encode_int, values))
        else:
            return values
    def get_nr_classes(self):
        return self.nr_classes
    def get_unknown(self, cat_type='float'):
        if cat_type == 'float':
            return self.unknown_float_code
        elif cat_type == 'int':
            return self.unknown_int_code
        else:
            return np.nan
    def append_to_file(self, encoding_file):
        for c in self.float_encoding.keys():
            encoding_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.feature, self.nr_classes, c, self.float_encoding[c]))
        encoding_file.write("{0}\t{1}\t{2}\t{3}\n".format(self.feature, self.nr_classes, 'unknown', self.unknown_float_code))
    def get_encoding(self, encoding_df):
        df = encoding_df[encoding_df['Feature'] == self.feature]
        self.float_encoding = {}
        for idx in df.index:
            cat = df.loc[idx, 'Category']
            if cat == 'unknown':
                self.unknown_float_code = df.loc[idx, 'Value']
            elif cat == 'False':
                self.float_encoding[False] = df.loc[idx, 'Value']
            elif cat == 'True':
                self.float_encoding[True] = df.loc[idx, 'Value']
            elif cat == 'nan':
                self.float_encoding[np.nan] = df.loc[idx, 'Value']
            else:
                self.float_encoding[cat] = df.loc[idx, 'Value']
        self.nr_classes = len(self.float_encoding.keys())
        self.encode_int_from_float()
    @staticmethod
    def read_cat_encodings(dataset: str, peptide_type: str):
        #encoding_file = GlobalParameters().get_cat_to_num_info_file(dataset, peptide_type)
        encoding_file = get_cat_to_num_info_file(dataset, peptide_type)
        encoding_df = pd.read_csv(encoding_file, header=0, sep="\t", comment='#')
        features = encoding_df['Feature'].unique()
        encoders = {}
        for f in features:
            encoder = CatEncoder(f)
            encoder.get_encoding(encoding_df)
            encoders[f] = encoder
        return encoders
    @staticmethod
    def get_file_header():
        return "Feature\tNr_Categories\tCategory\tValue"


#@staticmethod
def get_cat_to_num_info_file(dataset: str, peptide_type: str):
    if dataset in datasets_encoding:
        return cat_to_num_info_files[peptide_type][dataset]
    else:
        return None

def load_filter_data(data: pd.DataFrame, patient: str, dataset: str, response_types: list) -> list:
    idx = np.full(data.shape[0], True)
    if patient != "":
        idx = np.logical_and(idx, data['patient'] == patient)
    elif dataset != "":
        if dataset == 'NCI_train':
            idx = np.logical_and(idx, (data['dataset'] == 'NCI') & (data['train_test'] == 'train'))
        elif dataset == 'NCI_test':
            idx = np.logical_and(idx, (data['dataset'] == 'NCI') & (data['train_test'] == 'test'))
        else:
            idx = np.logical_and(idx, data['dataset'] == dataset)
    response_types = set(response_types)
    if 0 < len(response_types) < 3:
        idx = np.logical_and(idx, data.response_type.apply(lambda row: row in response_types))
    return data.loc[idx, :]

'''
parser = argparse.ArgumentParser(description='Preprocess of neopep and mutation data')
parser.add_argument('-pt', '--peptide_type', type=str, choices=peptide_types,
                    help='Peptide type (mutation  or neopep)')
parser.add_argument('-ds', '--dataset', type=str, choices=datasets_encoding,
                    help='Dataset used for encoding (NCI_train or NCI)')

if __name__ == "__main__":
    args = parser.parse_args()
    # GlobalParameters has predefined filenames for NCI and NCI_train datasets
    with open(get_cat_to_num_info_file(args.dataset, args.peptide_type), mode='w') as encoding_file:
        #for arg in vars(args):
            #encoding_file.write(f"#{arg}={getattr(args, arg)}\n")
            #print(f"{arg}={getattr(args, arg)}")
        #data_train = DataManager.load_filter_data(peptide_type=args.peptide_type, dataset=args.dataset,
                                                  #response_types=['CD8', 'negative'])
        encoding_file.write(CatEncoder.get_file_header() + "\n")
        data_train = args.dataset
        #ml_features = ml_features_neopep if args.peptide_type == 'neopep' else ml_features_mutation
        ml_features = feature_types_encode
        y = np.array(data_train.response_type.apply(lambda row: int(row == 'CD8')), dtype=int)
        for f in data_train.columns:
            if f in ml_features and (data_train[f].dtype.name == 'category' or data_train[f].dtype == bool):
                print("Encoding feature {0} ...".format(f))
                l_enc = CatEncoder(f)
                l_enc.fit(data_train[f].values, y)
                l_enc.append_to_file(encoding_file)
'''

feature_types_encode = {
    'Clonality': 'category',
    'mut_is_binding_pos': 'bool',
    'bestWTMatchType_I': 'category',
    'mutation_driver_statement_Intogen': 'category',
    'gene_driver_Intogen': 'category',
    'seq_len': 'category'
}

cat_to_num_info_files = {
        'neopep': {'NCI_train': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_neopep_NCI_train.txt'),
                   'NCI': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_neopep_NCI_all.txt')},
        'mutation': {'NCI_train': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_mutation_NCI_train.txt'),
                     'NCI': os.path.join(data_dir, 'cat_encoding', 'Cat_to_num_info_mutation_NCI_all.txt')}
    }

data_dir = '/data/NeoRanking_self/'
datasets_encoding = ['NCI', 'NCI_train']
data_file = '/data/NeoRanking_self/Neopep_data_org_bigmhc_mhcflurry_prime_mixmhcpred.txt'
data_df = pd.read_csv(data_file, header=0, sep="\t", comment='#')
#data_df[(data_df['dataset'] == 'NCI') & (data_df['train_test'] == 'train')]
#data = data_df
patient = ''
dataset = 'NCI_train'
peptide_type = 'neopep'
#response_types = ['CD8', 'negative', 'not_tested']
response_types = ['CD8', 'negative']
max_netmhc_rank = 20

data = load_filter_data(data_df, patient, dataset, response_types)
data = data.loc[data.mutant_rank_netMHCpan.apply(lambda r: r < max_netmhc_rank)]
#data.dtypes
#data_train[f].dtype.name == 'category'
data = data.astype(feature_types_encode)
#pd_epitope['HLA Allele'].astype('category').cat.categories.values
#df_patient_fill_norm_comb = df_patient_fill_norm_comb.astype({'seq_len': 'object'})
y = np.array(data.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)

with open(get_cat_to_num_info_file(dataset, peptide_type), mode='w') as encoding_file:
    #for arg in vars(args):
        #encoding_file.write(f"#{arg}={getattr(args, arg)}\n")
        #print(f"{arg}={getattr(args, arg)}")
    #data_train = DataManager.load_filter_data(peptide_type=args.peptide_type, dataset=args.dataset,
                                              #response_types=['CD8', 'negative'])
    encoding_file.write(CatEncoder.get_file_header() + "\n")
    data_train = data
    #ml_features = ml_features_neopep if args.peptide_type == 'neopep' else ml_features_mutation
    ml_features = feature_types_encode
    y = np.array(data_train.response_type.apply(lambda row: int(row == 'CD8')), dtype=int)
    for f in data_train.columns:
        if f in ml_features and (data_train[f].dtype.name == 'category' or data_train[f].dtype == bool):
            print("Encoding feature {0} ...".format(f))
            l_enc = CatEncoder(f)
            l_enc.fit(data_train[f].values, y)
            l_enc.append_to_file(encoding_file)