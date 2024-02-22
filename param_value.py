# coding=gbk

import numpy as np
import pandas as pd
import random
from collections import Counter
from typing import Final
from sklearn.preprocessing import QuantileTransformer, StandardScaler, PowerTransformer, MinMaxScaler, FunctionTransformer
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from sklearn.svm import SVC
from sklearn.metrics import make_scorer
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import f1_score
from sklearn.metrics import auc
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.model_selection import cross_val_score
from sklearn.ensemble import RandomForestClassifier
import joblib
from catboost import CatBoostClassifier
from xgboost import XGBClassifier
from sklearn.exceptions import UndefinedMetricWarning
from hyperopt import hp, fmin, tpe, rand, STATUS_OK, Trials, space_eval
from hyperopt.pyll import scope
from multiprocessing import Pool
import datetime
from os import path
import os
import pickle
import time
import re
import glob
from filelock import FileLock
import warnings
warnings.filterwarnings(action='ignore', category=UndefinedMetricWarning)
warnings.filterwarnings(action='ignore', category=UserWarning)
warnings.filterwarnings(action='ignore', category=RuntimeWarning)
warnings.filterwarnings("ignore")

from sklearn.utils.extmath import softmax
class RidgeClassifier(RidgeClassifier):
    def predict_proba(self, X):
        d = self.decision_function(X)
        d_2d = np.c_[-d, d]
        return softmax(d_2d)

"""
ml_features_neopep: Final[list] = \
    ['CCF', 'Clonality', 'rnaseq_TPM', 'rnaseq_alt_support', 'CSCAPE_score',
     'mutant_other_significant_alleles', 'mutant_rank', 'mutant_rank_PRIME',
     'mutant_rank_netMHCpan', 'Sample_Tissue_expression_GTEx',
     'GTEx_all_tissues_expression_mean', 'TCGA_Cancer_expression',
     'gene_driver_Intogen', 'nb_same_mutation_Intogen',
     'mutation_driver_statement_Intogen', 'bestWTMatchScore_I',
     'bestWTMatchOverlap_I', 'bestMutationScore_I', 'bestWTMatchType_I',
     'bestWTPeptideCount_I', 'mut_Rank_Stab', 'mut_netchop_score_ct',
     'TAP_score', 'mut_is_binding_pos', 'mut_binding_score', 'mut_aa_coeff',
     'seq_len', 'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab', 'DAI_MixMHC_mbp']
#len(ml_features_neopep)
#31
ml_normalize_features_neopep: Final[list] = \
    ['CCF', 'rnaseq_TPM', 'rnaseq_alt_support', 'CSCAPE_score',
     'mutant_rank', 'mutant_rank_PRIME',
     'mutant_rank_netMHCpan', 'Sample_Tissue_expression_GTEx',
     'GTEx_all_tissues_expression_mean', 'TCGA_Cancer_expression',
     'nb_same_mutation_Intogen',
     'bestWTMatchScore_I',
     'bestWTMatchOverlap_I', 'bestMutationScore_I',
     'bestWTPeptideCount_I', 'mut_Rank_Stab', 'mut_netchop_score_ct',
     'TAP_score', 'mut_binding_score', 'mut_aa_coeff',
     'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab', 'DAI_MixMHC_mbp']
#len(ml_normalize_features_neopep)
#24
"""

ml_features_neopep: Final[list] = \
    ['rnaseq_TPM', 'CSCAPE_score', 'mutant_rank', 'mutant_rank_PRIME','mutant_rank_netMHCpan', 'mut_Rank_Stab', 'bigmhc_score', 'mhcflurry_presentation_score',
     'TAP_score', 'DAI_NetMHC', 'DAI_MixMHC', 'DAI_NetStab', 'score_prime', 'score_mixmhcpred_mt']
ml_normalize_features_neopep = ml_features_neopep
#, 'foreignness_score'

#set(ml_features_neopep).difference(set(ml_normalize_features_neopep))
#{'seq_len', 'bestWTMatchType_I', 'mutation_driver_statement_Intogen', 'Clonality', 'mut_is_binding_pos', 'gene_driver_Intogen', 'mutant_other_significant_alleles'}
#'mutant_other_significant_alleles'特征没有标准化处理

feature_types_neopep: Final[dict] = {
    'foreignness_score' : 'float64',
    'score_prime' : 'float64',
    'score_mixmhcpred_mt': 'float64',
    'bigmhc_score' : 'float64',
    'mhcflurry_presentation_score': 'float64',
    'patient': 'str',
    'dataset': 'category',
    'train_test': 'category',
    'response_type': 'category',
    'Nb_Samples': 'str',
    'Sample_Tissue': 'str',
    'Cancer_Type': 'str',
    'chromosome': 'str',
    'genomic_coord': 'int64',
    'ref': 'str',
    'alt': 'str',
    'gene': 'str',
    'protein_coord': 'int32',
    'aa_mutant': 'category',
    'aa_wt': 'category',
    'mutant_seq': 'str',
    'wt_seq': 'str',
    'pep_mut_start': 'int8',
    'TumorContent': 'float64',
    'CCF': 'float64',
    'Clonality': 'category',
    'Zygosity': 'category',
    'mutation_type': 'category',
    'mutant_rank': 'float64',
    'mutant_rank_netMHCpan': 'float64',
    'mutant_rank_PRIME': 'float64',
    'mut_Rank_Stab': 'float64',
    'TAP_score': 'float64',
    'mut_netchop_score_ct': 'float64',
    'mut_binding_score': 'float64',
    'mut_is_binding_pos': 'bool',
    'mut_aa_coeff': 'float64',
    'DAI_NetMHC': 'float64',
    'DAI_MixMHC': 'float64',
    'DAI_NetStab': 'float64',
    'mutant_other_significant_alleles': 'int8',
    'DAI_MixMHC_mbp': 'float64',
    'rnaseq_TPM': 'float64',
    'rnaseq_alt_support': 'float64',
    'GTEx_all_tissues_expression_mean': 'float64',
    'Sample_Tissue_expression_GTEx': 'float64',
    'TCGA_Cancer_expression': 'float64',
    'bestWTMatchScore_I': 'float64',
    'bestWTMatchOverlap_I': 'float64',
    'bestMutationScore_I': 'float64',
    'bestWTPeptideCount_I': 'int32',
    'bestWTMatchType_I': 'category',
    'CSCAPE_score': 'float64',
    'nb_same_mutation_Intogen': 'int32',
    'mutation_driver_statement_Intogen': 'category',
    'gene_driver_Intogen': 'category',
    'seq_len': 'category'
}

ml_feature_mv_neopep: Final[dict] = {
    'foreignness_score' : 'min',
    'bigmhc_score' : 'min',
    'mhcflurry_presentation_score': 'min',
    'mutant_rank': 'max',
    'mutant_rank_netMHCpan': 'max',
    'mutant_rank_PRIME': 'max',
    'mut_Rank_Stab': 'max',
    'TAP_score': 'min',
    'mut_netchop_score_ct': 'min',
    'mut_binding_score': 'min',
    'mut_is_binding_pos': 'cnt',
    'mut_aa_coeff': 'cnt',
    'DAI_NetMHC': 'max',
    'DAI_MixMHC': 'max',
    'DAI_NetStab': 'max',
    'mutant_other_significant_alleles': 'min',
    'DAI_MixMHC_mbp': 'max',
    'rnaseq_TPM': 'min',
    'rnaseq_alt_support': 'min',
    'GTEx_all_tissues_expression_mean': 'min',
    'Sample_Tissue_expression_GTEx': 'min',
    'TCGA_Cancer_expression': 'min',
    'bestWTMatchScore_I': 'min',
    'bestWTMatchOverlap_I': 'min',
    'bestMutationScore_I': 'min',
    'bestWTPeptideCount_I': 'min',
    'CSCAPE_score': 'min',
    'CCF': 0.9,
    'nb_same_mutation_Intogen': 'min'
}

excluded_genes: Final[list] = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DRB1', 'HLA-DRB3', 'HLA-DRB4', 'HLA-DRB5',
                               'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DMA', 'TRBV3', 'TRBV5',
                               'TRBV6', 'TRBV6-1', 'TRBV10', 'TRBV10-1', 'TRBV11', 'TRAV12', 'KRT1', 'PRSS3']
"""
features_neopep: Final[list] = \
    ['patient', 'dataset', 'train_test', 'response_type', 'Nb_Samples', 'Sample_Tissue', 'Cancer_Type',
     'chromosome', 'genomic_coord', 'ref', 'alt', 'gene', 'protein_coord', 'aa_mutant', 'aa_wt',
     'pep_mut_start', 'TumorContent', 'Zygosity', 'mutation_type'] + ml_features_neopep
"""

peptide_type = 'neopep'
objective = 'ml'
main_dir = '/data/NeoRanking_self'

#preprocessing.py
encoding_file = '/data/NeoRanking_code/data/cat_encoding/Cat_to_num_info_neopep_NCI_train.txt'
max_mixmhc_rank: int = -1
max_netmhc_rank: int = 20
max_stab_rank: int = -1
max_prime_rank: int = 40
#org_file = '/data/NeoRanking-master/data/Neopep_data_org.txt'
#org_file = main_dir + '/Neopep_data_org_bigmhc_mhcflurry_prime_mixmhcpred.txt'
org_file = main_dir + '/Neopep_data_org_bigmhc_mhcflurry_prime_mixmhcpred_foreignness.txt'
preprocessed_data_file = main_dir + '/Neopep_preprocessed_data.txt'
normalizer_tag = 'q'

#training.py
train_data_reset = True
#文章提供的数据集中固定NCI_train_test标签，可能是因为这种划分可以使NCI_test中的样本尽可能含有阳性肽段，否则NCI_test样本中可能本身没有阳性肽段，也就无法检测并排序阳性肽段。
user_defined_split = False
patient = ''
dataset = 'NCI_train'
response_types = ['CD8', 'negative', 'not_tested']
#response_types = ['CD8', 'negative']
nr_non_immuno_neopeps: int = 100000
sample = True
shuffle = False
data_train_file = main_dir + '/Neopep_train_data.txt'
X_train_file = main_dir + '/Neopep_X_train_data.txt'
y_train_file = main_dir + '/Neopep_y_train_data.txt'

#classifier_tag = 'LR'
classifier_tag = 'XGBoost'
#classifier_tag = 'RF'
#classifier_tag = 'RC'
#classifier_tag = 'LRE'
#metric = 'sum_exp_rank'
#metric = 'precision'
best_loss = np.Inf
best_classifier = None
best_params = None
nr_hyperopt_iter = 200
nr_hyperopt_cv = 5
nr_hyperopt_rep = 10

neopep_alpha = 0.001
shuffle = False
random_state = None
run_tag = 'test'
classifier_model_dir = main_dir + '/classifier_models/'
scorer_name = 'sum_exp_rank'
#scorer_name = 'nr_correct_top30'
#scorer_name = 'sum_prob_top100'


#testing.py
user_defined_test_file = main_dir + '/Neopep_preprocessed_user_defined_test_data.txt'
#alpha = 0.05
#classifier_file_re = '*sav'  #LR,RC,LRE
classifier_file_re = '*xgbm'  #XGBoost
#classifier_file_re = '*joblib'  #RF
classifier_result_dir = main_dir + '/classifier_results/'
verbose: int = 0
write_header: bool = True

#testing_vc.py
weight = 0.5
classifier1_file_re = 'LR*sav'  #LR
classifier2_file_re = '*xgbm'  #XGBoost
#classifier2_file_re = 'RC*sav'  #RC