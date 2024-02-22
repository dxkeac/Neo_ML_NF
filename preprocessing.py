# coding=gbk

from param_value import *

'''
def load_patient_neopep(df):
    if df.shape[0] == 0:
        return None, None, None
    y = np.array(df.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
    if objective == 'ml':
        X = df.loc[:, ml_features_neopep]
    else:
        X = df.copy()
    if objective == 'ml':
        X = fill_missing_values(X)
    X = normalize(X)
    if objective == 'ml':
        X = encode_cat_features(X)
    X = X.astype(get_processed_types(peptide_type, objective))
    return df, X, y
'''

#fill_missing_values
def is_cat_type(type_):
    return type_ in ['category']

def fill_missing_values(x_):
    value_dict = {}
    types = feature_types_neopep if peptide_type == 'neopep' \
        else feature_types_mutation
    order_rels = ml_feature_mv_neopep if peptide_type == 'neopep' \
        else ml_feature_mv_mutation
    #for i, c in enumerate(x_.columns):
    for i, c in enumerate(list(ml_feature_mv_neopep.keys())):
        if c == 'rnaseq_alt_support':
            x_ = impute_rnaseq_cov(x_)
        if not is_cat_type(types[c]):
            order_rel = order_rels[c]
            if order_rel == 'min':
                fill_val = np.min(x_[c])
            elif order_rel == 'max':
                fill_val = np.max(x_[c])
            elif order_rel == 'cnt':
                counter = Counter(x_[c])
                fill_val = counter.most_common(1)[0][0]
            elif type(order_rel) == float:
                fill_val = order_rel
            else:
                fill_val = np.mean(x_[c])
            value_dict[c] = fill_val
    x_.fillna(value=value_dict, inplace=True)
    return x_

'''
原始的，有错
def impute_rnaseq_cov(x_):
    quartiles = np.quantile(a=x_['rnaseq_TPM'], q=[0.5, 0.75])
    mv = np.where(x_['rnaseq_TPM'] < quartiles[0], 0, 11)
    mv = np.where(x_['rnaseq_TPM'] > quartiles[1], mv, 23)
    x_['rnaseq_alt_support'] = \
        x_['rnaseq_alt_support'].fillna(data=pd.Series(mv), index=x_.index)
    return x_
#测试
a = list(range(10))
quartiles = np.quantile(a, q=[0.5, 0.75])
mv = np.where(a < quartiles[0], 0, 11)
mv
array([ 0,  0,  0,  0,  0, 11, 11, 11, 11, 11])
mv = np.where(a > quartiles[1], 23, mv)
mv
array([ 0,  0,  0,  0,  0, 11, 11, 23, 23, 23])
'''

#修改后的
def impute_rnaseq_cov(x_):
    quartiles = np.quantile(a=x_['rnaseq_TPM'], q=[0.5, 0.75])
    mv = np.where(x_['rnaseq_TPM'] < quartiles[0], 0, 11)
    mv = np.where(x_['rnaseq_TPM'] > quartiles[1], 23, mv)
    x_.loc[:,'rnaseq_alt_support'] = x_['rnaseq_alt_support'].fillna(pd.Series(mv))
    #x_['rnaseq_alt_support'] = x_['rnaseq_alt_support'].fillna(pd.Series(mv))
    return x_

#normalize
def is_cont_type(type_):
    '''
    Plot as continuous values, but natural order
    :param type_:
    :return: True if type_ any 'floatxx' or 'intxx' (xx>8)
    '''
    return type_.startswith('float') or type_ in ['int16', 'int32', 'int64']

def normalize(x_):
    if normalizer is None:
        return x_
    #for i, c in enumerate(x_.columns):
    for i, c in enumerate(ml_normalize_features_neopep):
        if is_cont_type(x_[c].dtype.name):
            if type(normalizer) is dict:
                if c in normalizer:
                    norm_transform = normalizer[c]
                else:
                    norm_transform = None
            else:
                norm_transform = normalizer
            if norm_transform:
                v = x_[c].to_numpy().reshape(-1, 1)
                if type(norm_transform) is FunctionTransformer and norm_transform.func.__name__ == 'log10':
                    if sum(v > 0) > 0:
                        v[v <= 0] = min(v[v > 0])/10
                    else:
                        v[v <= 0] = 1
                x_.loc[:, c] = norm_transform.fit_transform(v)
    return x_

def get_normalizer(normalizer_tag):
    if normalizer_tag == 'q':
        return QuantileTransformer()
    elif normalizer_tag == 'z':
        return StandardScaler()
    elif normalizer_tag == 'p':
        return PowerTransformer()
    elif normalizer_tag == 'i':
        return MinMaxScaler()
    elif normalizer_tag == 'l':
        return FunctionTransformer(np.log10, inverse_func=lambda x: np.power(10, x), validate=False, check_inverse=True)
    elif normalizer_tag == 'a':
        return FunctionTransformer(np.arcsinh, inverse_func=np.sinh, validate=False, check_inverse=True)
    elif normalizer_tag == 'n':
        return None
    elif normalizer_tag == 'plot':
        d = {}
        for k, v in plot_normalization.items():
            d[k] = get_normalizer(v)
        return d

def encode_cat_features(X, encoders: dict, encoding_df_any_na_row):
    for k,v in encoders.items():
        X.loc[:, k] = X[k].replace(encoders[k])
    for col in encoding_df_any_na_row['Feature']:
        #print(col)
        X[col].fillna(encoding_df_any_na_row[encoding_df_any_na_row['Feature'] == col]['Value'].tolist()[0], inplace=True)
    return X

def get_processed_types(peptide_type: str, objective: str):
    if objective == 'ml':
        if peptide_type == 'neopep':
            type_dict = {feature: 'float' for feature in ml_features_neopep}
            #type_dict['seq_len'] = 'int'
            return type_dict

def filter_rows_neopep(df):
    if df.shape[0] > 0:
        df = df.loc[df.mutation_type.apply(lambda r: r == 'SNV')]
    if df.shape[0] > 0 and max_mixmhc_rank > 0:
        df = df.loc[df.mutant_rank.apply(lambda r: r < max_mixmhc_rank)]
    if df.shape[0] > 0 and max_netmhc_rank > 0:
        df = df.loc[df.mutant_rank_netMHCpan.apply(lambda r: r < max_netmhc_rank)]
    if df.shape[0] > 0 and max_stab_rank > 0:
        df = df.loc[df.mut_Rank_Stab.apply(lambda r: r < max_stab_rank)]
    if df.shape[0] > 0 and max_prime_rank > 0:
        df = df.loc[df.mutant_rank_PRIME.apply(lambda r: r < max_prime_rank)]
    if df.shape[0] > 0 and len(excluded_genes) > 0:
        df = df.loc[df.gene.apply(lambda g: g not in excluded_genes)]
    df.reset_index(inplace=True, drop=True)
    return df


#encode_cat_features
encoding_df = pd.read_csv(encoding_file, header=0, sep="\t", comment='#')
#过滤出含有na的行
encoding_df_any_na_row = encoding_df[encoding_df.isna().any(axis=1)]
encoding_df = encoding_df[~encoding_df.isna().any(axis=1)]
#encoding_df.dropna()
#encoding_df.dropna(subset=['Value'])
features = encoding_df['Feature'].unique()
'''
array(['mut_is_binding_pos', 'bestWTMatchType_I', 'Clonality',
       'mutation_driver_statement_Intogen', 'gene_driver_Intogen',
       'seq_len'], dtype=object)
categories = encoding_df['Category'].unique()
'''
encoders = {}
for feature in features:
    #encoder = CatEncoder(f)
    #encoder.get_encoding(encoding_df)
    df_feature = encoding_df[encoding_df['Feature'] == feature]
    #print(df_feature)
    if feature not in encoders:
        encoders[feature] = {}
    categories = df_feature['Category'].unique()
    for category in categories:
        encoders[feature][category] = float(df_feature[df_feature['Category'] == category]['Value'].to_list()[0])
    #encoders[f] = encoder
    #return encoders


#preprocess running
df = pd.read_csv(org_file, sep='\t',low_memory=False)
'''
df.dtypes
data = pd.read_csv(data_file, sep="\t", header=0)
df = df.astype(dtype=feature_types_neopep, errors='ignore')
data.replace('', 'nan', inplace=True)
df['response_type'].value_counts()
response_type
not_tested    1364625
negative       422907
CD8               178
df.shape
(1787710, 57)
'''
df_filter_rows = filter_rows_neopep(df)
'''
df_filter_rows.shape
(852045, 57)
#和文章代码运行结果一致
wc -l *txt
    852046 Neopep_data_ml_norm.txt
    852046 Neopep_data_ml_sel.txt
   1787711 Neopep_data_org.txt
df_filter_rows[['dataset','train_test']].value_counts()
dataset  train_test
NCI      train         487670
         test          172451
TESLA    test          142546
HiTIDE   test           49378
'''

#每个样本填补缺失值，标准化，类转数
normalizer = get_normalizer(normalizer_tag)
patients = df_filter_rows['patient'].unique()
#DataManager._immunogenic_patients[peptide_type] = []

for i, p in enumerate(sorted(patients)):
    print("processing patient {0}".format(p))
    df_patient = df_filter_rows[df_filter_rows['patient'] == p]
    df_patient_fill = fill_missing_values(df_patient)
    df_patient_fill_norm = normalize(df_patient_fill)
    if i == 0:
        df_patient_fill_norm_comb = df_patient_fill_norm
    else:
        df_patient_fill_norm_comb = pd.concat([df_patient_fill_norm_comb, df_patient_fill_norm], ignore_index=True)

#def filter_selected_data(peptide_type: str, patient: str = "", dataset: str = "",
#                         response_types: list = response_types) -> pd.DataFrame:

if 'mutation_driver_statement_Intogen' in df_patient_fill_norm_comb.columns and 'seq_len' in df_patient_fill_norm_comb.columns:
    df_patient_fill_norm_comb['mutation_driver_statement_Intogen'].replace({'POLYMORPHISM':'unknown'},inplace = True)
    #Intogen['driver_statement'] = Intogen['driver_statement'].replace({'predicted passenger':'PREDICTED_PASSENGER','polymorphism':'POLYMORPHISM','not protein-affecting':'NOT_PROTEIN_AFFECTING'})
    #Intogen['driver_statement'] = Intogen['driver_statement'].replace('predicted driver.*','PREDICTED_DRIVER',regex=True)
    #df_patient_fill_norm_comb = df_patient_fill_norm_comb.astype({'seq_len': 'object'})
    #encoders['seq_len']
    seq_len_dict = {12: 0.0208332380835462, 8: 0.0392458784892738, 11: 0.0399536898050818, 10: 0.2189206723038542, 9: 0.6810465213182438, 'unknown': 0.0}
    #df_patient_fill_norm_comb_cat['seq_len'].astype("string").replace({'9':'0.1'})
    #df_patient_fill_norm_comb_cat['seq_len'].astype("float16").replace({9.0:0.1})
    df_patient_fill_norm_comb['seq_len'] = df_patient_fill_norm_comb['seq_len'].replace(seq_len_dict)
else:
    pass

df_patient_fill_norm_comb_cat = encode_cat_features(df_patient_fill_norm_comb, encoders, encoding_df_any_na_row)
df_patient_fill_norm_comb_cat = df_patient_fill_norm_comb_cat.astype(get_processed_types(peptide_type, objective))
df_patient_fill_norm_comb_cat.to_csv(preprocessed_data_file, header=True, index=False,sep='\t')
