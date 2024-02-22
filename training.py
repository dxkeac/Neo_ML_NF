# coding=gbk

from param_value import *

#def add(a, b):
#    return a + b

#过滤出NCI_train数据
def get_filtered_data_index(data: pd.DataFrame, patient: str, dataset: str, response_types: list) -> list:
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
    return idx

def sample_rows(data, X, y) -> tuple:
    if sum(y == 0) < nr_non_immuno_neopeps:
        return data, X, y
    idx = random.sample(range(sum(y == 0)), nr_non_immuno_neopeps)
    X_1 = X.loc[y == 1, :]
    #fraction = 0.8
    #sampled_data = X_1.sample(frac=fraction, random_state=42)
    X_0 = X.loc[y == 0, :]  #.loc按index或者条件提取,从0开始
    if X_0.shape[0] > nr_non_immuno_neopeps:
        X_0 = X_0.iloc[idx, :]  #.iloc按第几行/个提取,从0开始
    X_s = pd.concat([X_1, X_0])
    data_1 = data.loc[y == 1, :]
    data_0 = data.loc[y == 0, :]
    if data_0.shape[0] > nr_non_immuno_neopeps:
        data_0 = data_0.iloc[idx, :]
    data_s = pd.concat([data_1, data_0])
    y_0 = y[y == 0]
    y_1 = y[y == 1]
    y_0 = y_0[idx]
    y_s = np.append(y_1, y_0)
    return data_s, X_s, y_s

def shuffle_rows(data, X, y) -> tuple:
    idx = random.sample(range(len(y)), k=len(y))
    X_s = X.iloc[idx, :]
    data_s = data.iloc[idx, :]
    y_s = y[idx]
    return data_s, X_s, y_s

#data = df_patient_fill_norm_comb_cat
preprocessed_data_df = pd.read_csv(preprocessed_data_file, sep='\t',low_memory=False)
preprocessed_data_df_X = preprocessed_data_df[ml_features_neopep]
#X.isna().sum()
preprocessed_data_df_y = preprocessed_data_df.response_type
#preprocessed_data_np_y = np.array(preprocessed_data_df.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
#range(sum(y == 0))
#range(0, 851867)
if train_data_reset == True:
    if user_defined_split == True:
        preprocessed_data_df_NCI = preprocessed_data_df[preprocessed_data_df['dataset'] == 'NCI']
        NCI_index = preprocessed_data_df_NCI.index
        preprocessed_data_df_X_NCI = preprocessed_data_df_X.loc[NCI_index,:]
        preprocessed_data_df_y_NCI = preprocessed_data_df_y.loc[NCI_index]
        """
        #X = np.array([[1, 2], [3, 4], [1, 2], [3, 4], [1, 2], [3, 4]])
        #y = np.array([0, 0, 0, 1, 1, 1])
        preprocessed_data_df_y_NCI = np.array(preprocessed_data_df_y_NCI.apply(lambda rt: int(rt == 'CD8')))
        #按比例并分层分割数据，没有考虑到单个样本作为整体抽样，造成单个样本的CD8分散在训练和测试集
        sss = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=0)
        for i, (train_index, test_index) in enumerate(sss.split(preprocessed_data_df_X_NCI, preprocessed_data_df_y_NCI)):
            pass
            print(f"Fold {i}:")
            print(f"  Train: index={train_index}")
            print(f"  Test:  index={test_index}")
        pd.DataFrame(preprocessed_data_df_y_NCI[train_index]).value_counts()
        0    528014
        1        82
        pd.DataFrame(preprocessed_data_df_y_NCI[test_index]).value_counts()
        0    132004
        1        21
        data_train_all, data_test_all = preprocessed_data_df_NCI.iloc[train_index,:], preprocessed_data_df_NCI.iloc[test_index,:]
        X_train_all, X_test_all = preprocessed_data_df_X.iloc[train_index, :], preprocessed_data_df_X.iloc[test_index, :]
        y_train_all, y_test_all = preprocessed_data_df_y_NCI[train_index], preprocessed_data_df_y_NCI[test_index]
        data_train['response_type'].value_counts()
        response_type
        not_tested    366688
        negative      161326
        CD8               82
        user_defined_test = pd.concat([preprocessed_data_df[preprocessed_data_df['dataset'] != 'NCI'],preprocessed_data_df_NCI.iloc[test_index,:]])
        user_defined_test.to_csv(user_defined_test_file, header=True, index=False,sep='\t')
        """

        NCI_patients = preprocessed_data_df_NCI['patient'].unique()
        selected_patients = []
        selected_patients_CD8 = 0
        for i in range(preprocessed_data_df_NCI.shape[0]):
            selected_patient = np.random.choice(NCI_patients,replace=False)
            if selected_patient in selected_patients:
                continue
            selected_patient_CD8 = sum(preprocessed_data_df_NCI[preprocessed_data_df_NCI['patient'] == selected_patient]['response_type'] == 'CD8')
            selected_patients_CD8 += selected_patient_CD8
            selected_patients.append(selected_patient)
            #print(selected_patients_CD8)
            #if i == 0:
                #preprocessed_data_df_NCI_test = preprocessed_data_df_NCI[preprocessed_data_df_NCI['patient'] == selected_patient]
            #else:
                #preprocessed_data_df_NCI_test = pd.concat([preprocessed_data_df_NCI_test,preprocessed_data_df_NCI[preprocessed_data_df_NCI['patient'] == selected_patient]])
            if selected_patients_CD8 < round(sum(preprocessed_data_df_NCI['response_type'] == 'CD8')*0.2):
                continue
            elif selected_patients_CD8 == round(sum(preprocessed_data_df_NCI['response_type'] == 'CD8')*0.2):
                break
            else:
                selected_patients = []
                selected_patients_CD8 = 0
                #preprocessed_data_df_NCI_test = pd.DataFrame()
        #sum(preprocessed_data_df_NCI.loc[preprocessed_data_df_NCI.apply(lambda row : row['patient'] in selected_patients, axis=1),:]['response_type'] == 'CD8')
        #21
        filtered_data_test = preprocessed_data_df_NCI[preprocessed_data_df_NCI['patient'].astype(str).isin(selected_patients)]
        #filtered_data = filtered_data.drop_duplicates(ignore_index=True)
        #sum(filtered_data_test['response_type'] == 'CD8')
        #21
        filtered_data_train = preprocessed_data_df_NCI[~preprocessed_data_df_NCI['patient'].astype(str).isin(selected_patients)]
        #sum(filtered_data_train['response_type'] == 'CD8')
        #82
        train_index = filtered_data_train.index
        test_index = filtered_data_test.index
        data_train_all, data_test_all = filtered_data_train, filtered_data_test
        X_train_all, X_test_all = preprocessed_data_df_X_NCI.loc[train_index, :], preprocessed_data_df_X_NCI.loc[test_index, :]
        y_train_all, y_test_all = preprocessed_data_df_y_NCI.loc[train_index], preprocessed_data_df_y_NCI.loc[test_index]
        y_train_all = np.array(y_train_all.apply(lambda rt: int(rt == 'CD8')))
        y_test_all = np.array(y_test_all.apply(lambda rt: int(rt == 'CD8')))

        #y_test = np.array(data_test.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
        #y_test_pd = pd.DataFrame(y_test)
        #y_test_pd.index = data_test.index

        """
        data_train['response_type'].value_counts()
        response_type
        not_tested    366688
        negative      161326
        CD8               82
        """
        user_defined_test = pd.concat([preprocessed_data_df[preprocessed_data_df['dataset'] != 'NCI'],filtered_data_test])
        user_defined_test.to_csv(user_defined_test_file, header=True, index=False,sep='\t')
    else:
        idx = get_filtered_data_index(preprocessed_data_df, patient, dataset, response_types)
        if not all(idx):
        #if any(idx):
            data_train_all = preprocessed_data_df.loc[idx, :]
            X_train_all = preprocessed_data_df_X.loc[idx, :]
            y_train_all = preprocessed_data_df_y.loc[idx]
            y_train_all = np.array(y_train_all.apply(lambda rt: int(rt == 'CD8')))
        #[487670 rows x 57 columns]
        #和直接用NCI_train标签过滤效果一样
        #data[(data['dataset'] == 'NCI') & (data['train_test'] == 'train')]
        #[487670 rows x 31 columns]

    if sample:
        data_train, X_train, y_train = sample_rows(data_train_all, X_train_all, y_train_all)
    elif shuffle:
        data_train, X_train, y_train = shuffle_rows(data_train, X_train, y_train)

    data_train.to_csv(data_train_file, header=True, index=False,sep='\t')
    X_train.to_csv(X_train_file, header=True, index=False,sep='\t')
    pd.DataFrame(y_train).to_csv(y_train_file, header=True, index=False,sep='\t')


"""
data.shape
(100082, 57)
X.shape
(100082, 31)
y.shape
(100082,)
"""

#skf = StratifiedKFold(n_splits=5,random_state=0,shuffle=True)
skf = StratifiedKFold(n_splits=nr_hyperopt_cv, shuffle=True)

def sum_rank_correct(y_true: np.ndarray, y_pred: np.ndarray):
    """
    Rank_score optimization score used in the paper
    Args:
        y_true (np.ndarray): array with true immunogenicity indicators (0: non-immunogenic, 1: immunogenic)
        y_pred (np.ndarray): array with predicted probabilities that peptide is immunogenic
    Returns:
        rank_score (float): sum of rank_scores for all immunogenic peptides
    """
    idx = np.argsort(-y_pred)
    y_true = y_true[idx]
    #print(idx)
    r = np.where(y_true == 1)[0]
    return np.sum(np.exp(np.multiply(-neopep_alpha, r)))

def sum_rank_correct_pp(y_true, y_pred, patients):
    score = 0.0
    for p in set(patients):
        idx = p == patients
        y_true_p = y_true[idx]
        y_pred_p = y_pred[idx]
        idx = np.argsort(-y_pred_p)
        y_true_p = y_true_p[idx]
        r = np.where(y_true_p == 1)[0]
        score += np.sum(np.exp(np.multiply(-alpha, r)))
    return score

def nr_correct_top100(y_true, y_pred, max_rank=100):
    n = min(len(y_true), max_rank)
    idx = np.argsort(-y_pred)
    return np.sum(y_true[idx][:n] == 1)

def nr_correct_top30(y_true, y_pred, max_rank=30):
    n = min(len(y_true), max_rank)
    idx = np.argsort(-y_pred)
    return np.sum(y_true[idx][:n] == 1)

def sum_prob_correct(y_true, y_pred, max_rank=100):
    idx = np.argsort(-y_pred)
    n = min(len(y_true), max_rank)
    y_true = y_true[idx][:n]
    y_pred = y_pred[idx][:n]
    n = np.sum(y_true)
    if n == 0:
        s = 0
    else:
        s = np.sum(y_pred[y_true == 1])/n
    return s - np.sum(y_pred[y_true != 1])/(100-n)

def get_scorer(scorer_name: str, data: pd.DataFrame):
    """
    Defines a scorer object used as loss score in Hyperopt optimization
    Args:
        scorer_name (str): name of loss score function
        data (pd.DataFrame): dataframe (only used in scorer_name=='sum_exp_rank_pp')
    Returns:
        scorer (Callable): Callable object that returns a scalar score; greater is better.
    """
    if scorer_name == 'sum_exp_rank':
        return make_scorer(sum_rank_correct, needs_threshold=True)
    elif scorer_name == 'sum_exp_rank_pp':
        return make_scorer(sum_rank_correct_pp, patients=data['patient'].to_numpy(), needs_threshold=True)
    elif scorer_name == 'nr_correct_top100':
        return make_scorer(nr_correct_top100, needs_threshold=True)
    elif scorer_name == 'nr_correct_top30':
        return make_scorer(nr_correct_top30, needs_threshold=True)
    elif scorer_name == 'sum_prob_top100':
        return make_scorer(sum_prob_correct, needs_threshold=True)
    else:
        print('No scorer with name '+str(scorer_name)+' implemented. Abort')
        return None

def get_classifier(classifier_tag: str, params: dict):
    """
    Creates classifier object with hyperparameters corresponding to params
    Args:
        classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
        params (dict): dictionary with parameter values for classifier specified in classifier_tag
    Returns:
        classifier object (Object): classifier object with hyperparameters corresponding to params
    """
    if classifier_tag == "SVM":
        return SVC(probability=True, kernel='rbf', C=params['C'], gamma=params['gamma'],
                   class_weight=params['class_weight'])
    elif classifier_tag == "SVM-lin":
        return SVC(probability=True, kernel='linear', C=params['C'], class_weight=params['class_weight'])
    elif classifier_tag == "LR":
        return LogisticRegression(solver='saga',  penalty=params['penalty'], C=params['C'],
                                  class_weight=params['class_weight'])
    elif classifier_tag == "LRE":
        return LogisticRegression(solver='saga',  penalty="elasticnet", C=params['C'], l1_ratio=params['l1_ratio'],
                                  class_weight=params['class_weight'])
    elif classifier_tag == "CatBoost":
        return CatBoostClassifier(
            loss_function='Logloss',
            iterations=params['iterations'],
            subsample=params['subsample'],
            random_strength=params['random_strength'],
            learning_rate=params['learning_rate'],
            l2_leaf_reg=params['l2_leaf_reg'],
            leaf_estimation_iterations=params['leaf_estimation_iterations'],
            depth=params['depth'],
            bagging_temperature=params['bagging_temperature'],
            random_seed=42,
            use_best_model=False,
            cat_features=cat_idx,
            auto_class_weights=params['auto_class_weights'],
            silent=True)
    elif classifier_tag == "XGBoost":
        #SystemError
        #null argument to internal routine  系统资源不够，减小nr_hyperopt_rep数目即可
        return XGBClassifier(
            enable_categorical=True,
            max_depth=params['max_depth'],
            learning_rate=params['learning_rate'],
            n_estimators=params['n_estimators'],
            eval_metric='logloss',
            verbosity=0,
            silent=None,
            objective='binary:logistic',
            booster=params['booster'],
            tree_method='hist',
            n_jobs=int(os.cpu_count()/nr_hyperopt_rep - 1 ),
            nthread=None,
            gamma=0,
            min_child_weight=params['min_child_weight'],
            max_delta_step=0,
            subsample=params['subsample'],
            colsample_bytree=params['colsample_bytree'],
            colsample_bylevel=params['colsample_bylevel'],
            colsample_bynode=1,
            reg_alpha=params['reg_alpha'],
            reg_lambda=1,
            scale_pos_weight=params['scale_pos_weight'],
            base_score=0.5,
            random_state=0,
            seed=None)
    elif classifier_tag == "RF":
        #return RandomForestClassifier(**params)
        return RandomForestClassifier(
            n_jobs=int(os.cpu_count()/nr_hyperopt_rep-1),
            n_estimators=params['n_estimators'],
            criterion=params['criterion'],
            max_depth=params['max_depth'],
            max_features=params['max_features'],
            min_samples_leaf=params['min_samples_leaf'],
            min_samples_split=params['min_samples_split'],
            class_weight=params['class_weight']
            )
    elif classifier_tag == "RC":
        #return RandomForestClassifier(**params)
        return RidgeClassifier(
        #return RidgeClassifierWithProba(
            alpha=params['alpha'],
            fit_intercept=params['fit_intercept'],
            #solver=params['solver'],
            class_weight=params['class_weight']
        )

def get_class_weights( ):
    if not class_ratio or class_ratio == 0 or class_ratio >= 1:
        return 'balanced'
    else:
        cws = []
        v = int(2.0/class_ratio)
        for cw in range(1, int(2.0/class_ratio), round(v/20)):
            cws.append({1: cw})
        return hp.choice('class_weight', cws)

def get_xgb_pos_weights( ):
    if not class_ratio or class_ratio == 0 or class_ratio >= 1:
        return 1
    else:
        v = int(2.0/class_ratio)
        return hp.choice('scale_pos_weight', range(1, v, round(v/20)))

def get_param_space(classifier_tag: str):
    """
    Defines parameter space searched durich Hyperopt loop for each classifier
    Args:
        classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
    Returns:
        parameter_space (dict): dictionary with parameter space for each classifier hyperparameter used in the
        Hyperopt loop
    """
    #global parameter_space
    #parameter_space = {}
    if classifier_tag == "SVM":
        parameter_space = {
            'C': hp.uniform('C', 0.005, 1.0),
            'gamma': hp.uniform('gamma', 0, 2),
            'class_weight': get_class_weights()
        }
    elif classifier_tag == "SVM-lin":
        parameter_space = {
            'C': hp.uniform('C', 0.005, 1.0),
            'class_weight': get_class_weights()
        }
    elif classifier_tag == "LR":
        parameter_space = {
            'penalty': hp.choice('penalty', ['l1', 'l2']),
            'C': hp.uniform('C', 0.0, 5.0),
            'class_weight': get_class_weights()
        }
    elif classifier_tag == "LRE":
        parameter_space = {
            'C': hp.uniform('C', 0.0, 5.0),
            'l1_ratio': hp.uniform('l1_ratio', 0.0, 1.0),
            'class_weight': get_class_weights()
        }
    elif classifier_tag == "XGBoost":
        parameter_space = {
            'booster': hp.choice('booster', ['gbtree', 'gblinear']),
            'max_depth': hp.choice('max_depth', [3, 4, 5, 7, 9]),
            'min_child_weight': hp.choice('min_child_weight', np.round(np.arange(0.0,  0.2, 0.01), 5)),
            'learning_rate': hp.loguniform('learning_rate', np.log(0.01), np.log(0.1)),
            'subsample': hp.uniform('subsample', 0.3, 1.0),
            'colsample_bylevel': hp.uniform('colsample_bylevel', 0.4, 1.0),
            'colsample_bytree': hp.uniform('colsample_bytree', 0.4, 1.0),
            'n_estimators': hp.choice('n_estimators', np.arange(50, 1500, 50)),
            'reg_alpha': hp.loguniform('reg_alpha', np.log(0.0001), np.log(1)),
            'gamma':  hp.uniform('gamma', 0.0, 10.0),
            'scale_pos_weight': get_xgb_pos_weights(),
        }
    elif classifier_tag == "CatBoost":
        parameter_space = {
            'iterations': hp.choice('iterations', np.round(np.arange(100, 1500, 100))),
            'auto_class_weights': hp.choice('auto_class_weights', ['None', 'Balanced']),
            'subsample': hp.uniform('subsample', 0.3, 1.0),
            'random_strength': scope.int(hp.quniform('random_strength', 1, 20, 1)),
            'learning_rate': hp.loguniform('learning_rate', np.log(0.0001), np.log(1)),
            'l2_leaf_reg': hp.loguniform('l2_leaf_reg', np.log(1), np.log(10)),
            'leaf_estimation_iterations': scope.int(hp.quniform('leaf_estimation_iterations', 1, 20, 1)),
            'depth': scope.int(hp.quniform('depth', 5, 10, 1)),
            'bagging_temperature': hp.uniform('bagging_temperature', 0.0, 1.0)
        }
    elif classifier_tag == "RF":
        parameter_space = {
            'n_estimators': hp.choice('n_estimators', np.arange(50, 1500, 50)),
            'criterion': hp.choice('criterion', ['entropy','gini','log_loss']),
            'max_depth': hp.choice('max_depth', np.arange(10, 100, 5)),
            #'max_features': hp.choice('max_features',['auto','sqrt','log2','None']),
            'max_features': hp.choice('max_features',['sqrt','log2',None]),
            #'min_samples_split': hp.uniform('min_samples_split',0,1),
            'min_samples_split': hp.choice('min_samples_split', np.arange(2, 20, 1)),
            #'min_samples_leaf': hp.uniform('min_samples_leaf',0,0.5),
            'min_samples_leaf': hp.choice('min_samples_leaf', np.arange(1, 20, 1)),
            'class_weight': get_class_weights()
        }
    elif classifier_tag == "RC":
        parameter_space = {
            #'alpha': hp.choice('alpha', [0.01, 0.05, 0.1, 1, 5, 10]),
            'alpha': hp.uniform('alpha', 0.0, 5.0),
            'fit_intercept': hp.choice('fit_intercept', [True,False]),
            #'solver': hp.choice('solver',['auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg']),
            'class_weight': get_class_weights()
        }

    return parameter_space  #上边需要空一行，否则报错UnboundLocalError: cannot access local variable 'parameter_space' where it is not associated with a value

def save_classifier(classifier_tag: str, classifier, classifier_file: str):
    """
    Saves classifier to file.
    Args:
        classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
        classifier (object): classifier to be fitted. Classifier object must implement the sklearn fit
                             method
        classifier_file (str): file name of classifier model file
    """
    if classifier_tag in ['CatBoost', 'XGBoost']:
        classifier.save_model(classifier_file)
    elif classifier_tag == "RF":
        joblib.dump(classifier, classifier_file)
    else:
        pickle.dump(classifier, open(classifier_file, 'wb'))


def fit_classifier(x: pd.DataFrame, y: np.ndarray, classifier=None, params: dict = None) -> object:
    """
    Calls classifier.fit. If params is None, then classifier object is used. Otherwise, a classifier is
    constructed with hyperparameters defined in params
    Args:
        x (pd.DataFrame): processed dataframe with rows and columns selected for ML
        y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
        classifier (object): classifier to be fitted. Classifier object must implement the sklearn fit method
        params (dict): dictionary with classifiers hyperparameters
    Returns:
        fitted classifier object
    """
    assert classifier is not None or params is not None
    if params is None:
        clf = classifier
    else:
        clf = get_classifier(classifier_tag, params)
    if classifier_tag == 'CatBoost':
        clf.fit(x, y, plot=False)
    else:
        clf.fit(x, y)
    return clf

def get_classifier_file(clf_name, sub_dir, run_tag, run_idx, peptide_type):
    #file_dir = os.path.join(classifier_model_dir, sub_dir)
    file_dir = sub_dir
    date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
    if clf_name in ['LR', 'LRE', 'SVM', 'SVM-lin', 'RC']:
        ext = 'sav'
    elif clf_name == 'XGBoost':
        ext = 'xgbm'
    elif clf_name == 'CatBoost':
        ext = 'cbm'
    elif clf_name == 'RF':
        ext = 'joblib'
    file_name = '{0}_{1}_{2}_{3}_{4}_clf.{5}'.format(clf_name, run_tag, run_idx, peptide_type, date_time_str, ext)
    result_file = path.join(file_dir, file_name)
    # make sure file does not already exist
    while os.path.isfile(result_file):
        date_time_str = datetime.datetime.now().strftime("%m.%d.%Y-%H.%M.%S")
        file_name = '{0}_{1}_{2}_{3}_{4}_model.{5}'.format(clf_name, run_tag, run_idx, peptide_type, date_time_str, ext)
        result_file = path.join(file_dir, file_name)
    return result_file

def score(params):
    #classifier = SVC(probability=True, kernel='rbf', C=params['C'], gamma=params['gamma'], class_weight=params['class_weight'])
    #global best_loss, best_classifier, best_params
    """
    #classifier = SVC(**params)
    #classifier = clf_xgb
    classifier = get_classifier(classifier_tag, params)
    # cross_val_score calls the metric function with arguments metric(y, classifier.predict(X))
    #loss = 1 - cross_val_score(classifier, x, y, cv=stratifiedKFold, scoring='accuracy').mean()
    loss = 1 - cross_val_score(classifier, x, y, cv=skf, scoring=make_scorer(sum_rank_correct)).mean()
    if loss < best_loss:
        best_loss = loss
        best_classifier = classifier
        best_params = params
    return loss
    """
    classifier = get_classifier(classifier_tag, params)
    #if classifier_tag in ['SVM', 'SVM-lin', 'LR', 'XGBoost', 'RF']:
    if classifier_tag != 'CatBoost':
        # cross_val_score calls the metric function with arguments metric(y, classifier.predict(X))
        loss = 1 - cross_val_score(classifier, x, y, cv=skf, scoring=get_scorer(scorer_name, data)).mean()
        #loss = 1 - cross_val_score(classifier, x, y, cv=skf, scoring=make_scorer(sum_rank_correct)).mean()
        #loss = -cross_val_score(classifier, x, y, cv=skf, scoring=make_scorer(nr_correct_top30)).mean()
        #if loss < best_loss:
            #best_loss = loss
            #best_classifier = classifier
            #best_params = params
        return loss
    elif classifier_tag == 'CatBoost':
        loss = 0
        for train_index, valid_index in skf.split(x, y):
            X_train, X_valid = x.iloc[train_index, :], x.iloc[valid_index, :]
            y_train, y_valid = y[train_index], y[valid_index]
            classifier.fit(X_train, y_train, use_best_model=True, eval_set=(X_valid, y_valid))
            res = classifier.get_best_score()
            loss += res['validation']['Logloss']
        loss = loss/nr_hyperopt_cv
        #if loss < best_loss:
            #best_loss = loss
            #best_classifier = classifier
            #best_params = params
        return loss

#params = get_param_space(classifier_tag)
#classifier = get_classifier(classifier_tag, params)
def run_training_new(run_index):
    start = time.time()
    trials = Trials()
    best = fmin(fn=score,
                space=get_param_space(classifier_tag),
                algo=tpe.suggest,
                max_evals=nr_hyperopt_iter,
                trials=trials,
                rstate=np.random.RandomState(42+run_index*997),verbose=1)
    elapsed_time_hopt = time.time() - start
    classifier_file = get_classifier_file(classifier_tag, classifier_model_dir, run_tag, run_index, peptide_type)
    #'/data/NeoRanking_code/simple_code/classifier_model/LR_test_0_neopep_01.15.2024-11.14.58_clf.sav'
    best_params = space_eval(get_param_space(classifier_tag), best)
    #print("best_params = ",best_params)
    #best_params =  {'C': 1.8844593808692833, 'class_weight': {1: 123}, 'penalty': 'l2'}
    min_loss = min(trials.losses())
    best_clf = fit_classifier(x, y, classifier_tag, best_params)
    save_classifier(classifier_tag, best_clf, classifier_file)
    clf_param_file = re.sub("_clf\\.\\w+$", "_param.txt", classifier_file)
    with open(clf_param_file, mode='w') as param_file:
        verbose = 0
        if verbose > 0:
            print('Classifier = {0:s}, run index = {1:d}\nBest training params: {2:s}\nsum_exp_rank = {3:.3f}\nSaved to {4:s}'.
                  format(str(best_clf), run_index, str(best_params), min_loss, classifier_file))
            print("Hyperopt: Score={0:.3f}, Time={1:f}, Params={2:s}".
                  format(((1 - min_loss) * 100), elapsed_time_hopt, str(best_params)))
        if param_file is not None:
            param_file.write("Hyperopt: Score={0:.3f}; Time={1:f}; Params={2:s}\n".
                              format(((1 - min_loss) * 100), elapsed_time_hopt,
                                     str(best_params)))
            #param_file.write('Training dataset: {0}\n'.format(dataset_train))
            param_file.write('Saved to {0:s}\n'.format(classifier_file))
            param_file.flush()


print(classifier_tag)

if train_data_reset == True:
    data = data_train
    x = X_train
    y = y_train
    class_ratio = sum(y == 1)/sum(y == 0)
    with Pool(processes=nr_hyperopt_rep) as pool:
        pool.map(run_training_new, range(nr_hyperopt_rep))
else:
    data = pd.read_csv(data_train_file, sep='\t',low_memory=False)
    x = pd.read_csv(X_train_file, sep='\t',low_memory=False)
    y = pd.read_csv(y_train_file, sep='\t',low_memory=False)
    y = np.array(y)
    class_ratio = sum(y == 1)/sum(y == 0)
    with Pool(processes=nr_hyperopt_rep) as pool:
        pool.map(run_training_new, range(nr_hyperopt_rep))



"""
#split分割数据时，返回的是从0开始的行号(不是行索引)，因此应该用iloc. 
a_df
   a  b
0  A  B
1  A  B
2  C  D
3  E  F
aa_df = a_df[a_df['a'] != 'C']
x = pd.concat([aa_df,aa_df])
x
   a  b
0  A  B
1  A  B
3  E  F
0  A  B
1  A  B
3  E  F
y = np.array([0,1,1,0,1,1])
for i, (train_index, test_index) in enumerate(sss.split(x, y)):
    print(f"Fold {i}:")
    print(f"  Train: index={train_index}")
    print(f"  Test:  index={test_index}")
Fold 0:
  Train: index=[3 1 4 2]
  Test:  index=[0 5]
x.iloc[train_index,:]
   a  b
0  A  B
1  A  B
1  A  B
3  E  F

for train_index, test_index in skf.split(x, y):
    print('TRAIN:', train_index, "TEST:", test_index)
    x_train, x_test = x.iloc[train_index, :], x.iloc[test_index, :]
    y_train, y_test = y[train_index], y[test_index]
TRAIN: [0 1 4 5] TEST: [2 3]
TRAIN: [1 2 3 5] TEST: [0 4]
TRAIN: [0 2 3 4] TEST: [1 5]
"""


"""
try:
    if int(k)
except ValueError:
    print "Error: invalid literal for int() with base"
else:
    print "内容写入文件成功"
    fh.close()

num = 0
for i in ml_features_neopep:
    #print(i)
    num += 1
    cor = round((df_patient_fill_norm_comb_cat[str(i)].astype('float16')).corr(Neopep_data_ml_norm_df[str(i)].astype('float16')),6)
    print(f'{num}\t{i}\t{cor}')
"""