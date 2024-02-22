# coding=gbk

from param_value import *
#from training import sum_rank_correct
#training.py中的from param_value import *会影响from training import sum_rank_correct
#另外，sum_rank_correct中有np也会影响from training import sum_rank_correct


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

def nr_correct_top30(y_true, y_pred, max_rank=30):
    n = min(len(y_true), max_rank)
    idx = np.argsort(-y_pred)
    return np.sum(y_true[idx][:n] == 1)

def load_classifier(classifier_tag, classifier_file) -> object:
    """
    Saves classifier to file.
    Args:
        classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
        optimization_params (OptimizationParams): OptimizationParams object
        classifier_file (str): file name of classifier model file
    Returns:
        classifier object
    """
    if classifier_tag in ['CatBoost', 'XGBoost']:
        classifier = get_base_classifier(classifier_tag)
        classifier.load_model(classifier_file)
    elif classifier_tag == "RF":
        classifier = joblib.load(classifier_file)
    else:
        classifier = pickle.load(open(classifier_file, 'rb'))
    return classifier


def get_categorical_feature_idx(peptide_type: str, x_: pd.DataFrame) -> list:
    """
    Retrieves indexes of columns with categorical features
    Args:
        peptide_type (str):  either 'neopep' or 'mutation'
        x_ (pd.DataFrame): ML dataframe with feature columns in order
    Returns:
        list of indexes of columns with categorical features
    """
    if peptide_type == 'neopep':
        idx = [i for i, c in enumerate(x_.columns) if feature_types_neopep[c] == 'category']
    else:
        idx = [i for i, c in enumerate(x_.columns) if feature_types_mutation[c] == 'category']
    return idx

def get_base_classifier(classifier_tag):
    """
    Creates classifier object with default hyperparameters

    Args:
        classifier_tag (str): tag of classifier ('SVM', 'SVM-lin', 'LR', 'XGBoost', 'CatBoost')
    Returns:
        classifier object (Object): classifier object with default hyperparameters
    """
    if classifier_tag == "SVM":
        return SVC(probability=True, kernel='rbf')
    elif classifier_tag == "SVM-lin":
        return SVC(probability=True, kernel='linear')
    elif classifier_tag == "LR":
        return LogisticRegression(solver='saga')
    elif classifier_tag == "RF":
        return RandomForestClassifier( )
    elif classifier_tag == "CatBoost":
        cb = CatBoostClassifier(
            loss_function='Logloss',
            iterations=10,
            learning_rate=0.01,
            random_seed=42,
            logging_level='Silent',
            use_best_model=False,
            cat_features=cat_idx,
            auto_class_weights='None',
            silent=True
        )
        return cb
    elif classifier_tag == "XGBoost":
        clf_xgb = XGBClassifier(
            enable_categorical=True,
            max_depth=8,
            learning_rate=0.1,
            n_estimators=1000,
            verbosity=0,
            silent=None,
            objective='binary:logistic',
            booster='gbtree',
            tree_method='hist',
            n_jobs=int(os.cpu_count()/nr_hyperopt_rep),
            nthread=None,
            gamma=0,
            min_child_weight=1,
            max_delta_step=0,
            subsample=0.7,
            colsample_bytree=1,
            colsample_bylevel=1,
            colsample_bynode=1,
            reg_alpha=0,
            reg_lambda=1,
            scale_pos_weight=1,
            base_score=0.5,
            random_state=0,
            seed=None)
        return clf_xgb


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


def test_classifier(classifier, peptide_type: str, patient: str, data: pd.DataFrame, x, y, max_rank=20,
                    report_file=None, sort_columns=[]) -> tuple:
    """
    Tests classifier on given patient's data
    Args:
        classifier (object): classifier to be tested. Classifier object must implement the sklearn predict_proba
                             method
        data (pd.DataFrame): unprocessed dataframe with rows selected for ML
        x (pd.DataFrame): processed dataframe with rows and columns selected for ML
        y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
        max_rank (int): number of ranks taken into account to calculate topN counts. only needed dor report
                        output, but does not affect testing results
        report_file (str): file name to write results to. If None nothing is written
        sort_columns (list): additional sort columns to resolve ties in predict_proba sorting
    Returns:
        (predicted_proba, x_sorted_annot, nr_correct, nr_immuno, ranks, score)
            predicted_proba: predicted immunogenic probability
            x_sorted_annot: x sorted by predict_proba probability with additional columns:
                            ML_pred: predict_proba values
                            response: response (0/1), Corresponds to sorted vector y
                            mutant_seq: peptide with mutation
                            gene: mutated gene containing mutant_seq
            nr_correct: nr of immunogenic peptides in top max_rank ranks after sorting
            nr_immuno: nr immunogenic peptides for this patient
            ranks: ranks of the immunogenic peptide of this patients
            score: score for the ranking of this patient (as defined by scorer_name in __init__)
    """
    classifier_scorer = get_scorer(scorer_name, data)
    write_header = True
    if verbose > 1 and write_header:
        print("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
              "CD8_ranks\tCD8_mut_seqs\tCD8_genes".format(max_rank))
    if report_file is not None:
        lock = FileLock(report_file+".lock")
        with lock:
            if os.path.getsize(report_file) == 0:
                with open(report_file, mode='w') as file:
                    file.write("Patient\tNr_correct_top{0}\tNr_immunogenic\tMax_rank\tNr_peptides\tClf_score\t"
                               "CD8_ranks\tCD8_mut_seqs\tCD8_genes\n".format(max_rank))
    write_header = False
    y_pred = classifier.predict_proba(x)[:, 1]    #classifier.predict_proba(x).shape: (100082, 2): probability of class 0 and class 1
    X_r = x.copy()
    X_r['ML_pred'] = y_pred
    X_r['response'] = y
    X_r.loc[:, 'gene'] = data.loc[:, 'gene']
    X_r.loc[:, 'mutant_seq'] = data.loc[:, 'mutant_seq']
    data['ML_pred'] = y_pred
    data_sorted = data.sort_values(by=['ML_pred'], ascending=False)
    for c in sort_columns:
        if peptide_type == 'neopep':
            if ml_feature_mv_neopep[c] == 'max':
                X_r.loc[:, c] = -X_r.loc[:, c]
        elif peptide_type == 'mutation':
            if ml_feature_mv_mutation[c] == 'max':
                X_r.loc[:, c] = -X_r.loc[:, c]
    sort_columns = ['ML_pred'] + sort_columns
    X_r = X_r.sort_values(by=sort_columns, ascending=False)
    r = np.where(X_r['response'] == 1)[0]
    nr_correct = sum(r < max_rank)
    nr_immuno = sum(y == 1)
    score = classifier_scorer._score_func(y, y_pred)
    sort_idx = np.argsort(r)
    ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
    mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
    mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
    genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
    gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])
    if verbose > 1:
        print("%s\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s" %
              (patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), score, ranks_str,
               mut_seqs_str, gene_str))
    if report_file is not None:
        lock = FileLock(report_file+".lock")
        with lock:
            with open(report_file, mode='a') as file:
                file.write("{0:s}\t{1:d}\t{2:d}\t{3:d}\t{4:d}\t{5:.5f}\t{6:s}\t{7:s}\t{8:s}\n".
                           format(patient, nr_correct, nr_immuno, np.min((max_rank, len(y))), len(y), score,
                                  ranks_str, mut_seqs_str, gene_str))
    return data_sorted, X_r['ML_pred'], X_r, nr_correct, nr_immuno, r, score



if user_defined_split == True:
    data_test = pd.read_csv(user_defined_test_file, sep='\t',low_memory=False)
    X_test = data_test[ml_features_neopep]
    y_test = np.array(data_test.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
    y_test_pd = pd.DataFrame(y_test)
else:
    preprocessed_data = pd.read_csv(preprocessed_data_file, sep='\t',low_memory=False)
    data_test = preprocessed_data[preprocessed_data['train_test'] == 'test']
    X_test = data_test[ml_features_neopep]
    y_test = np.array(data_test.response_type.apply(lambda rt: int(rt == 'CD8')), dtype=int)
    y_test_pd = pd.DataFrame(y_test)
    y_test_pd.index = data_test.index


cat_idx = get_categorical_feature_idx(peptide_type, X_test)

classifier_model_files = glob.glob(os.path.join(classifier_model_dir, classifier_file_re))
classifier_result_files = []
#for wc in classifier_file_re:
    #classifier_model_files = \
        #classifier_model_files + glob.glob(os.path.join(classifier_model_dir, wc))
for classifier_model_file in classifier_model_files:
    #result_file_name = re.sub("_clf\\.\\w+$", "_test.txt", os.path.basename(classifier_model_file))
    result_file_name = re.sub(r"_clf\.\w+$", "_test.txt", os.path.basename(classifier_model_file))
    result_file_name = os.path.join(classifier_result_dir, result_file_name)
    classifier_result_files.append(result_file_name)
    open(result_file_name, mode='w').close()


def count_CD8_ranks(row, rank):
    count = 0
    for i in row['CD8_ranks'].split(','):
        if int(i) <= rank:
            count += 1
    return count
#data_sorted_all = pd.DataFrame()
#for patient in sorted(data_test.patient.unique()):
for model_file, result_file in zip(classifier_model_files, classifier_result_files):
    #if not DataManager.has_immunogenic_peptides(peptide_type, patient):
        #continue
    #data_test, X_test, y_test = \
        #DataManager.filter_processed_data(peptide_type=peptide_type, objective='ml', patient=patient,
                                          #sample=False)
    #data_test_patient, X_test, y_test三者都需要按照样本分割
    #for model_file, result_file in zip(classifier_model_files, classifier_result_files):
    for patient in sorted(data_test['patient'].unique()):
        if not sum(data_test[data_test['patient'] == patient]['response_type'] == 'CD8') > 0:
            continue
        #print(patient)
        data_test_patient = data_test[data_test['patient'] == patient]
        X_test_patient = X_test.loc[data_test_patient.index,:]
        y_test_pd_patient = y_test_pd.loc[data_test_patient.index,:]
        y_test_patient = np.array(y_test_pd_patient).reshape(-1)
        #clf_name = os.path.basename(model_file).split("_")[0]
        #clf_mgr = get_clf_mgr(peptide_type, dataset_train, clf_name, X_test)
        classifier = load_classifier(classifier_tag, model_file)
        data_sorted, y_pred_sorted, X_sorted, nr_correct, nr_immuno, r, score = \
            test_classifier(classifier, peptide_type, patient, data_test_patient, X_test_patient, y_test_patient,
                                    report_file=result_file)
        #data_sorted_all = pd.concat([data_sorted_all, data_sorted])

    result_file_pd = pd.read_csv(result_file, sep='\t',low_memory=False)
    result_file_pd['Dataset'] = result_file_pd['Patient'].apply(lambda i: "TESLA" if i.startswith('TESLA') else ("HiTIDE" if i.startswith('Patient') else "NCI"))
    print(result_file)
    result_file_pd['Nr_correct_top30'] = result_file_pd.apply(lambda row: count_CD8_ranks(row, 30), axis=1)
    result_file_pd['Nr_correct_top50'] = result_file_pd.apply(lambda row: count_CD8_ranks(row, 50), axis=1)
    result_file_pd['Nr_correct_top100'] = result_file_pd.apply(lambda row: count_CD8_ranks(row, 100), axis=1)
    print(result_file_pd.groupby(['Dataset'])[["Nr_correct_top20","Nr_correct_top30","Nr_correct_top50","Nr_correct_top100"]].sum())
