# coding=gbk
from param_value import *

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

def test_voting_classifier(classifiers: list, weights: list, patient: str, data: pd.DataFrame,
                           x: pd.DataFrame, y: np.ndarray, report_file: str = None, sort_columns: list = []) \
        -> tuple:
    """
    Tests classifier on given patient's data

    Args:
        classifiers (list): classifiers to be combined into voting classifier and tested. Classifier objects must
                            implement the sklearn predict_proba method
        weights (list(float)): weights of each classifier in classifiers. weights must be positive.
        data (pd.DataFrame): unprocessed dataframe with rows selected for ML
        x (pd.DataFrame): processed dataframe with rows and columns selected for ML
        y (np.ndarray): 0/1 array indicating immunogenicity (value == 1)
        report_file (str): file name to write results to. If None nothing is written
        sort_columns (list): additional sort columns to resolve ties in predict_proba sorting
    Returns:
        (predicted_proba, x_sorted_annot, nr_correct20, nr_tested20, nr_correct50, nr_tested50, nr_correct100,
         nr_tested100, nr_immuno, ranks, score)
            predicted_proba: predicted immunogenic probability
            x_sorted_annot: x sorted by predict_proba probability with additional columns:
                            ML_pred: predict_proba values
                            response: response (0/1), Corresponds to sorted vector y
                            mutant_seq: peptide with mutation
                            gene: mutated gene containing mutant_seq
            nr_correct20: nr of immunogenic peptides in top 20 ranks after sorting
            nr_tested20: nr of tested peptides in top 20 ranks after sorting
            nr_correct50: nr of immunogenic peptides in top 50 ranks after sorting
            nr_tested50: nr of tested peptides in top 50 ranks after sorting
            nr_correct100: nr of immunogenic peptides in top 100 ranks after sorting
            nr_tested100: nr of tested peptides in top 100 ranks after sorting
            nr_immuno: nr immunogenic peptides for this patient
            ranks: ranks of the immunogenic peptide of this patients
            score: score for the ranking of this patient (as defined by scorer_name in __init__)
    """

    classifier_scorer = get_scorer(scorer_name, data)
    write_header = True
    if verbose > 1 and write_header:
        print("Patient\tNr_correct_top20\tNr_tested_top20\tNr_correct_top30\tNr_tested_top30\tNr_correct_top50\tNr_tested_top50\t"
              "Nr_correct_top100\tNr_tested_top100\tNr_immunogenic\tNr_peptides\tClf_score\t"
              "CD8_ranks\tCD8_mut_seqs\tCD8_genes")
    if report_file and os.path.getsize(report_file.name) == 0:
        report_file.write("Patient\tNr_correct_top20\tNr_tested_top20\tNr_correct_top30\tNr_tested_top30\tNr_correct_top50\tNr_tested_top50\t"
                          "Nr_correct_top100\tNr_tested_top100\tNr_immunogenic\tNr_peptides\tClf_score\t"
                          "CD8_ranks\tCD8_mut_seqs\tCD8_genes\n")

    write_header = False
    y_pred = np.full(len(y), 0.0)
    for (w, clf) in zip(weights, classifiers):
        #print(clf)
        y_pred = np.add(y_pred, np.array(clf[1].predict_proba(x)[:, 1]) * w)

    X_r = x.copy()
    X_r['ML_pred'] = y_pred
    X_r['response'] = y
    X_r['response_type'] = data['response_type']
    X_r.loc[:, 'gene'] = data.loc[:, 'gene']
    X_r.loc[:, 'mutant_seq'] = data.loc[:, 'mutant_seq']
    #for c in sort_columns:
        #if GlobalParameters().get_order_relation(c) == '<':
        #if ml_feature_mv_neopep[c] == 'max':
            #X_r.loc[:, c] = -X_r.loc[:, c]
    sort_columns = ['ML_pred'] + sort_columns
    X_r = X_r.sort_values(by=sort_columns, ascending=False)

    r = np.where(X_r['response'] == 1)[0]
    rt = np.where(X_r['response_type'] == 'negative')[0]
    nr_correct20 = sum(r < 20)
    nr_tested20 = nr_correct20 + sum(rt < 20)
    nr_correct30 = sum(r < 30)
    nr_tested30 = nr_correct30 + sum(rt < 30)
    nr_correct50 = sum(r < 50)
    nr_tested50 = nr_correct50 + sum(rt < 50)
    nr_correct100 = sum(r < 100)
    nr_tested100 = nr_correct100 + sum(rt < 100)
    nr_immuno = sum(y == 1)
    score = classifier_scorer._score_func(y, y_pred)
    sort_idx = np.argsort(r)
    ranks_str = ",".join(["{0:.0f}".format(np.floor(r+1)) for r in r[sort_idx]])
    mut_seqs = X_r.loc[X_r['response'] == 1, 'mutant_seq'].to_numpy()
    mut_seqs_str = ",".join(["{0}".format(s) for s in mut_seqs[sort_idx]])
    genes = X_r.loc[X_r['response'] == 1, 'gene'].to_numpy()
    gene_str = ",".join(["{0}".format(s) for s in genes[sort_idx]])

    if verbose > 1:
        print("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s" %
              (patient, nr_correct20, nr_tested20, nr_correct30, nr_tested30, nr_correct50, nr_tested50, nr_correct100, nr_tested100,
               nr_immuno, len(y), score, ranks_str, mut_seqs_str, gene_str))

    if report_file:
        report_file.write("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s\n" %
                          (patient, nr_correct20, nr_tested20, nr_correct30, nr_tested30, nr_correct50, nr_tested50,
                           nr_correct100, nr_tested100, nr_immuno, len(y), score, ranks_str,
                           mut_seqs_str, gene_str))
        report_file.flush()

    return X_r['ML_pred'], X_r, nr_correct20, nr_tested20, nr_correct30, nr_tested30, nr_correct50, nr_tested50, nr_correct100, nr_tested100,\
           nr_immuno, r, score



#@staticmethod
def write_tesla_scores(patient: str, dataset: str, nr_correct20: int, nr_tested20: int, nr_correct100: int,
                       nr_immuno: int, x_sorted: pd.DataFrame, y_pred_sorted: np.ndarray, report_file: str):
    """
    Calculates TESLA FP, TTIF, and AUPRC scores for a patient and writes them to report file

    Args:
        patient (str): patient id.
        dataset (bool): dataset id. if not provided all datasets are considered
        nr_correct20 (int): nr immunogenic peptides ranked in top 20 for this patient
        nr_tested20 (int): nr tested peptides ranked in top 20 for this patient
        nr_correct100 (int): nr immunogenic peptides ranked in top 100 for this patient
        nr_immuno (int): nr immunogenic peptides for this patient
        x_sorted (pd.DataFrame): ml data matrix sorted by predict_proba
        y_pred_sorted (np.ndarray): sorted response types (1 = immunogenic)
        report_file (str): report file
    """
    idx = x_sorted['response_type'] != 'not_tested'
    y_pred_tesla = y_pred_sorted[idx].to_numpy()
    y_tesla = x_sorted.loc[idx, 'response'].to_numpy()
    ttif = nr_correct20/nr_tested20 if nr_tested20 > 0 else 0
    fr = nr_correct100/nr_immuno if nr_immuno > 0 else 0
    precision, recall, _ = precision_recall_curve(y_tesla, y_pred_tesla)
    auprc = auc(recall, precision)
    report_file.write("{0}\t{1}\t{2:.3f}\t{3:.3f}\t{4:.3}\n".format(dataset, patient, ttif, fr, auprc))

"""
immunogenic_patients: dict = {'mutation': None, 'neopep': None}
def has_immunogenic_peptides(peptide_type: str, patient: str) -> bool:
    Checks whether patient has immunogenic mutations or neo-peptides
    Args:
        peptide_type (str): either 'neopep' or 'mutation'
        patient (str, optional): patient id.
    Returns:
        True if patient has immunogenic mutations or neo-peptides, False otherwise
    patients = data_test['patient'].unique()
    immunogenic_patients[peptide_type] = []
    for p in patients:
        data_p = data_test[data_test['patient'] == p]
        if sum(data_p['response_type'] == 'CD8') > 0:
            immunogenic_patients[peptide_type].append(p)
    return patient in immunogenic_patients[peptide_type]
"""

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


classifier1_model_files = glob.glob(os.path.join(classifier_model_dir, classifier1_file_re))
classifier2_model_files = glob.glob(os.path.join(classifier_model_dir, classifier2_file_re))
vc_result_file = os.path.join(classifier_result_dir, 'Voting_classifier_{0}_{1:.2f}_test.txt'.format(peptide_type,weight))
open(vc_result_file, mode='w').close()

tesla_score_file = os.path.join(classifier_result_dir, 'Voting_classifier_{0}_{1:.2f}_tesla_scores.txt'.format(peptide_type,weight))
with open(tesla_score_file, mode='w') as score_file:
    score_file.write("Dataset\tPatient\tTTIF\tFR\tAUPRC\n")


def count_CD8_ranks(row, rank):
    count = 0
    for i in row['CD8_ranks'].split(','):
        if int(i) <= rank:
            count += 1
    return count
voting_clfs = []
for patient in sorted(data_test['patient'].unique()):
    if not sum(data_test[data_test['patient'] == patient]['response_type'] == 'CD8') > 0:
        continue
    #print(f'patient {patient}')
    #if not has_immunogenic_peptides(peptide_type, patient):
        #continue
    #data_test, X_test, y_test = \
        #DataManager.filter_processed_data(peptide_type=args.peptide_type, objective='ml', patient=patient,
                                          #sample=False)
    #data_test_patient = data_test[(data_test['patient'] == patient) & (data_test['dataset'] != 'NCI')]
    data_test_patient = data_test[data_test['patient'] == patient]
    X_test_patient = X_test.loc[data_test_patient.index,:]
    y_test_pd_patient = y_test_pd.loc[data_test_patient.index,:]
    y_test_patient = np.array(y_test_pd_patient).reshape(-1)
    dataset = data_test_patient['patient'].apply(lambda i: "TESLA" if i.startswith('TESLA') else ("HiTIDE" if i.startswith('Patient') else "NCI")).tolist()[0]
    #print(dataset)

    weights = []
    for model_file in classifier1_model_files:
        clf_name = os.path.basename(model_file).split("_")[0]
        #clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, clf_name, X_test)
        #classifier = clf_mgr.load_classifier(clf_name, clf_mgr.get_optimization_params(), model_file)
        classifier1 = load_classifier(clf_name, model_file)
        voting_clfs.append((clf_name, classifier1))
        weights.append(1.0-weight)
    for model_file in classifier2_model_files:
        clf_name = os.path.basename(model_file).split("_")[0]
        #clf_mgr = get_clf_mgr(args.peptide_type, args.dataset_train, clf_name, X_test)
        #classifier = clf_mgr.load_classifier(clf_name, clf_mgr.get_optimization_params(), model_file)
        classifier2 = load_classifier(clf_name, model_file)
        voting_clfs.append((clf_name, classifier2))
        weights.append(weight)
    with open(vc_result_file, mode='a') as result_file:
        y_pred_sorted, X_sorted, nr_correct20, nr_tested20, nr_correct30, nr_tested30, nr_correct50, nr_tested50, nr_correct100, \
        nr_tested100, nr_immuno, r, score = \
            test_voting_classifier(voting_clfs, weights, patient, data_test_patient, X_test_patient, y_test_patient,
                                           report_file=result_file)
        #if peptide_type == 'neopep':
    with open(tesla_score_file, mode='a') as score_file:
        write_tesla_scores(patient, dataset, nr_correct20, nr_tested20, nr_correct100,
                                             nr_immuno, X_sorted, y_pred_sorted, score_file)


result_file_pd = pd.read_csv(vc_result_file, sep='\t',low_memory=False)
result_file_pd['Dataset'] = result_file_pd['Patient'].apply(lambda i: "TESLA" if i.startswith('TESLA') else ("HiTIDE" if i.startswith('Patient') else "NCI"))
print(vc_result_file)
result_file_pd['Nr_correct_top30'] = result_file_pd.apply(lambda row: count_CD8_ranks(row, 30), axis=1)
result_file_pd['Nr_correct_top50'] = result_file_pd.apply(lambda row: count_CD8_ranks(row, 50), axis=1)
result_file_pd['Nr_correct_top100'] = result_file_pd.apply(lambda row: count_CD8_ranks(row, 100), axis=1)
print(result_file_pd.groupby(['Dataset'])[["Nr_correct_top20","Nr_correct_top30","Nr_correct_top50","Nr_correct_top100"]].sum())

