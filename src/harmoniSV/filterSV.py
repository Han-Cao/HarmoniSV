#!/usr/bin/env python3

# Random forest filter for SVs
# Created: 12/9/2022
# Author: Han Cao

import argparse
import logging
import pickle

import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RandomizedSearchCV
from sklearn.metrics import roc_auc_score, precision_score
from sklearn.utils import resample

# # SVM
# from sklearn.pipeline import Pipeline
# from sklearn.preprocessing import StandardScaler
# from sklearn.svm import SVC

from utils import read_vcf, vcf_to_df, parse_cmdargs

# parse arguments
parser = argparse.ArgumentParser(prog="harmonisv filter",
                                 description="Random forest filter for SVs", 
                                 add_help=False)
required = parser.add_argument_group('required arguments')
required.add_argument("-i", "--invcf", metavar="vcf", type=str, required=True,
                      help="input vcf file")
required.add_argument("-o", "--output", metavar="PREFIX", type=str, required=True,
                      help="ouput prefix")
required.add_argument("--sv-type", metavar="TYPE", type=str, required=True,
                      help="which SVTYPE stored in INFO/SVTYPE to work on")
required.add_argument("--feature", metavar="tags", type=str, required=True,
                      help="comma-separated INFO tags used as training features")

train_group_args = parser.add_argument_group('training SV group arguments')
train_group_args.add_argument("--merge", metavar="FILE", type=str, required=False,
                              help="SV merge results file, one line represent merged SVs separated by comma")
train_group_args.add_argument("--target-prefix", metavar="PREFIX", type=str, required=False,
                              help="target ID prefix for SV call set")
train_group_args.add_argument("--train-size", metavar="0.9", type=float, required=False, default=0.9,
                              help="proportion of SVs used for training (default: 0.9)")
train_group_args.add_argument("--max-train", metavar="100000", type=int, required=False,
                              help="maximum number of SVs used for training")
train_group_args.add_argument("--min-support-set", metavar="N", type=int, required=False, default=4,
                              help="minimum number of supporting SV call set for positive SVs (default: 4)")
train_group_args.add_argument("--min-re", metavar="N", type=int, required=False, default=None,
                              help="minimum value of MAX_RE for non-negative SVs (default: None)")
train_group_args.add_argument("--min-support-caller", metavar="N", type=int, required=False, default=None,
                              help="minimum number of SUPP_CALLER to be excluded from negative SV in training (default: None)")
train_group_args.add_argument("--min-support-method", metavar="N", type=int, required=False, default=None,
                              help="minimum number of SUPP_METHOD to be excluded from negative SV in training, overide --min-support-caller (default: None)")
train_group_args.add_argument("--train-sites", metavar="FILE", type=str, required=False,
                              help="user-specified positive and negative SV sites, tab-separated columns: SV_ID (without sample prefix), label (1: true, 0: false)")
train_group_args.add_argument("--bench-vcf", metavar="FILE", type=str, required=False,
                              help="benchmark SV calling vcf file, if not given, will search benchmark SV in --invcf")
train_group_args.add_argument("--bench-sites", metavar="FILE", type=str, required=False,
                              help="ground truth SV for benchmark, tab-separated columns: SV_ID (with sample prefix) and label (1: true, 0: false)")


model_args = parser.add_argument_group('model arguments')
model_args.add_argument("--n-estimator", metavar="N", type=int, required=False, default=300,
                        help="number of trees in the random forest (default: 300)")
model_args.add_argument("--k-fold", metavar="N", type=int, required=False, default=10,
                        help="number of folds for random search cross validation model selection (default: 10)")
model_args.add_argument("--n-iter", metavar="N", type=int, required=False, default=100,
                        help="number of iterations for random search cross validation model selection (default: 100)")
model_args.add_argument("--max-features", metavar="str", type=str, required=False, default="auto",
                        help="RandomForestClassifier parameter max_features, format: single value or comma-separated list  or from:to:step (default: 1/2 * sqrt(n_features) to 2/3 * n_features)")
model_args.add_argument("--max-depth", metavar="str", type=str, required=False, default="7:14:1",
                        help="RandomForestClassifier parameter max_depth, format: single value or comma-separated list  or from:to:step (default: 7:14:1)")
model_args.add_argument("--min-samples-leaf", metavar="str", type=str, required=False, default="1,2,5,10,20,50",
                        help="RandomForestClassifier parameter min_samples_leaf, format: single value or comma-separated list  or from:to:step (default: 1,2,5,10,20,50)")
model_args.add_argument("--apply-model", metavar="FILE", type=str, required=False,
                        help="skip training and apply pre-trained model to vcf")
model_args.add_argument("--model-thread", metavar="N", type=int, required=False, default=4,
                        help="number of threads for RandomForestClassifier (default: 4). Final requested number of cores is --model-thread * --cv-thread")
model_args.add_argument("--cv-thread", metavar="N", type=int, required=False, default=4,
                        help="number of threads for cross validation (default: 4). Final requested number of cores is --model-thread * --cv-thread")

optional_arg = parser.add_argument_group('optional arguments')
optional_arg.add_argument("--seed", metavar="N", type=int, required=False, default=42,
                          help="random seed (default: 42)")
optional_arg.add_argument("-h", "--help", action='help', help='show this help message and exit')


def parse_hyperparam(param_str):
    """Parse hyperparameter string"""
    if ',' in param_str:
        param_list = [int(x) for x in param_str.split(',')]
    elif ':' in param_str:
        start, end, step = [int(x) for x in param_str.split(':')]
        param_list = list(range(start, end + 1, step))
    else:
        param_list = [int(param_str)]
    
    return param_list


def hyperparam_range(n_features: int, max_features: str, max_depth: str, min_samples_leaf: str) -> dict:
    """Get hyperparameter range for random search"""
    if max_features == "auto":
        max_features = f"{int(np.sqrt(n_features) * 0.5)}:{int(n_features * 2/3)}:1"
    
    param_dict = {
        'max_features': parse_hyperparam(max_features),
        'max_depth': parse_hyperparam(max_depth),
        'min_samples_leaf': parse_hyperparam(min_samples_leaf)
    }

    return param_dict


def extract_uniq_param(param_dict: dict) -> dict:
    """If only one set of parameters in the dictionary, convert list to single value"""
    new_param_dict = {}
    for param in param_dict:
        if len(param_dict[param]) > 1:
            return None
        else:
            new_param_dict[param] = param_dict[param][0]
    
    return new_param_dict
    


def merge_list_to_df(merge_list: list, target_id: str) -> pd.DataFrame:
    """Summarise SV occurrences in merge list to df"""
    target_merge = [x for x in merge_list if target_id in x]
    merge_dict_list = []
    for merge in target_merge:
        new_record = {}
        for sv_id in merge.split(','):
            set_id = sv_id.split('.')[0]
            if target_id == set_id:
                new_record['ID'] = sv_id
            else:
                new_record[set_id] = 1
        merge_dict_list.append(new_record)

    df_merge = pd.DataFrame(merge_dict_list).set_index('ID')
    df_merge = df_merge.fillna(0)
    df_merge['support_set'] = df_merge.sum(axis=1)

    return df_merge


def read_known_sv(file: str) -> pd.DataFrame:
    """Read known positive and negative SVs"""
    logger = logging.getLogger(__name__)
    df = pd.read_csv(file, sep='\t', header=None)
    df.columns = ['ID', 'label']
    df = df.set_index('ID')
    if not df['label'].isin([0, 1]).all():
        logger.error(f"Invalid label in {file}, only accept 0 for negative and 1 for positive")
        raise SystemExit()

    return df


def get_df_vcf(vcf: pysam.VariantFile, 
               sv_type: str,
               feature_list: list, 
               merge_list: list=None, 
               target_id: str=None, 
               min_support_set: int=None, 
               min_support_caller: int=None,
               min_support_method: int=None,
               min_re: int=None,
               pos_set: set=None, 
               neg_set: set=None,
               feature_only: bool=False) -> pd.DataFrame:
    """Get dataframe of features and labels for training"""
    # convert vcf to df
    if isinstance(feature_list, list):
        col_list = feature_list.copy()
        if 'SVTYPE' not in col_list:
            col_list.append('SVTYPE')
        if 'AC' not in col_list:
            col_list.append('AC')
        if 'SUPP_CALLER' not in col_list:
            col_list.append('SUPP_CALLER')
        if 'SUPP_METHOD' not in col_list:
            col_list.append('SUPP_METHOD')
        if 'MAX_RE' not in col_list:
            col_list.append('MAX_RE')
    else:
        col_list = feature_list
        
    df_vcf = vcf_to_df(vcf, col_list)
    df_vcf = df_vcf[(df_vcf['SVTYPE'] == sv_type) & (df_vcf['AC'] > 0)]
    df_vcf = df_vcf.set_index('ID')
    df_vcf['ID_REPRESENT'] = df_vcf.index.str.replace(r'^.*?\.', '', regex=True)
    df_vcf['Sample'] = df_vcf.index.str.replace(r'\..*$', '', regex=True)

    if not feature_only:
        df_vcf['label'] = -9

        # get positive and negative SVs
        if merge_list is not None:
            df_merge = merge_list_to_df(merge_list, target_id)
            merge_pos = set(df_merge.index[df_merge['support_set'] >= min_support_set])
            merge_novel = set(df_merge.index[df_merge['support_set'] == 0])
            df_vcf = df_vcf.merge(df_merge['support_set'], left_on='ID_REPRESENT', right_index=True, how='left')
            df_vcf.loc[df_vcf['ID_REPRESENT'].isin(merge_pos), 'label'] = 1
            if min_support_method is None:
                df_vcf.loc[df_vcf['ID_REPRESENT'].isin(merge_novel) & (df_vcf['SUPP_CALLER'] < min_support_caller), 'label'] = 0
            else:
                df_vcf.loc[df_vcf['ID_REPRESENT'].isin(merge_novel) & (df_vcf['SUPP_METHOD'] < min_support_method), 'label'] = 0
        
        if pos_set is not None and neg_set is not None:
            df_vcf.loc[df_vcf['ID_REPRESENT'].isin(pos_set), 'label'] = 1
            df_vcf.loc[df_vcf['ID_REPRESENT'].isin(neg_set), 'label'] = 0
        
        if min_re is not None:
            df_vcf.loc[df_vcf['MAX_RE'] < min_re, 'label'] = 0

    return df_vcf


def even_train_test_split(X: np.ndarray, y: np.ndarray, train_size: float, max_train: int = 100000, random_state: int=42) -> tuple:
    """
    Split imbalanced data into balanced training data and test data
    Return: X_train, X_test, y_train, y_test
    """
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    min_class_n = min(len(pos_idx), len(neg_idx))
    if min_class_n > max_train:
        min_class_n = max_train

    if len(pos_idx) > min_class_n:
        pos_idx = resample(pos_idx, replace=False, n_samples=min_class_n, random_state=random_state)
    elif len(neg_idx) > min_class_n:
        neg_idx = resample(neg_idx, replace=False, n_samples=min_class_n, random_state=random_state)
    
    
    train_pos_idx = resample(pos_idx, replace=False, n_samples=int(min_class_n * train_size), random_state=random_state)
    train_neg_idx = resample(neg_idx, replace=False, n_samples=int(min_class_n * train_size), random_state=random_state)
    train_idx = np.concatenate((train_pos_idx, train_neg_idx))

    X_train = X[train_idx]
    y_train = y[train_idx]

    if train_size == 1:
        X_test = None
        y_test = None
    else:
        X_test = np.delete(X, train_idx, axis=0)
        y_test = np.delete(y, train_idx, axis=0)
        keep_idx = np.where(y_test != -9)[0]
        X_test = X_test[keep_idx]
        y_test = y_test[keep_idx]
    
    return X_train, X_test, y_train, y_test
    

def stat_by_cutoff(y_true: np.ndarray, y_pred: np.ndarray, cutoff_list: list, stat) -> pd.DataFrame:
    """Calculate statitics by y_pred cutoffs"""
    stat_list = []
    for cutoff in cutoff_list:
        y_pred_cutoff = np.where(y_pred >= cutoff, 1, 0)
        stat_list.append({'cutoff': cutoff,
                          'TP': sum(y_true[y_pred_cutoff == 1]),
                          'FP': sum(y_pred_cutoff[y_true == 0]),
                          'stat': stat(y_true, y_pred_cutoff)})
    return pd.DataFrame(stat_list)


def count_by_cutoff(df_vcf: pd.DataFrame, y_pred: np.ndarray, cutoff_list: list, prefix='N') -> pd.DataFrame:
    """Count per sample number of SVs by y_pred cutoffs"""
    stat_list = []
    for cutoff in cutoff_list:
        y_pred_cutoff = y_pred >= cutoff
        n_sv = df_vcf.loc[y_pred_cutoff].value_counts('Sample')
        stat_list.append({'cutoff': cutoff, 
                          f'{prefix}_mean_sv': n_sv.mean(),
                          f'{prefix}_min_sv': n_sv.min(),
                          f'{prefix}_max_sv': n_sv.max()})
    return pd.DataFrame(stat_list)


def density_by_group(score, label):
    df_plot = pd.DataFrame({'score': score, 'label': label})
    fig, ax = plt.subplots(figsize=(8,6))
    for group, df in df_plot.groupby('label'):
        df['score'].plot.kde(ax=ax, label=group)
    ax.legend()

    return fig, ax


def cv_table(cv_results: dict) -> pd.DataFrame:
    """Convert cv_results to table"""
    df = pd.DataFrame(cv_results)
    keep_col = [col for col in df.columns if col.startswith('param_')] + ['mean_test_score', 'std_test_score', 'rank_test_score']
    df = df[keep_col]
    df = df.sort_values('rank_test_score')
    df = df.reset_index(drop=True)
    return df


def select_model(model: RandomForestClassifier,
                 X: np.ndarray,
                 y: np.ndarray,
                 out_prefix: str,
                 k_fold: int,
                 n_iter: int,
                 param_dict:dict,
                 threads: int,
                 random_state: int) -> RandomForestClassifier:
    """Train model and save statistics"""
    logger = logging.getLogger(__name__)

    clf = RandomizedSearchCV(model, param_dict, n_iter=n_iter, cv=k_fold, 
                             n_jobs=threads, random_state=random_state, scoring='roc_auc', verbose=2)
    logger.info('Start model selection')
    search = clf.fit(X, y)
    logger.info(f'Best parameters: {search.best_params_}')
    
    df_cv = cv_table(search.cv_results_)
    df_cv.to_csv(f'{out_prefix}.cv_summary.txt', sep='\t', index=False)
    logger.info(f'Save model selection cross-validation summary to {out_prefix}.cv_summary.txt')
    
    return search.best_estimator_


def test_model(model: RandomForestClassifier, 
               X: np.ndarray, 
               y: np.ndarray, 
               out_prefix: str, 
               cutoff_list = np.arange(0, 1, 0.001)) -> np.ndarray:
    """Apply model and save statistics"""
    logger = logging.getLogger(__name__)

    y_pred = model.predict_proba(X)[:, 1]
    fig, ax = density_by_group(y_pred, y)
    fig.savefig(f'{out_prefix}.score_distribution.pdf')
    logger.info(f'Save score distribution to {out_prefix}.score_distribution.pdf')
    logger.info(f'Model roc auc statistic: {round(roc_auc_score(y, y_pred), 4)}')

    df_stat = stat_by_cutoff(y, y_pred, cutoff_list, precision_score)
    df_stat = df_stat.rename(columns={'stat': 'precision'})
    df_stat.to_csv(f'{out_prefix}.cutoff_summary.txt', sep='\t', index=False)
    logger.info(f'Save statistics by cutoff to {out_prefix}.cutoff_summary.txt')

    return y_pred


def model_worker(df_train: pd.DataFrame,
                 features: list,
                 out_prefix: str,
                 args: argparse.Namespace,
                 df_bench: pd.DataFrame = None) -> RandomForestClassifier:
    """Train and evaluate model"""
    logger = logging.getLogger(__name__)

    # train model
    hyperparam_dict = hyperparam_range(n_features=len(features), 
                                       max_features=args.max_features, 
                                       max_depth=args.max_depth,
                                       min_samples_leaf=args.min_samples_leaf)
    X_all = df_train[features].to_numpy()
    y_all = df_train['label'].to_numpy()

    X_train, X_test, y_train, y_test = even_train_test_split(X_all, y_all, train_size=args.train_size, random_state=args.seed)
    logger.info(f'Training set: {sum(y_train==1)} positive SV, {sum(y_train==0)} negative SV')
    if y_test is not None:
        logger.info(f'Testing set: {sum(y_test==1)} positive SV, {sum(y_test==0)} negatives SV')

    logger.info('Train random forest model')
    model = RandomForestClassifier(n_estimators=args.n_estimator, random_state=args.seed, n_jobs=args.model_thread)
    
    uniq_hyperparam_dict = extract_uniq_param(hyperparam_dict)
    if uniq_hyperparam_dict is not None:
        logger.info(f'Use user-specified hyperparameters: {uniq_hyperparam_dict}')
        model.set_params(**uniq_hyperparam_dict)
        model = model.fit(X_train, y_train)
    else:
        model = select_model(model=model, 
                            X=X_train, 
                            y=y_train, 
                            out_prefix=out_prefix,
                            k_fold=args.k_fold,
                            n_iter=args.n_iter,
                            param_dict=hyperparam_dict,
                            threads=args.cv_thread,
                            random_state=args.seed)
    
    # evaluate model
    if y_test is not None:
        logger.info('Evaluate model on test set')
        y_pred_test = test_model(model=model, 
                                 X=X_test, 
                                 y=y_test, 
                                 out_prefix=out_prefix + '.test')
    
    if df_bench is not None:
        logger.info('Evaluate model on benchmark set')
        y_pred_bench = test_model(model=model, 
                                  X=df_bench[features].to_numpy(), 
                                  y=df_bench['label'].to_numpy(), 
                                  out_prefix=out_prefix + '.bench')
    
    return model


def write_vcf_filter(invcf: pysam.VariantFile,
                     outvcf: pysam.VariantFile,
                     sv_type: str,
                     score: dict) -> None:
    """Write score and hard filter to vcf"""
    for variant in invcf.fetch():
        new_variant = variant.copy()
        if variant.info['SVTYPE'] == sv_type:
            new_variant.info['RF_SCORE'] = score.get(variant.id, None)
        
        outvcf.write(new_variant)
    
    return None


def apply_vcf(invcf: pysam.VariantFile,
              output: str,
              sv_type: str,
              model: RandomForestClassifier,
              df_vcf: pd.DataFrame,
              features: list) -> None:
    """Apply model to vcf"""

    # apply model to vcf
    df_apply = df_vcf.copy()
    df_apply['RF_SCORE'] = model.predict_proba(df_apply[features].to_numpy())[:,1]
    score = df_apply['RF_SCORE'].round(6).to_dict()

    # write to vcf
    new_header = invcf.header
    if 'RF_SCORE' not in new_header.info:
        new_header.info.add('RF_SCORE', 1, 'Float', 'Random forest probability score')
    outvcf = pysam.VariantFile(output, 'w', header=new_header)

    write_vcf_filter(invcf=invcf,
                     outvcf=outvcf,
                     sv_type=sv_type,
                     score=score)
    
    outvcf.close()

    return None

def filterSV_main(cmdargs):
    """Main function to filter SVs"""
    logger = logging.getLogger(__name__)

    args = parse_cmdargs(parser, cmdargs)

    # read input
    invcf = read_vcf(args.invcf)
    features = args.feature.split(',')

    # TODO: apply model only
    if args.apply_model is not None:
        model = pickle.load(open(args.apply_model, 'rb'))
        logger.info(f'Load model from {args.apply_model}')

        df_vcf = get_df_vcf(vcf=invcf, sv_type=args.sv_type, feature_list=features, feature_only=True)


    if args.merge is not None:
        merge_list = [line.strip() for line in open(args.merge)]
    else:
        merge_list = None
    
    if args.train_sites is not None:
        df_train_sites = read_known_sv(args.train_sites)
        pos_set = set(df_train_sites.index[df_train_sites['label'] == 1])
        neg_set = set(df_train_sites.index[df_train_sites['label'] == 0])
    else:
        pos_set = None
        neg_set = None
    

    df_vcf = get_df_vcf(vcf=invcf,
                        sv_type=args.sv_type,
                        feature_list=features,
                        merge_list=merge_list,
                        target_id=args.target_prefix,
                        min_support_set=args.min_support_set,
                        min_support_caller=args.min_support_caller,
                        min_support_method=args.min_support_method,
                        min_re=args.min_re,
                        pos_set=pos_set,
                        neg_set=neg_set)
    
    df_train = df_vcf.copy()
    if args.bench_sites is not None:
        bench_set = read_known_sv(args.bench_sites)
        if args.bench_vcf is not None:
            bench_vcf = read_vcf(args.bench_vcf)
            df_bench = get_df_vcf(vcf=bench_vcf, sv_type=args.sv_type, feature_list=features, feature_only=True)
        else:
            logger.info('Benchmark sites from input vcf are removed when training model')
            df_train = df_vcf.loc[df_vcf.index.difference(bench_set.index)]
        
        df_bench = df_bench.merge(bench_set, left_index=True, right_index=True, how='inner')

    else:
        df_bench = None
    
    # train and evaluate model
    logger.info(f'Train model for SVTYPE: {args.sv_type}')
    model = model_worker(df_train=df_train,
                         features=features,
                         out_prefix=args.output,
                         args=args,
                         df_bench=df_bench)
    

    # save model to file
    logger.info(f'Save model to {args.output}.model')
    pickle.dump(model, open(f'{args.output}.model', "wb"))

    # apply model to input vcf
    logger.info(f'Apply model on input vcf, write to {args.output}.vcf.gz')
    apply_vcf(invcf=invcf,
              output=args.output + '.vcf.gz',
              sv_type=args.sv_type,
              model=model,
              df_vcf=df_vcf,
              features=features)

    # apply model to bench vcf if provided
    if args.bench_vcf is not None:
        logger.info(f'Apply model on benchmark vcf, write to {args.output}.bench.vcf.gz')
        apply_vcf(invcf=bench_vcf,
                  output=args.output + '.bench.vcf.gz',
                  sv_type=args.sv_type,
                  model=model,
                  df_vcf=df_bench,
                  features=features)
    
    invcf.close()

    return None

