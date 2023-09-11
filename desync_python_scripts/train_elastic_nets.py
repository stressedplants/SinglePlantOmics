# This is changed to include a parallel loop and more verbosity of output!

import numpy as np
import pandas as pd
from sklearn.model_selection import GridSearchCV, cross_val_score, RepeatedKFold, LeaveOneOut, train_test_split
from sklearn.linear_model import ElasticNet, ElasticNetCV, LogisticRegressionCV
from sklearn.metrics import make_scorer, log_loss
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor

import matplotlib.pyplot as plt
import pickle

from typing import Callable

from multiprocessing import Pool
import os
import functools

import argparse

def custom_loocv_single(
    sample_pos: int,
    sample_names: list[str],
    input_pheno: np.array,
    input_expr: pd.DataFrame,
    model_fn,
    verbose = True
):
    """Truncates the dataset based on inputs then returns the output of a single
    fit from model_fn.
    """
    
    sample = sample_names[sample_pos]
    if verbose: print(f'Performing fitting without sample {sample}')

    # Remove the chosen sample from the list, then create the dataset
    remaining_positions = list(range(len(sample_names)))
    remaining_positions.remove(sample_pos)
        
    # make truncated datasets
    trunc_pheno = input_pheno[remaining_positions]
    trunc_expr = input_expr.iloc[:, remaining_positions]

    y_values = trunc_pheno
    x_values = np.transpose(trunc_expr.values)

    # get the output of the fitting (eg ---.fit(x_values, y_values))
    return(model_fn(x_values, y_values))

def run_model_with_parallel_loocv(
    input_expr: pd.DataFrame,
    input_pheno: np.array,
    model_fn,
    verbose = True
) -> dict:
    """Returns a dictionary with the model results, where the keys are IDs of 
    each sample and values are the outputs of model_fn.

    Args:
        input_expr (pd.DataFrame): represents the gene expression values.
        input_pheno (np.array): represents the phenotypes to be predicted - eg floats
            for regression or bools for classification.
        model_fn: A function which will take an array of expression values and an
            array of phenotypes, and return an arbitrary modelling output.
        verbose (bool, optional): whether to print steps. Defaults to True.
    """

    # validate whether there are the same number of samples as phenotypes
    sample_names = input_expr.columns.values.tolist()
    if not (len(input_pheno) == len(sample_names)):
        raise(ValueError("Length of samples and phenotypes must be equal"))

    results_dict = dict.fromkeys(sample_names)
    range_of_samples = range(len(sample_names))

    # # make a helper function that fixes the data inputs, allowing only sample_id to change
    # def specific_looc(sample_id): 
    #     custom_loocv_single(
    #         sample_id,
    #         sample_names,
    #         input_pheno,
    #         input_expr,
    #         model_fn,
    #         verbose = verbose
    #     )

    # given all the sample indices, run the loocv parallelised
    pool = Pool()
    result_async = pool.map_async(
        functools.partial(  # use this to keep data inputs fixed
            custom_loocv_single,
            sample_names = sample_names,
            input_pheno = input_pheno,
            input_expr = input_expr,
            model_fn = model_fn,
            verbose = verbose
        ), 
        range_of_samples
    )  # this will take time
    result_async.wait() # until completion
    
    list_of_results = result_async.get()
    
    # Cast this into a results dictionary
    for new_key, new_val in zip(results_dict.keys(), list_of_results):
        results_dict[new_key] = new_val
    
    return results_dict

# regularised linear regressions with enet
def elastic_net_model_fn(
    x_values: np.array,
    y_values: np.array,
    l1_list =  [.1, .5, .7, .9, .95, .99, 1]
) -> ElasticNetCV:
    
    regr = ElasticNetCV(
        l1_ratio = l1_list,
        cv = 5,  # automatically chooses set to 5-fold validation
        random_state = 0,  # ensures a repeatable outcome
        verbose = 1,
        n_jobs = 1 # do NOT parallelise the 5-fold CV
    )

    # fit with CV
    return regr.fit(x_values, y_values)


# # regularised logistic regressions with enet
# def logistic_model_fn(
#     x_values: np.array,
#     y_values: np.array  # dtype must be bool
# ) -> LogisticRegressionCV:

#     l1_list =  [.1, .5, .7, .9, .95, .99, 1]

#     regr = LogisticRegressionCV(
#         Cs = 10,
#         l1_ratios = l1_list,
#         penalty = 'elasticnet',
#         cv = 5,  # automatically chooses set to 5-fold validation
#         random_state = 0,
#         solver = 'saga',
#         scoring = 'neg_log_loss',
#         max_iter = 10**3,
#         verbose = 1,
#         n_jobs = 1 # do NOT parallelise the 5-fold CV
#     )

#     return regr.fit(x_values, y_values)


if __name__ == "__main__":

    # load the datasets
    expr_df = pd.read_csv("data/TPM_z_scored.csv", index_col = 0)
    regs_df = pd.read_csv("data/TPM_z_scored_only_regs.csv", index_col = 0)
    pheno_df = pd.read_csv("data/phenos_to_predict.csv", index_col = 0)

    # # DURING TESTING - MAKE A BABY DATASET!
    # nrow, ncol = expr_df.shape
    # choose_samples = np.random.choice(np.arange(0, ncol), size = 20, replace = False)
    # choose_genes = np.random.choice(np.arange(0, nrow), size = 500, replace = False)

    # baby_expr = expr_df.iloc[choose_genes, choose_samples]
    # baby_pheno = pheno_df.iloc[choose_samples, :]

    # # test linear regression
    # baby_leaf_output = run_model_with_parallel_loocv(
    #     baby_expr,
    #     baby_pheno["leaf_avg"].values,
    #     elastic_net_model_fn
    # )

    # # test logistic regression
    # baby_bolting_output = run_model_with_parallel_loocv(
    #     baby_expr,
    #     (baby_pheno["bolting"].values == "Y"),
    #     logistic_model_fn
    # )


    # Use argparse to decide whether to run all possible l1_ratio or just fixed ...

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", action="store_true")
    args = parser.parse_args()

    if (args.f):
        # run linear regression with only fixed l1_ratio
        print("Running elastic nets with fixed l1_ratio ...")

        print("Leaf size with all genes, l1_ratio = 0.5")
        leaf_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["leaf_avg"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(leaf_enet_output, open("outputs/leaf_enet_05.sav", "wb"))

        print("Leaf size with regulators, l1_ratio = 0.1")
        leaf_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            pheno_df["leaf_avg"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.1]
            )
        )
        pickle.dump(leaf_regs_enet_output, open("outputs/leaf_regs_enet_01.sav", "wb"))

        print("Biomass with all genes, l1_ratio = 0.1")
        biomass_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["biomass"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.1]
            )
        )
        pickle.dump(biomass_enet_output, open("outputs/biomass_enet_01.sav", "wb"))

        print("Biomass with regulators, l1_ratio = 0.1")
        biomass_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            pheno_df["biomass"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.1]
            )
        )

        pickle.dump(biomass_regs_enet_output, open("outputs/biomass_regs_enet_01.sav", "wb"))
    
    else: 
        # run linear regression in all possible combos
        print("Running elastic nets with varying l1_ratio ...")

        print("Leaf size with all genes")
        leaf_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["leaf_avg"],
            elastic_net_model_fn
        )
        pickle.dump(leaf_enet_output, open("outputs/leaf_enet.sav", "wb"))

        print("Leaf size with regulators")
        leaf_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            pheno_df["leaf_avg"],
            elastic_net_model_fn
        )
        pickle.dump(leaf_regs_enet_output, open("outputs/leaf_regs_enet.sav", "wb"))

        print("Biomass with all genes")
        biomass_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["biomass"],
            elastic_net_model_fn
        )
        pickle.dump(biomass_enet_output, open("outputs/biomass_enet.sav", "wb"))

        print("Biomass with regulators")
        biomass_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            pheno_df["biomass"],
            elastic_net_model_fn
        )

        pickle.dump(biomass_regs_enet_output, open("outputs/biomass_regs_enet.sav", "wb"))

