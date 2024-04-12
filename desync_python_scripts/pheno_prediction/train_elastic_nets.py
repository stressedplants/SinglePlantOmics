import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV

import matplotlib.pyplot as plt
import pickle

from multiprocessing import Pool
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
        random_state = 123,  # ensures a repeatable outcome
        verbose = 1,
        n_jobs = 1 # do NOT parallelise the 5-fold CV
    )

    # fit with CV
    return regr.fit(x_values, y_values)


# regularised logistic regressions with enet
def logistic_model_fn(
    x_values: np.array,
    y_values: np.array,  # dtype must be bool
    l1_list =  [.1, .5, .7, .9, .95, .99, 1]
) -> LogisticRegressionCV:

    regr = LogisticRegressionCV(
        Cs = 10,
        l1_ratios = l1_list,
        penalty = 'elasticnet',
        cv = 5,  # automatically chooses set to 5-fold validation
        random_state = 123,
        solver = 'saga',
        # scoring = 'neg_log_loss',  # use accuracy instead
        max_iter = 10**3,
        verbose = 0, # too many output lines to be useful if verbose = 1!
        n_jobs = 1 # do NOT parallelise the 5-fold CV
    )

    return regr.fit(x_values, y_values)


if __name__ == "__main__":

    # load the datasets
    expr_df = pd.read_csv("data/TPM_z_scored.csv", index_col = 0)
    regs_df = pd.read_csv("data/TPM_z_scored_only_regs.csv", index_col = 0)
    pheno_df = pd.read_csv("data/phenos_to_predict.csv", index_col = 0)

    # Use argparse to decide whether to run all possible l1_ratio or just fixed ...

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", action="store_true")
    args = parser.parse_args()

    if (args.f):
        # run linear regression with only fixed l1_ratio
        print("Running elastic nets with fixed l1_ratio ...")

        print("Bolting (logistic regression) with all genes. l1_ratio = 0.5")
        bolting_enet_output = run_model_with_parallel_loocv(
            expr_df,
            (pheno_df["bolting"].values == "Y"),
            functools.partial(
                logistic_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(bolting_enet_output, open("outputs/bolting_enet_05.sav", "wb"))

        print("Bolting (logistic regression) with regulators. l1_ratio = 0.1")
        bolting_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            (pheno_df["bolting"].values == "Y"),
            functools.partial(
                logistic_model_fn,
                l1_list = [.1]
            )
        )
        pickle.dump(bolting_regs_enet_output, open("outputs/bolting_regs_enet_01.sav", "wb"))


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

        print("Leaf size with regulators, l1_ratio = 0.5")
        leaf_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            pheno_df["leaf_avg"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(leaf_regs_enet_output, open("outputs/leaf_regs_enet_05.sav", "wb"))

        print("Biomass with all genes, l1_ratio = 0.5")
        biomass_enet_output = run_model_with_parallel_loocv(
            expr_df,
            pheno_df["biomass"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )
        pickle.dump(biomass_enet_output, open("outputs/biomass_enet_05.sav", "wb"))

        print("Biomass with regulators, l1_ratio = 0.5")
        biomass_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            pheno_df["biomass"],
            functools.partial(
                elastic_net_model_fn,
                l1_list = [.5]
            )
        )

        pickle.dump(biomass_regs_enet_output, open("outputs/biomass_regs_enet_05.sav", "wb"))

    else: 
        # run linear regression in all possible combos
        print("Running elastic nets with varying l1_ratio ...")

        print("Bolting (logistic regression) with all genes")
        bolting_enet_output = run_model_with_parallel_loocv(
            expr_df,
            (pheno_df["bolting"].values == "Y"),
            logistic_model_fn
        )
        pickle.dump(bolting_enet_output, open("outputs/bolting_enet.sav", "wb"))

        print("Bolting (logistic regression) with regulators")
        bolting_regs_enet_output = run_model_with_parallel_loocv(
            regs_df,
            (pheno_df["bolting"].values == "Y"),
            logistic_model_fn
        )
        pickle.dump(bolting_regs_enet_output, open("outputs/bolting_regs_enet.sav", "wb"))

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

