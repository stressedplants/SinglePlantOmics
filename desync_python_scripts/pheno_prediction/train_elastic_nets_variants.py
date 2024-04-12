import pandas as pd
import pickle
 
# import the code from train_elastic_nets.py
from train_elastic_nets import run_model_with_parallel_loocv, elastic_net_model_fn, logistic_model_fn

if __name__ == "__main__":

    # load the datasets
    variant_df = pd.read_csv("data/removed_all_heterozygous.csv", index_col = 0)
    pheno_df = pd.read_csv("data/phenos_to_predict.csv", index_col = 0)

    # run linear regression in all possible combos
    print("Running elastic nets with varying l1_ratio ...")

    print("Bolting (logistic regression) with variants")
    bolting_enet_output = run_model_with_parallel_loocv(
        variant_df,
        (pheno_df["bolting"].values == "Y"),
        logistic_model_fn
    )
    pickle.dump(bolting_enet_output, open("outputs/bolting_variant_enet.sav", "wb"))


    print("Leaf size with variants")
    leaf_enet_output = run_model_with_parallel_loocv(
        variant_df,
        pheno_df["leaf_avg"],
        elastic_net_model_fn
    )
    pickle.dump(leaf_enet_output, open("outputs/leaf_variant_enet.sav", "wb"))

    print("Biomass with variants")
    biomass_enet_output = run_model_with_parallel_loocv(
        variant_df,
        pheno_df["biomass"],
        elastic_net_model_fn
    )
    pickle.dump(biomass_enet_output, open("outputs/biomass_variant_enet.sav", "wb"))

