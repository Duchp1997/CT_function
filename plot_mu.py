"""
-*- coding: utf-8 -*-

@author: Du Changping
@time: 2021/10/16 10:55
@file name: plot_mu
@software：PyCharm

Do not smash your computer!

"""
import os
import pandas as pd
import numpy as np
from scipy import interpolate, integrate
import matplotlib.pyplot as plt

# data file path set
cur_dir = os.getcwd()
csv_data_path = os.path.join(cur_dir, "csv_data")

element_look_up_path = os.path.join(cur_dir, "./csv_data/element.csv")
elem_energy_dir_path = os.path.join(cur_dir, "./csv_data")

comp_data_path = os.path.join(csv_data_path, "compound.csv")
compound_data = pd.read_csv(comp_data_path)
compound_data.round(9)
compound_data.drop('Unnamed: 0', axis=1, inplace=True)


def get_material_density(mat_name_, compound_data_):
    """
    look up the standard density according to data in website:
    https://physics.nist.gov/PhysRefData/XrayMassCoef/tab2.html

    Parameters
    ----------
    mat_name_: tatget materail name, please notice there is an additional space at
                the end of the string
    compound_data_: look-table, here use user-made csv file

    Returns
    -------
    density data
    """
    assert type(mat_name_) is str, "please input the right parameters"
    try:
        result_ = compound_data_[compound_data_["Material"] == mat_name_]["Density"].dropna().to_numpy()
        result_ = np.max(result_)
        print("get the density data of input material:{}".format(mat_name_))
        return result_
    except LookupError:
        print("please check the input name\n."
              "please notice there is a space at"
              "at the end of the words")
    finally:
        print("the input name is :{}".format(mat_name_))


def get_target_tested_data(target_mat_):
    """
    get tested data based on nist website
    https://physics.nist.gov/PhysRefData/XrayMassCoef/tab4.html

    Parameters
    ----------
    target_mat_: name string

    Returns
    -------
    pd.DataFrame
    """
    tested_data_path_ = os.path.join(cur_dir, "tmp_data")
    tested_data_path_ = tested_data_path_ + "/" + target_mat_ + ".txt"
    tested_data_ = pd.read_table(tested_data_path_, sep=" ", header=None, skiprows=[0, 1, 2, 3])
    tested_data_.dropna(axis=1, inplace=True)
    tested_data_.columns = ["energy", "mu_rou", "muen_rou"]
    return tested_data_


def make_the_mu_plot(mat_name_, compound_data_):
    density_ = get_material_density(mat_name_, compound_data_)
    target_material_tested_data_ = get_target_tested_data(mat_name_)
    energy_ = target_material_tested_data_["energy"].to_numpy() * 1000
    mu_ = target_material_tested_data_["mu_rou"].to_numpy() * density_
    fig_, ax_ = plt.subplots(1, 1)
    ax_.set_title(mat_name_)
    # ax_.set_yscale("log")
    ax_.set_ylabel(r"$\mu$: $cm^{-1}$")
    ax_.set_ylim(0, 20)
    ax_.set_xlim(1, 200)
    ax_.set_xlabel(r"energy: $KeV$")
    ax_.plot(energy_, mu_)
    fig_.tight_layout()
    plt.show()


def get_random_energy(mean_, std_, x_min_, x_max_, x_num_=100):
    assert (mean_ < x_max_) & (mean_ > x_min_), "please notice the relative number"
    x_np = np.linspace(x_min_, x_max_, num=x_num_)
    result_energy_ = 1 / (std_ * pow(2 * np.pi, 0.5)) * np.exp(-((x_np - mean_) ** 2) / (2 * std_ ** 2))
    # plt.plot(x_np, result_energy_)
    # plt.show()
    return result_energy_


if __name__ == "__main__":
    target_mat = "Adipose Tissue (ICRU-44) "
    # result = get_material_density(target_mat, compound_data)
    make_the_mu_plot(target_mat, compound_data)
    # print(result)
