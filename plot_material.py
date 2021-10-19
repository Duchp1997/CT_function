"""
-*- coding: utf-8 -*-

@author: Du Changping
@time: 2021/10/10 10:09
@file name: plot_material
@software：PyCharm

Do not smash your computer!

******Done******
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate

# data file path set
cur_dir = os.getcwd()
csv_data_path = os.path.join(cur_dir, "csv_data")

element_look_up_path = os.path.join(cur_dir, "./csv_data/element.csv")
elem_energy_dir_path = os.path.join(cur_dir, "./csv_data")

comp_data_path = os.path.join(csv_data_path, "compound.csv")
compound_data = pd.read_csv(comp_data_path)
compound_data.round(9)
compound_data.drop('Unnamed: 0', axis=1, inplace=True)

# 注意这里面结尾是有一个空格的
target_mat = "Adipose Tissue (ICRU-44) "
# target_mat_ = "Water, Liquid "
# target_mat_ = "A-150 Tissue-Equivalent Plastic"
comp_coeff = compound_data[compound_data["Material"] == target_mat]['Composition']


def get_element_mac(ele_no, element_look_up_path_, elem_energy_dir_path_):
    elem_lookup = pd.read_csv(element_look_up_path_)
    elem_lookup.round(9)
    elem_syb = elem_lookup[elem_lookup["Z"] == ele_no]['syb'].to_numpy()[0]
    print(elem_syb)
    assert type(elem_syb) is str, "here we assume that the type of symbol"
    elem_energy_file_path = elem_syb + '.csv'
    elem_energy_file_path = os.path.join(elem_energy_dir_path_, elem_energy_file_path)
    if not os.path.exists(elem_energy_file_path):
        if not os.path.exists(elem_energy_dir_path_):
            print("{} does not exist, please check the dir path".format(elem_energy_dir_path_))
        else:
            print("{} does not exist".format(elem_energy_file_path))
            exit()
    elem_energy_data = pd.read_csv(elem_energy_file_path)
    elem_energy_data.round(9)
    energy = elem_energy_data["Energy"]
    mu_rou = elem_energy_data["mu/rou"]
    elem_pd = pd.DataFrame({"energy": energy, "mu_rou" + str(ele_no): mu_rou})
    elem_pd.round(9)
    return elem_pd


def get_compound_coeffi(target_mat_, compound_data_):
    comp_coeff_ = compound_data_[compound_data_["Material"] == target_mat_]['Composition'].to_numpy()
    elem_weight_ = []
    total_mac_ = pd.DataFrame({"energy": [0.001]})
    total_mac_.round(9)
    for elem_ in comp_coeff_:
        elem_mac = get_element_mac(int(elem_.split(":")[0]),
                                   element_look_up_path,
                                   elem_energy_dir_path)
        total_mac_ = total_mac_.merge(elem_mac, how='outer', on="energy")

        elem_weight_.append(float(elem_.split(":")[1]))
    total_mac_.sort_values(by="energy", inplace=True)

    return total_mac_, elem_weight_


def costum_interpolate(input_ndarry, if_inverse=True):
    keys = input_ndarry.keys()
    if if_inverse:
        energy_to_interp = input_ndarry[keys[0]].apply(lambda x: 1/x).to_numpy()
    else:
        energy_to_interp = input_ndarry[keys[0]].to_numpy()
    for key in keys[1:]:
        two_array = input_ndarry[[keys[0], key]].dropna(axis=0, how='any')
        if if_inverse:
            energy_data_without_nan = two_array[keys[0]].apply(lambda x: 1 / x).to_numpy()
        else:
            energy_data_without_nan = two_array[keys[0]].to_numpy()
        to_interp_data_without_nan = two_array[key].to_numpy()
        f_interp = interpolate.interp1d(energy_data_without_nan,
                                        to_interp_data_without_nan,
                                        kind="cubic",
                                        fill_value="extrapolate")
        tmp_pred = f_interp(energy_to_interp)
        input_ndarry[key] = tmp_pred
    return input_ndarry


def get_target_materail_compound(target_mat_, compound_data_):
    total_mac_, elem_weight_ = get_compound_coeffi(target_mat_, compound_data_)
    tmp_ = total_mac_.reset_index()
    tmp_.drop("index", axis=1, inplace=True)
    tmp_ = costum_interpolate(tmp_)
    # tmp_ = tmp_.interpolate(method="polynomial",
    #                         order=3,
    #                         limit_direction="both",
    #                         index=tmp_.index, values=tmp_['energy'])
    tmp_['result'] = 0
    tmp_keys_ = tmp_.keys()
    for idx_ in range(len(elem_weight_)):
        tmp_['result'] = tmp_['result'] + tmp_[tmp_keys_[idx_ + 1]] * elem_weight_[idx_]
    return tmp_


def get_tested_data():
    # get tested data based on nist website
    tested_data_path_ = os.path.join(cur_dir, "tmp_data")
    tested_data_path_ = tested_data_path_ + "/" + target_mat + ".txt"
    tested_data_ = pd.read_table(tested_data_path_, sep=" ", header=None, skiprows=[0, 1, 2, 3])
    tested_data_.dropna(axis=1, inplace=True)
    tested_data_.columns = ["energy", "mu_rou", "muen_rou"]
    return tested_data_


def plotting(tested_data_, tmp_data_):
    fig_, ax_ = plt.subplots(1, 1)
    ax_.set_xscale("log")
    ax_.set_yscale("log")
    ax_.set(title=target_mat)
    ax_.set_ylim(bottom=0.1)
    ax_.set_xlim(1, 200)
    ax_.set_xlabel("energy: Kev")
    ax_.set_ylabel(r"$\mu / \rho$: $cm^{2}/g$")

    # ax.set_xlim(bottom=0.001)
    ax_.plot(tmp_data_['energy'].to_numpy()*1000, tmp_data_['result'].to_numpy(), "*b")
    ax_.plot(tested_data_['energy'].to_numpy()*1000, tested_data_['mu_rou'].to_numpy(), 'r')
    fig_.tight_layout()
    plt.show()


if __name__ == "__main__":
    # 注意这里面结尾是有一个空格的 give the name of the target material
    target_mat = "Adipose Tissue (ICRU-44) "
    # target_mat_ = "Water, Liquid "
    # target_mat_ = "A-150 Tissue-Equivalent Plastic"

    # get the target materaial compound based on published file
    comp_coeff = compound_data[compound_data["Material"] == target_mat]['Composition']
    # sample
    # 1: 0.114000
    # 6: 0.598000
    # 7: 0.007000
    # 8: 0.278000
    # 11: 0.001000
    # 16: 0.001000
    # 17: 0.001000

    tested_data = get_tested_data()
    target_data = get_target_materail_compound(target_mat_=target_mat,
                                               compound_data_=compound_data)
    plotting(tested_data_=tested_data,
             tmp_data_=target_data)


