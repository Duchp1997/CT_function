"""
-*- coding: utf-8 -*-

@author: Du Changping
@time: 2021/10/17 11:18
@file name: simulate_mu_and_depth_final
@software：PyCharm

Do not smash your computer!
this file is used to simulate the relationship between depth of 
compound material and mu
Done
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import integrate, interpolate
from plot_material import get_target_materail_compound, get_compound_coeffi


def set_cur_dir_and_compound_file():
    """
    data file path set
    Returns
    -------

    """
    cur_dir_ = os.getcwd()
    csv_data_path = os.path.join(cur_dir_, "csv_data")

    element_look_up_path = os.path.join(cur_dir_, "./csv_data/element.csv")
    elem_energy_dir_path = os.path.join(cur_dir_, "./csv_data")

    comp_data_path = os.path.join(csv_data_path, "compound.csv")
    compound_data_ = pd.read_csv(comp_data_path)
    compound_data_.round(9)
    compound_data_.drop('Unnamed: 0', axis=1, inplace=True)
    return cur_dir_, compound_data_


def prob_density_func(x_, std_, mean_, unit_="Mev"):
    """
    here we omit the physical unit, but keep in mind
    when encountering the look-up table, change the
    physical unit between kev and mev
    Parameters
    ----------
    x_
    std_
    mean_
    unit_

    Returns
    -------

    """
    result_ = 1 / (std_ * pow(2 * np.pi, 0.5)) * np.exp(-((x_ - mean_) ** 2) / (2 * std_ ** 2))

    return result_


def get_random_energy(mean_, std_, x_min_, x_max_, x_num_=101, unit_="kev"):
    """
    we use this function to get some energy spectrum of some distribution
    in future work, we can just customize any distribution
    Parameters
    ----------
    mean_
    std_
    x_min_
    x_max_
    x_num_
    unit_

    Returns
    -------

    """
    print("the physical unit is {}, keep this in mind!!!".format(unit_))
    assert (mean_ < x_max_) & (mean_ > x_min_), "please notice the relative number"
    x_np = np.linspace(x_min_, x_max_, num=x_num_)
    result_energy_ = 1 / (std_ * pow(2 * np.pi, 0.5)) * np.exp(-((x_np - mean_) ** 2) / (2 * std_ ** 2))
    # plt.plot(x_np, result_energy_)
    # plt.show()
    return result_energy_


def customized_intepolate(target_material_compound_,
                          to_intepolate_energy_,
                          compound_data_,
                          target_material_,
                          unit_="keV"):
    """
    here we use cubic interpolation to fill those wanted data
    Parameters
    ----------
    target_material_compound_
    to_intepolate_energy_
    compound_data_
    target_material_
    unit_

    Returns
    -------

    """
    print("physical unit used is {}, but the unit in target material compound "
          "is Mev. \n Keep this in mind!!!".format(unit_))
    keys_ = target_material_compound_.keys()
    assert len(keys_) > 2, "there must be sth wrong with data"
    # here we times 0.001 to convert kev to mev
    cus_to_interpolate_ = pd.DataFrame({"energy": to_intepolate_energy_ * 0.001})
    cus_to_interpolate_np_ = cus_to_interpolate_["energy"].apply(lambda x: 1 / x).to_numpy()
    tested_energy_ = target_material_compound_[keys_[0]].apply(lambda x: 1 / x)
    for key_ in keys_[1:-1]:
        tested_key_np_ = target_material_compound_[key_]
        f_interp_ = interpolate.interp1d(tested_energy_,
                                         tested_key_np_,
                                         kind="cubic",
                                         fill_value="extrapolate")
        key_interpolated_np_ = f_interp_(cus_to_interpolate_np_)
        cus_to_interpolate_[key_] = key_interpolated_np_

    cus_to_interpolate_['result'] = 0
    _, elem_weights_ = get_compound_coeffi(target_material_, compound_data_)
    for idx_ in range(len(elem_weights_)):
        cus_to_interpolate_['result'] = cus_to_interpolate_['result'] + \
                                        cus_to_interpolate_[keys_[idx_ + 1]] * elem_weights_[idx_]
    # cos_to_interpolate_result_np = cus_to_interpolate['result'].to_numpy()
    return cus_to_interpolate_


def get_mu(energy_, cus_to_interpolate_, unit_="kev"):
    """
    check the  $\mu$  based on customized interpolation
    Parameters
    ----------
    energy_
    cus_to_interpolate_
    unit_

    Returns
    -------

    """
    if unit_ == "kev" or unit_ == "Kev":
        result_ = cus_to_interpolate_[cus_to_interpolate_["energy"] == energy_ * 0.001]["result"].to_numpy()
    elif unit_ == "mev" or unit_ == "Mev":
        result_ = cus_to_interpolate_[cus_to_interpolate_["energy"] == energy_ * 0.001]["result"].to_numpy()
    else:
        raise TypeError("please notice the physical unit: kev or mev")
    return result_


def calc_compound_mu_with_depth(depth_,
                                energy_distribution_,
                                to_intepolate_energy_,
                                cus_interpolation_, unit_="kev"):
    """
    compute the prob function of the scattered photon based on Beer-Lambert law
    in coding, we just the previous energy distribution:
              $\int \Omega(E) * \exp(-depth * \mu(E) $
    Parameters
    ----------
    depth_
    energy_distribution_
    to_intepolate_energy_
    cus_interpolation_
    unit_

    Returns
    -------

    """
    exp_func_part_ = np.exp(-1 * depth_ * get_mu(to_intepolate_energy_, cus_interpolation_, unit_=unit_))
    assert exp_func_part_.shape == energy_distribution_.shape, "something wrong"
    integ_func = energy_distribution_ * exp_func_part_
    results = integrate.trapz(integ_func, to_intepolate_energy_)
    return results


def calculation(target_energy_prob_distribution_,
                energy_interval_,
                thickness_):
    initial_energy_ = integrate.trapz(target_energy_prob_distribution_,
                                      energy_interval_)
    print("initial photon number is {}".format(initial_energy_))

    # get the target materaial compound based on published file
    cur_dir_, compound_data_ = set_cur_dir_and_compound_file()
    # comp_coeff = compound_data_[compound_data_["Material"] == target_mat]['Composition']
    target_material_compound_ = get_target_materail_compound(target_mat_=target_mat,
                                                             compound_data_=compound_data_)

    cus_to_interpolate_ = customized_intepolate(target_material_compound_=target_material_compound_,
                                                to_intepolate_energy_=np.linspace(X_MIN, X_MAX, ENERGY_NUM),
                                                compound_data_=compound_data_,
                                                target_material_=target_mat)

    mu_au_ = [calc_compound_mu_with_depth(depth_=depth_,
                                          energy_distribution_=target_energy_prob_distribution,
                                          to_intepolate_energy_=energy_interval,
                                          cus_interpolation_=cus_to_interpolate_,
                                          ) for depth_ in L]
    mu_result_ = np.log(initial_energy_) - np.log(np.array(mu_au_))
    mu_result_ = mu_result_ / thickness_
    return mu_result_


def plotting(fig_1, ax_1, mu_result_, x_axis_, label_, ):
    # ax_.set_ylim(bottom=0.1)
    # ax_.plot(np.linspace(X_MIN, X_MAX, num=100), mu_result)
    ax_1.plot(x_axis_, mu_result_, label=label_)
    ax_1.set_xlabel(r'depth :$cm$')
    ax_1.set_ylabel(r'$\mu$: $cm^{-1}$')
    # fig_1.tight_layout()
    # plt.show()


if __name__ == "__main__":
    # tested material 注意这里面结尾是有一个空格的
    target_mat = "Adipose Tissue (ICRU-44) "
    # target_mat_ = "Water, Liquid "
    # target_mat_ = "A-150 Tissue-Equivalent Plastic"

    # energy interval
    X_MIN = 1
    X_MAX = 200

    # get target energy (unit: Kev)
    # here we set the incident photon number is unit 1
    # STD = 20
    # STD = 0.5
    MEAN = 100
    ENERGY_NUM = 100
    L = np.linspace(0.1, 30, 50)
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_title(target_mat)

    for std_ids in [5, 10, 20]:
        target_energy_prob_distribution = get_random_energy(mean_=MEAN, std_=std_ids, x_min_=X_MIN, x_max_=X_MAX,
                                                            x_num_=ENERGY_NUM, unit_="kev")
        energy_interval = np.linspace(X_MIN, X_MAX, ENERGY_NUM)
        mu_result = calculation(target_energy_prob_distribution_=target_energy_prob_distribution,
                                energy_interval_=energy_interval,
                                thickness_=L)
        ax.plot(L, mu_result, label="std:{}".format(std_ids))

    fig.tight_layout()
    plt.legend(loc="best")
    plt.show()
    print(2)


