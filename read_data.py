"""
-*- coding: utf-8 -*-

@author: Du Changping
@time: 2021/10/8 11:13
@file name: read_data
@software：PyCharm

Do not smash your computer!

"""

import os
import numpy as np
import pandas as pd

cur_dir = os.getcwd()
source_dir = os.path.join(cur_dir, "elem_data")


# element
element_data_path = os.path.join(source_dir, "element.txt")
element_data = pd.read_table(element_data_path, sep="\t")
element_data.round(9)
element_data.rename(columns={'Unnamed: 1': 'syb'}, inplace=True)
element_data.to_csv("./csv_data/element.csv")


# compound
compound_data_path = os.path.join(cur_dir, "elem_data", "compound.txt")
compound_data = pd.read_table(compound_data_path, 
                              sep="\t", 
                              header=None,
                              skiprows=[0, 1, 2])
compound_data.round(9)
data_np = compound_data.to_numpy()

tmp_index = -1
for i in range(data_np.shape[0]):
    if data_np[i][4] is not np.nan:
        tmp_index = i
        continue
    else:
        data_np[i][4] = data_np[i][0]
        data_np[i][1] = data_np[i][2] = data_np[i][3] = None
        data_np[i][0] = data_np[tmp_index][0]

data_pd = pd.DataFrame(data_np)
data_pd.round(9)
data_pd.columns = ["Material", "ZA", "I", "Density", "Composition"]
# data_pd = data_pd.reindex(columns=["Material", "ZA", "I", "Density", "Composition"])
data_pd.head()
data_pd.to_csv("./csv_data/compound.csv")   


# element mass atteunation coeffecient
def get_element_energy(path, save_path, save=True, is_print=False):
    data = pd.read_table(path, sep=" ", header=None, skiprows=[0, 1, 2, 3])
    data.round(9)
    data.dropna(axis=1, inplace=True)
    data.columns = ["Energy", "mu/rou", "mu_en/rou"]
    if is_print:
        print(data.head(10))
    if save:
        data.to_csv(save_path)

# H
H_data_path = os.path.join(cur_dir, "elem_data", "H.txt")
save_path = "./csv_data/H.csv"
get_element_energy(H_data_path, save_path, save=True, is_print=False)

# C
C_data_path = os.path.join(cur_dir, "elem_data", "C.txt")
save_path = "./csv_data/C.csv"
get_element_energy(C_data_path, save_path, save=True, is_print=False)

# N
N_data_path = os.path.join(cur_dir, "elem_data", "N.txt")
save_path = "./csv_data/N.csv"
get_element_energy(N_data_path, save_path, save=True, is_print=False)

# O
O_data_path = os.path.join(cur_dir, "elem_data", "O.txt")
save_path = "./csv_data/O.csv"
get_element_energy(O_data_path, save_path, save=True, is_print=False)

# Na
Na_data_path = os.path.join(cur_dir, "elem_data", "Na.txt")
save_path = "./csv_data/Na.csv"
get_element_energy(Na_data_path, save_path, save=True, is_print=False)


# S
S_data_path = os.path.join(cur_dir, "elem_data", "S.txt")
save_path = "./csv_data/S.csv"
get_element_energy(S_data_path, save_path, save=True, is_print=False)

# Ca
Ca_data_path = os.path.join(cur_dir, "elem_data", "Ca.txt")
save_path = "./csv_data/Ca.csv"
get_element_energy(Ca_data_path, save_path, save=True, is_print=False)

# Cl
Cl_data_path = os.path.join(cur_dir, "elem_data", "Cl.txt")
save_path = "./csv_data/Cl.csv"
get_element_energy(Cl_data_path, save_path, save=True, is_print=False)
