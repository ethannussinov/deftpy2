""" Get crystal features for structures in Yu Kumagai's Physical Review Materials Paper """
import sys
import tarfile
from glob import glob

import matplotlib.pyplot as plt
import pandas as pd
from pymatgen.core import Composition


def main():
    # Li2O
    with tarfile.open(glob("playground/kumagai_data/oxygen_vacancies_db_data/Li2O/*0.tar.gz")[0], "r:gz") as tar:
        # Read the member with the name "CONTCAR-finish"
        f = tar.extractfile("CONTCAR-finish")
        print(f.read().decode("utf-8"))

    df = pd.read_csv("playground/kumagai_data/vacancy_formation_energy_ml/charge0.csv")
    df["is_binary"] = df.formula.apply(lambda x: len(Composition(x)) == 2)
    df_binary = df.loc[df.is_binary == True].reset_index(drop=True)

    # plot
    df_plot = df_binary[["formula", "full_name", "band_gap", "formation_energy", "nn_ave_eleneg", "o2p_center_from_vbm",
                         "vacancy_formation_energy"]].reset_index(drop=True)
    n_atoms = []
    metals = []
    n_metals = []
    oxi_states = []
    for _, row in df_plot.iterrows():
        formula = row["formula"]
        composition = Composition(formula)
        metal = [el.symbol for el in composition.elements if el.symbol != "O"][0]
        n_metal = composition.get_el_amt_dict()[metal]
        n_oxygen = composition.get_el_amt_dict()["O"]
        oxi_state = 2 * n_oxygen / n_metal
        n_atoms.append(composition.num_atoms)
        metals.append(metal)
        n_metals.append(n_metal)
        oxi_states.append(oxi_state)
    df_plot["n_atoms"] = n_atoms
    df_plot["metal"] = metals
    df_plot["n_metal"] = n_metals
    df_plot["oxi_state"] = oxi_states
    df_plot["vr"] = df_plot["n_atoms"] * df_plot["formation_energy"] / df_plot["n_metal"] / df_plot["oxi_state"]
    plt.plot(df_plot["formation_energy"], df_plot["vacancy_formation_energy"], "o")
    plt.plot(df_plot["vr"], df_plot["vacancy_formation_energy"], "o")
    plt.show()


if __name__ == "__main__":
    main()
