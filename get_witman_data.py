""" Get crystal features for structures in Matt Witman's Nature Computational Science Paper """
from glob import glob

import numpy as np
import pandas as pd
from pymatgen.core import Structure, Composition
from tqdm import tqdm

from deftpy.crystal_analysis import Crystal


def main():
    csv_paths = sorted(glob("playground/witman_data/data_01_03_22/csvs/*.csv"))
    poscar_path = "playground/witman_data/data_01_03_22/poscars"
    oxstate_paths = sorted(glob("playground/witman_data/data_01_03_22/oxstate/*_oxstate"))

    # Read in all the data and concatenate it into one dataframe
    # Create a column for the csv file name, not including the path and .csv extension
    df = pd.concat(
        [pd.read_csv(csv_path).assign(filename=csv_path.split("/")[-1][:-4]) for csv_path in csv_paths]
    )

    # Create a column for the defectid, which is the string after the last period in the filename
    df["defectid"] = df["filename"].str.split(".").str[-1]

    # Drop rows with defectname != "V_O"
    df = df[df["defectname"] == "V_O"]

    # Reset the index
    df = df.reset_index(drop=True)

    # Create a column for the poscar
    # The poscars are in poscar_path, and the filename is filename + "_POSCAR_wyck"
    df["poscar"] = df["filename"].apply(lambda x: open(f"{poscar_path}/{x}_POSCAR_wyck").read())

    # Create a column for the pymatgen structure
    df["structure"] = df["poscar"].apply(lambda x: Structure.from_str(x, fmt="poscar"))

    # Add oxidation states to the structure
    for _, row in df.iterrows():
        oxstate_path = f"playground/witman_data/data_01_03_22/oxstate/{row['filename']}_oxstate"
        oxstate = []
        for x in open(oxstate_path).read().split():
            oxstate += int(x.split("*")[0]) * [float(x.split("*")[1])]
        structure = row["structure"]
        structure.add_oxidation_state_by_site(oxstate)

    # Binary?
    df["is_binary"] = df["formula"].apply(lambda x: len(Composition(x).elements) == 2)

    # Sort by defectid and then by site
    df = df.sort_values(["defectid", "site"])
    # df.to_csv("oxygen_vacancies.csv", index=False)

    # Calculate crystal features for binary structures
    df = df[df["is_binary"]]
    df_cf = pd.DataFrame()
    for defectid in tqdm(df["defectid"].unique()):
        df_defectid = df[df["defectid"] == defectid]
        structure = df_defectid["structure"].iloc[0]
        crystal = Crystal(pymatgen_structure=structure)

        CN = crystal.cn_dicts
        Eb = crystal.bond_dissociation_enthalpies
        Vr = crystal.reduction_potentials

        # Calculate CN-weighted Eb sum
        Eb_sum = []
        for CN_dict, Eb_dict in zip(CN, Eb):
            CN_array = np.array(list(CN_dict.values()))
            Eb_array = np.array(list(Eb_dict.values()))
            Eb_sum.append(np.sum(CN_array * Eb_array))

        # Calculate maximum Vr
        Vr_max = []
        for Vr_dict in Vr:
            try:
                Vr_max.append(max(Vr_dict.values()))
            except ValueError:
                Vr_max.append(np.nan)

        # Make a dataframe
        formula = df_defectid["formula"].values
        defectid = df_defectid["defectid"].values
        site = df_defectid["site"].values
        Eg = df_defectid["bandgap_eV"].values
        Ev = df_defectid["dH_eV"].values
        try:
            df_cf = pd.concat(
                [
                    df_cf,
                    pd.DataFrame(
                        {
                            "formula": formula,
                            "defectid": defectid,
                            "site": site,
                            "Eb_sum": Eb_sum,
                            "Vr_max": Vr_max,
                            "Eg": Eg,
                            "Ev": Ev,
                        }
                    ),
                ]
            )
        except ValueError:
            pass
    df_cf = df_cf.reset_index(drop=True)
    df_cf.to_csv("witman_data.csv", index=False)


if __name__ == "__main__":
    main()
