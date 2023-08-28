""" Get crystal features for structures in Matt Witman's Nature Computational Science Paper """
from glob import glob

import pandas as pd
from pymatgen.core import Structure
import crystal


def main():
    csv_paths = sorted(glob("oxygen_vacancies/structures_production/data_01_03_22/csvs/*.csv"))
    poscar_path = "oxygen_vacancies/structures_production/data_01_03_22/poscars"

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

    for structure in df["structure"]:
        structure.add_oxidation_state_by_guess()
        c = crystal.Crystal(structure)
        print(crystal.crystal_data(c))
        break


if __name__ == "__main__":
    main()
