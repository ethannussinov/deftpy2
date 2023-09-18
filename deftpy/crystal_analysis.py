import re
from enum import Enum
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Callable

import pandas as pd
from ase.visualize import view
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Species, Structure
from pymatgen.io.ase import AseAtomsAdaptor

SPECIES_SYMBOL = "O"


class DataFiles(Enum):
    EB = ("Eb", Path("../data/Eb.csv"))
    VR = ("Vr", Path("../data/Vr.csv"))

    @property
    def column_name(self):
        return self.value[0]

    @property
    def file_path(self):
        return self.value[1]


class Crystal:
    """
    A class for analyzing crystal structures.

    Attributes:
        structure: The pymatgen Structure object.
        nn_finder: The CrystalNN object.
        eb: The bond dissociation enthalpies.
        vr: The reduction potentials.
        cn_dicts: A list of coordination number dictionaries.
        bond_dissociation_enthalpies: A list of bond dissociation enthalpies.
        reduction_potentials: A list of reduction potentials.


    Methods:
        visualize: Visualizes the crystal structure using ASE's view function.


    Examples:
        # TODO: Add examples
    """

    @staticmethod
    def _split_before_first_number(s: str) -> List[str]:
        """
        Splits a string before the first number.

        Args:
            s: The string to split.

        Returns:
            A list of strings.

        Examples:
            # TODO: Add examples
        """
        return re.split(r"(?=\d)", s, maxsplit=1)

    @staticmethod
    def _parse_species_string(species_string: str) -> Tuple[Optional[Species], str, int]:
        """
        Parses a species string.

        Args:
            species_string: The species string to parse.

        Returns:
            A tuple of the species, the symbol, and the oxidation state.

        Examples:
            # TODO: Add examples
        """
        # Check if the string is of valid format before trying to parse
        if not re.match(r"[A-Za-z]+\d+\+", species_string):
            split_str = Crystal._split_before_first_number(species_string)
            return None, split_str[0], round(float(split_str[1][:-1]))

        species = Species(species_string)
        return species, species.symbol, species.oxi_state

    def __init__(self, filepath: Optional[str] = None, poscar_string: Optional[str] = None,
                 pymatgen_structure: Optional[Structure] = None, nn_finder: Optional[CrystalNN] = None,
                 use_weights: Optional[bool] = False):
        """
        Initializes the Crystal object.

        Args:
            filepath: The filepath to the POSCAR file.
            poscar_string: The POSCAR string.
            pymatgen_structure: The pymatgen Structure object.
            nn_finder: The CrystalNN object.

        Raises:
            ValueError: If neither filepath, poscar_string, nor pymatgen_structure is specified.

        Examples:
            # TODO: Add examples
        """
        if filepath:
            self.structure = Structure.from_file(filepath)
        elif poscar_string:
            self.structure = Structure.from_str(poscar_string, fmt="poscar")
        elif pymatgen_structure:
            self.structure = pymatgen_structure
        else:
            raise ValueError("Specify either filepath, poscar_string, or pymatgen_structure.")

        self.nn_finder = nn_finder or CrystalNN()
        self.use_weights = use_weights

        package_dir = Path(__file__).parent
        self.eb = pd.read_csv(package_dir / DataFiles.EB.file_path)
        self.vr = pd.read_csv(package_dir / DataFiles.VR.file_path)

        self._cn_dicts_initialized = False
        self.cn_dicts = []
        self.bond_dissociation_enthalpies = self._get_values(self.eb, DataFiles.EB.column_name,
                                                             lambda df, os: df.os == os)
        self.reduction_potentials = self._get_values(self.vr, DataFiles.VR.column_name, lambda df, os: df.n == os)

    def _initialize_structure_analysis(self) -> List[Dict[str, int]]:
        """
        Initializes the structure analysis.

        Returns:
            A list of coordination number dictionaries.

        Examples:
            # TODO: Add examples
        """
        if self._cn_dicts_initialized:
            return self.cn_dicts
        self.structure.add_oxidation_state_by_guess()
        vacancy_generator = VacancyGenerator()
        vacancies = vacancy_generator.get_defects(self.structure)
        indices = [v.defect_site_index for v in vacancies if v.site.specie.symbol == SPECIES_SYMBOL]
        self.cn_dicts = [self.nn_finder.get_cn_dict(self.structure, i, use_weights=self.use_weights) for i in indices]
        self._cn_dicts_initialized = True
        return self.cn_dicts

    def _get_values(self, dataframe: pd.DataFrame, column_name: str, comparison: Callable[[pd.DataFrame, int], bool]) \
            -> List[Dict[str, float]]:
        """
        Gets the values from the dataframe.

        Args:
            dataframe: The dataframe.
            column_name: The column name.
            comparison: The comparison function.

        Returns:
            A list of dictionaries.

        Examples:
            # TODO: Add examples
        """
        return [
            {
                species_string: dataframe.loc[
                    (dataframe.elem == symbol) & comparison(dataframe, oxidation_state), column_name].values[0]
                for species_string, (species, symbol, oxidation_state) in
                ((ss, self._parse_species_string(ss)) for ss in cn_dict.keys())
                if not dataframe.loc[(dataframe.elem == symbol)
                                     & comparison(dataframe, oxidation_state), column_name].empty
            }
            for cn_dict in self._initialize_structure_analysis()
        ]

    def visualize(self):
        """
        Visualizes the crystal structure using ASE's view function.

        Examples:
            # TODO: Add examples
        """
        atoms = AseAtomsAdaptor.get_atoms(self.structure)
        view(atoms)

    def __repr__(self) -> str:
        """
        Returns the string representation of the Crystal object.

        Returns:
            The string representation of the Crystal object.

        Examples:
            # TODO: Add examples
        """
        return f"Crystal({self.structure})"
