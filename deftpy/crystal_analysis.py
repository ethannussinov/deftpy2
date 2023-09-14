# Standard Library Imports
import os

# Third-Party Imports
import pandas as pd
from pymatgen.analysis.defects.generators import VacancyGenerator
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core import Species
from pymatgen.core.structure import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view

class Crystal:
    """
    A class representing a crystal structure.

    Parameters:
        filepath (str): The filepath to the crystal file.

    Attributes:
        structure (Structure): The pymatgen Structure object representing the crystal structure.
        eb (DataFrame): A pandas DataFrame of bond dissociation enthalpies, read from Eb.csv.
        vr (DataFrame): A pandas DataFrame of reduction potentials, read from Vr.csv.
        bond_dissociation_enthalpies (dict): A dictionary of bond dissociation enthalpies for the crystal structure.
        reduction_potentials (dict): A dictionary of reduction potentials for the crystal structure.
    
    Methods:
        __init__(self, filepath: str): Initializes a new Crystal object with the given filepath.
        _initialize_structure_analysis(self): Initializes the structure analysis.
        visualize(self): Opens a visualization window of the crystal structure using ASE's view function.
        __repr__(self): Returns a string representation of the Crystal object.
    """

    # Define class level variables for important file paths
    eb_path = "../data/Eb.csv"
    vr_path = "../data/Vr.csv"

    def __init__(self, filepath):
        self.structure = Structure.from_file(filepath)
        package_dir = os.path.dirname(os.path.abspath(__file__))
        self.eb = pd.read_csv(os.path.join(package_dir, self.eb_path))
        self.vr = pd.read_csv(os.path.join(package_dir, self.vr_path))
        self._cn_dicts_initialized = False 
        self.bond_dissociation_enthalpies = self.get_bond_dissociation_enthalpies()
        self.reduction_potentials = self.get_reduction_potentials()

    def _initialize_structure_analysis(self):
        cn_dicts = []
        
        if not self._cn_dicts_initialized:
            structure = self.structure

            # Add oxidation states to structure
            structure.add_oxidation_state_by_guess()

            # Get oxygen vacancies
            vacancy_generator = VacancyGenerator()
            vacancies = vacancy_generator.get_defects(structure)
            indices = [vacancy.defect_site_index for vacancy in vacancies if vacancy.site.specie.symbol == "O"]

            # Get coordination numbers
            crystal_nn = CrystalNN()
            cn_dicts = [crystal_nn.get_cn_dict(structure, i) for i in indices]

            self._cn_dicts_initialized = True

        return cn_dicts

    def get_bond_dissociation_enthalpies(self):
        """
        Get the crystal bond dissociation enthalpies for a given structure.

        Returns:
            dict: A dictionary of crystal bond dissociation enthalpies.
        """
        #initalize dictionaries
        cn_dicts = self._initialize_structure_analysis()

        # Get crystal bond dissociation enthalpies
        crystal_bond_dissociation_enthalpies = {}
        for cn_dict in cn_dicts:
            species = Species(list(cn_dict.keys())[0])
            symbol = species.symbol
            oxidation_state = species.oxi_state

            try:
                crystal_bond_dissociation_enthalpy = self.eb.loc[(self.eb.elem == symbol) & (self.eb.os == oxidation_state), "Eb"].values[0]
                crystal_bond_dissociation_enthalpies[str(species)] = crystal_bond_dissociation_enthalpy
            except IndexError:
                # error handling needed here (TODO)
                pass


        return crystal_bond_dissociation_enthalpies
    
    def get_reduction_potentials(self):
        """
        
        Get the crystal reduction potentials for a given structure.

        Returns:
            dict: A dictionary of crystal reduction potentials.
        """
        #initalize dictionaries
        cn_dicts = self._initialize_structure_analysis()

        # Get crystal reduction potentials
        crystal_reduction_potentials = {}
        for cn_dict in cn_dicts:
            species = Species(list(cn_dict.keys())[0])
            symbol = species.symbol
            oxidation_state = species.oxi_state

            try:
                crystal_reduction_potential = self.vr.loc[(self.vr.elem == symbol) & (self.vr.os == oxidation_state), "Vr"].values[0]
                crystal_reduction_potentials[str(species)] = crystal_reduction_potential
            except IndexError:
                # error handling needed here (TODO)
                pass


        return crystal_reduction_potentials

    def visualize(self):
        """
        Opens a visualization window of the crystal structure using ASE's view function.
        """
        atoms = AseAtomsAdaptor.get_atoms(self.structure)
        view(atoms)

    def __repr__(self):
        return f"Crystal({self.structure})"
    