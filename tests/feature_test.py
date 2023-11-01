import unittest
from deftpy.crystal_analysis import Crystal
from pathlib import Path
import re

class TestCrystalVacancies(unittest.TestCase):
    def load_poscar_data(self, compound_name):
        poscar_file_path = Path(__file__).parent / 'data' / 'test_files' / f'OQMD_{compound_name}_POSCAR.txt'
        with open(poscar_file_path, 'r') as poscar_file:
            poscar_data = poscar_file.read()
        return poscar_data

    def determine_compound_type(self, compound_name):
        # Using a simple heuristic based on the number of unique elements in the compound
        unique_elements = set(re.findall(r'[A-Z][a-z]*', compound_name))
        if len(unique_elements) == 2:
            return 'binary'
        elif len(unique_elements) == 3:
            return 'ternary'
        else:
            return 'other'

    def run_test(self, compound_name, element_symbol, expected_eb_value, expected_vr_value):
        poscar_data = self.load_poscar_data(compound_name)
        compound_type = self.determine_compound_type(compound_name)

        print(f"Testing {compound_name} ({compound_type} compound) for element {element_symbol}...")

        crystal = Crystal(poscar_string=poscar_data)
        eb_value = crystal.bond_dissociation_enthalpies[0][element_symbol]
        vr_value = crystal.reduction_potentials[0][element_symbol]

        self.assertAlmostEqual(eb_value, expected_eb_value, msg=f"Test failed for {compound_name} on element {element_symbol} for Eb value")
        self.assertAlmostEqual(vr_value, expected_vr_value, msg=f"Test failed for {compound_name} on element {element_symbol} for Vr value")

    def test_generic_compound(self):
        # Example usage:
        # self.run_test("CaTiO3", "Ca", expected_eb_value_for_Ca_in_CaTiO3, expected_vr_value_for_Ca_in_CaTiO3)
        # self.run_test("CaAlO3", "Al", expected_eb_value_for_Al_in_CaAlO3, expected_vr_value_for_Al_in_CaAlO3)
        pass

if __name__ == '__main__':
    unittest.main()
